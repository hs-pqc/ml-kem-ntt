# ML-KEM NTT Implementation (FIPS 203)

FIPS 203 표준을 따르는 ML-KEM(CRYSTALS-Kyber)의
Number Theoretic Transform(NTT) C 구현체.

Barrett reduction과 Montgomery reduction 두 가지 방식을 모두 구현하고
성능을 비교한다.

---

## Parameters

| Parameter | Value | 근거 |
|---|---|---|
| N (polynomial degree) | 256 | FIPS 203 §4.1 |
| Q (modulus) | 3329 | NTT-friendly 소수 (2¹² < Q < 2¹³) |
| 128⁻¹ mod Q | 3303 | INTT 스케일링 상수 |
| Primitive root ω | 17 | ord(17) = 256 mod 3329 |
| Barrett v | 20159 | ⌊2²⁶ / Q⌋ |
| Montgomery R | 2¹⁶ | R > Q, gcd(R, Q) = 1 |
| Montgomery QINV | 62209 | Q⁻¹ mod 2¹⁶ |

---

## 구현 구조

ntt_mlkem.c      — Barrett reduction 기반 NTT (int32_t)

ntt_mlkem_mont.c — Montgomery reduction 기반 NTT (int16_t)

benchmark.c      — 두 방식 성능 비교

### 함수 목록

| 함수 | 설명 | FIPS 203 대응 |
|---|---|---|
| `ntt()` | 7-layer forward NTT (len=128→2) | Algorithm 9 |
| `inv_ntt()` | 7-layer INTT + 128⁻¹ 스케일링 | Algorithm 10 |
| `base_mul()` | degree-1 다항식 곱 (최종 레이어) | Algorithm 11 |
| `poly_mul()` | 전체 다항식 곱 mod (x²⁵⁶+1, Q) | — |

---

## 핵심 구현 세부사항

### NTT가 7레이어인 이유

ML-KEM NTT는 len=2에서 멈춤. 마지막 레이어(len=1)는
`base_mul()`로 별도 처리:
(a₀ + a₁x)(b₀ + b₁x) mod (x² - ζ) = (a₀b₀ + a₁b₁ζ) + (a₀b₁ + a₁b₀)x

이유: x²⁵⁶ + 1 = ∏(x² - ζ²ⁱ⁺¹) mod Q (i = 0..127)
→ 256차 환을 128개의 2차 환의 곱으로 분해.

### Barrett Reduction

나눗셈을 곱셈 + 비트시프트로 대체:

```c
// v = floor(2^26 / Q) = 20159
static int32_t barrett_reduce(int32_t a) {
    int32_t t = (int32_t)(((int64_t)20159 * a) >> 26);
    return a - t * 3329;
}
```

### Montgomery Reduction

R = 2¹⁶, QINV = Q⁻¹ mod R = 62209:

```c
static int16_t montgomery_reduce(int32_t a) {
    int16_t t = (int16_t)((uint16_t)(a & 0xFFFF) * (uint16_t)62209);
    return (int16_t)((a - (int32_t)t * 3329) >> 16);
}
```

---

## Barrett vs Montgomery 성능 비교

gcc -O2, Docker(Linux) 환경 기준 (단위: ns/op, 100만 회 평균):

| 연산 | Barrett | Montgomery | 비율 |
|---|---|---|---|
| `ntt()` | 4,039 ns | 1,928 ns | **2.10x** |
| `inv_ntt()` | 3,833 ns | 1,874 ns | **2.05x** |

**Montgomery가 약 2배 빠름.**

이유: Barrett는 64-bit 정수 나눗셈 근사(v=⌊2²⁶/Q⌋)를 매 계수마다 수행.
Montgomery는 16-bit 곱셈으로 reduction을 대체하고 int16_t 연산으로
메모리 대역폭도 절반으로 줄어듦.
단, Montgomery는 사전에 Montgomery form 변환 오버헤드가 있어서
단발성 연산에서는 Barrett가 유리할 수 있음.

* 측정 환경: Docker gcc:latest, -O2, CLOCK_MONOTONIC

---

## Build & Test

```bash
# Barrett 버전
gcc -O2 -o ntt_mlkem ntt_mlkem.c && ./ntt_mlkem

# Montgomery 버전
gcc -O2 -o ntt_mlkem_mont ntt_mlkem_mont.c && ./ntt_mlkem_mont

# 벤치마크
gcc -O2 -o benchmark benchmark.c && ./benchmark
```

### 테스트 결과

[Test 1] NTT -> INTT 복원: PASS

[Test 2] (1+x)^2 = 1+2x+x^2: PASS

[Test 3] 3 * 7 = 21: PASS

[Test 4] x^255 * x = -1 (mod x^256+1): PASS

---

## 보안 파라미터와의 연결

ML-KEM 세 파라미터셋 모두 동일한 NTT를 사용:

| 파라미터셋 | k | 보안 레벨 | NTT 호출 수 (KeyGen 기준) |
|---|---|---|---|
| ML-KEM-512 | 2 | Level 1 (128-bit) | k²+k = 6회 |
| ML-KEM-768 | 3 | Level 3 (192-bit) | k²+k = 12회 |
| ML-KEM-1024 | 4 | Level 5 (256-bit) | k²+k = 20회 |

→ NTT 성능이 전체 ML-KEM 성능의 핵심 병목.

---

## 향후 작업

- [ ] AVX2 SIMD 최적화
- [ ] Constant-time 검증
- [ ] ML-KEM 전체 KeyGen/Encaps/Decaps 연결
- [ ] lattice-estimator 보안 추정과 연계 분석

---

## References

- [NIST FIPS 203](https://csrc.nist.gov/pubs/fips/203/final) — ML-KEM Standard
- [pq-crystals/kyber](https://github.com/pq-crystals/kyber) — Reference implementation
- 관련 실험: [lattice-estimator-mlkem](https://github.com/hs-pqc/lattice-estimator-mlkem)
