# ML-KEM NTT Implementation (FIPS 203)

Implementation of the Number Theoretic Transform (NTT) 
used in ML-KEM (CRYSTALS-Kyber), following the FIPS 203 standard.

## Parameters

| Parameter | Value |
|-----------|-------|
| N (polynomial degree) | 256 |
| Q (modulus) | 3329 |
| 128^(-1) mod Q | 3303 |
| Primitive root | 17 |

## Structure

- `ntt()` — 7-layer NTT (len=128 to len=2)
- `inv_ntt()` — 7-layer INTT + 128^(-1) scaling
- `base_mul()` — degree-1 polynomial multiplication (final layer)
- `poly_mul()` — full polynomial multiplication mod (x^256+1, Q)

## Key Implementation Details

**Barrett Reduction**  
Replaces expensive division with multiplication + bit shift.  
`v = floor(2^26 / Q) = 20159`

**Why 7 layers, not 8?**  
ML-KEM NTT stops at len=2. The final layer (len=1) 
is handled separately via `base_mul()`,  
which computes degree-1 polynomial products mod (x^2 - zeta).

**base_mul**  
Each pair of 4 coefficients is multiplied as:  
`(a0 + a1*x) * (b0 + b1*x) mod (x^2 - zeta)`

## Build & Run

```bash
gcc -O2 -o ntt_mlkem ntt_mlkem.c
./ntt_mlkem
```

## Test Results
```
[Test 1] NTT -> INTT 복원: PASS

[Test 2] (1+x)^2 = 1+2x+x^2: PASS

[Test 3] 3 * 7 = 21: PASS

[Test 4] x^255 * x = -1 (mod x^256+1): PASS
```
## References

- [NIST FIPS 203](https://csrc.nist.gov/pubs/fips/203/final) — ML-KEM Standard
- [pq-crystals/kyber](https://github.com/pq-crystals/kyber) — Reference implementation
