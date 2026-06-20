/*
 * ML-KEM NTT 구현 (FIPS 203 기준)
 * 파라미터: N=256, Q=3329
 *
 * 구조:
 *  - 7레이어 NTT (len=128~2)
 *  - base_mul: 마지막 레이어 degree-1 다항식 곱셈
 *  - INTT: 7레이어 역변환 + 128^(-1) 스케일링
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define KYBER_Q    3329
#define KYBER_N    256
#define KYBER_NINV 3303   /* 128^(-1) mod 3329 */

/* FIPS 203 zetas: 17^(BitRev7(k)) mod 3329 */
static const int16_t zetas[128] = {
       1, 1729, 2580, 3289, 2642,  630, 1897,  848,
    1062, 1919,  193,  797, 2786, 3260,  569, 1746,
     296, 2447, 1339, 1476, 3046,   56, 2240, 1333,
    1426, 2094,  535, 2882, 2393, 2879, 1974,  821,
     289,  331, 3253, 1756, 1197, 2304, 2277, 2055,
     650, 1977, 2513,  632, 2865,   33, 1320, 1915,
    2319, 1435,  807,  452, 1438, 2868, 1534, 2402,
    2647, 2617, 1481,  648, 2474, 3110, 1227,  910,
      17, 2761,  583, 2649, 1637,  723, 2288, 1100,
    1409, 2662, 3281,  233,  756, 2156, 3015, 3050,
    1703, 1651, 2789, 1789, 1847,  952, 1461, 2687,
     939, 2308, 2437, 2388,  733, 2337,  268,  641,
    1584, 2298, 2037, 3220,  375, 2549, 2090, 1645,
    1063,  319, 2773,  757, 2099,  561, 2466, 2594,
    2804, 1092,  403, 1026, 1143, 2150, 2775,  886,
    1722, 1212, 1874, 1029, 2110, 2935,  885, 2154,
};

/* ── 기본 연산 ──────────────────────────── */

static int32_t fqmul(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * b) % KYBER_Q);
}

static int32_t normalize(int32_t a) {
    a %= KYBER_Q;
    if (a < 0) a += KYBER_Q;
    return a;
}

/* ── NTT (7레이어, len=128~2) ───────────── */

void ntt(int32_t a[KYBER_N]) {
    int32_t k = 1;
    for (int32_t length = 128; length >= 2; length >>= 1) {
        for (int32_t start = 0; start < KYBER_N; start += 2 * length) {
            int32_t zeta = zetas[k++];
            for (int32_t j = start; j < start + length; j++) {
                int32_t t      = fqmul(zeta, a[j + length]);
                a[j + length]  = normalize(a[j] - t);
                a[j]           = normalize(a[j] + t);
            }
        }
    }
}

/* ── INTT (7레이어, len=2~128) ──────────── */

void inv_ntt(int32_t a[KYBER_N]) {
    int32_t k = 127;
    for (int32_t length = 2; length <= 128; length <<= 1) {
        for (int32_t start = 0; start < KYBER_N; start += 2 * length) {
            int32_t zeta = zetas[k--];
            for (int32_t j = start; j < start + length; j++) {
                int32_t t      = a[j];
                a[j]           = normalize(t + a[j + length]);
                a[j + length]  = fqmul(zeta, normalize(a[j + length] - t));
            }
        }
    }
    /* n^(-1) 스케일링 */
    for (int32_t j = 0; j < KYBER_N; j++)
        a[j] = fqmul(KYBER_NINV, a[j]);
}

/* ── base_mul: (a0+a1*x)*(b0+b1*x) mod (x^2 - zeta) ── */

static void base_mul(int32_t a0, int32_t a1,
                     int32_t b0, int32_t b1,
                     int32_t zeta,
                     int32_t *r0, int32_t *r1) {
    *r0 = normalize(fqmul(a0, b0) + fqmul(fqmul(a1, b1), zeta));
    *r1 = normalize(fqmul(a0, b1) + fqmul(a1, b0));
}

/* ── 다항식 곱셈: result = f * g mod (x^256+1, Q) ── */

void poly_mul(int32_t f[KYBER_N],
              int32_t g[KYBER_N],
              int32_t result[KYBER_N]) {
    int32_t fa[KYBER_N], ga[KYBER_N];

    for (int32_t i = 0; i < KYBER_N; i++) { fa[i] = f[i]; ga[i] = g[i]; }

    ntt(fa);
    ntt(ga);

    /* base_mul: 4원소씩 처리 (128쌍) */
    for (int32_t i = 0; i < KYBER_N; i += 4) {
        int32_t zeta = zetas[64 + i / 4];
        base_mul(fa[i],   fa[i+1], ga[i],   ga[i+1],
                 zeta, &result[i], &result[i+1]);
        base_mul(fa[i+2], fa[i+3], ga[i+2], ga[i+3],
                 normalize(-zeta), &result[i+2], &result[i+3]);
    }

    inv_ntt(result);
}

/* ── 테스트 ─────────────────────────────── */

static void test_ntt_intt(void) {
    int32_t a[KYBER_N], orig[KYBER_N];
    int pass = 1;

    for (int i = 0; i < KYBER_N; i++)
        a[i] = orig[i] = (i + 1) % KYBER_Q;

    ntt(a);
    inv_ntt(a);

    for (int i = 0; i < KYBER_N; i++) {
        if (a[i] != orig[i]) {
            pass = 0;
            printf("  FAIL [%d]: got %d expected %d\n", i, a[i], orig[i]);
            break;
        }
    }
    printf("[Test 1] NTT -> INTT 복원: %s\n", pass ? "PASS" : "FAIL");
}

static void test_poly_mul_simple(void) {
    int32_t f[KYBER_N], g[KYBER_N], r[KYBER_N];
    int pass = 1;

    memset(f, 0, sizeof f);
    memset(g, 0, sizeof g);
    f[0] = 1; f[1] = 1;   /* 1 + x */
    g[0] = 1; g[1] = 1;   /* 1 + x */

    poly_mul(f, g, r);

    /* 기댓값: 1 + 2x + x^2 */
    if (r[0] != 1 || r[1] != 2 || r[2] != 1) pass = 0;
    for (int i = 3; i < KYBER_N; i++) if (r[i]) { pass = 0; break; }

    printf("[Test 2] (1+x)^2 = 1+2x+x^2: %s\n", pass ? "PASS" : "FAIL");
    printf("  [0]=%d [1]=%d [2]=%d [3]=%d\n", r[0], r[1], r[2], r[3]);
}

static void test_const_mul(void) {
    int32_t f[KYBER_N], g[KYBER_N], r[KYBER_N];

    memset(f, 0, sizeof f);
    memset(g, 0, sizeof g);
    f[0] = 3; g[0] = 7;

    poly_mul(f, g, r);

    printf("[Test 3] 3 * 7 = 21: %s\n", r[0] == 21 ? "PASS" : "FAIL");
    printf("  result[0] = %d\n", r[0]);
}

/* x^255 * x = x^256 = -1 mod (x^256+1)  →  result[0] = Q-1 */
static void test_wrap_around(void) {
    int32_t f[KYBER_N], g[KYBER_N], r[KYBER_N];
    int pass = 1;

    memset(f, 0, sizeof f);
    memset(g, 0, sizeof g);
    f[255] = 1;   /* x^255 */
    g[1]   = 1;   /* x     */

    poly_mul(f, g, r);

    if (r[0] != KYBER_Q - 1) pass = 0;
    for (int i = 1; i < KYBER_N; i++) if (r[i]) { pass = 0; break; }

    printf("[Test 4] x^255 * x = -1 (mod x^256+1): %s\n", pass ? "PASS" : "FAIL");
    printf("  result[0] = %d (expected %d)\n", r[0], KYBER_Q - 1);
}

int main(void) {
    printf("=== ML-KEM NTT 구현 (FIPS 203) ===\n");
    printf("N=%d  Q=%d  128^-1 mod Q=%d\n\n", KYBER_N, KYBER_Q, KYBER_NINV);

    test_ntt_intt();
    test_poly_mul_simple();
    test_const_mul();
    test_wrap_around();

    printf("\n=== 완료 ===\n");
    return 0;
}
