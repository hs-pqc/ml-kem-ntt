/*
 * ML-KEM NTT - Montgomery Reduction 버전 (FIPS 203)
 *
 * Barrett vs Montgomery:
 *   Barrett:    v = floor(2^26/Q), 나눗셈 근사
 *   Montgomery: QINV = Q^-1 mod R, 나눗셈 정확 대체
 *
 * 핵심 상수:
 *   R      = 2^16 = 65536
 *   QINV   = 3329^(-1) mod 2^16 = 62209
 *   zetas  = zeta_real * R mod Q  (Montgomery form)
 *   F_NTT  = R * 128^(-1) mod Q = 512   (NTT->INTT 복원)
 *   F_MUL  = R^2 * 128^(-1) mod Q = 1441 (poly_mul)
 *
 * 주의: INTT 덧셈 후 계수가 int16_t 범위를 초과하므로
 *       barrett_reduce로 중간 reduction 필요
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define KYBER_Q   3329
#define KYBER_N   256
#define MONT_QINV 62209
#define BARRETT_V 20159   /* floor(2^26 / Q) */
#define F_NTT     512     /* R * 128^(-1) mod Q */
#define F_MUL     1441    /* R^2 * 128^(-1) mod Q */

/* zetas: zeta_real * R mod Q */
static const int16_t zetas[128] = {
    2285, 2571, 2970, 1812, 1493, 1422,  287,  202,
    3158,  622, 1577,  182,  962, 2127, 1855, 1468,
     573, 2004,  264,  383, 2500, 1458, 1727, 3199,
    2648, 1017,  732,  608, 1787,  411, 3124, 1758,
    1223,  652, 2777, 1015, 2036, 1491, 3047, 1785,
     516, 3321, 3009, 2663, 1711, 2167,  126, 1469,
    2476, 3239, 3058,  830,  107, 1908, 3082, 2378,
    2931,  961, 1821, 2604,  448, 2264,  677, 2054,
    2226,  430,  555,  843, 2078,  871, 1550,  105,
     422,  587,  177, 3094, 3038, 2869, 1574, 1653,
    3083,  778, 1159, 3182, 2552, 1483, 2727, 1119,
    1739,  644, 2457,  349,  418,  329, 3173, 3254,
     817, 1097,  603,  610, 1322, 2044, 1864,  384,
    2114, 3193, 1218, 1994, 2455,  220, 2142, 1670,
    2144, 1799, 2051,  794, 1819, 2475, 2459,  478,
    3221, 3021,  996,  991,  958, 1869, 1522, 1628,
};

/* Montgomery reduction: a * R^(-1) mod Q */
static int16_t montgomery_reduce(int32_t a) {
    int16_t t = (int16_t)((uint16_t)(a & 0xFFFF) * (uint16_t)MONT_QINV);
    return (int16_t)((a - (int32_t)t * KYBER_Q) >> 16);
}

/* fqmul: a * b * R^(-1) mod Q */
static int16_t fqmul(int16_t a, int16_t b) {
    return montgomery_reduce((int32_t)a * (int32_t)b);
}

/* Barrett reduction: 중간 계수 범위 제어용 */
static int16_t barrett_reduce(int32_t a) {
    int32_t t = ((int64_t)BARRETT_V * a) >> 26;
    return (int16_t)(a - t * KYBER_Q);
}

/* [0, Q) 정규화 */
static int16_t normalize(int32_t a) {
    a %= KYBER_Q;
    if (a < 0) a += KYBER_Q;
    return (int16_t)a;
}

/* NTT (7레이어) */
void ntt(int16_t a[KYBER_N]) {
    int32_t k = 1;
    for (int32_t len = 128; len >= 2; len >>= 1) {
        for (int32_t start = 0; start < KYBER_N; start += 2*len) {
            int16_t zeta = zetas[k++];
            for (int32_t j = start; j < start+len; j++) {
                int16_t t = fqmul(zeta, a[j+len]);
                a[j+len]  = a[j] - t;
                a[j]      = a[j] + t;
            }
        }
    }
}

/* INTT 내부 (f: 스케일 상수) */
static void _intt(int16_t a[KYBER_N], int16_t f) {
    int32_t k = 127;
    for (int32_t len = 2; len <= 128; len <<= 1) {
        for (int32_t start = 0; start < KYBER_N; start += 2*len) {
            int16_t zeta = zetas[k--];
            for (int32_t j = start; j < start+len; j++) {
                int16_t t = a[j];
                /* 덧셈 후 범위 초과 방지: barrett_reduce 적용 */
                a[j]     = barrett_reduce((int32_t)t + a[j+len]);
                a[j+len] = fqmul(zeta, a[j+len] - t);
            }
        }
    }
    for (int32_t j = 0; j < KYBER_N; j++)
        a[j] = normalize(fqmul(f, a[j]));
}

/* NTT -> INTT 복원 */
void inv_ntt(int16_t a[KYBER_N]) {
    _intt(a, (int16_t)F_NTT);
}

/* poly_mul 내부용 */
static void inv_ntt_polymul(int16_t a[KYBER_N]) {
    _intt(a, (int16_t)F_MUL);
}

/* base_mul: (a0+a1*x)*(b0+b1*x) mod (x^2 - zeta) */
static void base_mul(int16_t a0, int16_t a1,
                     int16_t b0, int16_t b1,
                     int16_t zeta,
                     int16_t *r0, int16_t *r1) {
    *r0 = normalize(fqmul(a0,b0) + fqmul(fqmul(a1,b1), zeta));
    *r1 = normalize(fqmul(a0,b1) + fqmul(a1,b0));
}

/* 다항식 곱셈: result = f * g mod (x^256+1, Q) */
void poly_mul(int16_t f[KYBER_N],
              int16_t g[KYBER_N],
              int16_t result[KYBER_N]) {
    int16_t fa[KYBER_N], ga[KYBER_N];
    for (int i = 0; i < KYBER_N; i++) { fa[i]=f[i]; ga[i]=g[i]; }
    ntt(fa); ntt(ga);
    for (int i = 0; i < KYBER_N; i += 4) {
        int16_t zeta = zetas[64 + i/4];
        base_mul(fa[i],  fa[i+1], ga[i],  ga[i+1],
                 zeta, &result[i], &result[i+1]);
        base_mul(fa[i+2],fa[i+3], ga[i+2],ga[i+3],
                 normalize(-zeta), &result[i+2], &result[i+3]);
    }
    inv_ntt_polymul(result);
}

/* ── 테스트 ──────────────────────────────── */

static void test_ntt_intt(void) {
    int16_t a[KYBER_N], orig[KYBER_N];
    int pass = 1;
    for (int i = 0; i < KYBER_N; i++) a[i] = orig[i] = (i+1) % KYBER_Q;
    ntt(a); inv_ntt(a);
    for (int i = 0; i < KYBER_N; i++) {
        if (normalize(a[i]) != normalize(orig[i])) {
            pass = 0;
            printf("  FAIL [%d]: got %d expected %d\n",
                   i, normalize(a[i]), normalize(orig[i]));
            break;
        }
    }
    printf("[Test 1] NTT -> INTT 복원: %s\n", pass?"PASS":"FAIL");
}

static void test_poly_mul_simple(void) {
    int16_t f[KYBER_N], g[KYBER_N], r[KYBER_N];
    int pass = 1;
    memset(f,0,sizeof f); memset(g,0,sizeof g);
    f[0]=1; f[1]=1; g[0]=1; g[1]=1;
    poly_mul(f,g,r);
    if (normalize(r[0])!=1||normalize(r[1])!=2||normalize(r[2])!=1) pass=0;
    for (int i=3;i<KYBER_N;i++) if (normalize(r[i])) {pass=0;break;}
    printf("[Test 2] (1+x)^2 = 1+2x+x^2: %s\n", pass?"PASS":"FAIL");
    printf("  [0]=%d [1]=%d [2]=%d [3]=%d\n",
           normalize(r[0]),normalize(r[1]),normalize(r[2]),normalize(r[3]));
}

static void test_const_mul(void) {
    int16_t f[KYBER_N], g[KYBER_N], r[KYBER_N];
    memset(f,0,sizeof f); memset(g,0,sizeof g);
    f[0]=3; g[0]=7;
    poly_mul(f,g,r);
    printf("[Test 3] 3*7=21: %s\n", normalize(r[0])==21?"PASS":"FAIL");
    printf("  result[0]=%d\n", normalize(r[0]));
}

static void test_wrap_around(void) {
    int16_t f[KYBER_N], g[KYBER_N], r[KYBER_N];
    int pass=1;
    memset(f,0,sizeof f); memset(g,0,sizeof g);
    f[255]=1; g[1]=1;
    poly_mul(f,g,r);
    if (normalize(r[0])!=KYBER_Q-1) pass=0;
    for (int i=1;i<KYBER_N;i++) if (normalize(r[i])) {pass=0;break;}
    printf("[Test 4] x^255*x=-1 (wrap): %s\n", pass?"PASS":"FAIL");
    printf("  result[0]=%d (expected %d)\n", normalize(r[0]), KYBER_Q-1);
}

int main(void) {
    printf("=== ML-KEM NTT (Montgomery Reduction) ===\n");
    printf("N=%d  Q=%d  QINV=%d  F_NTT=%d  F_MUL=%d\n\n",
           KYBER_N, KYBER_Q, MONT_QINV, F_NTT, F_MUL);
    test_ntt_intt();
    test_poly_mul_simple();
    test_const_mul();
    test_wrap_around();
    printf("\n=== 완료 ===\n");
    return 0;
}
