#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define KYBER_N    256
#define KYBER_Q    3329
#define KYBER_NINV 3303
#define ITERATIONS 1000000

/* ── zetas (Barrett 버전) ── */
static const int16_t zetas_b[128] = {
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

/* ── zetas (Montgomery 버전) ── */
static const int16_t zetas_m[128] = {
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

/* ════════════════════════════════════════
   Barrett 구현 (int32_t)
   ════════════════════════════════════════ */
static int32_t fqmul_b(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * b) % KYBER_Q);
}
static int32_t norm_b(int32_t a) {
    a %= KYBER_Q; if (a < 0) a += KYBER_Q; return a;
}
static void ntt_b(int32_t a[KYBER_N]) {
    int32_t k = 1;
    for (int32_t len = 128; len >= 2; len >>= 1)
        for (int32_t s = 0; s < KYBER_N; s += 2*len) {
            int32_t z = zetas_b[k++];
            for (int32_t j = s; j < s+len; j++) {
                int32_t t = fqmul_b(z, a[j+len]);
                a[j+len]  = norm_b(a[j] - t);
                a[j]      = norm_b(a[j] + t);
            }
        }
}
static void inv_ntt_b(int32_t a[KYBER_N]) {
    int32_t k = 127;
    for (int32_t len = 2; len <= 128; len <<= 1)
        for (int32_t s = 0; s < KYBER_N; s += 2*len) {
            int32_t z = zetas_b[k--];
            for (int32_t j = s; j < s+len; j++) {
                int32_t t = a[j];
                a[j]      = norm_b(t + a[j+len]);
                a[j+len]  = fqmul_b(z, norm_b(a[j+len] - t));
            }
        }
    for (int32_t j = 0; j < KYBER_N; j++)
        a[j] = fqmul_b(KYBER_NINV, a[j]);
}

/* ════════════════════════════════════════
   Montgomery 구현 (int16_t)
   ════════════════════════════════════════ */
#define MONT_QINV 62209
#define F_NTT     512
static int16_t mont_reduce(int32_t a) {
    int16_t t = (int16_t)((uint16_t)(a & 0xFFFF) * (uint16_t)MONT_QINV);
    return (int16_t)((a - (int32_t)t * KYBER_Q) >> 16);
}
static int16_t fqmul_m(int16_t a, int16_t b) {
    return mont_reduce((int32_t)a * b);
}
static int16_t barrett_reduce(int32_t a) {
    int32_t t = ((int64_t)20159 * a) >> 26;
    return (int16_t)(a - t * KYBER_Q);
}
static int16_t norm_m(int32_t a) {
    a %= KYBER_Q; if (a < 0) a += KYBER_Q; return (int16_t)a;
}
static void ntt_m(int16_t a[KYBER_N]) {
    int32_t k = 1;
    for (int32_t len = 128; len >= 2; len >>= 1)
        for (int32_t s = 0; s < KYBER_N; s += 2*len) {
            int16_t z = zetas_m[k++];
            for (int32_t j = s; j < s+len; j++) {
                int16_t t = fqmul_m(z, a[j+len]);
                a[j+len]  = a[j] - t;
                a[j]      = a[j] + t;
            }
        }
}
static void inv_ntt_m(int16_t a[KYBER_N]) {
    int32_t k = 127;
    for (int32_t len = 2; len <= 128; len <<= 1)
        for (int32_t s = 0; s < KYBER_N; s += 2*len) {
            int16_t z = zetas_m[k--];
            for (int32_t j = s; j < s+len; j++) {
                int16_t t = a[j];
                a[j]      = barrett_reduce((int32_t)t + a[j+len]);
                a[j+len]  = fqmul_m(z, a[j+len] - t);
            }
        }
    for (int32_t j = 0; j < KYBER_N; j++)
        a[j] = norm_m(fqmul_m((int16_t)F_NTT, a[j]));
}

/* ════════════════════════════════════════
   타이머 & 벤치마크
   ════════════════════════════════════════ */
static double now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e9 + ts.tv_nsec;
}

static double bench_ntt_b(void) {
    int32_t a[KYBER_N];
    for (int i = 0; i < KYBER_N; i++) a[i] = (i+1) % KYBER_Q;
    /* 워밍업 */
    for (int i = 0; i < 10000; i++) { int32_t t[KYBER_N]; memcpy(t,a,sizeof t); ntt_b(t); }
    double s = now_ns();
    for (int i = 0; i < ITERATIONS; i++) { int32_t t[KYBER_N]; memcpy(t,a,sizeof t); ntt_b(t); }
    return (now_ns()-s)/ITERATIONS;
}
static double bench_inv_ntt_b(void) {
    int32_t a[KYBER_N];
    for (int i = 0; i < KYBER_N; i++) a[i] = (i+1) % KYBER_Q;
    ntt_b(a);
    double s = now_ns();
    for (int i = 0; i < ITERATIONS; i++) { int32_t t[KYBER_N]; memcpy(t,a,sizeof t); inv_ntt_b(t); }
    return (now_ns()-s)/ITERATIONS;
}
static double bench_ntt_m(void) {
    int16_t a[KYBER_N];
    for (int i = 0; i < KYBER_N; i++) a[i] = (i+1) % KYBER_Q;
    for (int i = 0; i < 10000; i++) { int16_t t[KYBER_N]; memcpy(t,a,sizeof t); ntt_m(t); }
    double s = now_ns();
    for (int i = 0; i < ITERATIONS; i++) { int16_t t[KYBER_N]; memcpy(t,a,sizeof t); ntt_m(t); }
    return (now_ns()-s)/ITERATIONS;
}
static double bench_inv_ntt_m(void) {
    int16_t a[KYBER_N];
    for (int i = 0; i < KYBER_N; i++) a[i] = (i+1) % KYBER_Q;
    ntt_m(a);
    double s = now_ns();
    for (int i = 0; i < ITERATIONS; i++) { int16_t t[KYBER_N]; memcpy(t,a,sizeof t); inv_ntt_m(t); }
    return (now_ns()-s)/ITERATIONS;
}

int main(void) {
    double b_ntt     = bench_ntt_b();
    double b_inv_ntt = bench_inv_ntt_b();
    double m_ntt     = bench_ntt_m();
    double m_inv_ntt = bench_inv_ntt_m();

    printf("=== ML-KEM NTT Benchmark (gcc -O2, %d iterations) ===\n\n", ITERATIONS);
    printf("%-15s %14s %14s %10s\n", "op", "Barrett(ns)", "Montgomery(ns)", "speedup");
    printf("%-15s %14.1f %14.1f %9.2fx\n", "ntt()",     b_ntt,     m_ntt,     b_ntt/m_ntt);
    printf("%-15s %14.1f %14.1f %9.2fx\n", "inv_ntt()", b_inv_ntt, m_inv_ntt, b_inv_ntt/m_inv_ntt);
    return 0;
}
