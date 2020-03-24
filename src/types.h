#include <math.h>
#include <cstdlib>
#include <map>
#include <vector>


extern "C"
{

#ifndef __FOUND_ROOT
#define ___FOUND_ROOT
#define FOUND_ROOT ((roots_found[i_root] == 1) || (roots_found[i_root] == -1))
#endif 

#ifndef __CONSTANTS
#define __CONSTANTS

// Default constants //
extern double CONST_G;
extern double CONST_G_P2;
extern double CONST_G_P3;
extern double CONST_C_LIGHT;
extern double CONST_C_LIGHT_P2;
extern double CONST_C_LIGHT_P4;
extern double CONST_C_LIGHT_P5;
extern double CONST_MSUN;
extern double CONST_R_SUN;
extern double CONST_L_SUN;
extern double CONST_KM_PER_S;
extern double CONST_PER_PC3;

// Default parameters //
extern double relative_tolerance;
extern double absolute_tolerance_eccentricity_vectors;
extern bool include_quadrupole_order_terms;
extern bool include_octupole_order_binary_pair_terms;
extern bool include_octupole_order_binary_triplet_terms;
extern bool include_hexadecupole_order_binary_pair_terms;
extern bool include_dotriacontupole_order_binary_pair_terms;
extern bool include_double_averaging_corrections;
extern int orbital_phases_random_seed;
extern double epsilon;

#define c_1div2             (double)    1.0/2.0
#define c_1div3             (double)    1.0/3.0
#define c_1div4             (double)    1.0/4.0
#define c_1div5             (double)    1.0/5.0
#define c_1div6             (double)    1.0/6.0
#define c_1div7             (double)    1.0/7.0
#define c_1div8             (double)    1.0/8.0
#define c_1div10            (double)    1.0/10.0
#define c_1div15            (double)    1.0/15.0
#define c_1div16            (double)    1.0/16.0
#define c_1div30            (double)    1.0/30.0
#define c_2div3             (double)    2.0/3.0
#define c_3div2             (double)    3.0/2.0
#define c_3div4             (double)    3.0/4.0
#define c_3div5             (double)    3.0/5.0
#define c_3div8             (double)    3.0/8.0
#define c_3div32            (double)    3.0/32.0
#define c_3div1024          (double)    3.0/1024.0
#define c_4div15            (double)    4.0/15.0
#define c_5div2             (double)    5.0/2.0
#define c_5div8             (double)    5.0/8.0
#define c_5div16            (double)    5.0/16.0
#define c_5div64            (double)    5.0/64.0
#define c_7div8             (double)    7.0/8.0
#define c_8div5             (double)    8.0/5.0
#define c_8div7             (double)    8.0/7.0
#define c_9div2             (double)    9.0/2.0
#define c_9div8             (double)    9.0/8.0
#define c_9div16            (double)    9.0/16.0
#define c_9div32            (double)    9.0/32.0
#define c_9div64            (double)    9.0/64.0
#define c_10div3            (double)    10.0/3.0
#define c_11div18           (double)    11.0/18.0
#define c_15div2            (double)    15.0/2.0
#define c_15div4            (double)    15.0/4.0
#define c_15div8            (double)    15.0/8.0
#define c_15div16           (double)    15.0/16.0
#define c_16div5            (double)    16.0/5.0
#define c_25div16           (double)    25.0/16.0
#define c_25div64           (double)    25.0/64.0
#define c_27div4            (double)    27.0/4.0
#define c_27div64           (double)    27.0/64.0
#define c_31div2            (double)    31.0/2.0
#define c_32div5            (double)    32.0/5.0
#define c_35div3            (double)    35.0/3.0
#define c_37div96           (double)    37.0/96.0
#define c_45div8            (double)    45.0/8.0
#define c_64div3            (double)    64.0/3.0
#define c_64div5            (double)    64.0/5.0
#define c_65div3            (double)    65.0/3.0
#define c_73div24           (double)    73.0/24.0
#define c_105div4096        (double)    105.0/4096.0
#define c_121div304         (double)    121.0/304.0
#define c_185div16          (double)    185.0/16.0
#define c_255div8           (double)    255.0/8.0
#define c_304div15          (double)    304.0/15.0
#endif


#ifndef __TABLES
#define __TABLES
#define MAX_ORDER (int) 5

#define TABLEWIDTH_A (int) 3
#define TABLELENGTH_A (int) 10
const double A_TABLE[TABLELENGTH_A][TABLEWIDTH_A] =
{{2, 0, -0.500000000000000}, {2, 2, 1.50000000000000}, {3, 
  1, -1.50000000000000}, {3, 3, 2.50000000000000}, {4, 0, 
  0.375000000000000}, {4, 2, -3.75000000000000}, {4, 4, 
  4.37500000000000}, {5, 1, 1.87500000000000}, {5, 
  3, -8.75000000000000}, {5, 5, 7.87500000000000}};

#define A_mn_table_max_order (int) 15
const double A_mn_table[15][15+1] = {{-0.500000000000000, 0, 1.50000000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0}, {0, -1.50000000000000, 0, 2.50000000000000, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0}, {0.375000000000000, 0, -3.75000000000000, 
  0, 4.37500000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 
  1.87500000000000, 0, -8.75000000000000, 0, 7.87500000000000, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0}, {-0.312500000000000, 0, 6.56250000000000, 
  0, -19.6875000000000, 0, 14.4375000000000, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, -2.18750000000000, 0, 19.6875000000000, 
  0, -43.3125000000000, 0, 26.8125000000000, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0.273437500000000, 0, -9.84375000000000, 0, 54.1406250000000, 
  0, -93.8437500000000, 0, 50.2734375000000, 0, 0, 0, 0, 0, 0, 0}, {0,
   2.46093750000000, 0, -36.0937500000000, 0, 140.765625000000, 
  0, -201.093750000000, 0, 94.9609375000000, 0, 0, 0, 0, 0, 
  0}, {-0.246093750000000, 0, 13.5351562500000, 0, -117.304687500000, 
  0, 351.914062500000, 0, -427.324218750000, 0, 180.425781250000, 0, 
  0, 0, 0, 0}, {0, -2.70703125000000, 0, 58.6523437500000, 
  0, -351.914062500000, 0, 854.648437500000, 0, -902.128906250000, 0, 
  344.449218750000, 0, 0, 0, 0}, {0.225585937500000, 
  0, -17.5957031250000, 0, 219.946289062500, 0, -997.089843750000, 0, 
  2029.79003906250, 0, -1894.47070312500, 0, 660.194335937500, 0, 0, 
  0}, {0, 2.93261718750000, 0, -87.9785156250000, 0, 747.817382812500,
   0, -2706.38671875000, 0, 4736.17675781250, 0, -3961.16601562500, 0,
   1269.60449218750, 0, 0}, {-0.209472656250000, 0, 21.9946289062500, 
  0, -373.908691406250, 0, 2368.08837890625, 0, -7104.26513671875, 0, 
  10893.2065429688, 0, -8252.42919921875, 0, 2448.52294921875, 
  0}, {0, -3.14208984375000, 0, 124.636230468750, 
  0, -1420.85302734375, 0, 7104.26513671875, 0, -18155.3442382813, 0, 
  24757.2875976563, 0, -17139.6606445313, 0, 4733.81103515625}};

#define TABLEWIDTH_B (int) 8
#define TABLELENGTH_B (int) 47
#define HIGHEST_POWER_ECCP2_IN_B_TABLE (int) 3
const double B_TABLE[TABLELENGTH_B][TABLEWIDTH_B] =
{{2, 0, 0, 0, 1.00000000000000, 1.50000000000000, 0, 0}, {2, 1, 1, 
  0, -2.00000000000000, -0.500000000000000, 0, 0}, {2, 2, 0, 0, 
  0.500000000000000, -0.500000000000000, 0, 0}, {2, 2, 0, 
  2, -0.500000000000000, 0, 0, 0}, {2, 2, 2, 0, 2.50000000000000, 0, 
  0, 0}, {3, 0, 0, 0, 1.00000000000000, 3.00000000000000, 
  0.375000000000000, 0}, {3, 1, 1, 
  0, -2.50000000000000, -1.87500000000000, 0, 0}, {3, 2, 0, 0, 
  0.500000000000000, -0.375000000000000, -0.125000000000000, 0}, {3, 
  2, 0, 2, -0.500000000000000, -0.125000000000000, 0, 0}, {3, 2, 2, 0,
   3.75000000000000, 0.625000000000000, 0, 0}, {3, 3, 1, 
  0, -1.87500000000000, 1.87500000000000, 0, 0}, {3, 3, 1, 2, 
  1.87500000000000, 0, 0, 0}, {3, 3, 3, 0, -4.37500000000000, 0, 0, 
  0}, {4, 0, 0, 0, 1.00000000000000, 5.00000000000000, 
  1.87500000000000, 0}, {4, 1, 1, 
  0, -3.00000000000000, -4.50000000000000, -0.375000000000000, 0}, {4,
   2, 0, 0, 0.500000000000000, -0.125000000000000, -0.375000000000000,
   0}, {4, 2, 0, 2, -0.500000000000000, -0.375000000000000, 0, 0}, {4,
   2, 2, 0, 5.25000000000000, 2.62500000000000, 0, 0}, {4, 3, 1, 
  0, -2.25000000000000, 1.87500000000000, 0.375000000000000, 0}, {4, 
  3, 1, 2, 2.25000000000000, 0.375000000000000, 0, 0}, {4, 3, 3, 
  0, -7.00000000000000, -0.875000000000000, 0, 0}, {4, 4, 0, 0, 
  0.375000000000000, -0.750000000000000, 0.375000000000000, 0}, {4, 4,
   0, 2, -0.750000000000000, 0.750000000000000, 0, 0}, {4, 4, 0, 4, 
  0.375000000000000, 0, 0, 0}, {4, 4, 2, 0, 
  5.25000000000000, -5.25000000000000, 0, 0}, {4, 4, 2, 
  2, -5.25000000000000, 0, 0, 0}, {4, 4, 4, 0, 7.87500000000000, 0, 0,
   0}, {5, 0, 0, 0, 1.00000000000000, 7.50000000000000, 
  5.62500000000000, 0.312500000000000}, {5, 1, 1, 
  0, -3.50000000000000, -8.75000000000000, -2.18750000000000, 0}, {5, 
  2, 0, 0, 0.500000000000000, 
  0.250000000000000, -0.687500000000000, -0.0625000000000000}, {5, 2, 
  0, 2, -0.500000000000000, -0.750000000000000, -0.0625000000000000, 
  0}, {5, 2, 2, 0, 7.00000000000000, 7.00000000000000, 
  0.437500000000000, 0}, {5, 3, 1, 0, -2.62500000000000, 
  1.31250000000000, 1.31250000000000, 0}, {5, 3, 1, 2, 
  2.62500000000000, 1.31250000000000, 0, 0}, {5, 3, 3, 
  0, -10.5000000000000, -3.93750000000000, 0, 0}, {5, 4, 0, 0, 
  0.375000000000000, -0.687500000000000, 0.250000000000000, 
  0.0625000000000000}, {5, 4, 0, 2, -0.750000000000000, 
  0.625000000000000, 0.125000000000000, 0}, {5, 4, 0, 4, 
  0.375000000000000, 0.0625000000000000, 0, 0}, {5, 4, 2, 0, 
  7.00000000000000, -6.12500000000000, -0.875000000000000, 0}, {5, 4, 
  2, 2, -7.00000000000000, -0.875000000000000, 0, 0}, {5, 4, 4, 0, 
  13.1250000000000, 1.31250000000000, 0, 0}, {5, 5, 1, 
  0, -2.18750000000000, 4.37500000000000, -2.18750000000000, 0}, {5, 
  5, 1, 2, 4.37500000000000, -4.37500000000000, 0, 0}, {5, 5, 1, 
  4, -2.18750000000000, 0, 0, 0}, {5, 5, 3, 0, -13.1250000000000, 
  13.1250000000000, 0, 0}, {5, 5, 3, 2, 13.1250000000000, 0, 0, 
  0}, {5, 5, 5, 0, -14.4375000000000, 0, 0, 0}};

#define TABLEWIDTH_D (int) 8
#define TABLELENGTH_D (int) 15
const double D_TABLE[TABLELENGTH_D][TABLEWIDTH_D] = 
{
    {2, 0, 0, 0, 0, 0, 0, 1}, \
    {2, 0, 2, 0, 0, 0, 0, 2}, \
    {2, 0, 2, 0, 2, 0, 0, 3}, \
    {2, 0, 2, 0, 0, 0, 2, 4}, \
    {2, 2, 0, 0, 0, 0, 0, 5}, \
    {2, 2, 0, 2, 0, 0, 0, 6}, \
    {2, 2, 0, 0, 0, 2, 0, 7}, \
    {3, 1, 0, 1, 0, 0, 0, 8}, \
    {3, 1, 2, 0, 1, 1, 1, 9}, \
    {3, 1, 2, 1, 0, 0, 0, 10}, \
    {3, 1, 2, 1, 0, 0, 2, 11}, \
    {3, 1, 2, 1, 2, 0, 0, 12}, \
    {3, 3, 0, 1, 0, 0, 0, 13}, \
    {3, 3, 0, 1, 0, 2, 0, 14}, \
    {3, 3, 0, 3, 0, 0, 0, 15}
};

/* the numbers in each last entry of the D table refer to the functions defined below */

inline double D_TABLE_FUNC1(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*(sqrt_ef_p2_minus_one + asec_minus_ef);
}
inline double D_TABLE_FUNC1_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC2(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return (1.0 - ep_p2)*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}
inline double D_TABLE_FUNC2_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -2.0*ep*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}

inline double D_TABLE_FUNC3(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) - asec_minus_ef;
}
inline double D_TABLE_FUNC3_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC4(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_2div3*one_div_ef_p2*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC4_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC5(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return ep_p2*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}
inline double D_TABLE_FUNC5_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*ep*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}

inline double D_TABLE_FUNC6(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_2div3*one_div_ef_p2*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC6_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC7(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) - asec_minus_ef;
}
inline double D_TABLE_FUNC7_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC8(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*(c_1div3*one_div_ef_p1*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + ef*asec_minus_ef);
}
inline double D_TABLE_FUNC8_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC9(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div15*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - ef*asec_minus_ef;
}
inline double D_TABLE_FUNC9_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC10(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return (1.0 - ep_p2)*( c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_1div2*ef*asec_minus_ef);
}
inline double D_TABLE_FUNC10_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -2.0*ep*( c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_1div2*ef*asec_minus_ef);
}
inline double D_TABLE_FUNC11(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_4div15*one_div_ef_p3*ef_p2_minus_one*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC11_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC12(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - c_1div2*ef*asec_minus_ef;
}
inline double D_TABLE_FUNC12_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC13(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return ep_p2*( c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_3div2*ef*asec_minus_ef );
}
inline double D_TABLE_FUNC13_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*ep*( c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_3div2*ef*asec_minus_ef );
}
inline double D_TABLE_FUNC14(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - c_3div2*ef*asec_minus_ef;
}
inline double D_TABLE_FUNC14_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC15(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_4div15*one_div_ef_p3*ef_p2_minus_one*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC15_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}



#endif

/*	ODE solver macros	*/
#ifndef __ODE_MACROS
#define __ODE_MACROS
    #define Ith(v,i)    NV_Ith_S(v,i-1)       		/* Ith numbers components 1..NEQ */
    #define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 		/* IJth numbers rows,cols 1..NEQ */
    #ifndef max
        #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
    #endif
    #ifndef min
        #define min(X,Y) ((X) < (Y) ? (X) : (Y))
    #endif
#endif

/* vector operators */
#ifndef __VECTOR_OPERATORS
#define __VECTOR_OPERATORS
inline void cross3(double a[3], double b[3], double result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}
inline double norm3(double v[3])
{
    double result = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return result;
}
inline double norm3_squared(double v[3])
{
    double result = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return result;
}
inline double dot3(double a[3], double b[3])
{
    double result = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    return result;
}
inline void matrix_on_vector( int N_x, int N_y, double M[][3], double v[], double result[] )
{
    
    int i,j;

    for (j=0; j<N_y; j++)
    {
        result[j] = 0.0;
    }

    for (j=0; j<N_y; j++)
    {
        for (i=0; i<N_x; i++)
        {
            result[j] += M[i][j]*v[i];
        }
    }
}

inline void scalar_times_matrix(int N_x, int N_y, double M[][3], double common_factor)
{
    int i,j;
    for (i=0; i<N_x; i++)
    {
        for (j=0; j<N_y; j++)
        {
            M[i][j] *= common_factor;
        }
    }
}

/* KS regularization */
inline double dot4(double a[4], double b[4])
{
    double result = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]);
    return result;
}
inline double norm4(double v[4])
{
    double result = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
    return result;
}
inline double norm4_squared(double v[4])
{
    double result = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
    return result;
}   

inline void transform_r_to_u_d0(double r[3], double u_d0[4])
{
  if (r[0] >= 0.0)
  {
    u_d0[0] = sqrt(0.5*(r[0] + norm3(r)));
    u_d0[1] = 0.5*r[1]/u_d0[0];
    u_d0[2] = 0.5*r[2]/u_d0[0];
    u_d0[3] = 0.0;
  }
  else if (r[0] < 0.0)
  {
    u_d0[1] = sqrt(0.5*(norm3(r) - r[0]));
    u_d0[0] = 0.5*r[1]/u_d0[1];
    u_d0[3] = 0.5*r[2]/u_d0[1];
    u_d0[2] = 0.0;
  }
}

inline void transform_v_to_u_d1(double u_d0[4], double v[3], double u_d1[4])
{
  u_d1[0] = 0.5*(u_d0[0]*v[0] + u_d0[1]*v[1] + u_d0[2]*v[2]);
  u_d1[1] = 0.5*(-u_d0[1]*v[0] + u_d0[0]*v[1] + u_d0[3]*v[2]);
  u_d1[2] = 0.5*(-u_d0[2]*v[0] - u_d0[3]*v[1] + u_d0[0]*v[2]);
  u_d1[3] = 0.5*(u_d0[3]*v[0] - u_d0[2]*v[1] + u_d0[1]*v[2]);
}

inline void transform_u_d0_to_r(double u_d0[4], double r[3])
{
  r[0] = u_d0[0]*u_d0[0] - u_d0[1]*u_d0[1] - u_d0[2]*u_d0[2] + u_d0[3]*u_d0[3];
  r[1] = 2.0*(u_d0[0]*u_d0[1] - u_d0[2]*u_d0[3]);
  r[2] = 2.0*(u_d0[0]*u_d0[2] + u_d0[1]*u_d0[3]);
}

inline void transform_u_d1_to_v(double u_d0[4], double u_d1[4], double v[3])
{
	double r_norm = norm4_squared(u_d0);
  v[0] = 2.0*(u_d0[0]*u_d1[0] - u_d0[1]*u_d1[1] - u_d0[2]*u_d1[2] + u_d0[3]*u_d1[3])/r_norm;
  v[1] = 2.0*(u_d0[1]*u_d1[0] + u_d0[0]*u_d1[1] - u_d0[3]*u_d1[2] - u_d0[2]*u_d1[3])/r_norm;
  v[2] = 2.0*(u_d0[2]*u_d1[0] + u_d0[3]*u_d1[1] + u_d0[0]*u_d1[2] + u_d0[1]*u_d1[3])/r_norm;
}

inline void transform_u_star_to_v(double u[4], double u_star[4], double omega, double v[3])
{
    double r_norm = norm4_squared(u);
    double factor = 4.0*omega/r_norm;
    
    v[0] = factor*(u[0]*u_star[0] - u[1]*u_star[1] - u[2]*u_star[2] + u[3]*u_star[3]);
    v[1] = factor*(u[1]*u_star[0] + u[0]*u_star[1] - u[3]*u_star[2] - u[2]*u_star[3]);
    v[2] = factor*(u[2]*u_star[0] + u[3]*u_star[1] + u[0]*u_star[2] + u[1]*u_star[3]);
}
    
inline void transform_v_to_u_star(double u[4], double v[3], double omega, double u_star[4])
{
    double factor = 0.25/omega;
    u_star[0] = factor*(u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
    u_star[1] = factor*(-u[1]*v[0] + u[0]*v[1] + u[3]*v[2]);
    u_star[2] = factor*(-u[2]*v[0] - u[3]*v[1] + u[0]*v[2]);
    u_star[3] = factor*(u[3]*v[0] - u[2]*v[1] + u[1]*v[2]);
}

inline void LT_u_on_vec3(double u_d0[4], double vec3[3], double result[4])
{
  result[0] = u_d0[0]*vec3[0] + u_d0[1]*vec3[1] + u_d0[2]*vec3[2];
  result[1] = -u_d0[1]*vec3[0] + u_d0[0]*vec3[1] + u_d0[3]*vec3[2];
  result[2] = -u_d0[2]*vec3[0] - u_d0[3]*vec3[1] + u_d0[0]*vec3[2];
  result[3] = u_d0[3]*vec3[0] - u_d0[2]*vec3[1] + u_d0[1]*vec3[2];
}

inline void transform_alpha_beta_to_u_u_star(double alpha[4], double beta[4], double E, double u[4], double u_star[4])
{
    double cos_E_div2 = cos(c_1div2*E);
    double sin_E_div2 = sin(c_1div2*E);
    for (int i=0; i<4; i++)
    {
        u[i] = alpha[i]*cos_E_div2 + beta[i]*sin_E_div2;
        u_star[i] = -c_1div2*alpha[i]*sin_E_div2 + c_1div2*beta[i]*cos_E_div2;
    }
}
#endif

/* classes */
#ifndef __Particle
#define __Particle
class Particle
{
    public:
    /* generic properties */
    int index,child1,child2;
    int parent;
    int sibling;
    std::vector<int> parents;
    std::vector<int> connecting_child_in_parents;
    int level,highest_level;
    bool is_binary;
    double mass,mass_dot,child1_mass,child2_mass,total_system_mass;
    double mu; /* reduced mass */
    int integration_method;

    /*******************
    /* body properties *
     * ****************/
    /* general */
    double radius,radius_dot,radius_ddot;
    double spin_vec_x_dot,spin_vec_y_dot,spin_vec_z_dot;
    int stellar_type;
    
    /* Absolute position/velocities (relative to an arbitrary inertial reference frame) */
    double X,Y,Z;
    double VX,VY,VZ;

    /* used in ODE solver only */
    double spin_vec[3],dspin_vec_dt[3]; 
    double spin_vec_norm;
    double dmass_dt,dradius_dt;    
    
    //void set_ODE_quantities(double delta_time);
    //void reset_ODE_quantities();
    
    /*********************
    /* binary properties *
     * ******************/
    /* general */
   
    /* Relative position/velocities (applies to binaries only) */
    double r,v;
    double r_vec[3];
    double v_vec[3];
    double a_vec[3];
    double r_p2,r_p3;
    double r_pm1,r_pm2,r_pm3;

    /* KS-related */
    double KS_alpha_vec[4],KS_beta_vec[4];
    double KS_u_vec[4],KS_u_star_vec[4];
    double KS_omega,KS_E,KS_a_vec[3];
    double KS_V;
        
    double KS_dalpha_vec_dt[4],KS_dbeta_vec_dt[4];
    double KS_domega_dt,KS_dE_dt;
    
    bool KS_use_perturbing_potential;
    
    /* phases */
    double true_anomaly;
    double initial_mean_anomaly; /* used to track phases of averaged orbits */
    
    /* PN terms */
    bool include_pairwise_1PN_terms,include_pairwise_25PN_terms;
        
    /* tidal friction */
    int include_tidal_friction_terms,tides_method,include_tidal_bulges_precession_terms,include_rotation_precession_terms;
    double minimum_eccentricity_for_tidal_precession;
    double tides_Q_prime; /* depricated */
    double tides_apsidal_motion_constant, tides_time_lag, tides_gyration_radius;
    double tides_viscous_time_scale;
    int tides_viscous_time_scale_prescription;
    double convective_envelope_mass,convective_envelope_radius,luminosity;

    /* root finding */
    bool check_for_secular_breakdown,secular_breakdown_has_occurred;
    bool check_for_dynamical_instability,dynamical_instability_has_occurred;
    int dynamical_instability_criterion;
    int dynamical_instability_central_particle;
    double dynamical_instability_K_parameter;
    bool check_for_physical_collision_or_orbit_crossing,physical_collision_or_orbit_crossing_has_occurred;
    bool check_for_minimum_periapse_distance,minimum_periapse_distance_has_occurred;
    double check_for_minimum_periapse_distance_value;
    bool check_for_RLOF_at_pericentre,check_for_RLOF_at_pericentre_use_sepinsky_fit,RLOF_at_pericentre_has_occurred;
    bool check_for_GW_condition,GW_condition_has_occurred;    

    /* used in ODE solver only */
    double e_vec[3],h_vec[3];
    double e_vec_unit[3],h_vec_unit[3],q_vec_unit[3];
    double de_vec_dt[3],dh_vec_dt[3];
    double child1_mass_plus_child2_mass,child1_mass_minus_child2_mass,child1_mass_times_child2_mass;
    double child1_mass_dot,child2_mass_dot;
    double e,e_p2;
    double j,j_p2,j_p3,j_p4,j_p5; // j=sqrt(1-e^2)
    double h,a;
    
    /* user-specified instantaneous perturbations */
    bool sample_orbital_phase_randomly;
    double instantaneous_perturbation_delta_mass;
    double instantaneous_perturbation_delta_X,instantaneous_perturbation_delta_Y,instantaneous_perturbation_delta_Z;
    double instantaneous_perturbation_delta_VX,instantaneous_perturbation_delta_VY,instantaneous_perturbation_delta_VZ;

    bool is_external;
    double external_t_ref,external_e,external_r_p;

    /* VRR */
    int VRR_model;
    //int VRR_include_GR_precession;
    int VRR_include_mass_precession;
    //int VRR_include_frame_dragging;
    double VRR_mass_precession_rate;
    double VRR_initial_time, VRR_final_time;
    
    /* Dipole model */
    double VRR_Omega_vec_x,VRR_Omega_vec_y,VRR_Omega_vec_z;
    
    /* Distortion model */
    //double VRR_dipole_fraction;
    //double VRR_initial_distortion_axis_vec_x,VRR_initial_distortion_axis_vec_y,VRR_initial_distortion_axis_vec_z;
    //double VRR_final_distortion_axis_vec_x,VRR_final_distortion_axis_vec_y,VRR_final_distortion_axis_vec_z;
    //double VRR_AD_parameter,VRR_AM_parameter,VRR_AJ_parameter;
    
    /* Bar-Or model */
    double VRR_eta_20_init,VRR_eta_a_22_init,VRR_eta_b_22_init,VRR_eta_a_21_init,VRR_eta_b_21_init;
    double VRR_eta_20_final,VRR_eta_a_22_final,VRR_eta_b_22_final,VRR_eta_a_21_final,VRR_eta_b_21_final;
    
    Particle(int index, bool is_binary) : index(index), is_binary(is_binary)
    {
        /* default values */
        integration_method = 0; /* assume full averaging by default */
        KS_omega = KS_E = 0.0;
        KS_use_perturbing_potential = true;
        
        check_for_secular_breakdown = false;
        check_for_dynamical_instability = false;
        check_for_physical_collision_or_orbit_crossing = false;
        check_for_minimum_periapse_distance = false;
        check_for_RLOF_at_pericentre = false;
        check_for_GW_condition = false;
        
        secular_breakdown_has_occurred = false;
        dynamical_instability_has_occurred = false;
        physical_collision_or_orbit_crossing_has_occurred = false;
        minimum_periapse_distance_has_occurred = false;
        RLOF_at_pericentre_has_occurred = false;
        GW_condition_has_occurred = false;
        
        stellar_type = 0;
        
        include_pairwise_1PN_terms = false;
        include_pairwise_25PN_terms = false;
        include_tidal_friction_terms = false;
        include_tidal_bulges_precession_terms = false;
        include_rotation_precession_terms = false;

        radius = 1.0e-10; /* this must be set (to nonzero), otherwise ODE solver will have invalid ewt values */
        tides_method = 1; /* `full' tides equations of motion including spin-orbit terms for all orientations */
        tides_viscous_time_scale_prescription = 0; /* constant, user-specified t_V */
        minimum_eccentricity_for_tidal_precession = 1.0e-3;
        spin_vec_x_dot = spin_vec_y_dot = spin_vec_z_dot = 0.0;
        mass_dot = radius_dot = radius_ddot = 0.0;
        
        sample_orbital_phase_randomly = false; /* false: not sampled randomly */
        instantaneous_perturbation_delta_mass = 0.0;
        instantaneous_perturbation_delta_X = instantaneous_perturbation_delta_Y = instantaneous_perturbation_delta_Z = 0.0;
        instantaneous_perturbation_delta_VX = instantaneous_perturbation_delta_VY = instantaneous_perturbation_delta_VZ = 0.0;
        
        is_external = false;
        external_t_ref = 0.0;
        external_r_p = 1.0;
        external_e = 10.0;
        
        /* VRR */
        VRR_model = 0;
        //VRR_include_GR_precession = 0;
        VRR_include_mass_precession = 0;
        //VRR_include_frame_dragging = 0;
        VRR_mass_precession_rate = 0.0;
        VRR_initial_time = 0.0;
        VRR_final_time = 1.0;
        
        VRR_Omega_vec_x = VRR_Omega_vec_y = VRR_Omega_vec_z = 0.0;
        //VRR_dipole_fraction = 0.5;
        //VRR_initial_distortion_axis_vec_x = 0.0;
        //VRR_initial_distortion_axis_vec_y = 0.0;
        //VRR_initial_distortion_axis_vec_z = 0.0;
        //VRR_final_distortion_axis_vec_x = 0.0;
        //VRR_final_distortion_axis_vec_y = 0.0;
        //VRR_final_distortion_axis_vec_z = 0.0;
        //VRR_AD_parameter = 0.0;
        //VRR_AM_parameter = 0.0;
        //VRR_AJ_parameter = 0.0;
        
        VRR_eta_20_init = 0.0;
        VRR_eta_a_22_init = 0.0;
        VRR_eta_b_22_init = 0.0;
        VRR_eta_a_21_init = 0.0;
        VRR_eta_b_21_init = 0.0;

        VRR_eta_20_final = 0.0;
        VRR_eta_a_22_final = 0.0;
        VRR_eta_b_22_final = 0.0;
        VRR_eta_a_21_final = 0.0;
        VRR_eta_b_21_final = 0.0;

    }
};
#endif

typedef std::map<int, Particle *> ParticlesMap;
typedef std::map<int, Particle *>::iterator ParticlesMapIterator;


extern int highest_particle_index;
extern ParticlesMap particlesMap;


/* CVODE UserData */
#ifndef __UserData
#define __UserData
typedef struct {
	ParticlesMap *particlesMap;
    double hamiltonian;
    int N_root_finding;
    double start_time;
} *UserData;
#endif

}
