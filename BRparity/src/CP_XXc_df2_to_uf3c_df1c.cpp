#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XXc_df2_to_uf3c_df1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_df2_to_uf3c_df1c(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_d = param.m_d;
    const real_t m_s = param.m_s;
    const real_t m_t = param.m_t;
    const real_t GbRt = param.GbRt;
    const real_t GcRt = param.GcRt;
    const real_t GdRt = param.GdRt;
    const real_t GsRt = param.GsRt;
    const real_t GtRt = param.GtRt;
    const real_t GuRt = param.GuRt;
    const real_t gwdR = param.gwdR;
    const real_t gwuR = param.gwuR;
    const real_t m_XX = param.m_XX;
    const real_t s_12 = param.s_12;
    const real_t s_13 = param.s_13;
    const real_t s_14 = param.s_14;
    const real_t s_23 = param.s_23;
    const real_t s_24 = param.s_24;
    const real_t s_34 = param.s_34;
    const real_t m_dtR1 = param.m_dtR1;
    const real_t m_dtR2 = param.m_dtR2;
    const real_t m_dtR3 = param.m_dtR3;
    const real_t m_utR1 = param.m_utR1;
    const real_t m_utR2 = param.m_utR2;
    const real_t m_utR3 = param.m_utR3;
    const real_t reg_prop = param.reg_prop;
    const complex_t lpp_001 = param.lpp_001;
    const complex_t lpp_010 = param.lpp_010;
    const complex_t lpp_101 = param.lpp_101;
    const complex_t lpp_110 = param.lpp_110;
    const complex_t lpp_201 = param.lpp_201;
    const complex_t lpp_202 = param.lpp_202;
    const complex_t lpp_210 = param.lpp_210;
    const complex_t lpp_212 = param.lpp_212;
    const complex_t lpp_220 = param.lpp_220;
    const complex_t lpp_221 = param.lpp_221;
    const complex_t U_dtR_00 = param.U_dtR_00;
    const complex_t U_dtR_01 = param.U_dtR_01;
    const complex_t U_dtR_02 = param.U_dtR_02;
    const complex_t U_dtR_10 = param.U_dtR_10;
    const complex_t U_dtR_11 = param.U_dtR_11;
    const complex_t U_dtR_12 = param.U_dtR_12;
    const complex_t U_dtR_20 = param.U_dtR_20;
    const complex_t U_dtR_21 = param.U_dtR_21;
    const complex_t U_dtR_22 = param.U_dtR_22;
    const complex_t U_utR_00 = param.U_utR_00;
    const complex_t U_utR_01 = param.U_utR_01;
    const complex_t U_utR_02 = param.U_utR_02;
    const complex_t U_utR_10 = param.U_utR_10;
    const complex_t U_utR_11 = param.U_utR_11;
    const complex_t U_utR_12 = param.U_utR_12;
    const complex_t U_utR_20 = param.U_utR_20;
    const complex_t U_utR_21 = param.U_utR_21;
    const complex_t U_utR_22 = param.U_utR_22;
    const complex_t IT_0000 = std::pow(m_t, 2);
    const complex_t IT_0001 = std::pow(m_XX, 2);
    const complex_t IT_0002 = std::pow(m_utR2, 2);
    const complex_t IT_0003 = std::pow((-2)*s_13 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0000 + IT_0001 + -IT_0002 + reg_prop, -1);
    const complex_t IT_0004 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_21);
    const complex_t IT_0005 = (complex_t{0, 1})*IT_0003*IT_0004;
    const complex_t IT_0006 = lpp_010*U_utR_01;
    const complex_t IT_0007 = lpp_110*U_utR_11;
    const complex_t IT_0008 = lpp_210*U_utR_21;
    const complex_t IT_0009 = (complex_t{0, 1})*(IT_0006 + IT_0007 + IT_0008);
    const complex_t IT_0010 = IT_0005*IT_0009;
    const complex_t IT_0011 = lpp_001*U_utR_01;
    const complex_t IT_0012 = lpp_101*U_utR_11;
    const complex_t IT_0013 = lpp_201*U_utR_21;
    const complex_t IT_0014 = (complex_t{0, 1})*(IT_0011 + IT_0012 + IT_0013);
    const complex_t IT_0015 = IT_0005*IT_0014;
    const complex_t IT_0016 = std::pow(m_utR3, 2);
    const complex_t IT_0017 = std::pow((-2)*s_13 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0000 + IT_0001 + -IT_0016 + reg_prop, -1);
    const complex_t IT_0018 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_22);
    const complex_t IT_0019 = (complex_t{0, 1})*IT_0017*IT_0018;
    const complex_t IT_0020 = lpp_010*U_utR_02;
    const complex_t IT_0021 = lpp_110*U_utR_12;
    const complex_t IT_0022 = lpp_210*U_utR_22;
    const complex_t IT_0023 = (complex_t{0, 1})*(IT_0020 + IT_0021 + IT_0022);
    const complex_t IT_0024 = IT_0019*IT_0023;
    const complex_t IT_0025 = lpp_001*U_utR_02;
    const complex_t IT_0026 = lpp_101*U_utR_12;
    const complex_t IT_0027 = lpp_201*U_utR_22;
    const complex_t IT_0028 = (complex_t{0, 1})*(IT_0025 + IT_0026 + IT_0027);
    const complex_t IT_0029 = IT_0019*IT_0028;
    const complex_t IT_0030 = std::pow(m_utR1, 2);
    const complex_t IT_0031 = std::pow((-2)*s_13 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0000 + IT_0001 + -IT_0030 + reg_prop, -1);
    const complex_t IT_0032 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_20);
    const complex_t IT_0033 = (complex_t{0, 1})*IT_0031*IT_0032;
    const complex_t IT_0034 = lpp_010*U_utR_00;
    const complex_t IT_0035 = lpp_110*U_utR_10;
    const complex_t IT_0036 = lpp_210*U_utR_20;
    const complex_t IT_0037 = (complex_t{0, 1})*(IT_0034 + IT_0035 + IT_0036);
    const complex_t IT_0038 = lpp_001*U_utR_00;
    const complex_t IT_0039 = lpp_101*U_utR_10;
    const complex_t IT_0040 = lpp_201*U_utR_20;
    const complex_t IT_0041 = (complex_t{0, 1})*(IT_0038 + IT_0039 + IT_0040);
    const complex_t IT_0042 = IT_0037 + -IT_0041;
    const complex_t IT_0043 = IT_0033*IT_0042;
    const complex_t IT_0044 = (-0.25)*IT_0010 + 0.25*IT_0015 + (-0.25)*IT_0024
       + 0.25*IT_0029 + (-0.25)*IT_0043;
    const complex_t IT_0045 = 0.25*IT_0010 + (-0.25)*IT_0015 + 0.25*IT_0024 + 
      (-0.25)*IT_0029 + 0.25*IT_0043;
    const complex_t IT_0046 = std::conj(IT_0044) + std::conj(IT_0045);
    const complex_t IT_0047 = lpp_201*U_dtR_00;
    const complex_t IT_0048 = lpp_210*U_dtR_00;
    const complex_t IT_0049 = lpp_212*U_dtR_20;
    const complex_t IT_0050 = lpp_221*U_dtR_20;
    const complex_t IT_0051 = (complex_t{0, 1})*(IT_0047 + -IT_0048 + -IT_0049
       + IT_0050);
    const complex_t IT_0052 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_00);
    const complex_t IT_0053 = IT_0051*IT_0052;
    const complex_t IT_0054 = std::pow(m_s, 2);
    const complex_t IT_0055 = std::pow(m_dtR1, 2);
    const complex_t IT_0056 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0000 + IT_0054 + -IT_0055 + reg_prop, -1);
    const complex_t IT_0057 = (complex_t{0, 1})*IT_0053*IT_0056;
    const complex_t IT_0058 = lpp_201*U_dtR_01;
    const complex_t IT_0059 = lpp_210*U_dtR_01;
    const complex_t IT_0060 = lpp_212*U_dtR_21;
    const complex_t IT_0061 = lpp_221*U_dtR_21;
    const complex_t IT_0062 = (complex_t{0, 1})*(IT_0058 + -IT_0059 + -IT_0060
       + IT_0061);
    const complex_t IT_0063 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_01);
    const complex_t IT_0064 = IT_0062*IT_0063;
    const complex_t IT_0065 = std::pow(m_dtR2, 2);
    const complex_t IT_0066 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0000 + IT_0054 + -IT_0065 + reg_prop, -1);
    const complex_t IT_0067 = (complex_t{0, 1})*IT_0064*IT_0066;
    const complex_t IT_0068 = lpp_201*U_dtR_02;
    const complex_t IT_0069 = lpp_210*U_dtR_02;
    const complex_t IT_0070 = lpp_212*U_dtR_22;
    const complex_t IT_0071 = lpp_221*U_dtR_22;
    const complex_t IT_0072 = (complex_t{0, 1})*(IT_0068 + -IT_0069 + -IT_0070
       + IT_0071);
    const complex_t IT_0073 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_02);
    const complex_t IT_0074 = IT_0072*IT_0073;
    const complex_t IT_0075 = std::pow(m_dtR3, 2);
    const complex_t IT_0076 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0000 + IT_0054 + -IT_0075 + reg_prop, -1);
    const complex_t IT_0077 = (complex_t{0, 1})*IT_0074*IT_0076;
    const complex_t IT_0078 = (-0.5)*IT_0057 + (-0.5)*IT_0067 + (-0.5)*IT_0077;
    const complex_t IT_0079 = s_14*s_23;
    const complex_t IT_0080 = s_13*s_24;
    const complex_t IT_0081 = -IT_0080;
    const complex_t IT_0082 = IT_0079 + IT_0081;
    const complex_t IT_0083 = 8*IT_0082;
    const complex_t IT_0084 = IT_0078*IT_0083;
    const complex_t IT_0085 = (-0.25)*IT_0057 + (-0.25)*IT_0067 + (-0.25)
      *IT_0077;
    const complex_t IT_0086 = m_d*m_t*s_12;
    const complex_t IT_0087 = m_s*m_XX*s_34;
    const complex_t IT_0088 = 96*IT_0087;
    const complex_t IT_0089 = m_d*m_s*m_t*m_XX;
    const complex_t IT_0090 = 192*IT_0089;
    const complex_t IT_0091 = (-64)*IT_0080;
    const complex_t IT_0092 = s_12*s_34;
    const complex_t IT_0093 = 64*IT_0092;
    const complex_t IT_0094 = 96*IT_0086 + IT_0088 + IT_0090 + IT_0091 +
       IT_0093;
    const complex_t IT_0095 = 128*IT_0080;
    const complex_t IT_0096 = 128*IT_0079;
    const complex_t IT_0097 = (-64)*IT_0092;
    const complex_t IT_0098 = IT_0095 + IT_0096 + IT_0097;
    const complex_t IT_0099 = (-192)*IT_0089;
    const complex_t IT_0100 = IT_0098 + IT_0099;
    const complex_t IT_0101 = std::conj(IT_0078)*IT_0083;
    const complex_t IT_0102 = lpp_201*U_dtR_10;
    const complex_t IT_0103 = lpp_202*U_dtR_20;
    const complex_t IT_0104 = lpp_210*U_dtR_10;
    const complex_t IT_0105 = lpp_220*U_dtR_20;
    const complex_t IT_0106 = (complex_t{0, -1})*(IT_0102 + IT_0103 + -IT_0104
       + -IT_0105);
    const complex_t IT_0107 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_10);
    const complex_t IT_0108 = IT_0106*IT_0107;
    const complex_t IT_0109 = std::pow(2*s_12 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0001 + IT_0054 + -IT_0055 + reg_prop, -1);
    const complex_t IT_0110 = (complex_t{0, 1})*IT_0108*IT_0109;
    const complex_t IT_0111 = lpp_201*U_dtR_11;
    const complex_t IT_0112 = lpp_202*U_dtR_21;
    const complex_t IT_0113 = lpp_210*U_dtR_11;
    const complex_t IT_0114 = lpp_220*U_dtR_21;
    const complex_t IT_0115 = (complex_t{0, -1})*(IT_0111 + IT_0112 + -IT_0113
       + -IT_0114);
    const complex_t IT_0116 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_11);
    const complex_t IT_0117 = IT_0115*IT_0116;
    const complex_t IT_0118 = std::pow(2*s_12 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0001 + IT_0054 + -IT_0065 + reg_prop, -1);
    const complex_t IT_0119 = (complex_t{0, 1})*IT_0117*IT_0118;
    const complex_t IT_0120 = lpp_201*U_dtR_12;
    const complex_t IT_0121 = lpp_202*U_dtR_22;
    const complex_t IT_0122 = lpp_210*U_dtR_12;
    const complex_t IT_0123 = lpp_220*U_dtR_22;
    const complex_t IT_0124 = (complex_t{0, -1})*(IT_0120 + IT_0121 + -IT_0122
       + -IT_0123);
    const complex_t IT_0125 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_12);
    const complex_t IT_0126 = IT_0124*IT_0125;
    const complex_t IT_0127 = std::pow(2*s_12 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0001 + IT_0054 + -IT_0075 + reg_prop, -1);
    const complex_t IT_0128 = (complex_t{0, 1})*IT_0126*IT_0127;
    const complex_t IT_0129 = IT_0033*IT_0037;
    const complex_t IT_0130 = IT_0033*IT_0041;
    const complex_t IT_0131 = (-0.5)*IT_0010 + 0.5*IT_0015 + (-0.5)*IT_0024 +
       0.5*IT_0029 + IT_0110 + IT_0119 + IT_0128 + (-0.5)*IT_0129 + 0.5*IT_0130;
    const complex_t IT_0132 = (-8)*IT_0082;
    const complex_t IT_0133 = std::conj(IT_0131)*IT_0132;
    const complex_t IT_0134 = IT_0090 + IT_0098;
    const complex_t IT_0135 = 0.25*IT_0057 + 0.25*IT_0067 + 0.25*IT_0077;
    const complex_t IT_0136 = (-128)*IT_0080;
    const complex_t IT_0137 = IT_0083*std::conj(IT_0131);
    const complex_t IT_0138 = std::conj(IT_0078)*IT_0132;
    const complex_t IT_0139 = (-2)*IT_0092;
    const complex_t IT_0140 = std::conj(IT_0085) + std::conj(IT_0135);
    const complex_t IT_0141 = 0.5*IT_0140;
    const complex_t IT_0142 = IT_0046*IT_0084 + IT_0045*(std::conj(IT_0085)
      *IT_0094 + std::conj(IT_0045)*IT_0100 + IT_0101 + IT_0133) + IT_0044*
      (IT_0101 + IT_0133 + std::conj(IT_0044)*IT_0134 + std::conj(IT_0135)*((
      -128)*IT_0079 + IT_0093 + IT_0099 + IT_0136)) + IT_0085*(std::conj(IT_0045
      )*IT_0094 + std::conj(IT_0085)*IT_0100 + IT_0137 + IT_0138) + IT_0135*
      (IT_0134*std::conj(IT_0135) + std::conj(IT_0044)*((-128)*IT_0079 + IT_0093
       + IT_0099 + IT_0136) + IT_0137 + IT_0138) + 2*IT_0131*(IT_0092*std::conj
      (IT_0131) + 0.5*IT_0046*IT_0132 + 0.5*std::conj(IT_0078)*IT_0139 + IT_0083
      *IT_0141) + 2*IT_0078*(std::conj(IT_0078)*IT_0092 + 0.5*std::conj(IT_0131)
      *IT_0139 + IT_0132*IT_0141);
    return IT_0142;
}
} // End of namespace brparity
