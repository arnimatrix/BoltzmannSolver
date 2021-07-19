#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XX_uf1_to_df2c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_uf1_to_df2c_df3c(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_b = param.m_b;
    const real_t m_s = param.m_s;
    const real_t m_u = param.m_u;
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
    const complex_t lpp_002 = param.lpp_002;
    const complex_t lpp_010 = param.lpp_010;
    const complex_t lpp_012 = param.lpp_012;
    const complex_t lpp_020 = param.lpp_020;
    const complex_t lpp_021 = param.lpp_021;
    const complex_t lpp_112 = param.lpp_112;
    const complex_t lpp_121 = param.lpp_121;
    const complex_t lpp_212 = param.lpp_212;
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
    const complex_t IT_0000 = lpp_002*U_dtR_00;
    const complex_t IT_0001 = lpp_012*U_dtR_10;
    const complex_t IT_0002 = lpp_020*U_dtR_00;
    const complex_t IT_0003 = lpp_021*U_dtR_10;
    const complex_t IT_0004 = (complex_t{0, 1})*(IT_0000 + IT_0001 + -IT_0002 
      + -IT_0003);
    const complex_t IT_0005 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_10);
    const complex_t IT_0006 = IT_0004*IT_0005;
    const complex_t IT_0007 = std::pow(m_s, 2);
    const complex_t IT_0008 = std::pow(m_XX, 2);
    const complex_t IT_0009 = std::pow(m_dtR1, 2);
    const complex_t IT_0010 = std::pow((-2)*s_13 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0007 + IT_0008 + -IT_0009 + reg_prop, -1);
    const complex_t IT_0011 = (complex_t{0, 1})*IT_0006*IT_0010;
    const complex_t IT_0012 = lpp_002*U_dtR_01;
    const complex_t IT_0013 = lpp_012*U_dtR_11;
    const complex_t IT_0014 = lpp_020*U_dtR_01;
    const complex_t IT_0015 = lpp_021*U_dtR_11;
    const complex_t IT_0016 = (complex_t{0, 1})*(IT_0012 + IT_0013 + -IT_0014 
      + -IT_0015);
    const complex_t IT_0017 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_11);
    const complex_t IT_0018 = IT_0016*IT_0017;
    const complex_t IT_0019 = std::pow(m_dtR2, 2);
    const complex_t IT_0020 = std::pow((-2)*s_13 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0007 + IT_0008 + -IT_0019 + reg_prop, -1);
    const complex_t IT_0021 = (complex_t{0, 1})*IT_0018*IT_0020;
    const complex_t IT_0022 = lpp_002*U_dtR_02;
    const complex_t IT_0023 = lpp_012*U_dtR_12;
    const complex_t IT_0024 = lpp_020*U_dtR_02;
    const complex_t IT_0025 = lpp_021*U_dtR_12;
    const complex_t IT_0026 = (complex_t{0, 1})*(IT_0022 + IT_0023 + -IT_0024 
      + -IT_0025);
    const complex_t IT_0027 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_12);
    const complex_t IT_0028 = IT_0026*IT_0027;
    const complex_t IT_0029 = std::pow(m_dtR3, 2);
    const complex_t IT_0030 = std::pow((-2)*s_13 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + IT_0008 + -IT_0029 + reg_prop, -1);
    const complex_t IT_0031 = (complex_t{0, 1})*IT_0028*IT_0030;
    const complex_t IT_0032 = (-0.25)*IT_0011 + (-0.25)*IT_0021 + (-0.25)
      *IT_0031;
    const complex_t IT_0033 = 0.25*IT_0011 + 0.25*IT_0021 + 0.25*IT_0031;
    const complex_t IT_0034 = std::conj(IT_0032) + std::conj(IT_0033);
    const complex_t IT_0035 = 0.5*IT_0011 + 0.5*IT_0021 + 0.5*IT_0031;
    const complex_t IT_0036 = s_14*s_23;
    const complex_t IT_0037 = s_13*s_24;
    const complex_t IT_0038 = -IT_0037;
    const complex_t IT_0039 = IT_0036 + IT_0038;
    const complex_t IT_0040 = 8*IT_0039;
    const complex_t IT_0041 = IT_0035*IT_0040;
    const complex_t IT_0042 = 128*IT_0037;
    const complex_t IT_0043 = 128*IT_0036;
    const complex_t IT_0044 = s_12*s_34;
    const complex_t IT_0045 = (-64)*IT_0044;
    const complex_t IT_0046 = IT_0042 + IT_0043 + IT_0045;
    const complex_t IT_0047 = m_b*m_s*m_u*m_XX;
    const complex_t IT_0048 = 192*IT_0047;
    const complex_t IT_0049 = IT_0046 + IT_0048;
    const complex_t IT_0050 = lpp_001*U_dtR_00;
    const complex_t IT_0051 = lpp_010*U_dtR_00;
    const complex_t IT_0052 = lpp_012*U_dtR_20;
    const complex_t IT_0053 = lpp_021*U_dtR_20;
    const complex_t IT_0054 = (complex_t{0, 1})*(IT_0050 + -IT_0051 + -IT_0052
       + IT_0053);
    const complex_t IT_0055 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_20);
    const complex_t IT_0056 = IT_0054*IT_0055;
    const complex_t IT_0057 = std::pow(m_u, 2);
    const complex_t IT_0058 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0007 + -IT_0009 + IT_0057 + reg_prop, -1);
    const complex_t IT_0059 = (complex_t{0, 1})*IT_0056*IT_0058;
    const complex_t IT_0060 = lpp_001*U_dtR_01;
    const complex_t IT_0061 = lpp_010*U_dtR_01;
    const complex_t IT_0062 = lpp_012*U_dtR_21;
    const complex_t IT_0063 = lpp_021*U_dtR_21;
    const complex_t IT_0064 = (complex_t{0, 1})*(IT_0060 + -IT_0061 + -IT_0062
       + IT_0063);
    const complex_t IT_0065 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_21);
    const complex_t IT_0066 = IT_0064*IT_0065;
    const complex_t IT_0067 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0007 + -IT_0019 + IT_0057 + reg_prop, -1);
    const complex_t IT_0068 = (complex_t{0, 1})*IT_0066*IT_0067;
    const complex_t IT_0069 = lpp_001*U_dtR_02;
    const complex_t IT_0070 = lpp_010*U_dtR_02;
    const complex_t IT_0071 = lpp_012*U_dtR_22;
    const complex_t IT_0072 = lpp_021*U_dtR_22;
    const complex_t IT_0073 = (complex_t{0, 1})*(IT_0069 + -IT_0070 + -IT_0071
       + IT_0072);
    const complex_t IT_0074 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_22);
    const complex_t IT_0075 = IT_0073*IT_0074;
    const complex_t IT_0076 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + -IT_0029 + IT_0057 + reg_prop, -1);
    const complex_t IT_0077 = (complex_t{0, 1})*IT_0075*IT_0076;
    const complex_t IT_0078 = 0.25*IT_0059 + 0.25*IT_0068 + 0.25*IT_0077;
    const complex_t IT_0079 = (-128)*IT_0037;
    const complex_t IT_0080 = (-128)*IT_0036;
    const complex_t IT_0081 = 64*IT_0044;
    const complex_t IT_0082 = IT_0079 + IT_0080 + IT_0081;
    const complex_t IT_0083 = (-192)*IT_0047;
    const complex_t IT_0084 = IT_0082 + IT_0083;
    const complex_t IT_0085 = std::conj(IT_0035)*IT_0040;
    const complex_t IT_0086 = lpp_012*U_utR_00;
    const complex_t IT_0087 = lpp_112*U_utR_10;
    const complex_t IT_0088 = lpp_212*U_utR_20;
    const complex_t IT_0089 = (complex_t{0, 1})*(IT_0086 + IT_0087 + IT_0088);
    const complex_t IT_0090 = -IT_0089;
    const complex_t IT_0091 = lpp_021*U_utR_00;
    const complex_t IT_0092 = lpp_121*U_utR_10;
    const complex_t IT_0093 = lpp_221*U_utR_20;
    const complex_t IT_0094 = (complex_t{0, 1})*(IT_0091 + IT_0092 + IT_0093);
    const complex_t IT_0095 = IT_0090 + IT_0094;
    const complex_t IT_0096 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_00);
    const complex_t IT_0097 = IT_0095*IT_0096;
    const complex_t IT_0098 = std::pow(m_utR1, 2);
    const complex_t IT_0099 = std::pow(2*s_12 + (complex_t{0, 1})*GuRt*m_utR1 
      + IT_0008 + IT_0057 + -IT_0098 + reg_prop, -1);
    const complex_t IT_0100 = (complex_t{0, 1})*IT_0097*IT_0099;
    const complex_t IT_0101 = lpp_012*U_utR_01;
    const complex_t IT_0102 = lpp_112*U_utR_11;
    const complex_t IT_0103 = lpp_212*U_utR_21;
    const complex_t IT_0104 = (complex_t{0, 1})*(IT_0101 + IT_0102 + IT_0103);
    const complex_t IT_0105 = -IT_0104;
    const complex_t IT_0106 = lpp_021*U_utR_01;
    const complex_t IT_0107 = lpp_121*U_utR_11;
    const complex_t IT_0108 = lpp_221*U_utR_21;
    const complex_t IT_0109 = (complex_t{0, 1})*(IT_0106 + IT_0107 + IT_0108);
    const complex_t IT_0110 = IT_0105 + IT_0109;
    const complex_t IT_0111 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_01);
    const complex_t IT_0112 = IT_0110*IT_0111;
    const complex_t IT_0113 = std::pow(m_utR2, 2);
    const complex_t IT_0114 = std::pow(2*s_12 + (complex_t{0, 1})*GcRt*m_utR2 
      + IT_0008 + IT_0057 + -IT_0113 + reg_prop, -1);
    const complex_t IT_0115 = (complex_t{0, 1})*IT_0112*IT_0114;
    const complex_t IT_0116 = lpp_012*U_utR_02;
    const complex_t IT_0117 = lpp_112*U_utR_12;
    const complex_t IT_0118 = lpp_212*U_utR_22;
    const complex_t IT_0119 = (complex_t{0, 1})*(IT_0116 + IT_0117 + IT_0118);
    const complex_t IT_0120 = -IT_0119;
    const complex_t IT_0121 = lpp_021*U_utR_02;
    const complex_t IT_0122 = lpp_121*U_utR_12;
    const complex_t IT_0123 = lpp_221*U_utR_22;
    const complex_t IT_0124 = (complex_t{0, 1})*(IT_0121 + IT_0122 + IT_0123);
    const complex_t IT_0125 = IT_0120 + IT_0124;
    const complex_t IT_0126 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_02);
    const complex_t IT_0127 = IT_0125*IT_0126;
    const complex_t IT_0128 = std::pow(m_utR3, 2);
    const complex_t IT_0129 = std::pow(2*s_12 + (complex_t{0, 1})*GtRt*m_utR3 
      + IT_0008 + IT_0057 + -IT_0128 + reg_prop, -1);
    const complex_t IT_0130 = (complex_t{0, 1})*IT_0127*IT_0129;
    const complex_t IT_0131 = 0.5*IT_0059 + 0.5*IT_0068 + 0.5*IT_0077 + 
      -IT_0100 + -IT_0115 + -IT_0130;
    const complex_t IT_0132 = (-8)*IT_0039;
    const complex_t IT_0133 = std::conj(IT_0131)*IT_0132;
    const complex_t IT_0134 = IT_0046 + IT_0083;
    const complex_t IT_0135 = (-0.25)*IT_0059 + (-0.25)*IT_0068 + (-0.25)
      *IT_0077;
    const complex_t IT_0136 = (-2)*IT_0044;
    const complex_t IT_0137 = IT_0040*std::conj(IT_0131);
    const complex_t IT_0138 = std::conj(IT_0035)*IT_0132;
    const complex_t IT_0139 = IT_0034*IT_0041 + IT_0032*(std::conj(IT_0032)
      *IT_0049 + std::conj(IT_0078)*IT_0084 + IT_0085 + IT_0133) + IT_0033*
      (IT_0085 + IT_0133 + std::conj(IT_0033)*IT_0134 + (IT_0048 + IT_0082)
      *std::conj(IT_0135)) + 2*IT_0131*(IT_0044*std::conj(IT_0131) + 0.5*IT_0034
      *IT_0132 + 0.5*IT_0040*(std::conj(IT_0078) + std::conj(IT_0135)) + 0.5
      *std::conj(IT_0035)*IT_0136) + 2*IT_0035*(std::conj(IT_0035)*IT_0044 + 0.5
      *IT_0132*(std::conj(IT_0078) + std::conj(IT_0135)) + 0.5*std::conj(IT_0131
      )*IT_0136) + IT_0078*(IT_0049*std::conj(IT_0078) + std::conj(IT_0032)
      *IT_0084 + IT_0137 + IT_0138) + IT_0135*(std::conj(IT_0033)*(IT_0048 +
       IT_0082) + IT_0134*std::conj(IT_0135) + IT_0137 + IT_0138);
    return IT_0139;
}
} // End of namespace brparity
