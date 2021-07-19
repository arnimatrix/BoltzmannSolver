#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XXc_uf1c_to_df2_df3.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_uf1c_to_df2_df3(
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
    const complex_t IT_0000 = std::conj(lpp_001)*std::conj(U_dtR_00);
    const complex_t IT_0001 = std::conj(lpp_010)*std::conj(U_dtR_00);
    const complex_t IT_0002 = std::conj(lpp_012)*std::conj(U_dtR_20);
    const complex_t IT_0003 = std::conj(lpp_021)*std::conj(U_dtR_20);
    const complex_t IT_0004 = (complex_t{0, 1})*(IT_0000 + -IT_0001 + -IT_0002
       + IT_0003);
    const complex_t IT_0005 = (complex_t{0, -1})*gw*gwdR*U_dtR_20;
    const complex_t IT_0006 = IT_0004*IT_0005;
    const complex_t IT_0007 = std::pow(m_s, 2);
    const complex_t IT_0008 = std::pow(m_u, 2);
    const complex_t IT_0009 = std::pow(m_dtR1, 2);
    const complex_t IT_0010 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0007 + IT_0008 + -IT_0009 + reg_prop, -1);
    const complex_t IT_0011 = (complex_t{0, 1})*IT_0006*IT_0010;
    const complex_t IT_0012 = std::conj(lpp_001)*std::conj(U_dtR_01);
    const complex_t IT_0013 = std::conj(lpp_010)*std::conj(U_dtR_01);
    const complex_t IT_0014 = std::conj(lpp_012)*std::conj(U_dtR_21);
    const complex_t IT_0015 = std::conj(lpp_021)*std::conj(U_dtR_21);
    const complex_t IT_0016 = (complex_t{0, 1})*(IT_0012 + -IT_0013 + -IT_0014
       + IT_0015);
    const complex_t IT_0017 = (complex_t{0, -1})*gw*gwdR*U_dtR_21;
    const complex_t IT_0018 = IT_0016*IT_0017;
    const complex_t IT_0019 = std::pow(m_dtR2, 2);
    const complex_t IT_0020 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0007 + IT_0008 + -IT_0019 + reg_prop, -1);
    const complex_t IT_0021 = (complex_t{0, 1})*IT_0018*IT_0020;
    const complex_t IT_0022 = std::conj(lpp_001)*std::conj(U_dtR_02);
    const complex_t IT_0023 = std::conj(lpp_010)*std::conj(U_dtR_02);
    const complex_t IT_0024 = std::conj(lpp_012)*std::conj(U_dtR_22);
    const complex_t IT_0025 = std::conj(lpp_021)*std::conj(U_dtR_22);
    const complex_t IT_0026 = (complex_t{0, 1})*(IT_0022 + -IT_0023 + -IT_0024
       + IT_0025);
    const complex_t IT_0027 = (complex_t{0, -1})*gw*gwdR*U_dtR_22;
    const complex_t IT_0028 = IT_0026*IT_0027;
    const complex_t IT_0029 = std::pow(m_dtR3, 2);
    const complex_t IT_0030 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + IT_0008 + -IT_0029 + reg_prop, -1);
    const complex_t IT_0031 = (complex_t{0, 1})*IT_0028*IT_0030;
    const complex_t IT_0032 = (-0.25)*IT_0011 + (-0.25)*IT_0021 + (-0.25)
      *IT_0031;
    const complex_t IT_0033 = std::conj(lpp_002)*std::conj(U_dtR_01);
    const complex_t IT_0034 = std::conj(lpp_012)*std::conj(U_dtR_11);
    const complex_t IT_0035 = std::conj(lpp_020)*std::conj(U_dtR_01);
    const complex_t IT_0036 = std::conj(lpp_021)*std::conj(U_dtR_11);
    const complex_t IT_0037 = (complex_t{0, 1})*(IT_0033 + IT_0034 + -IT_0035 
      + -IT_0036);
    const complex_t IT_0038 = (complex_t{0, -1})*gw*gwdR*U_dtR_11;
    const complex_t IT_0039 = IT_0037*IT_0038;
    const complex_t IT_0040 = std::pow(m_XX, 2);
    const complex_t IT_0041 = std::pow((-2)*s_13 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0007 + -IT_0019 + IT_0040 + reg_prop, -1);
    const complex_t IT_0042 = (complex_t{0, 1})*IT_0039*IT_0041;
    const complex_t IT_0043 = std::conj(lpp_002)*std::conj(U_dtR_02);
    const complex_t IT_0044 = std::conj(lpp_012)*std::conj(U_dtR_12);
    const complex_t IT_0045 = std::conj(lpp_020)*std::conj(U_dtR_02);
    const complex_t IT_0046 = std::conj(lpp_021)*std::conj(U_dtR_12);
    const complex_t IT_0047 = (complex_t{0, 1})*(IT_0043 + IT_0044 + -IT_0045 
      + -IT_0046);
    const complex_t IT_0048 = (complex_t{0, -1})*gw*gwdR*U_dtR_12;
    const complex_t IT_0049 = IT_0047*IT_0048;
    const complex_t IT_0050 = std::pow((-2)*s_13 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + -IT_0029 + IT_0040 + reg_prop, -1);
    const complex_t IT_0051 = (complex_t{0, 1})*IT_0049*IT_0050;
    const complex_t IT_0052 = 0.25*IT_0042 + 0.25*IT_0051;
    const complex_t IT_0053 = m_b*m_s*m_u*m_XX;
    const complex_t IT_0054 = 192*IT_0053;
    const complex_t IT_0055 = s_13*s_24;
    const complex_t IT_0056 = (-128)*IT_0055;
    const complex_t IT_0057 = s_14*s_23;
    const complex_t IT_0058 = (-128)*IT_0057;
    const complex_t IT_0059 = s_12*s_34;
    const complex_t IT_0060 = 64*IT_0059;
    const complex_t IT_0061 = IT_0054 + IT_0056 + IT_0058 + IT_0060;
    const complex_t IT_0062 = std::conj(lpp_002)*std::conj(U_dtR_00);
    const complex_t IT_0063 = std::conj(lpp_012)*std::conj(U_dtR_10);
    const complex_t IT_0064 = std::conj(lpp_020)*std::conj(U_dtR_00);
    const complex_t IT_0065 = std::conj(lpp_021)*std::conj(U_dtR_10);
    const complex_t IT_0066 = (complex_t{0, 1})*(IT_0062 + IT_0063 + -IT_0064 
      + -IT_0065);
    const complex_t IT_0067 = (complex_t{0, -1})*gw*gwdR*U_dtR_10;
    const complex_t IT_0068 = IT_0066*IT_0067;
    const complex_t IT_0069 = std::pow((-2)*s_13 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0007 + -IT_0009 + IT_0040 + reg_prop, -1);
    const complex_t IT_0070 = (complex_t{0, 1})*IT_0068*IT_0069;
    const complex_t IT_0071 = 0.25*IT_0070;
    const complex_t IT_0072 = m_b*m_s*s_12;
    const complex_t IT_0073 = m_u*m_XX*s_34;
    const complex_t IT_0074 = 96*IT_0073;
    const complex_t IT_0075 = (-64)*IT_0057;
    const complex_t IT_0076 = IT_0054 + IT_0060 + 96*IT_0072 + IT_0074 +
       IT_0075;
    const complex_t IT_0077 = 128*IT_0055;
    const complex_t IT_0078 = 128*IT_0057;
    const complex_t IT_0079 = (-64)*IT_0059;
    const complex_t IT_0080 = IT_0077 + IT_0078 + IT_0079;
    const complex_t IT_0081 = IT_0054 + IT_0080;
    const complex_t IT_0082 = (-192)*IT_0053;
    const complex_t IT_0083 = IT_0080 + IT_0082;
    const complex_t IT_0084 = 0.25*IT_0042 + 0.25*IT_0051 + 0.25*IT_0070;
    const complex_t IT_0085 = IT_0056 + IT_0058 + IT_0060 + IT_0082;
    const complex_t IT_0086 = std::conj(lpp_012)*std::conj(U_utR_00);
    const complex_t IT_0087 = std::conj(lpp_112)*std::conj(U_utR_10);
    const complex_t IT_0088 = std::conj(lpp_212)*std::conj(U_utR_20);
    const complex_t IT_0089 = std::conj(lpp_021)*std::conj(U_utR_00);
    const complex_t IT_0090 = std::conj(lpp_121)*std::conj(U_utR_10);
    const complex_t IT_0091 = std::conj(lpp_221)*std::conj(U_utR_20);
    const complex_t IT_0092 = (complex_t{0, 1})*(IT_0086 + IT_0087 + IT_0088 +
       -IT_0089 + -IT_0090 + -IT_0091);
    const complex_t IT_0093 = (complex_t{0, -1})*gw*gwuR*U_utR_00;
    const complex_t IT_0094 = IT_0092*IT_0093;
    const complex_t IT_0095 = std::pow(m_utR1, 2);
    const complex_t IT_0096 = std::pow(2*s_12 + (complex_t{0, 1})*GuRt*m_utR1 
      + IT_0008 + IT_0040 + -IT_0095 + reg_prop, -1);
    const complex_t IT_0097 = (complex_t{0, 1})*IT_0094*IT_0096;
    const complex_t IT_0098 = std::conj(lpp_012)*std::conj(U_utR_01);
    const complex_t IT_0099 = std::conj(lpp_112)*std::conj(U_utR_11);
    const complex_t IT_0100 = std::conj(lpp_212)*std::conj(U_utR_21);
    const complex_t IT_0101 = std::conj(lpp_021)*std::conj(U_utR_01);
    const complex_t IT_0102 = std::conj(lpp_121)*std::conj(U_utR_11);
    const complex_t IT_0103 = std::conj(lpp_221)*std::conj(U_utR_21);
    const complex_t IT_0104 = (complex_t{0, 1})*(IT_0098 + IT_0099 + IT_0100 +
       -IT_0101 + -IT_0102 + -IT_0103);
    const complex_t IT_0105 = (complex_t{0, -1})*gw*gwuR*U_utR_01;
    const complex_t IT_0106 = IT_0104*IT_0105;
    const complex_t IT_0107 = std::pow(m_utR2, 2);
    const complex_t IT_0108 = std::pow(2*s_12 + (complex_t{0, 1})*GcRt*m_utR2 
      + IT_0008 + IT_0040 + -IT_0107 + reg_prop, -1);
    const complex_t IT_0109 = (complex_t{0, 1})*IT_0106*IT_0108;
    const complex_t IT_0110 = std::conj(lpp_012)*std::conj(U_utR_02);
    const complex_t IT_0111 = std::conj(lpp_112)*std::conj(U_utR_12);
    const complex_t IT_0112 = std::conj(lpp_212)*std::conj(U_utR_22);
    const complex_t IT_0113 = std::conj(lpp_021)*std::conj(U_utR_02);
    const complex_t IT_0114 = std::conj(lpp_121)*std::conj(U_utR_12);
    const complex_t IT_0115 = std::conj(lpp_221)*std::conj(U_utR_22);
    const complex_t IT_0116 = (complex_t{0, 1})*(IT_0110 + IT_0111 + IT_0112 +
       -IT_0113 + -IT_0114 + -IT_0115);
    const complex_t IT_0117 = (complex_t{0, -1})*gw*gwuR*U_utR_02;
    const complex_t IT_0118 = IT_0116*IT_0117;
    const complex_t IT_0119 = std::pow(m_utR3, 2);
    const complex_t IT_0120 = std::pow(2*s_12 + (complex_t{0, 1})*GtRt*m_utR3 
      + IT_0008 + IT_0040 + -IT_0119 + reg_prop, -1);
    const complex_t IT_0121 = (complex_t{0, 1})*IT_0118*IT_0120;
    const complex_t IT_0122 = -IT_0097 + -IT_0109 + -IT_0121;
    const complex_t IT_0123 = -IT_0055;
    const complex_t IT_0124 = IT_0057 + IT_0123;
    const complex_t IT_0125 = 8*IT_0124;
    const complex_t IT_0126 = std::conj(IT_0122)*IT_0125;
    const complex_t IT_0127 = (-0.5)*IT_0042 + (-0.5)*IT_0051 + (-0.5)*IT_0070;
    const complex_t IT_0128 = (-8)*IT_0124;
    const complex_t IT_0129 = std::conj(IT_0127)*IT_0128;
    const complex_t IT_0130 = IT_0126 + IT_0129;
    const complex_t IT_0131 = (-0.5)*IT_0011 + (-0.5)*IT_0021 + (-0.5)*IT_0031;
    const complex_t IT_0132 = 2*IT_0059;
    const complex_t IT_0133 = std::conj(IT_0122) + std::conj(IT_0131);
    const complex_t IT_0134 = 64*IT_0057;
    const complex_t IT_0135 = (-96)*IT_0073;
    const complex_t IT_0136 = IT_0127*IT_0128;
    const complex_t IT_0137 = (-2)*IT_0059;
    const complex_t IT_0138 = std::conj(IT_0032)*(IT_0052*IT_0061 + IT_0071
      *IT_0076 + IT_0032*(IT_0081 + IT_0083) + IT_0084*IT_0085) + (IT_0052 +
       IT_0071)*IT_0130 + std::conj(IT_0084)*(IT_0081*IT_0084 + IT_0032*IT_0085 
      + IT_0128*(IT_0122 + IT_0131)) + IT_0125*(std::conj(IT_0084)*IT_0127 +
       IT_0084*std::conj(IT_0127) + (std::conj(IT_0052) + std::conj(IT_0071))*
      (IT_0122 + IT_0131) + (IT_0052 + IT_0071)*std::conj(IT_0131)) + (IT_0127
      *std::conj(IT_0127) + IT_0122*(std::conj(IT_0122) + std::conj(IT_0131)) +
       IT_0131*(std::conj(IT_0122) + std::conj(IT_0131)))*IT_0132 + IT_0084
      *IT_0128*IT_0133 + std::conj(IT_0071)*(IT_0032*IT_0076 + IT_0071*IT_0083 +
       IT_0052*((-96)*IT_0072 + IT_0079 + IT_0082 + IT_0134 + IT_0135) + IT_0136
      ) + std::conj(IT_0052)*(IT_0032*IT_0061 + IT_0052*IT_0083 + IT_0071*((-96)
      *IT_0072 + IT_0079 + IT_0082 + IT_0134 + IT_0135) + IT_0136) + (std::conj
      (IT_0127)*(IT_0122 + IT_0131) + IT_0127*IT_0133)*IT_0137;
    return IT_0138;
}
} // End of namespace brparity
