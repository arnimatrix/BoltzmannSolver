#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XXc_df1_to_uf3c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_df1_to_uf3c_df3c(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_b = param.m_b;
    const real_t m_d = param.m_d;
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
    const complex_t lpp_002 = param.lpp_002;
    const complex_t lpp_020 = param.lpp_020;
    const complex_t lpp_102 = param.lpp_102;
    const complex_t lpp_120 = param.lpp_120;
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
    const complex_t IT_0002 = std::pow(m_utR1, 2);
    const complex_t IT_0003 = std::pow((-2)*s_13 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0000 + IT_0001 + -IT_0002 + reg_prop, -1);
    const complex_t IT_0004 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_20);
    const complex_t IT_0005 = (complex_t{0, 1})*IT_0003*IT_0004;
    const complex_t IT_0006 = lpp_020*U_utR_00;
    const complex_t IT_0007 = lpp_120*U_utR_10;
    const complex_t IT_0008 = lpp_220*U_utR_20;
    const complex_t IT_0009 = (complex_t{0, 1})*(IT_0006 + IT_0007 + IT_0008);
    const complex_t IT_0010 = lpp_002*U_utR_00;
    const complex_t IT_0011 = lpp_102*U_utR_10;
    const complex_t IT_0012 = lpp_202*U_utR_20;
    const complex_t IT_0013 = (complex_t{0, 1})*(IT_0010 + IT_0011 + IT_0012);
    const complex_t IT_0014 = IT_0009 + -IT_0013;
    const complex_t IT_0015 = IT_0005*IT_0014;
    const complex_t IT_0016 = std::pow(m_utR2, 2);
    const complex_t IT_0017 = std::pow((-2)*s_13 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0000 + IT_0001 + -IT_0016 + reg_prop, -1);
    const complex_t IT_0018 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_21);
    const complex_t IT_0019 = (complex_t{0, 1})*IT_0017*IT_0018;
    const complex_t IT_0020 = lpp_020*U_utR_01;
    const complex_t IT_0021 = lpp_120*U_utR_11;
    const complex_t IT_0022 = lpp_220*U_utR_21;
    const complex_t IT_0023 = (complex_t{0, 1})*(IT_0020 + IT_0021 + IT_0022);
    const complex_t IT_0024 = IT_0019*IT_0023;
    const complex_t IT_0025 = lpp_002*U_utR_01;
    const complex_t IT_0026 = lpp_102*U_utR_11;
    const complex_t IT_0027 = lpp_202*U_utR_21;
    const complex_t IT_0028 = (complex_t{0, 1})*(IT_0025 + IT_0026 + IT_0027);
    const complex_t IT_0029 = IT_0019*IT_0028;
    const complex_t IT_0030 = std::pow(m_utR3, 2);
    const complex_t IT_0031 = std::pow((-2)*s_13 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0000 + IT_0001 + -IT_0030 + reg_prop, -1);
    const complex_t IT_0032 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_22);
    const complex_t IT_0033 = (complex_t{0, 1})*IT_0031*IT_0032;
    const complex_t IT_0034 = lpp_020*U_utR_02;
    const complex_t IT_0035 = lpp_120*U_utR_12;
    const complex_t IT_0036 = lpp_220*U_utR_22;
    const complex_t IT_0037 = (complex_t{0, 1})*(IT_0034 + IT_0035 + IT_0036);
    const complex_t IT_0038 = IT_0033*IT_0037;
    const complex_t IT_0039 = lpp_002*U_utR_02;
    const complex_t IT_0040 = lpp_102*U_utR_12;
    const complex_t IT_0041 = lpp_202*U_utR_22;
    const complex_t IT_0042 = (complex_t{0, 1})*(IT_0039 + IT_0040 + IT_0041);
    const complex_t IT_0043 = IT_0033*IT_0042;
    const complex_t IT_0044 = (-0.25)*IT_0015 + (-0.25)*IT_0024 + 0.25*IT_0029
       + (-0.25)*IT_0038 + 0.25*IT_0043;
    const complex_t IT_0045 = 0.25*IT_0015 + 0.25*IT_0024 + (-0.25)*IT_0029 +
       0.25*IT_0038 + (-0.25)*IT_0043;
    const complex_t IT_0046 = lpp_202*U_dtR_00;
    const complex_t IT_0047 = lpp_212*U_dtR_10;
    const complex_t IT_0048 = lpp_220*U_dtR_00;
    const complex_t IT_0049 = lpp_221*U_dtR_10;
    const complex_t IT_0050 = (complex_t{0, 1})*(IT_0046 + IT_0047 + -IT_0048 
      + -IT_0049);
    const complex_t IT_0051 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_00);
    const complex_t IT_0052 = IT_0050*IT_0051;
    const complex_t IT_0053 = std::pow(m_d, 2);
    const complex_t IT_0054 = std::pow(m_dtR1, 2);
    const complex_t IT_0055 = std::pow(2*s_12 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0001 + IT_0053 + -IT_0054 + reg_prop, -1);
    const complex_t IT_0056 = (complex_t{0, 1})*IT_0052*IT_0055;
    const complex_t IT_0057 = lpp_202*U_dtR_01;
    const complex_t IT_0058 = lpp_212*U_dtR_11;
    const complex_t IT_0059 = lpp_220*U_dtR_01;
    const complex_t IT_0060 = lpp_221*U_dtR_11;
    const complex_t IT_0061 = (complex_t{0, 1})*(IT_0057 + IT_0058 + -IT_0059 
      + -IT_0060);
    const complex_t IT_0062 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_01);
    const complex_t IT_0063 = IT_0061*IT_0062;
    const complex_t IT_0064 = std::pow(m_dtR2, 2);
    const complex_t IT_0065 = std::pow(2*s_12 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0001 + IT_0053 + -IT_0064 + reg_prop, -1);
    const complex_t IT_0066 = (complex_t{0, 1})*IT_0063*IT_0065;
    const complex_t IT_0067 = lpp_202*U_dtR_02;
    const complex_t IT_0068 = lpp_212*U_dtR_12;
    const complex_t IT_0069 = lpp_220*U_dtR_02;
    const complex_t IT_0070 = lpp_221*U_dtR_12;
    const complex_t IT_0071 = (complex_t{0, 1})*(IT_0067 + IT_0068 + -IT_0069 
      + -IT_0070);
    const complex_t IT_0072 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_02);
    const complex_t IT_0073 = IT_0071*IT_0072;
    const complex_t IT_0074 = std::pow(m_dtR3, 2);
    const complex_t IT_0075 = std::pow(2*s_12 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0001 + IT_0053 + -IT_0074 + reg_prop, -1);
    const complex_t IT_0076 = (complex_t{0, 1})*IT_0073*IT_0075;
    const complex_t IT_0077 = -IT_0056 + -IT_0066 + -IT_0076;
    const complex_t IT_0078 = s_14*s_23;
    const complex_t IT_0079 = s_13*s_24;
    const complex_t IT_0080 = -IT_0079;
    const complex_t IT_0081 = IT_0078 + IT_0080;
    const complex_t IT_0082 = (-8)*IT_0081;
    const complex_t IT_0083 = IT_0077*IT_0082;
    const complex_t IT_0084 = (-0.5)*IT_0015 + (-0.5)*IT_0024 + 0.5*IT_0029 + 
      (-0.5)*IT_0038 + 0.5*IT_0043;
    const complex_t IT_0085 = IT_0082*IT_0084;
    const complex_t IT_0086 = IT_0083 + IT_0085;
    const complex_t IT_0087 = std::conj(IT_0077) + std::conj(IT_0084);
    const complex_t IT_0088 = lpp_201*U_dtR_10;
    const complex_t IT_0089 = lpp_202*U_dtR_20;
    const complex_t IT_0090 = lpp_210*U_dtR_10;
    const complex_t IT_0091 = lpp_220*U_dtR_20;
    const complex_t IT_0092 = (complex_t{0, -1})*(IT_0088 + IT_0089 + -IT_0090
       + -IT_0091);
    const complex_t IT_0093 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_20);
    const complex_t IT_0094 = IT_0092*IT_0093;
    const complex_t IT_0095 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0000 + IT_0053 + -IT_0054 + reg_prop, -1);
    const complex_t IT_0096 = (complex_t{0, 1})*IT_0094*IT_0095;
    const complex_t IT_0097 = lpp_201*U_dtR_11;
    const complex_t IT_0098 = lpp_202*U_dtR_21;
    const complex_t IT_0099 = lpp_210*U_dtR_11;
    const complex_t IT_0100 = lpp_220*U_dtR_21;
    const complex_t IT_0101 = (complex_t{0, -1})*(IT_0097 + IT_0098 + -IT_0099
       + -IT_0100);
    const complex_t IT_0102 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_21);
    const complex_t IT_0103 = IT_0101*IT_0102;
    const complex_t IT_0104 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0000 + IT_0053 + -IT_0064 + reg_prop, -1);
    const complex_t IT_0105 = (complex_t{0, 1})*IT_0103*IT_0104;
    const complex_t IT_0106 = lpp_201*U_dtR_12;
    const complex_t IT_0107 = lpp_202*U_dtR_22;
    const complex_t IT_0108 = lpp_210*U_dtR_12;
    const complex_t IT_0109 = lpp_220*U_dtR_22;
    const complex_t IT_0110 = (complex_t{0, -1})*(IT_0106 + IT_0107 + -IT_0108
       + -IT_0109);
    const complex_t IT_0111 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_22);
    const complex_t IT_0112 = IT_0110*IT_0111;
    const complex_t IT_0113 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0000 + IT_0053 + -IT_0074 + reg_prop, -1);
    const complex_t IT_0114 = (complex_t{0, 1})*IT_0112*IT_0113;
    const complex_t IT_0115 = 0.25*IT_0096 + 0.25*IT_0105 + 0.25*IT_0114;
    const complex_t IT_0116 = 8*IT_0081;
    const complex_t IT_0117 = IT_0115*IT_0116;
    const complex_t IT_0118 = 128*IT_0079;
    const complex_t IT_0119 = 128*IT_0078;
    const complex_t IT_0120 = s_12*s_34;
    const complex_t IT_0121 = (-64)*IT_0120;
    const complex_t IT_0122 = IT_0118 + IT_0119 + IT_0121;
    const complex_t IT_0123 = m_b*m_d*m_t*m_XX;
    const complex_t IT_0124 = (-192)*IT_0123;
    const complex_t IT_0125 = IT_0122 + IT_0124;
    const complex_t IT_0126 = (-128)*IT_0079;
    const complex_t IT_0127 = (-128)*IT_0078;
    const complex_t IT_0128 = 64*IT_0120;
    const complex_t IT_0129 = IT_0126 + IT_0127 + IT_0128;
    const complex_t IT_0130 = 192*IT_0123;
    const complex_t IT_0131 = 0.5*IT_0096 + 0.5*IT_0105 + 0.5*IT_0114;
    const complex_t IT_0132 = IT_0082*std::conj(IT_0131);
    const complex_t IT_0133 = (-0.25)*IT_0096 + (-0.25)*IT_0105 + (-0.25)
      *IT_0114;
    const complex_t IT_0134 = IT_0122 + IT_0130;
    const complex_t IT_0135 = IT_0124 + IT_0129;
    const complex_t IT_0136 = IT_0116*std::conj(IT_0131);
    const complex_t IT_0137 = std::conj(IT_0077)*IT_0082;
    const complex_t IT_0138 = IT_0082*std::conj(IT_0084);
    const complex_t IT_0139 = 2*IT_0120;
    const complex_t IT_0140 = (-2)*IT_0120;
    const complex_t IT_0141 = IT_0116*std::conj(IT_0133);
    const complex_t IT_0142 = std::conj(IT_0115)*IT_0116;
    const complex_t IT_0143 = std::conj(IT_0131)*IT_0140;
    const complex_t IT_0144 = (std::conj(IT_0044) + std::conj(IT_0045))
      *IT_0086 + IT_0087*IT_0117 + IT_0115*(std::conj(IT_0115)*IT_0125 +
       std::conj(IT_0045)*(IT_0129 + IT_0130) + IT_0132) + IT_0133*(IT_0087
      *IT_0116 + IT_0132 + std::conj(IT_0133)*IT_0134 + std::conj(IT_0044)
      *IT_0135) + IT_0045*(std::conj(IT_0045)*IT_0125 + std::conj(IT_0115)*
      (IT_0129 + IT_0130) + IT_0136 + IT_0137 + IT_0138) + IT_0044*(std::conj
      (IT_0044)*IT_0134 + std::conj(IT_0133)*IT_0135 + IT_0136 + IT_0137 +
       IT_0138) + IT_0131*((std::conj(IT_0044) + std::conj(IT_0045))*IT_0116 +
       IT_0082*(std::conj(IT_0115) + std::conj(IT_0133)) + std::conj(IT_0131)
      *IT_0139 + IT_0087*IT_0140) + IT_0077*((std::conj(IT_0077) + std::conj
      (IT_0084))*IT_0139 + IT_0141 + IT_0142 + IT_0143) + IT_0084*(IT_0087
      *IT_0139 + IT_0141 + IT_0142 + IT_0143);
    return IT_0144;
}
} // End of namespace brparity
