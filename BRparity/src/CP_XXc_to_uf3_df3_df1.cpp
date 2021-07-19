#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XXc_to_uf3_df3_df1.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_to_uf3_df3_df1(
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
    const complex_t IT_0000 = std::conj(lpp_202)*std::conj(U_dtR_00);
    const complex_t IT_0001 = std::conj(lpp_212)*std::conj(U_dtR_10);
    const complex_t IT_0002 = std::conj(lpp_220)*std::conj(U_dtR_00);
    const complex_t IT_0003 = std::conj(lpp_221)*std::conj(U_dtR_10);
    const complex_t IT_0004 = (complex_t{0, 1})*(IT_0000 + IT_0001 + -IT_0002 
      + -IT_0003);
    const complex_t IT_0005 = (complex_t{0, -1})*gw*gwdR*U_dtR_00;
    const complex_t IT_0006 = IT_0004*IT_0005;
    const complex_t IT_0007 = std::pow(m_b, 2);
    const complex_t IT_0008 = std::pow(m_t, 2);
    const complex_t IT_0009 = std::pow(m_dtR1, 2);
    const complex_t IT_0010 = std::pow(2*s_23 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0007 + IT_0008 + -IT_0009 + reg_prop, -1);
    const complex_t IT_0011 = (complex_t{0, 1})*IT_0006*IT_0010;
    const complex_t IT_0012 = std::conj(lpp_202)*std::conj(U_dtR_01);
    const complex_t IT_0013 = std::conj(lpp_212)*std::conj(U_dtR_11);
    const complex_t IT_0014 = std::conj(lpp_220)*std::conj(U_dtR_01);
    const complex_t IT_0015 = std::conj(lpp_221)*std::conj(U_dtR_11);
    const complex_t IT_0016 = (complex_t{0, 1})*(IT_0012 + IT_0013 + -IT_0014 
      + -IT_0015);
    const complex_t IT_0017 = (complex_t{0, -1})*gw*gwdR*U_dtR_01;
    const complex_t IT_0018 = IT_0016*IT_0017;
    const complex_t IT_0019 = std::pow(m_dtR2, 2);
    const complex_t IT_0020 = std::pow(2*s_23 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0007 + IT_0008 + -IT_0019 + reg_prop, -1);
    const complex_t IT_0021 = (complex_t{0, 1})*IT_0018*IT_0020;
    const complex_t IT_0022 = std::conj(lpp_202)*std::conj(U_dtR_02);
    const complex_t IT_0023 = std::conj(lpp_212)*std::conj(U_dtR_12);
    const complex_t IT_0024 = std::conj(lpp_220)*std::conj(U_dtR_02);
    const complex_t IT_0025 = std::conj(lpp_221)*std::conj(U_dtR_12);
    const complex_t IT_0026 = (complex_t{0, 1})*(IT_0022 + IT_0023 + -IT_0024 
      + -IT_0025);
    const complex_t IT_0027 = (complex_t{0, -1})*gw*gwdR*U_dtR_02;
    const complex_t IT_0028 = IT_0026*IT_0027;
    const complex_t IT_0029 = std::pow(m_dtR3, 2);
    const complex_t IT_0030 = std::pow(2*s_23 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0007 + IT_0008 + -IT_0029 + reg_prop, -1);
    const complex_t IT_0031 = (complex_t{0, 1})*IT_0028*IT_0030;
    const complex_t IT_0032 = 0.5*IT_0011 + 0.5*IT_0021 + 0.5*IT_0031;
    const complex_t IT_0033 = s_12*s_34;
    const complex_t IT_0034 = 12*IT_0033;
    const complex_t IT_0035 = std::conj(lpp_002)*std::conj(U_utR_00);
    const complex_t IT_0036 = std::conj(lpp_102)*std::conj(U_utR_10);
    const complex_t IT_0037 = std::conj(lpp_202)*std::conj(U_utR_20);
    const complex_t IT_0038 = std::conj(lpp_020)*std::conj(U_utR_00);
    const complex_t IT_0039 = std::conj(lpp_120)*std::conj(U_utR_10);
    const complex_t IT_0040 = std::conj(lpp_220)*std::conj(U_utR_20);
    const complex_t IT_0041 = (complex_t{0, 1})*(IT_0035 + IT_0036 + IT_0037 +
       -IT_0038 + -IT_0039 + -IT_0040);
    const complex_t IT_0042 = (complex_t{0, -1})*gw*gwuR*U_utR_20;
    const complex_t IT_0043 = IT_0041*IT_0042;
    const complex_t IT_0044 = std::pow(m_XX, 2);
    const complex_t IT_0045 = std::pow(m_utR1, 2);
    const complex_t IT_0046 = std::pow((-2)*s_12 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0008 + IT_0044 + -IT_0045 + reg_prop, -1);
    const complex_t IT_0047 = (complex_t{0, 1})*IT_0043*IT_0046;
    const complex_t IT_0048 = std::conj(lpp_002)*std::conj(U_utR_01);
    const complex_t IT_0049 = std::conj(lpp_102)*std::conj(U_utR_11);
    const complex_t IT_0050 = std::conj(lpp_202)*std::conj(U_utR_21);
    const complex_t IT_0051 = std::conj(lpp_020)*std::conj(U_utR_01);
    const complex_t IT_0052 = std::conj(lpp_120)*std::conj(U_utR_11);
    const complex_t IT_0053 = std::conj(lpp_220)*std::conj(U_utR_21);
    const complex_t IT_0054 = (complex_t{0, 1})*(IT_0048 + IT_0049 + IT_0050 +
       -IT_0051 + -IT_0052 + -IT_0053);
    const complex_t IT_0055 = (complex_t{0, -1})*gw*gwuR*U_utR_21;
    const complex_t IT_0056 = IT_0054*IT_0055;
    const complex_t IT_0057 = std::pow(m_utR2, 2);
    const complex_t IT_0058 = std::pow((-2)*s_12 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0008 + IT_0044 + -IT_0057 + reg_prop, -1);
    const complex_t IT_0059 = (complex_t{0, 1})*IT_0056*IT_0058;
    const complex_t IT_0060 = std::conj(lpp_002)*std::conj(U_utR_02);
    const complex_t IT_0061 = std::conj(lpp_102)*std::conj(U_utR_12);
    const complex_t IT_0062 = std::conj(lpp_202)*std::conj(U_utR_22);
    const complex_t IT_0063 = std::conj(lpp_020)*std::conj(U_utR_02);
    const complex_t IT_0064 = std::conj(lpp_120)*std::conj(U_utR_12);
    const complex_t IT_0065 = std::conj(lpp_220)*std::conj(U_utR_22);
    const complex_t IT_0066 = (complex_t{0, 1})*(IT_0060 + IT_0061 + IT_0062 +
       -IT_0063 + -IT_0064 + -IT_0065);
    const complex_t IT_0067 = (complex_t{0, -1})*gw*gwuR*U_utR_22;
    const complex_t IT_0068 = IT_0066*IT_0067;
    const complex_t IT_0069 = std::pow(m_utR3, 2);
    const complex_t IT_0070 = std::pow((-2)*s_12 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0008 + IT_0044 + -IT_0069 + reg_prop, -1);
    const complex_t IT_0071 = (complex_t{0, 1})*IT_0068*IT_0070;
    const complex_t IT_0072 = IT_0047 + IT_0059 + IT_0071;
    const complex_t IT_0073 = std::conj(lpp_201)*std::conj(U_dtR_10);
    const complex_t IT_0074 = std::conj(lpp_202)*std::conj(U_dtR_20);
    const complex_t IT_0075 = std::conj(lpp_210)*std::conj(U_dtR_10);
    const complex_t IT_0076 = std::conj(lpp_220)*std::conj(U_dtR_20);
    const complex_t IT_0077 = (complex_t{0, -1})*(IT_0073 + IT_0074 + -IT_0075
       + -IT_0076);
    const complex_t IT_0078 = (complex_t{0, -1})*gw*gwdR*U_dtR_20;
    const complex_t IT_0079 = IT_0077*IT_0078;
    const complex_t IT_0080 = std::pow((-2)*s_13 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0007 + -IT_0009 + IT_0044 + reg_prop, -1);
    const complex_t IT_0081 = (complex_t{0, 1})*IT_0079*IT_0080;
    const complex_t IT_0082 = std::conj(lpp_201)*std::conj(U_dtR_11);
    const complex_t IT_0083 = std::conj(lpp_202)*std::conj(U_dtR_21);
    const complex_t IT_0084 = std::conj(lpp_210)*std::conj(U_dtR_11);
    const complex_t IT_0085 = std::conj(lpp_220)*std::conj(U_dtR_21);
    const complex_t IT_0086 = (complex_t{0, -1})*(IT_0082 + IT_0083 + -IT_0084
       + -IT_0085);
    const complex_t IT_0087 = (complex_t{0, -1})*gw*gwdR*U_dtR_21;
    const complex_t IT_0088 = IT_0086*IT_0087;
    const complex_t IT_0089 = std::pow((-2)*s_13 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0007 + -IT_0019 + IT_0044 + reg_prop, -1);
    const complex_t IT_0090 = (complex_t{0, 1})*IT_0088*IT_0089;
    const complex_t IT_0091 = std::conj(lpp_201)*std::conj(U_dtR_12);
    const complex_t IT_0092 = std::conj(lpp_202)*std::conj(U_dtR_22);
    const complex_t IT_0093 = std::conj(lpp_210)*std::conj(U_dtR_12);
    const complex_t IT_0094 = std::conj(lpp_220)*std::conj(U_dtR_22);
    const complex_t IT_0095 = (complex_t{0, -1})*(IT_0091 + IT_0092 + -IT_0093
       + -IT_0094);
    const complex_t IT_0096 = (complex_t{0, -1})*gw*gwdR*U_dtR_22;
    const complex_t IT_0097 = IT_0095*IT_0096;
    const complex_t IT_0098 = std::pow((-2)*s_13 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + -IT_0029 + IT_0044 + reg_prop, -1);
    const complex_t IT_0099 = (complex_t{0, 1})*IT_0097*IT_0098;
    const complex_t IT_0100 = 0.5*IT_0081 + 0.5*IT_0090 + 0.5*IT_0099;
    const complex_t IT_0101 = IT_0072 + IT_0100;
    const complex_t IT_0102 = (-12)*IT_0033;
    const complex_t IT_0103 = 0.25*IT_0011 + 0.25*IT_0021 + 0.25*IT_0031;
    const complex_t IT_0104 = (-0.25)*IT_0081 + (-0.25)*IT_0090 + (-0.25)
      *IT_0099;
    const complex_t IT_0105 = IT_0103 + IT_0104;
    const complex_t IT_0106 = s_13*s_24;
    const complex_t IT_0107 = 48*IT_0106;
    const complex_t IT_0108 = s_14*s_23;
    const complex_t IT_0109 = (-48)*IT_0108;
    const complex_t IT_0110 = IT_0107 + IT_0109;
    const complex_t IT_0111 = std::conj(IT_0072)*IT_0110;
    const complex_t IT_0112 = std::conj(IT_0032)*IT_0110;
    const complex_t IT_0113 = (-48)*IT_0106;
    const complex_t IT_0114 = 48*IT_0108;
    const complex_t IT_0115 = IT_0113 + IT_0114;
    const complex_t IT_0116 = std::conj(IT_0072)*IT_0115;
    const complex_t IT_0117 = std::conj(IT_0032)*IT_0115;
    const complex_t IT_0118 = (-768)*IT_0106;
    const complex_t IT_0119 = (-768)*IT_0108;
    const complex_t IT_0120 = 384*IT_0033;
    const complex_t IT_0121 = IT_0118 + IT_0119 + IT_0120;
    const complex_t IT_0122 = m_b*m_d*m_t*m_XX;
    const complex_t IT_0123 = (-1152)*IT_0122;
    const complex_t IT_0124 = IT_0121 + IT_0123;
    const complex_t IT_0125 = 1152*IT_0122;
    const complex_t IT_0126 = IT_0121 + IT_0125;
    const complex_t IT_0127 = 768*IT_0106;
    const complex_t IT_0128 = 768*IT_0108;
    const complex_t IT_0129 = (-384)*IT_0033;
    const complex_t IT_0130 = IT_0127 + IT_0128 + IT_0129;
    const complex_t IT_0131 = IT_0123 + IT_0130;
    const complex_t IT_0132 = IT_0125 + IT_0130;
    const complex_t IT_0133 = IT_0032*IT_0110;
    const complex_t IT_0134 = IT_0032*IT_0115;
    const complex_t IT_0135 = IT_0072*IT_0110;
    const complex_t IT_0136 = IT_0100*IT_0110;
    const complex_t IT_0137 = IT_0032*IT_0102;
    const complex_t IT_0138 = IT_0034*IT_0101;
    const complex_t IT_0139 = std::conj(IT_0032)*(IT_0032*IT_0034 + IT_0101
      *IT_0102) + IT_0105*(IT_0111 + IT_0112 + IT_0116 + IT_0117) + std::conj
      (IT_0103)*((IT_0072 + IT_0100)*IT_0115 + IT_0104*(IT_0124 + IT_0126) +
       IT_0103*(IT_0131 + IT_0132) + IT_0133 + IT_0134 + IT_0135 + IT_0136) +
       std::conj(IT_0104)*((IT_0072 + IT_0100)*IT_0115 + IT_0103*(IT_0124 +
       IT_0126) + IT_0104*(IT_0131 + IT_0132) + IT_0133 + IT_0134 + IT_0135 +
       IT_0136) + std::conj(IT_0072)*(IT_0137 + IT_0138) + std::conj(IT_0100)*
      (IT_0105*(IT_0110 + IT_0115) + IT_0137 + IT_0138);
    return IT_0139;
}
} // End of namespace brparity
