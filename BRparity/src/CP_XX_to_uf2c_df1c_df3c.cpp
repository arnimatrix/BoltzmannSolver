#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XX_to_uf2c_df1c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_to_uf2c_df1c_df3c(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_b = param.m_b;
    const real_t m_c = param.m_c;
    const real_t m_d = param.m_d;
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
    const complex_t lpp_101 = param.lpp_101;
    const complex_t lpp_102 = param.lpp_102;
    const complex_t lpp_110 = param.lpp_110;
    const complex_t lpp_112 = param.lpp_112;
    const complex_t lpp_120 = param.lpp_120;
    const complex_t lpp_121 = param.lpp_121;
    const complex_t lpp_202 = param.lpp_202;
    const complex_t lpp_220 = param.lpp_220;
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
    const complex_t IT_0000 = lpp_101*U_dtR_10;
    const complex_t IT_0001 = lpp_102*U_dtR_20;
    const complex_t IT_0002 = lpp_110*U_dtR_10;
    const complex_t IT_0003 = lpp_120*U_dtR_20;
    const complex_t IT_0004 = (complex_t{0, -1})*(IT_0000 + IT_0001 + -IT_0002
       + -IT_0003);
    const complex_t IT_0005 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_20);
    const complex_t IT_0006 = IT_0004*IT_0005;
    const complex_t IT_0007 = std::pow(m_c, 2);
    const complex_t IT_0008 = std::pow(m_d, 2);
    const complex_t IT_0009 = std::pow(m_dtR1, 2);
    const complex_t IT_0010 = std::pow(2*s_23 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0007 + IT_0008 + -IT_0009 + reg_prop, -1);
    const complex_t IT_0011 = (complex_t{0, 1})*IT_0006*IT_0010;
    const complex_t IT_0012 = lpp_101*U_dtR_11;
    const complex_t IT_0013 = lpp_102*U_dtR_21;
    const complex_t IT_0014 = lpp_110*U_dtR_11;
    const complex_t IT_0015 = lpp_120*U_dtR_21;
    const complex_t IT_0016 = (complex_t{0, -1})*(IT_0012 + IT_0013 + -IT_0014
       + -IT_0015);
    const complex_t IT_0017 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_21);
    const complex_t IT_0018 = IT_0016*IT_0017;
    const complex_t IT_0019 = std::pow(m_dtR2, 2);
    const complex_t IT_0020 = std::pow(2*s_23 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0007 + IT_0008 + -IT_0019 + reg_prop, -1);
    const complex_t IT_0021 = (complex_t{0, 1})*IT_0018*IT_0020;
    const complex_t IT_0022 = lpp_101*U_dtR_12;
    const complex_t IT_0023 = lpp_102*U_dtR_22;
    const complex_t IT_0024 = lpp_110*U_dtR_12;
    const complex_t IT_0025 = lpp_120*U_dtR_22;
    const complex_t IT_0026 = (complex_t{0, -1})*(IT_0022 + IT_0023 + -IT_0024
       + -IT_0025);
    const complex_t IT_0027 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_22);
    const complex_t IT_0028 = IT_0026*IT_0027;
    const complex_t IT_0029 = std::pow(m_dtR3, 2);
    const complex_t IT_0030 = std::pow(2*s_23 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0007 + IT_0008 + -IT_0029 + reg_prop, -1);
    const complex_t IT_0031 = (complex_t{0, 1})*IT_0028*IT_0030;
    const complex_t IT_0032 = (-0.5)*IT_0011 + (-0.5)*IT_0021 + (-0.5)*IT_0031;
    const complex_t IT_0033 = (-0.25)*IT_0011 + (-0.25)*IT_0021 + (-0.25)
      *IT_0031;
    const complex_t IT_0034 = 0.25*IT_0011 + 0.25*IT_0021 + 0.25*IT_0031;
    const complex_t IT_0035 = s_13*s_24;
    const complex_t IT_0036 = (-48)*IT_0035;
    const complex_t IT_0037 = s_14*s_23;
    const complex_t IT_0038 = 48*IT_0037;
    const complex_t IT_0039 = IT_0036 + IT_0038;
    const complex_t IT_0040 = 48*IT_0035;
    const complex_t IT_0041 = (-48)*IT_0037;
    const complex_t IT_0042 = IT_0040 + IT_0041;
    const complex_t IT_0043 = lpp_102*U_dtR_00;
    const complex_t IT_0044 = lpp_112*U_dtR_10;
    const complex_t IT_0045 = lpp_120*U_dtR_00;
    const complex_t IT_0046 = lpp_121*U_dtR_10;
    const complex_t IT_0047 = (complex_t{0, 1})*(IT_0043 + IT_0044 + -IT_0045 
      + -IT_0046);
    const complex_t IT_0048 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_00);
    const complex_t IT_0049 = IT_0047*IT_0048;
    const complex_t IT_0050 = std::pow(m_XX, 2);
    const complex_t IT_0051 = std::pow((-2)*s_13 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0008 + -IT_0009 + IT_0050 + reg_prop, -1);
    const complex_t IT_0052 = (complex_t{0, 1})*IT_0049*IT_0051;
    const complex_t IT_0053 = lpp_102*U_dtR_01;
    const complex_t IT_0054 = lpp_112*U_dtR_11;
    const complex_t IT_0055 = lpp_120*U_dtR_01;
    const complex_t IT_0056 = lpp_121*U_dtR_11;
    const complex_t IT_0057 = (complex_t{0, 1})*(IT_0053 + IT_0054 + -IT_0055 
      + -IT_0056);
    const complex_t IT_0058 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_01);
    const complex_t IT_0059 = IT_0057*IT_0058;
    const complex_t IT_0060 = std::pow((-2)*s_13 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0008 + -IT_0019 + IT_0050 + reg_prop, -1);
    const complex_t IT_0061 = (complex_t{0, 1})*IT_0059*IT_0060;
    const complex_t IT_0062 = lpp_102*U_dtR_02;
    const complex_t IT_0063 = lpp_112*U_dtR_12;
    const complex_t IT_0064 = lpp_120*U_dtR_02;
    const complex_t IT_0065 = lpp_121*U_dtR_12;
    const complex_t IT_0066 = (complex_t{0, 1})*(IT_0062 + IT_0063 + -IT_0064 
      + -IT_0065);
    const complex_t IT_0067 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_02);
    const complex_t IT_0068 = IT_0066*IT_0067;
    const complex_t IT_0069 = std::pow((-2)*s_13 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0008 + -IT_0029 + IT_0050 + reg_prop, -1);
    const complex_t IT_0070 = (complex_t{0, 1})*IT_0068*IT_0069;
    const complex_t IT_0071 = 0.25*IT_0052 + 0.25*IT_0061 + 0.25*IT_0070;
    const complex_t IT_0072 = (-0.25)*IT_0052 + (-0.25)*IT_0061 + (-0.25)
      *IT_0070;
    const complex_t IT_0073 = std::conj(IT_0071) + std::conj(IT_0072);
    const complex_t IT_0074 = s_12*s_34;
    const complex_t IT_0075 = 12*IT_0074;
    const complex_t IT_0076 = (-0.5)*IT_0052 + (-0.5)*IT_0061 + (-0.5)*IT_0070;
    const complex_t IT_0077 = lpp_002*U_utR_00;
    const complex_t IT_0078 = lpp_102*U_utR_10;
    const complex_t IT_0079 = lpp_202*U_utR_20;
    const complex_t IT_0080 = (complex_t{0, 1})*(IT_0077 + IT_0078 + IT_0079);
    const complex_t IT_0081 = -IT_0080;
    const complex_t IT_0082 = lpp_020*U_utR_00;
    const complex_t IT_0083 = lpp_120*U_utR_10;
    const complex_t IT_0084 = lpp_220*U_utR_20;
    const complex_t IT_0085 = (complex_t{0, 1})*(IT_0082 + IT_0083 + IT_0084);
    const complex_t IT_0086 = IT_0081 + IT_0085;
    const complex_t IT_0087 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_10);
    const complex_t IT_0088 = IT_0086*IT_0087;
    const complex_t IT_0089 = std::pow(m_utR1, 2);
    const complex_t IT_0090 = std::pow((-2)*s_12 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0007 + IT_0050 + -IT_0089 + reg_prop, -1);
    const complex_t IT_0091 = (complex_t{0, 1})*IT_0088*IT_0090;
    const complex_t IT_0092 = lpp_002*U_utR_01;
    const complex_t IT_0093 = lpp_102*U_utR_11;
    const complex_t IT_0094 = lpp_202*U_utR_21;
    const complex_t IT_0095 = (complex_t{0, 1})*(IT_0092 + IT_0093 + IT_0094);
    const complex_t IT_0096 = -IT_0095;
    const complex_t IT_0097 = lpp_020*U_utR_01;
    const complex_t IT_0098 = lpp_120*U_utR_11;
    const complex_t IT_0099 = lpp_220*U_utR_21;
    const complex_t IT_0100 = (complex_t{0, 1})*(IT_0097 + IT_0098 + IT_0099);
    const complex_t IT_0101 = IT_0096 + IT_0100;
    const complex_t IT_0102 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_11);
    const complex_t IT_0103 = IT_0101*IT_0102;
    const complex_t IT_0104 = std::pow(m_utR2, 2);
    const complex_t IT_0105 = std::pow((-2)*s_12 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0007 + IT_0050 + -IT_0104 + reg_prop, -1);
    const complex_t IT_0106 = (complex_t{0, 1})*IT_0103*IT_0105;
    const complex_t IT_0107 = lpp_002*U_utR_02;
    const complex_t IT_0108 = lpp_102*U_utR_12;
    const complex_t IT_0109 = lpp_202*U_utR_22;
    const complex_t IT_0110 = (complex_t{0, 1})*(IT_0107 + IT_0108 + IT_0109);
    const complex_t IT_0111 = -IT_0110;
    const complex_t IT_0112 = lpp_020*U_utR_02;
    const complex_t IT_0113 = lpp_120*U_utR_12;
    const complex_t IT_0114 = lpp_220*U_utR_22;
    const complex_t IT_0115 = (complex_t{0, 1})*(IT_0112 + IT_0113 + IT_0114);
    const complex_t IT_0116 = IT_0111 + IT_0115;
    const complex_t IT_0117 = (complex_t{0, -1})*gw*gwuR*std::conj(U_utR_12);
    const complex_t IT_0118 = IT_0116*IT_0117;
    const complex_t IT_0119 = std::pow(m_utR3, 2);
    const complex_t IT_0120 = std::pow((-2)*s_12 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0007 + IT_0050 + -IT_0119 + reg_prop, -1);
    const complex_t IT_0121 = (complex_t{0, 1})*IT_0118*IT_0120;
    const complex_t IT_0122 = -IT_0091 + -IT_0106 + -IT_0121;
    const complex_t IT_0123 = (-12)*IT_0074;
    const complex_t IT_0124 = IT_0033*IT_0042;
    const complex_t IT_0125 = IT_0034*IT_0042;
    const complex_t IT_0126 = std::conj(IT_0032)*IT_0123;
    const complex_t IT_0127 = std::conj(IT_0033)*IT_0042;
    const complex_t IT_0128 = std::conj(IT_0034)*IT_0042;
    const complex_t IT_0129 = (-768)*IT_0035;
    const complex_t IT_0130 = (-768)*IT_0037;
    const complex_t IT_0131 = 384*IT_0074;
    const complex_t IT_0132 = IT_0129 + IT_0130 + IT_0131;
    const complex_t IT_0133 = m_b*m_c*m_d*m_XX;
    const complex_t IT_0134 = (-1152)*IT_0133;
    const complex_t IT_0135 = IT_0132 + IT_0134;
    const complex_t IT_0136 = 768*IT_0035;
    const complex_t IT_0137 = 768*IT_0037;
    const complex_t IT_0138 = (-384)*IT_0074;
    const complex_t IT_0139 = IT_0136 + IT_0137 + IT_0138;
    const complex_t IT_0140 = 1152*IT_0133;
    const complex_t IT_0141 = IT_0139 + IT_0140;
    const complex_t IT_0142 = std::conj(IT_0032)*IT_0039;
    const complex_t IT_0143 = IT_0132 + IT_0140;
    const complex_t IT_0144 = IT_0039*IT_0076;
    const complex_t IT_0145 = std::conj(IT_0032)*IT_0042;
    const complex_t IT_0146 = IT_0039*std::conj(IT_0076);
    const complex_t IT_0147 = IT_0032*((std::conj(IT_0033) + std::conj(IT_0034
      ))*IT_0039 + IT_0042*IT_0073 + std::conj(IT_0032)*IT_0075 + (std::conj
      (IT_0076) + std::conj(IT_0122))*IT_0123) + std::conj(IT_0076)*(IT_0124 +
       IT_0125) + std::conj(IT_0122)*(IT_0039*(IT_0071 + IT_0072) + IT_0075
      *IT_0076 + IT_0124 + IT_0125) + IT_0076*(IT_0075*std::conj(IT_0076) +
       IT_0126 + IT_0127 + IT_0128) + IT_0122*(IT_0039*(std::conj(IT_0071) +
       std::conj(IT_0072)) + IT_0075*(std::conj(IT_0076) + std::conj(IT_0122)) +
       IT_0126 + IT_0127 + IT_0128) + IT_0034*(std::conj(IT_0072)*IT_0135 +
       std::conj(IT_0034)*IT_0141 + IT_0142) + IT_0033*(std::conj(IT_0033)*
      (IT_0134 + IT_0139) + IT_0142 + std::conj(IT_0071)*IT_0143) + IT_0073
      *IT_0144 + IT_0072*(std::conj(IT_0034)*IT_0135 + std::conj(IT_0072)
      *IT_0141 + IT_0145 + IT_0146) + IT_0071*(std::conj(IT_0071)*(IT_0134 +
       IT_0139) + std::conj(IT_0033)*IT_0143 + IT_0145 + IT_0146);
    return IT_0147;
}
} // End of namespace brparity
