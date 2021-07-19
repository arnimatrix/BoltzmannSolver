#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XX_df1c_to_uf3_df2.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_df1c_to_uf3_df2(
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
    const complex_t IT_0000 = std::conj(lpp_001)*std::conj(U_utR_00);
    const complex_t IT_0001 = std::conj(lpp_101)*std::conj(U_utR_10);
    const complex_t IT_0002 = std::conj(lpp_201)*std::conj(U_utR_20);
    const complex_t IT_0003 = std::conj(lpp_010)*std::conj(U_utR_00);
    const complex_t IT_0004 = std::conj(lpp_110)*std::conj(U_utR_10);
    const complex_t IT_0005 = std::conj(lpp_210)*std::conj(U_utR_20);
    const complex_t IT_0006 = (complex_t{0, 1})*(IT_0000 + IT_0001 + IT_0002 +
       -IT_0003 + -IT_0004 + -IT_0005);
    const complex_t IT_0007 = (complex_t{0, -1})*gw*gwuR*U_utR_20;
    const complex_t IT_0008 = IT_0006*IT_0007;
    const complex_t IT_0009 = std::pow(m_t, 2);
    const complex_t IT_0010 = std::pow(m_XX, 2);
    const complex_t IT_0011 = std::pow(m_utR1, 2);
    const complex_t IT_0012 = std::pow((-2)*s_13 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0009 + IT_0010 + -IT_0011 + reg_prop, -1);
    const complex_t IT_0013 = (complex_t{0, 1})*IT_0008*IT_0012;
    const complex_t IT_0014 = std::conj(lpp_001)*std::conj(U_utR_01);
    const complex_t IT_0015 = std::conj(lpp_101)*std::conj(U_utR_11);
    const complex_t IT_0016 = std::conj(lpp_201)*std::conj(U_utR_21);
    const complex_t IT_0017 = std::conj(lpp_010)*std::conj(U_utR_01);
    const complex_t IT_0018 = std::conj(lpp_110)*std::conj(U_utR_11);
    const complex_t IT_0019 = std::conj(lpp_210)*std::conj(U_utR_21);
    const complex_t IT_0020 = (complex_t{0, 1})*(IT_0014 + IT_0015 + IT_0016 +
       -IT_0017 + -IT_0018 + -IT_0019);
    const complex_t IT_0021 = (complex_t{0, -1})*gw*gwuR*U_utR_21;
    const complex_t IT_0022 = IT_0020*IT_0021;
    const complex_t IT_0023 = std::pow(m_utR2, 2);
    const complex_t IT_0024 = std::pow((-2)*s_13 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0009 + IT_0010 + -IT_0023 + reg_prop, -1);
    const complex_t IT_0025 = (complex_t{0, 1})*IT_0022*IT_0024;
    const complex_t IT_0026 = std::conj(lpp_001)*std::conj(U_utR_02);
    const complex_t IT_0027 = std::conj(lpp_101)*std::conj(U_utR_12);
    const complex_t IT_0028 = std::conj(lpp_201)*std::conj(U_utR_22);
    const complex_t IT_0029 = std::conj(lpp_010)*std::conj(U_utR_02);
    const complex_t IT_0030 = std::conj(lpp_110)*std::conj(U_utR_12);
    const complex_t IT_0031 = std::conj(lpp_210)*std::conj(U_utR_22);
    const complex_t IT_0032 = (complex_t{0, 1})*(IT_0026 + IT_0027 + IT_0028 +
       -IT_0029 + -IT_0030 + -IT_0031);
    const complex_t IT_0033 = (complex_t{0, -1})*gw*gwuR*U_utR_22;
    const complex_t IT_0034 = IT_0032*IT_0033;
    const complex_t IT_0035 = std::pow(m_utR3, 2);
    const complex_t IT_0036 = std::pow((-2)*s_13 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0009 + IT_0010 + -IT_0035 + reg_prop, -1);
    const complex_t IT_0037 = (complex_t{0, 1})*IT_0034*IT_0036;
    const complex_t IT_0038 = 0.25*IT_0013 + 0.25*IT_0025 + 0.25*IT_0037;
    const complex_t IT_0039 = s_12*s_34;
    const complex_t IT_0040 = (-64)*IT_0039;
    const complex_t IT_0041 = s_13*s_24;
    const complex_t IT_0042 = 128*IT_0041;
    const complex_t IT_0043 = s_14*s_23;
    const complex_t IT_0044 = 128*IT_0043;
    const complex_t IT_0045 = IT_0040 + IT_0042 + IT_0044;
    const complex_t IT_0046 = m_d*m_s*m_t*m_XX;
    const complex_t IT_0047 = 192*IT_0046;
    const complex_t IT_0048 = IT_0045 + IT_0047;
    const complex_t IT_0049 = std::conj(lpp_201)*std::conj(U_dtR_00);
    const complex_t IT_0050 = std::conj(lpp_210)*std::conj(U_dtR_00);
    const complex_t IT_0051 = std::conj(lpp_212)*std::conj(U_dtR_20);
    const complex_t IT_0052 = std::conj(lpp_221)*std::conj(U_dtR_20);
    const complex_t IT_0053 = (complex_t{0, 1})*(IT_0049 + -IT_0050 + -IT_0051
       + IT_0052);
    const complex_t IT_0054 = (complex_t{0, -1})*gw*gwdR*U_dtR_00;
    const complex_t IT_0055 = IT_0053*IT_0054;
    const complex_t IT_0056 = std::pow(m_d, 2);
    const complex_t IT_0057 = std::pow(m_dtR1, 2);
    const complex_t IT_0058 = std::pow(2*s_12 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0010 + IT_0056 + -IT_0057 + reg_prop, -1);
    const complex_t IT_0059 = (complex_t{0, 1})*IT_0055*IT_0058;
    const complex_t IT_0060 = std::conj(lpp_201)*std::conj(U_dtR_01);
    const complex_t IT_0061 = std::conj(lpp_210)*std::conj(U_dtR_01);
    const complex_t IT_0062 = std::conj(lpp_212)*std::conj(U_dtR_21);
    const complex_t IT_0063 = std::conj(lpp_221)*std::conj(U_dtR_21);
    const complex_t IT_0064 = (complex_t{0, 1})*(IT_0060 + -IT_0061 + -IT_0062
       + IT_0063);
    const complex_t IT_0065 = (complex_t{0, -1})*gw*gwdR*U_dtR_01;
    const complex_t IT_0066 = IT_0064*IT_0065;
    const complex_t IT_0067 = std::pow(m_dtR2, 2);
    const complex_t IT_0068 = std::pow(2*s_12 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0010 + IT_0056 + -IT_0067 + reg_prop, -1);
    const complex_t IT_0069 = (complex_t{0, 1})*IT_0066*IT_0068;
    const complex_t IT_0070 = std::conj(lpp_201)*std::conj(U_dtR_02);
    const complex_t IT_0071 = std::conj(lpp_210)*std::conj(U_dtR_02);
    const complex_t IT_0072 = std::conj(lpp_212)*std::conj(U_dtR_22);
    const complex_t IT_0073 = std::conj(lpp_221)*std::conj(U_dtR_22);
    const complex_t IT_0074 = (complex_t{0, 1})*(IT_0070 + -IT_0071 + -IT_0072
       + IT_0073);
    const complex_t IT_0075 = (complex_t{0, -1})*gw*gwdR*U_dtR_02;
    const complex_t IT_0076 = IT_0074*IT_0075;
    const complex_t IT_0077 = std::pow(m_dtR3, 2);
    const complex_t IT_0078 = std::pow(2*s_12 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0010 + IT_0056 + -IT_0077 + reg_prop, -1);
    const complex_t IT_0079 = (complex_t{0, 1})*IT_0076*IT_0078;
    const complex_t IT_0080 = -IT_0059 + -IT_0069 + -IT_0079;
    const complex_t IT_0081 = (-0.5)*IT_0013 + (-0.5)*IT_0025 + (-0.5)*IT_0037;
    const complex_t IT_0082 = std::conj(lpp_201)*std::conj(U_dtR_10);
    const complex_t IT_0083 = std::conj(lpp_202)*std::conj(U_dtR_20);
    const complex_t IT_0084 = std::conj(lpp_210)*std::conj(U_dtR_10);
    const complex_t IT_0085 = std::conj(lpp_220)*std::conj(U_dtR_20);
    const complex_t IT_0086 = (complex_t{0, -1})*(IT_0082 + IT_0083 + -IT_0084
       + -IT_0085);
    const complex_t IT_0087 = (complex_t{0, -1})*gw*gwdR*U_dtR_10;
    const complex_t IT_0088 = IT_0086*IT_0087;
    const complex_t IT_0089 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0009 + IT_0056 + -IT_0057 + reg_prop, -1);
    const complex_t IT_0090 = (complex_t{0, 1})*IT_0088*IT_0089;
    const complex_t IT_0091 = std::conj(lpp_201)*std::conj(U_dtR_11);
    const complex_t IT_0092 = std::conj(lpp_202)*std::conj(U_dtR_21);
    const complex_t IT_0093 = std::conj(lpp_210)*std::conj(U_dtR_11);
    const complex_t IT_0094 = std::conj(lpp_220)*std::conj(U_dtR_21);
    const complex_t IT_0095 = (complex_t{0, -1})*(IT_0091 + IT_0092 + -IT_0093
       + -IT_0094);
    const complex_t IT_0096 = (complex_t{0, -1})*gw*gwdR*U_dtR_11;
    const complex_t IT_0097 = IT_0095*IT_0096;
    const complex_t IT_0098 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0009 + IT_0056 + -IT_0067 + reg_prop, -1);
    const complex_t IT_0099 = (complex_t{0, 1})*IT_0097*IT_0098;
    const complex_t IT_0100 = std::conj(lpp_201)*std::conj(U_dtR_12);
    const complex_t IT_0101 = std::conj(lpp_202)*std::conj(U_dtR_22);
    const complex_t IT_0102 = std::conj(lpp_210)*std::conj(U_dtR_12);
    const complex_t IT_0103 = std::conj(lpp_220)*std::conj(U_dtR_22);
    const complex_t IT_0104 = (complex_t{0, -1})*(IT_0100 + IT_0101 + -IT_0102
       + -IT_0103);
    const complex_t IT_0105 = (complex_t{0, -1})*gw*gwdR*U_dtR_12;
    const complex_t IT_0106 = IT_0104*IT_0105;
    const complex_t IT_0107 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0009 + IT_0056 + -IT_0077 + reg_prop, -1);
    const complex_t IT_0108 = (complex_t{0, 1})*IT_0106*IT_0107;
    const complex_t IT_0109 = 0.5*IT_0090 + 0.5*IT_0099 + 0.5*IT_0108;
    const complex_t IT_0110 = std::conj(IT_0081) + std::conj(IT_0109);
    const complex_t IT_0111 = 2*IT_0039;
    const complex_t IT_0112 = 0.25*IT_0013 + 0.25*IT_0025;
    const complex_t IT_0113 = 0.25*IT_0037;
    const complex_t IT_0114 = -IT_0041;
    const complex_t IT_0115 = IT_0043 + IT_0114;
    const complex_t IT_0116 = 8*IT_0115;
    const complex_t IT_0117 = (-8)*IT_0115;
    const complex_t IT_0118 = m_d*m_XX*s_34;
    const complex_t IT_0119 = (-96)*IT_0118;
    const complex_t IT_0120 = m_s*m_t*s_12;
    const complex_t IT_0121 = (-96)*IT_0120;
    const complex_t IT_0122 = (-192)*IT_0046;
    const complex_t IT_0123 = 64*IT_0043;
    const complex_t IT_0124 = IT_0040 + IT_0119 + IT_0121 + IT_0122 + IT_0123;
    const complex_t IT_0125 = IT_0045 + IT_0122;
    const complex_t IT_0126 = std::conj(IT_0081)*IT_0117;
    const complex_t IT_0127 = (-2)*IT_0039;
    const complex_t IT_0128 = 0.25*IT_0090 + 0.25*IT_0099 + 0.25*IT_0108;
    const complex_t IT_0129 = IT_0113*IT_0125;
    const complex_t IT_0130 = std::conj(IT_0112)*IT_0124;
    const complex_t IT_0131 = IT_0038*std::conj(IT_0038)*IT_0048 + (IT_0080
      *std::conj(IT_0080) + (IT_0081 + IT_0109)*IT_0110)*IT_0111 + (std::conj
      (IT_0038)*(IT_0081 + IT_0109) + IT_0038*(std::conj(IT_0081) + std::conj
      (IT_0109)) + std::conj(IT_0080)*(IT_0112 + IT_0113) + IT_0080*(std::conj
      (IT_0112) + std::conj(IT_0113)))*IT_0116 + (std::conj(IT_0038)*IT_0080 +
       IT_0038*std::conj(IT_0080) + std::conj(IT_0109)*(IT_0112 + IT_0113) + 
      (IT_0081 + IT_0109)*(std::conj(IT_0112) + std::conj(IT_0113)))*IT_0117 +
       IT_0112*(std::conj(IT_0113)*IT_0124 + std::conj(IT_0112)*IT_0125 +
       IT_0126) + (std::conj(IT_0080)*(IT_0081 + IT_0109) + IT_0080*IT_0110)
      *IT_0127 + std::conj(IT_0128)*(IT_0112*IT_0124 + IT_0125*IT_0128 + IT_0048
      *(IT_0038 + IT_0128) + IT_0129) + std::conj(IT_0113)*IT_0129 + IT_0128*
      (std::conj(IT_0038)*IT_0048 + std::conj(IT_0113)*IT_0125 + IT_0130) +
       IT_0113*(IT_0126 + IT_0130);
    return IT_0131;
}
} // End of namespace brparity
