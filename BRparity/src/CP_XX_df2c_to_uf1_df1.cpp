#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XX_df2c_to_uf1_df1.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_df2c_to_uf1_df1(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_d = param.m_d;
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
    const complex_t lpp_101 = param.lpp_101;
    const complex_t lpp_110 = param.lpp_110;
    const complex_t lpp_201 = param.lpp_201;
    const complex_t lpp_210 = param.lpp_210;
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
    const complex_t IT_0005 = (complex_t{0, -1})*gw*gwdR*U_dtR_00;
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
    const complex_t IT_0017 = (complex_t{0, -1})*gw*gwdR*U_dtR_01;
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
    const complex_t IT_0027 = (complex_t{0, -1})*gw*gwdR*U_dtR_02;
    const complex_t IT_0028 = IT_0026*IT_0027;
    const complex_t IT_0029 = std::pow(m_dtR3, 2);
    const complex_t IT_0030 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + IT_0008 + -IT_0029 + reg_prop, -1);
    const complex_t IT_0031 = (complex_t{0, 1})*IT_0028*IT_0030;
    const complex_t IT_0032 = (-0.5)*IT_0011 + (-0.5)*IT_0021 + (-0.5)*IT_0031;
    const complex_t IT_0033 = std::conj(lpp_001)*std::conj(U_utR_00);
    const complex_t IT_0034 = std::conj(lpp_101)*std::conj(U_utR_10);
    const complex_t IT_0035 = std::conj(lpp_201)*std::conj(U_utR_20);
    const complex_t IT_0036 = std::conj(lpp_010)*std::conj(U_utR_00);
    const complex_t IT_0037 = std::conj(lpp_110)*std::conj(U_utR_10);
    const complex_t IT_0038 = std::conj(lpp_210)*std::conj(U_utR_20);
    const complex_t IT_0039 = (complex_t{0, 1})*(IT_0033 + IT_0034 + IT_0035 +
       -IT_0036 + -IT_0037 + -IT_0038);
    const complex_t IT_0040 = (complex_t{0, -1})*gw*gwuR*U_utR_00;
    const complex_t IT_0041 = IT_0039*IT_0040;
    const complex_t IT_0042 = std::pow(m_XX, 2);
    const complex_t IT_0043 = std::pow(m_utR1, 2);
    const complex_t IT_0044 = std::pow((-2)*s_13 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0008 + IT_0042 + -IT_0043 + reg_prop, -1);
    const complex_t IT_0045 = (complex_t{0, 1})*IT_0041*IT_0044;
    const complex_t IT_0046 = std::conj(lpp_001)*std::conj(U_utR_01);
    const complex_t IT_0047 = std::conj(lpp_101)*std::conj(U_utR_11);
    const complex_t IT_0048 = std::conj(lpp_201)*std::conj(U_utR_21);
    const complex_t IT_0049 = std::conj(lpp_010)*std::conj(U_utR_01);
    const complex_t IT_0050 = std::conj(lpp_110)*std::conj(U_utR_11);
    const complex_t IT_0051 = std::conj(lpp_210)*std::conj(U_utR_21);
    const complex_t IT_0052 = (complex_t{0, 1})*(IT_0046 + IT_0047 + IT_0048 +
       -IT_0049 + -IT_0050 + -IT_0051);
    const complex_t IT_0053 = (complex_t{0, -1})*gw*gwuR*U_utR_01;
    const complex_t IT_0054 = IT_0052*IT_0053;
    const complex_t IT_0055 = std::pow(m_utR2, 2);
    const complex_t IT_0056 = std::pow((-2)*s_13 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0008 + IT_0042 + -IT_0055 + reg_prop, -1);
    const complex_t IT_0057 = (complex_t{0, 1})*IT_0054*IT_0056;
    const complex_t IT_0058 = std::conj(lpp_001)*std::conj(U_utR_02);
    const complex_t IT_0059 = std::conj(lpp_101)*std::conj(U_utR_12);
    const complex_t IT_0060 = std::conj(lpp_201)*std::conj(U_utR_22);
    const complex_t IT_0061 = std::conj(lpp_010)*std::conj(U_utR_02);
    const complex_t IT_0062 = std::conj(lpp_110)*std::conj(U_utR_12);
    const complex_t IT_0063 = std::conj(lpp_210)*std::conj(U_utR_22);
    const complex_t IT_0064 = (complex_t{0, 1})*(IT_0058 + IT_0059 + IT_0060 +
       -IT_0061 + -IT_0062 + -IT_0063);
    const complex_t IT_0065 = (complex_t{0, -1})*gw*gwuR*U_utR_02;
    const complex_t IT_0066 = IT_0064*IT_0065;
    const complex_t IT_0067 = std::pow(m_utR3, 2);
    const complex_t IT_0068 = std::pow((-2)*s_13 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0008 + IT_0042 + -IT_0067 + reg_prop, -1);
    const complex_t IT_0069 = (complex_t{0, 1})*IT_0066*IT_0068;
    const complex_t IT_0070 = 0.5*IT_0045 + 0.5*IT_0057 + 0.5*IT_0069;
    const complex_t IT_0071 = std::conj(lpp_001)*std::conj(U_dtR_10);
    const complex_t IT_0072 = std::conj(lpp_002)*std::conj(U_dtR_20);
    const complex_t IT_0073 = std::conj(lpp_010)*std::conj(U_dtR_10);
    const complex_t IT_0074 = std::conj(lpp_020)*std::conj(U_dtR_20);
    const complex_t IT_0075 = (complex_t{0, -1})*(IT_0071 + IT_0072 + -IT_0073
       + -IT_0074);
    const complex_t IT_0076 = (complex_t{0, -1})*gw*gwdR*U_dtR_10;
    const complex_t IT_0077 = IT_0075*IT_0076;
    const complex_t IT_0078 = std::pow(2*s_12 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0007 + -IT_0009 + IT_0042 + reg_prop, -1);
    const complex_t IT_0079 = (complex_t{0, 1})*IT_0077*IT_0078;
    const complex_t IT_0080 = std::conj(lpp_001)*std::conj(U_dtR_11);
    const complex_t IT_0081 = std::conj(lpp_002)*std::conj(U_dtR_21);
    const complex_t IT_0082 = std::conj(lpp_010)*std::conj(U_dtR_11);
    const complex_t IT_0083 = std::conj(lpp_020)*std::conj(U_dtR_21);
    const complex_t IT_0084 = (complex_t{0, -1})*(IT_0080 + IT_0081 + -IT_0082
       + -IT_0083);
    const complex_t IT_0085 = (complex_t{0, -1})*gw*gwdR*U_dtR_11;
    const complex_t IT_0086 = IT_0084*IT_0085;
    const complex_t IT_0087 = std::pow(2*s_12 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0007 + -IT_0019 + IT_0042 + reg_prop, -1);
    const complex_t IT_0088 = (complex_t{0, 1})*IT_0086*IT_0087;
    const complex_t IT_0089 = std::conj(lpp_001)*std::conj(U_dtR_12);
    const complex_t IT_0090 = std::conj(lpp_002)*std::conj(U_dtR_22);
    const complex_t IT_0091 = std::conj(lpp_010)*std::conj(U_dtR_12);
    const complex_t IT_0092 = std::conj(lpp_020)*std::conj(U_dtR_22);
    const complex_t IT_0093 = (complex_t{0, -1})*(IT_0089 + IT_0090 + -IT_0091
       + -IT_0092);
    const complex_t IT_0094 = (complex_t{0, -1})*gw*gwdR*U_dtR_12;
    const complex_t IT_0095 = IT_0093*IT_0094;
    const complex_t IT_0096 = std::pow(2*s_12 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0007 + -IT_0029 + IT_0042 + reg_prop, -1);
    const complex_t IT_0097 = (complex_t{0, 1})*IT_0095*IT_0096;
    const complex_t IT_0098 = IT_0079 + IT_0088 + IT_0097;
    const complex_t IT_0099 = std::conj(IT_0070) + std::conj(IT_0098);
    const complex_t IT_0100 = s_12*s_34;
    const complex_t IT_0101 = 2*IT_0100;
    const complex_t IT_0102 = (-0.25)*IT_0045 + (-0.25)*IT_0057 + (-0.25)
      *IT_0069;
    const complex_t IT_0103 = (-0.25)*IT_0011 + (-0.25)*IT_0021 + (-0.25)
      *IT_0031;
    const complex_t IT_0104 = s_13*s_24;
    const complex_t IT_0105 = (-128)*IT_0104;
    const complex_t IT_0106 = s_14*s_23;
    const complex_t IT_0107 = (-128)*IT_0106;
    const complex_t IT_0108 = 64*IT_0100;
    const complex_t IT_0109 = IT_0105 + IT_0107 + IT_0108;
    const complex_t IT_0110 = m_d*m_s*m_u*m_XX;
    const complex_t IT_0111 = (-192)*IT_0110;
    const complex_t IT_0112 = IT_0109 + IT_0111;
    const complex_t IT_0113 = 192*IT_0110;
    const complex_t IT_0114 = 128*IT_0104;
    const complex_t IT_0115 = 128*IT_0106;
    const complex_t IT_0116 = (-64)*IT_0100;
    const complex_t IT_0117 = IT_0114 + IT_0115 + IT_0116;
    const complex_t IT_0118 = IT_0111 + IT_0117;
    const complex_t IT_0119 = (-2)*IT_0100;
    const complex_t IT_0120 = (IT_0032*std::conj(IT_0032) + (IT_0070 + IT_0098
      )*IT_0099)*IT_0101 + std::conj(IT_0102)*(IT_0103*(IT_0109 + IT_0112 +
       IT_0113) + IT_0102*(IT_0113 + IT_0117 + IT_0118)) + std::conj(IT_0103)*
      (IT_0102*(IT_0109 + IT_0112 + IT_0113) + IT_0103*(IT_0113 + IT_0117 +
       IT_0118)) + (std::conj(IT_0032)*(IT_0070 + IT_0098) + IT_0032*IT_0099)
      *IT_0119;
    return IT_0120;
}
} // End of namespace brparity
