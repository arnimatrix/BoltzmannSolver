#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XXc_to_uf1_df2_df2.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_to_uf1_df2_df2(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_s = param.m_s;
    const real_t m_u = param.m_u;
    const real_t GbRt = param.GbRt;
    const real_t GdRt = param.GdRt;
    const real_t GsRt = param.GsRt;
    const real_t gwdR = param.gwdR;
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
    const real_t reg_prop = param.reg_prop;
    const complex_t lpp_001 = param.lpp_001;
    const complex_t lpp_010 = param.lpp_010;
    const complex_t lpp_012 = param.lpp_012;
    const complex_t lpp_021 = param.lpp_021;
    const complex_t U_dtR_00 = param.U_dtR_00;
    const complex_t U_dtR_01 = param.U_dtR_01;
    const complex_t U_dtR_02 = param.U_dtR_02;
    const complex_t U_dtR_10 = param.U_dtR_10;
    const complex_t U_dtR_11 = param.U_dtR_11;
    const complex_t U_dtR_12 = param.U_dtR_12;
    const complex_t U_dtR_20 = param.U_dtR_20;
    const complex_t U_dtR_21 = param.U_dtR_21;
    const complex_t U_dtR_22 = param.U_dtR_22;
    const complex_t IT_0000 = std::conj(lpp_001)*std::conj(U_dtR_00);
    const complex_t IT_0001 = std::conj(lpp_010)*std::conj(U_dtR_00);
    const complex_t IT_0002 = std::conj(lpp_012)*std::conj(U_dtR_20);
    const complex_t IT_0003 = std::conj(lpp_021)*std::conj(U_dtR_20);
    const complex_t IT_0004 = (complex_t{0, 1})*(IT_0000 + -IT_0001 + -IT_0002
       + IT_0003);
    const complex_t IT_0005 = (complex_t{0, -1})*gw*gwdR*U_dtR_10;
    const complex_t IT_0006 = IT_0004*IT_0005;
    const complex_t IT_0007 = std::pow(m_s, 2);
    const complex_t IT_0008 = std::pow(m_u, 2);
    const complex_t IT_0009 = std::pow(m_dtR1, 2);
    const complex_t IT_0010 = std::pow(2*s_23 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0007 + IT_0008 + -IT_0009 + reg_prop, -1);
    const complex_t IT_0011 = (complex_t{0, 1})*IT_0006*IT_0010;
    const complex_t IT_0012 = (-0.25)*IT_0011;
    const complex_t IT_0013 = std::conj(lpp_001)*std::conj(U_dtR_01);
    const complex_t IT_0014 = std::conj(lpp_010)*std::conj(U_dtR_01);
    const complex_t IT_0015 = std::conj(lpp_012)*std::conj(U_dtR_21);
    const complex_t IT_0016 = std::conj(lpp_021)*std::conj(U_dtR_21);
    const complex_t IT_0017 = (complex_t{0, 1})*(IT_0013 + -IT_0014 + -IT_0015
       + IT_0016);
    const complex_t IT_0018 = (complex_t{0, -1})*gw*gwdR*U_dtR_11;
    const complex_t IT_0019 = IT_0017*IT_0018;
    const complex_t IT_0020 = std::pow(m_dtR2, 2);
    const complex_t IT_0021 = std::pow(2*s_23 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0007 + IT_0008 + -IT_0020 + reg_prop, -1);
    const complex_t IT_0022 = (complex_t{0, 1})*IT_0019*IT_0021;
    const complex_t IT_0023 = std::conj(lpp_001)*std::conj(U_dtR_02);
    const complex_t IT_0024 = std::conj(lpp_010)*std::conj(U_dtR_02);
    const complex_t IT_0025 = std::conj(lpp_012)*std::conj(U_dtR_22);
    const complex_t IT_0026 = std::conj(lpp_021)*std::conj(U_dtR_22);
    const complex_t IT_0027 = (complex_t{0, 1})*(IT_0023 + -IT_0024 + -IT_0025
       + IT_0026);
    const complex_t IT_0028 = (complex_t{0, -1})*gw*gwdR*U_dtR_12;
    const complex_t IT_0029 = IT_0027*IT_0028;
    const complex_t IT_0030 = std::pow(m_dtR3, 2);
    const complex_t IT_0031 = std::pow(2*s_23 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0007 + IT_0008 + -IT_0030 + reg_prop, -1);
    const complex_t IT_0032 = (complex_t{0, 1})*IT_0029*IT_0031;
    const complex_t IT_0033 = (-0.25)*IT_0022 + (-0.25)*IT_0032;
    const complex_t IT_0034 = (-0.5)*IT_0011 + (-0.5)*IT_0022 + (-0.5)*IT_0032;
    const complex_t IT_0035 = s_13*s_24;
    const complex_t IT_0036 = 24*IT_0035;
    const complex_t IT_0037 = s_14*s_23;
    const complex_t IT_0038 = (-24)*IT_0037;
    const complex_t IT_0039 = IT_0036 + IT_0038;
    const complex_t IT_0040 = IT_0034*IT_0039;
    const complex_t IT_0041 = std::pow(m_XX, 2);
    const complex_t IT_0042 = std::pow((-2)*s_13 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0007 + -IT_0009 + IT_0041 + reg_prop, -1);
    const complex_t IT_0043 = (complex_t{0, 1})*IT_0006*IT_0042;
    const complex_t IT_0044 = std::pow((-2)*s_13 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0007 + -IT_0020 + IT_0041 + reg_prop, -1);
    const complex_t IT_0045 = (complex_t{0, 1})*IT_0019*IT_0044;
    const complex_t IT_0046 = std::pow((-2)*s_13 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + -IT_0030 + IT_0041 + reg_prop, -1);
    const complex_t IT_0047 = (complex_t{0, 1})*IT_0029*IT_0046;
    const complex_t IT_0048 = (-0.5)*IT_0043 + (-0.5)*IT_0045 + (-0.5)*IT_0047;
    const complex_t IT_0049 = (-24)*IT_0035;
    const complex_t IT_0050 = 24*IT_0037;
    const complex_t IT_0051 = IT_0049 + IT_0050;
    const complex_t IT_0052 = IT_0048*IT_0051;
    const complex_t IT_0053 = IT_0040 + IT_0052;
    const complex_t IT_0054 = s_12*s_34;
    const complex_t IT_0055 = (-0.25)*IT_0011 + (-0.25)*IT_0022 + (-0.25)
      *IT_0032;
    const complex_t IT_0056 = (-6)*IT_0054;
    const complex_t IT_0057 = 0.25*IT_0043 + 0.25*IT_0045 + 0.25*IT_0047;
    const complex_t IT_0058 = m_u*m_XX*IT_0007;
    const complex_t IT_0059 = (-576)*IT_0058;
    const complex_t IT_0060 = (-384)*IT_0035;
    const complex_t IT_0061 = (-384)*IT_0037;
    const complex_t IT_0062 = 192*IT_0054;
    const complex_t IT_0063 = IT_0059 + IT_0060 + IT_0061 + IT_0062;
    const complex_t IT_0064 = m_u*m_XX*s_34;
    const complex_t IT_0065 = (-288)*IT_0064;
    const complex_t IT_0066 = s_12*IT_0007;
    const complex_t IT_0067 = 288*IT_0066;
    const complex_t IT_0068 = (-192)*IT_0037;
    const complex_t IT_0069 = IT_0059 + IT_0062 + IT_0065 + IT_0067 + IT_0068;
    const complex_t IT_0070 = 384*IT_0035;
    const complex_t IT_0071 = 384*IT_0037;
    const complex_t IT_0072 = (-192)*IT_0054;
    const complex_t IT_0073 = IT_0070 + IT_0071 + IT_0072;
    const complex_t IT_0074 = IT_0059 + IT_0073;
    const complex_t IT_0075 = 576*IT_0058;
    const complex_t IT_0076 = IT_0073 + IT_0075;
    const complex_t IT_0077 = IT_0060 + IT_0061 + IT_0062 + IT_0075;
    const complex_t IT_0078 = std::conj(IT_0034)*IT_0039;
    const complex_t IT_0079 = std::conj(IT_0048)*IT_0051;
    const complex_t IT_0080 = 192*IT_0037;
    const complex_t IT_0081 = 288*IT_0064;
    const complex_t IT_0082 = (-288)*IT_0066;
    const complex_t IT_0083 = IT_0039*std::conj(IT_0048);
    const complex_t IT_0084 = std::conj(IT_0034)*IT_0051;
    const complex_t IT_0085 = (std::conj(IT_0012) + std::conj(IT_0033))
      *IT_0053 + 6*IT_0048*(std::conj(IT_0048)*IT_0054 + 0.166666666666667
      *IT_0039*std::conj(IT_0055) + 0.166666666666667*std::conj(IT_0034)*IT_0056
      ) + 6*IT_0034*(std::conj(IT_0034)*IT_0054 + 0.166666666666667*IT_0051
      *std::conj(IT_0055) + 0.166666666666667*std::conj(IT_0048)*IT_0056) +
       std::conj(IT_0057)*(IT_0039*IT_0048 + IT_0034*(IT_0039 + IT_0051) +
       IT_0052 + IT_0033*IT_0063 + IT_0012*IT_0069 + IT_0057*(IT_0074 + IT_0076)
       + IT_0055*IT_0077) + IT_0033*(std::conj(IT_0033)*IT_0076 + IT_0078 +
       IT_0079 + std::conj(IT_0012)*(IT_0072 + IT_0075 + IT_0080 + IT_0081 +
       IT_0082)) + IT_0012*(std::conj(IT_0012)*IT_0076 + IT_0078 + IT_0079 +
       std::conj(IT_0033)*(IT_0072 + IT_0075 + IT_0080 + IT_0081 + IT_0082)) +
       IT_0055*(std::conj(IT_0055)*IT_0074 + IT_0083 + IT_0084) + IT_0057*
      (std::conj(IT_0034)*IT_0039 + std::conj(IT_0033)*IT_0063 + std::conj
      (IT_0012)*IT_0069 + std::conj(IT_0055)*IT_0077 + IT_0079 + IT_0083 +
       IT_0084);
    return IT_0085;
}
} // End of namespace brparity