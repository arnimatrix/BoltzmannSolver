#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XXc_uf2c_to_df1_df1.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_uf2c_to_df1_df1(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_c = param.m_c;
    const real_t m_d = param.m_d;
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
    const complex_t lpp_101 = param.lpp_101;
    const complex_t lpp_102 = param.lpp_102;
    const complex_t lpp_110 = param.lpp_110;
    const complex_t lpp_120 = param.lpp_120;
    const complex_t U_dtR_00 = param.U_dtR_00;
    const complex_t U_dtR_01 = param.U_dtR_01;
    const complex_t U_dtR_02 = param.U_dtR_02;
    const complex_t U_dtR_10 = param.U_dtR_10;
    const complex_t U_dtR_11 = param.U_dtR_11;
    const complex_t U_dtR_12 = param.U_dtR_12;
    const complex_t U_dtR_20 = param.U_dtR_20;
    const complex_t U_dtR_21 = param.U_dtR_21;
    const complex_t U_dtR_22 = param.U_dtR_22;
    const complex_t IT_0000 = s_12*s_34;
    const complex_t IT_0001 = std::conj(lpp_101)*std::conj(U_dtR_10);
    const complex_t IT_0002 = std::conj(lpp_102)*std::conj(U_dtR_20);
    const complex_t IT_0003 = std::conj(lpp_110)*std::conj(U_dtR_10);
    const complex_t IT_0004 = std::conj(lpp_120)*std::conj(U_dtR_20);
    const complex_t IT_0005 = (complex_t{0, -1})*(IT_0001 + IT_0002 + -IT_0003
       + -IT_0004);
    const complex_t IT_0006 = (complex_t{0, -1})*gw*gwdR*U_dtR_00;
    const complex_t IT_0007 = IT_0005*IT_0006;
    const complex_t IT_0008 = std::pow(m_d, 2);
    const complex_t IT_0009 = std::pow(m_XX, 2);
    const complex_t IT_0010 = std::pow(m_dtR1, 2);
    const complex_t IT_0011 = std::pow((-2)*s_13 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0008 + IT_0009 + -IT_0010 + reg_prop, -1);
    const complex_t IT_0012 = (complex_t{0, 1})*IT_0007*IT_0011;
    const complex_t IT_0013 = std::conj(lpp_101)*std::conj(U_dtR_11);
    const complex_t IT_0014 = std::conj(lpp_102)*std::conj(U_dtR_21);
    const complex_t IT_0015 = std::conj(lpp_110)*std::conj(U_dtR_11);
    const complex_t IT_0016 = std::conj(lpp_120)*std::conj(U_dtR_21);
    const complex_t IT_0017 = (complex_t{0, -1})*(IT_0013 + IT_0014 + -IT_0015
       + -IT_0016);
    const complex_t IT_0018 = (complex_t{0, -1})*gw*gwdR*U_dtR_01;
    const complex_t IT_0019 = IT_0017*IT_0018;
    const complex_t IT_0020 = std::pow(m_dtR2, 2);
    const complex_t IT_0021 = std::pow((-2)*s_13 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0008 + IT_0009 + -IT_0020 + reg_prop, -1);
    const complex_t IT_0022 = (complex_t{0, 1})*IT_0019*IT_0021;
    const complex_t IT_0023 = std::conj(lpp_101)*std::conj(U_dtR_12);
    const complex_t IT_0024 = std::conj(lpp_102)*std::conj(U_dtR_22);
    const complex_t IT_0025 = std::conj(lpp_110)*std::conj(U_dtR_12);
    const complex_t IT_0026 = std::conj(lpp_120)*std::conj(U_dtR_22);
    const complex_t IT_0027 = (complex_t{0, -1})*(IT_0023 + IT_0024 + -IT_0025
       + -IT_0026);
    const complex_t IT_0028 = (complex_t{0, -1})*gw*gwdR*U_dtR_02;
    const complex_t IT_0029 = IT_0027*IT_0028;
    const complex_t IT_0030 = std::pow(m_dtR3, 2);
    const complex_t IT_0031 = std::pow((-2)*s_13 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0008 + IT_0009 + -IT_0030 + reg_prop, -1);
    const complex_t IT_0032 = (complex_t{0, 1})*IT_0029*IT_0031;
    const complex_t IT_0033 = (-0.5)*IT_0012 + (-0.5)*IT_0022 + (-0.5)*IT_0032;
    const complex_t IT_0034 = std::pow(m_c, 2);
    const complex_t IT_0035 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0008 + -IT_0010 + IT_0034 + reg_prop, -1);
    const complex_t IT_0036 = (complex_t{0, 1})*IT_0007*IT_0035;
    const complex_t IT_0037 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0008 + -IT_0020 + IT_0034 + reg_prop, -1);
    const complex_t IT_0038 = (complex_t{0, 1})*IT_0019*IT_0037;
    const complex_t IT_0039 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0008 + -IT_0030 + IT_0034 + reg_prop, -1);
    const complex_t IT_0040 = (complex_t{0, 1})*IT_0029*IT_0039;
    const complex_t IT_0041 = (-0.5)*IT_0036 + (-0.5)*IT_0038 + (-0.5)*IT_0040;
    const complex_t IT_0053 = 0.25*IT_0012 + 0.25*IT_0022 + 0.25*IT_0032;
    const complex_t IT_0043 = s_13*s_24;
    const complex_t IT_0044 = 64*IT_0043;
    const complex_t IT_0045 = s_14*s_23;
    const complex_t IT_0046 = 64*IT_0045;
    const complex_t IT_0047 = (-32)*IT_0000;
    const complex_t IT_0048 = IT_0044 + IT_0046 + IT_0047;
    const complex_t IT_0049 = m_c*m_XX*IT_0008;
    const complex_t IT_0050 = (-96)*IT_0049;
    const complex_t IT_0051 = IT_0048 + IT_0050;
    const complex_t IT_0052 = 96*IT_0049;
    const complex_t IT_0042 = (-0.25)*IT_0036 + (-0.25)*IT_0038 + (-0.25)
      *IT_0040;
    const complex_t IT_0058 = 32*IT_0000;
    const complex_t IT_0054 = s_12*IT_0008;
    const complex_t IT_0055 = m_c*m_XX*s_34;
    const complex_t IT_0056 = 48*IT_0055;
    const complex_t IT_0057 = (-32)*IT_0045;
    const complex_t IT_0059 = IT_0052 + 48*IT_0054 + IT_0056 + IT_0057 +
       IT_0058;
    const complex_t IT_0060 = (-64)*IT_0043;
    const complex_t IT_0061 = -IT_0000;
    const complex_t IT_0062 = IT_0000*(IT_0033*std::conj(IT_0033) + IT_0041
      *std::conj(IT_0041)) + std::conj(IT_0053)*((IT_0048 + IT_0051 + IT_0052)
      *IT_0053 + IT_0042*((-64)*IT_0045 + IT_0050 + IT_0058 + IT_0059 + IT_0060)
      ) + std::conj(IT_0042)*(IT_0042*(IT_0048 + IT_0051 + IT_0052) + IT_0053*((
      -64)*IT_0045 + IT_0050 + IT_0058 + IT_0059 + IT_0060)) + (std::conj
      (IT_0033)*IT_0041 + IT_0033*std::conj(IT_0041))*IT_0061;
    return IT_0062;
}
} // End of namespace brparity
