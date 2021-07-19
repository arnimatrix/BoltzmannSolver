#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XXc_df2_to_uf2c_df2c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_df2_to_uf2c_df2c(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_c = param.m_c;
    const real_t m_s = param.m_s;
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
    const complex_t lpp_110 = param.lpp_110;
    const complex_t lpp_112 = param.lpp_112;
    const complex_t lpp_121 = param.lpp_121;
    const complex_t U_dtR_00 = param.U_dtR_00;
    const complex_t U_dtR_01 = param.U_dtR_01;
    const complex_t U_dtR_02 = param.U_dtR_02;
    const complex_t U_dtR_10 = param.U_dtR_10;
    const complex_t U_dtR_11 = param.U_dtR_11;
    const complex_t U_dtR_12 = param.U_dtR_12;
    const complex_t U_dtR_20 = param.U_dtR_20;
    const complex_t U_dtR_21 = param.U_dtR_21;
    const complex_t U_dtR_22 = param.U_dtR_22;
    const complex_t IT_0000 = lpp_101*U_dtR_00;
    const complex_t IT_0001 = lpp_110*U_dtR_00;
    const complex_t IT_0002 = lpp_112*U_dtR_20;
    const complex_t IT_0003 = lpp_121*U_dtR_20;
    const complex_t IT_0004 = (complex_t{0, 1})*(IT_0000 + -IT_0001 + -IT_0002
       + IT_0003);
    const complex_t IT_0005 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_10);
    const complex_t IT_0006 = IT_0004*IT_0005;
    const complex_t IT_0007 = std::pow(m_s, 2);
    const complex_t IT_0008 = std::pow(m_XX, 2);
    const complex_t IT_0009 = std::pow(m_dtR1, 2);
    const complex_t IT_0010 = std::pow(2*s_12 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0007 + IT_0008 + -IT_0009 + reg_prop, -1);
    const complex_t IT_0011 = (complex_t{0, 1})*IT_0006*IT_0010;
    const complex_t IT_0012 = lpp_101*U_dtR_01;
    const complex_t IT_0013 = lpp_110*U_dtR_01;
    const complex_t IT_0014 = lpp_112*U_dtR_21;
    const complex_t IT_0015 = lpp_121*U_dtR_21;
    const complex_t IT_0016 = (complex_t{0, 1})*(IT_0012 + -IT_0013 + -IT_0014
       + IT_0015);
    const complex_t IT_0017 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_11);
    const complex_t IT_0018 = IT_0016*IT_0017;
    const complex_t IT_0019 = std::pow(m_dtR2, 2);
    const complex_t IT_0020 = std::pow(2*s_12 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0007 + IT_0008 + -IT_0019 + reg_prop, -1);
    const complex_t IT_0021 = (complex_t{0, 1})*IT_0018*IT_0020;
    const complex_t IT_0022 = lpp_101*U_dtR_02;
    const complex_t IT_0023 = lpp_110*U_dtR_02;
    const complex_t IT_0024 = lpp_112*U_dtR_22;
    const complex_t IT_0025 = lpp_121*U_dtR_22;
    const complex_t IT_0026 = (complex_t{0, 1})*(IT_0022 + -IT_0023 + -IT_0024
       + IT_0025);
    const complex_t IT_0027 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_12);
    const complex_t IT_0028 = IT_0026*IT_0027;
    const complex_t IT_0029 = std::pow(m_dtR3, 2);
    const complex_t IT_0030 = std::pow(2*s_12 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0007 + IT_0008 + -IT_0029 + reg_prop, -1);
    const complex_t IT_0031 = (complex_t{0, 1})*IT_0028*IT_0030;
    const complex_t IT_0032 = -IT_0011 + -IT_0021 + -IT_0031;
    const complex_t IT_0033 = s_12*s_34;
    const complex_t IT_0034 = std::pow(m_c, 2);
    const complex_t IT_0035 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0007 + -IT_0009 + IT_0034 + reg_prop, -1);
    const complex_t IT_0036 = (complex_t{0, 1})*IT_0006*IT_0035;
    const complex_t IT_0037 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0007 + -IT_0019 + IT_0034 + reg_prop, -1);
    const complex_t IT_0038 = (complex_t{0, 1})*IT_0018*IT_0037;
    const complex_t IT_0039 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + -IT_0029 + IT_0034 + reg_prop, -1);
    const complex_t IT_0040 = (complex_t{0, 1})*IT_0028*IT_0039;
    const complex_t IT_0041 = (-0.25)*IT_0036 + (-0.25)*IT_0038 + (-0.25)
      *IT_0040;
    const complex_t IT_0042 = 0.25*IT_0036 + 0.25*IT_0038 + 0.25*IT_0040;
    const complex_t IT_0043 = s_14*s_23;
    const complex_t IT_0044 = s_13*s_24;
    const complex_t IT_0045 = -IT_0044;
    const complex_t IT_0046 = IT_0043 + IT_0045;
    const complex_t IT_0047 = 8*IT_0046;
    const complex_t IT_0048 = 0.5*IT_0036 + 0.5*IT_0038 + 0.5*IT_0040;
    const complex_t IT_0049 = (-2)*IT_0033;
    const complex_t IT_0050 = (-8)*IT_0046;
    const complex_t IT_0051 = m_c*m_XX*IT_0007;
    const complex_t IT_0052 = 128*IT_0044;
    const complex_t IT_0053 = (-64)*IT_0033;
    const complex_t IT_0054 = std::conj(IT_0048)*IT_0050;
    const complex_t IT_0055 = std::conj(IT_0032)*IT_0047;
    const complex_t IT_0056 = 2*IT_0032*(std::conj(IT_0032)*IT_0033 + 0.5*
      (std::conj(IT_0041) + std::conj(IT_0042))*IT_0047 + 0.5*std::conj(IT_0048)
      *IT_0049) + 2*IT_0048*(IT_0033*std::conj(IT_0048) + 0.5*std::conj(IT_0032)
      *IT_0049 + 0.5*(std::conj(IT_0041) + std::conj(IT_0042))*IT_0050) +
       IT_0041*(std::conj(IT_0041)*(128*IT_0043 + 192*IT_0051 + IT_0052 +
       IT_0053) + IT_0054 + IT_0055) + IT_0042*(std::conj(IT_0042)*(128*IT_0043 
      + (-192)*IT_0051 + IT_0052 + IT_0053) + IT_0054 + IT_0055);
    return IT_0056;
}
} // End of namespace brparity
