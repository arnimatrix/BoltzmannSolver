#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XXc_df3c_to_uf1_df3.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_df3c_to_uf1_df3(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_b = param.m_b;
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
    const complex_t lpp_002 = param.lpp_002;
    const complex_t lpp_012 = param.lpp_012;
    const complex_t lpp_020 = param.lpp_020;
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
    const complex_t IT_0000 = s_12*s_34;
    const complex_t IT_0001 = std::conj(lpp_002)*std::conj(U_dtR_00);
    const complex_t IT_0002 = std::conj(lpp_012)*std::conj(U_dtR_10);
    const complex_t IT_0003 = std::conj(lpp_020)*std::conj(U_dtR_00);
    const complex_t IT_0004 = std::conj(lpp_021)*std::conj(U_dtR_10);
    const complex_t IT_0005 = (complex_t{0, 1})*(IT_0001 + IT_0002 + -IT_0003 
      + -IT_0004);
    const complex_t IT_0006 = (complex_t{0, -1})*gw*gwdR*U_dtR_20;
    const complex_t IT_0007 = IT_0005*IT_0006;
    const complex_t IT_0008 = std::pow(m_b, 2);
    const complex_t IT_0009 = std::pow(m_XX, 2);
    const complex_t IT_0010 = std::pow(m_dtR1, 2);
    const complex_t IT_0011 = std::pow(2*s_12 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0008 + IT_0009 + -IT_0010 + reg_prop, -1);
    const complex_t IT_0012 = (complex_t{0, 1})*IT_0007*IT_0011;
    const complex_t IT_0013 = std::conj(lpp_002)*std::conj(U_dtR_01);
    const complex_t IT_0014 = std::conj(lpp_012)*std::conj(U_dtR_11);
    const complex_t IT_0015 = std::conj(lpp_020)*std::conj(U_dtR_01);
    const complex_t IT_0016 = std::conj(lpp_021)*std::conj(U_dtR_11);
    const complex_t IT_0017 = (complex_t{0, 1})*(IT_0013 + IT_0014 + -IT_0015 
      + -IT_0016);
    const complex_t IT_0018 = (complex_t{0, -1})*gw*gwdR*U_dtR_21;
    const complex_t IT_0019 = IT_0017*IT_0018;
    const complex_t IT_0020 = std::pow(m_dtR2, 2);
    const complex_t IT_0021 = std::pow(2*s_12 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0008 + IT_0009 + -IT_0020 + reg_prop, -1);
    const complex_t IT_0022 = (complex_t{0, 1})*IT_0019*IT_0021;
    const complex_t IT_0023 = std::conj(lpp_002)*std::conj(U_dtR_02);
    const complex_t IT_0024 = std::conj(lpp_012)*std::conj(U_dtR_12);
    const complex_t IT_0025 = std::conj(lpp_020)*std::conj(U_dtR_02);
    const complex_t IT_0026 = std::conj(lpp_021)*std::conj(U_dtR_12);
    const complex_t IT_0027 = (complex_t{0, 1})*(IT_0023 + IT_0024 + -IT_0025 
      + -IT_0026);
    const complex_t IT_0028 = (complex_t{0, -1})*gw*gwdR*U_dtR_22;
    const complex_t IT_0029 = IT_0027*IT_0028;
    const complex_t IT_0030 = std::pow(m_dtR3, 2);
    const complex_t IT_0031 = std::pow(2*s_12 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0008 + IT_0009 + -IT_0030 + reg_prop, -1);
    const complex_t IT_0032 = (complex_t{0, 1})*IT_0029*IT_0031;
    const complex_t IT_0033 = -IT_0012 + -IT_0022 + -IT_0032;
    const complex_t IT_0034 = std::pow(m_u, 2);
    const complex_t IT_0035 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0008 + -IT_0010 + IT_0034 + reg_prop, -1);
    const complex_t IT_0036 = (complex_t{0, 1})*IT_0007*IT_0035;
    const complex_t IT_0037 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0008 + -IT_0020 + IT_0034 + reg_prop, -1);
    const complex_t IT_0038 = (complex_t{0, 1})*IT_0019*IT_0037;
    const complex_t IT_0039 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0008 + -IT_0030 + IT_0034 + reg_prop, -1);
    const complex_t IT_0040 = (complex_t{0, 1})*IT_0029*IT_0039;
    const complex_t IT_0041 = 0.5*IT_0036 + 0.5*IT_0038 + 0.5*IT_0040;
    const complex_t IT_0042 = 0.25*IT_0036 + 0.25*IT_0038 + 0.25*IT_0040;
    const complex_t IT_0043 = s_14*s_23;
    const complex_t IT_0044 = s_13*s_24;
    const complex_t IT_0045 = 128*IT_0044;
    const complex_t IT_0046 = (-64)*IT_0000;
    const complex_t IT_0047 = 128*IT_0043 + IT_0045 + IT_0046;
    const complex_t IT_0048 = (-2)*IT_0000;
    return 2*IT_0000*(IT_0033*std::conj(IT_0033) + IT_0041*std::conj(IT_0041))
       + 2*IT_0042*std::conj(IT_0042)*IT_0047 + (std::conj(IT_0033)*IT_0041 +
       IT_0033*std::conj(IT_0041))*IT_0048;
}
} // End of namespace brparity
