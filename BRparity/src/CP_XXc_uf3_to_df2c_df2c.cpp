#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XXc_uf3_to_df2c_df2c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_uf3_to_df2c_df2c(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_s = param.m_s;
    const real_t m_t = param.m_t;
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
    const complex_t lpp_201 = param.lpp_201;
    const complex_t lpp_210 = param.lpp_210;
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
    const complex_t IT_0000 = lpp_201*U_dtR_00;
    const complex_t IT_0001 = lpp_210*U_dtR_00;
    const complex_t IT_0002 = lpp_212*U_dtR_20;
    const complex_t IT_0003 = lpp_221*U_dtR_20;
    const complex_t IT_0004 = (complex_t{0, 1})*(IT_0000 + -IT_0001 + -IT_0002
       + IT_0003);
    const complex_t IT_0005 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_10);
    const complex_t IT_0006 = IT_0004*IT_0005;
    const complex_t IT_0007 = std::pow(m_s, 2);
    const complex_t IT_0008 = std::pow(m_XX, 2);
    const complex_t IT_0009 = std::pow(m_dtR1, 2);
    const complex_t IT_0010 = std::pow((-2)*s_13 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0007 + IT_0008 + -IT_0009 + reg_prop, -1);
    const complex_t IT_0011 = (complex_t{0, 1})*IT_0006*IT_0010;
    const complex_t IT_0012 = lpp_201*U_dtR_01;
    const complex_t IT_0013 = lpp_210*U_dtR_01;
    const complex_t IT_0014 = lpp_212*U_dtR_21;
    const complex_t IT_0015 = lpp_221*U_dtR_21;
    const complex_t IT_0016 = (complex_t{0, 1})*(IT_0012 + -IT_0013 + -IT_0014
       + IT_0015);
    const complex_t IT_0017 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_11);
    const complex_t IT_0018 = IT_0016*IT_0017;
    const complex_t IT_0019 = std::pow(m_dtR2, 2);
    const complex_t IT_0020 = std::pow((-2)*s_13 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0007 + IT_0008 + -IT_0019 + reg_prop, -1);
    const complex_t IT_0021 = (complex_t{0, 1})*IT_0018*IT_0020;
    const complex_t IT_0022 = lpp_201*U_dtR_02;
    const complex_t IT_0023 = lpp_210*U_dtR_02;
    const complex_t IT_0024 = lpp_212*U_dtR_22;
    const complex_t IT_0025 = lpp_221*U_dtR_22;
    const complex_t IT_0026 = (complex_t{0, 1})*(IT_0022 + -IT_0023 + -IT_0024
       + IT_0025);
    const complex_t IT_0027 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_12);
    const complex_t IT_0028 = IT_0026*IT_0027;
    const complex_t IT_0029 = std::pow(m_dtR3, 2);
    const complex_t IT_0030 = std::pow((-2)*s_13 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + IT_0008 + -IT_0029 + reg_prop, -1);
    const complex_t IT_0031 = (complex_t{0, 1})*IT_0028*IT_0030;
    const complex_t IT_0032 = 0.25*IT_0011 + 0.25*IT_0021 + 0.25*IT_0031;
    const complex_t IT_0033 = s_13*s_24;
    const complex_t IT_0034 = 64*IT_0033;
    const complex_t IT_0035 = s_14*s_23;
    const complex_t IT_0036 = 64*IT_0035;
    const complex_t IT_0037 = s_12*s_34;
    const complex_t IT_0038 = (-32)*IT_0037;
    const complex_t IT_0039 = IT_0034 + IT_0036 + IT_0038;
    const complex_t IT_0040 = m_t*m_XX*IT_0007;
    const complex_t IT_0041 = 96*IT_0040;
    const complex_t IT_0042 = IT_0039 + IT_0041;
    const complex_t IT_0043 = std::pow(m_t, 2);
    const complex_t IT_0044 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0007 + -IT_0009 + IT_0043 + reg_prop, -1);
    const complex_t IT_0045 = (complex_t{0, 1})*IT_0006*IT_0044;
    const complex_t IT_0046 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0007 + -IT_0019 + IT_0043 + reg_prop, -1);
    const complex_t IT_0047 = (complex_t{0, 1})*IT_0018*IT_0046;
    const complex_t IT_0048 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0007 + -IT_0029 + IT_0043 + reg_prop, -1);
    const complex_t IT_0049 = (complex_t{0, 1})*IT_0028*IT_0048;
    const complex_t IT_0050 = (-0.25)*IT_0045 + (-0.25)*IT_0047 + (-0.25)
      *IT_0049;
    const complex_t IT_0051 = (-64)*IT_0033;
    const complex_t IT_0052 = (-64)*IT_0035;
    const complex_t IT_0053 = 32*IT_0037;
    const complex_t IT_0054 = IT_0051 + IT_0052 + IT_0053;
    const complex_t IT_0055 = (-96)*IT_0040;
    const complex_t IT_0056 = IT_0054 + IT_0055;
    const complex_t IT_0057 = 0.5*IT_0045 + 0.5*IT_0047 + 0.5*IT_0049;
    const complex_t IT_0058 = -IT_0033;
    const complex_t IT_0059 = IT_0035 + IT_0058;
    const complex_t IT_0060 = 4*IT_0059;
    const complex_t IT_0061 = std::conj(IT_0057)*IT_0060;
    const complex_t IT_0062 = (-0.25)*IT_0011 + (-0.25)*IT_0021 + (-0.25)
      *IT_0031;
    const complex_t IT_0063 = IT_0039 + IT_0055;
    const complex_t IT_0064 = 0.25*IT_0045 + 0.25*IT_0047 + 0.25*IT_0049;
    const complex_t IT_0065 = 0.5*IT_0011 + 0.5*IT_0021 + 0.5*IT_0031;
    const complex_t IT_0066 = (-4)*IT_0059;
    const complex_t IT_0067 = std::conj(IT_0057)*IT_0066;
    const complex_t IT_0068 = IT_0032*(std::conj(IT_0032)*IT_0042 + std::conj
      (IT_0050)*IT_0056 + IT_0061) + IT_0062*(IT_0061 + std::conj(IT_0062)
      *IT_0063 + (IT_0041 + IT_0054)*std::conj(IT_0064)) + std::conj(IT_0065)*
      (IT_0060*(IT_0050 + IT_0064) + (IT_0032 + IT_0062)*IT_0066) + IT_0065*
      (IT_0060*(std::conj(IT_0050) + std::conj(IT_0064)) + -IT_0037*(std::conj
      (IT_0057) + -std::conj(IT_0065)) + (std::conj(IT_0032) + std::conj(IT_0062
      ))*IT_0066) + IT_0057*(IT_0060*(std::conj(IT_0032) + std::conj(IT_0062)) +
       IT_0037*(std::conj(IT_0057) + -std::conj(IT_0065)) + (std::conj(IT_0050) 
      + std::conj(IT_0064))*IT_0066) + IT_0050*(IT_0042*std::conj(IT_0050) +
       std::conj(IT_0032)*IT_0056 + IT_0067) + IT_0064*((IT_0041 + IT_0054)
      *std::conj(IT_0062) + IT_0063*std::conj(IT_0064) + IT_0067);
    return IT_0068;
}
} // End of namespace brparity
