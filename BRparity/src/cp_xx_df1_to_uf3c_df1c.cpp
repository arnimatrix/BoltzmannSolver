#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_xx_df1_to_uf3c_df1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_df1_to_uf3c_df1c(
        param_t const &param
        )
{
    clearcache();
    auto const &gw = param.gw;
    auto const &mX = param.mX;
    auto const &m_d = param.m_d;
    auto const &GbRt = param.GbRt;
    auto const &GdRt = param.GdRt;
    auto const &GsRt = param.GsRt;
    auto const &gwdR = param.gwdR;
    auto const &s_23 = param.s_23;
    auto const &s_24 = param.s_24;
    auto const &s_25 = param.s_25;
    auto const &s_34 = param.s_34;
    auto const &s_35 = param.s_35;
    auto const &s_45 = param.s_45;
    auto const &m_dtR1 = param.m_dtR1;
    auto const &m_dtR2 = param.m_dtR2;
    auto const &m_dtR3 = param.m_dtR3;
    auto const &lpp_201 = param.lpp_201;
    auto const &lpp_202 = param.lpp_202;
    auto const &lpp_210 = param.lpp_210;
    auto const &lpp_220 = param.lpp_220;
    auto const &U_dtR_00 = param.U_dtR_00;
    auto const &U_dtR_01 = param.U_dtR_01;
    auto const &U_dtR_02 = param.U_dtR_02;
    auto const &U_dtR_10 = param.U_dtR_10;
    auto const &U_dtR_11 = param.U_dtR_11;
    auto const &U_dtR_12 = param.U_dtR_12;
    auto const &U_dtR_20 = param.U_dtR_20;
    auto const &U_dtR_21 = param.U_dtR_21;
    auto const &U_dtR_22 = param.U_dtR_22;
    const complex_t IT_0000 = lpp_201*U_dtR_11;
    const complex_t IT_0001 = lpp_202*U_dtR_21;
    const complex_t IT_0002 = IT_0000 + IT_0001;
    const complex_t IT_0003 = lpp_210*U_dtR_11;
    const complex_t IT_0004 = lpp_220*U_dtR_21;
    const complex_t IT_0005 = -IT_0003 + -IT_0004;
    const complex_t IT_0006 = IT_0002 + IT_0005;
    const complex_t IT_0007 = gw*gwdR*std::conj(U_dtR_01);
    const complex_t IT_0008 = IT_0006*IT_0007;
    const complex_t IT_0009 = std::pow(m_d, 2);
    const complex_t IT_0010 = 2*mX;
    const complex_t IT_0011 = std::pow(IT_0010, 2);
    const complex_t IT_0012 = std::pow(m_dtR2, 2);
    const complex_t IT_0013 = std::pow(2*s_23 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0009 + IT_0011 + -IT_0012, -1);
    const complex_t IT_0014 = (complex_t{0, 1})*IT_0008*IT_0013;
    const complex_t IT_0015 = -IT_0014;
    const complex_t IT_0016 = s_23*s_45;
    const complex_t IT_0017 = lpp_201*U_dtR_10;
    const complex_t IT_0018 = lpp_202*U_dtR_20;
    const complex_t IT_0019 = IT_0017 + IT_0018;
    const complex_t IT_0020 = lpp_210*U_dtR_10;
    const complex_t IT_0021 = lpp_220*U_dtR_20;
    const complex_t IT_0022 = -IT_0020 + -IT_0021;
    const complex_t IT_0023 = IT_0019 + IT_0022;
    const complex_t IT_0024 = gw*gwdR*std::conj(U_dtR_00);
    const complex_t IT_0025 = IT_0023*IT_0024;
    const complex_t IT_0026 = std::pow(m_dtR1, 2);
    const complex_t IT_0027 = std::pow((-2)*s_25 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0009 + IT_0011 + -IT_0026, -1);
    const complex_t IT_0028 = (complex_t{0, 1})*IT_0025*IT_0027;
    const complex_t IT_0029 = std::pow((-2)*s_25 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0009 + IT_0011 + -IT_0012, -1);
    const complex_t IT_0030 = (complex_t{0, 1})*IT_0008*IT_0029;
    const complex_t IT_0031 = lpp_201*U_dtR_12;
    const complex_t IT_0032 = lpp_202*U_dtR_22;
    const complex_t IT_0033 = IT_0031 + IT_0032;
    const complex_t IT_0034 = lpp_210*U_dtR_12;
    const complex_t IT_0035 = lpp_220*U_dtR_22;
    const complex_t IT_0036 = -IT_0034 + -IT_0035;
    const complex_t IT_0037 = IT_0033 + IT_0036;
    const complex_t IT_0038 = gw*gwdR*std::conj(U_dtR_02);
    const complex_t IT_0039 = IT_0037*IT_0038;
    const complex_t IT_0040 = std::pow(m_dtR3, 2);
    const complex_t IT_0041 = std::pow((-2)*s_25 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0009 + IT_0011 + -IT_0040, -1);
    const complex_t IT_0042 = (complex_t{0, 1})*IT_0039*IT_0041;
    const complex_t IT_0043 = -IT_0028 + -IT_0030 + -IT_0042;
    const complex_t IT_0044 = s_25*s_34;
    const complex_t IT_0045 = s_24*s_35;
    const complex_t IT_0046 = -IT_0045;
    const complex_t IT_0047 = IT_0016 + IT_0044 + IT_0046;
    const complex_t IT_0048 = std::pow(2*s_23 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0009 + IT_0011 + -IT_0026, -1);
    const complex_t IT_0049 = (complex_t{0, 1})*IT_0025*IT_0048;
    const complex_t IT_0050 = std::pow(2*s_23 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0009 + IT_0011 + -IT_0040, -1);
    const complex_t IT_0051 = (complex_t{0, 1})*IT_0039*IT_0050;
    const complex_t IT_0052 = IT_0049 + IT_0051;
    const complex_t IT_0053 = (-2)*IT_0016;
    const complex_t IT_0054 = -IT_0047;
    const complex_t IT_0055 = 0.5*std::conj(IT_0015);
    const complex_t IT_0056 = 2*IT_0015*(std::conj(IT_0015)*IT_0016 + 0.5
      *std::conj(IT_0043)*IT_0047 + 0.5*std::conj(IT_0052)*IT_0053) + 2*IT_0043*
      (std::conj(IT_0043)*IT_0044 + 0.5*std::conj(IT_0052)*IT_0054 + IT_0047
      *IT_0055) + 2*IT_0052*(IT_0016*std::conj(IT_0052) + 0.5*std::conj(IT_0043)
      *IT_0054 + IT_0053*IT_0055);
    return IT_0056;
}
} // End of namespace brparity
