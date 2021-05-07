#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_xx_uf3_to_df1c_df1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_uf3_to_df1c_df1c(
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
    const complex_t IT_0000 = s_24*s_35;
    const complex_t IT_0001 = lpp_201*U_dtR_10;
    const complex_t IT_0002 = lpp_202*U_dtR_20;
    const complex_t IT_0003 = IT_0001 + IT_0002;
    const complex_t IT_0004 = lpp_210*U_dtR_10;
    const complex_t IT_0005 = lpp_220*U_dtR_20;
    const complex_t IT_0006 = -IT_0004 + -IT_0005;
    const complex_t IT_0007 = IT_0003 + IT_0006;
    const complex_t IT_0008 = gw*gwdR*std::conj(U_dtR_00);
    const complex_t IT_0009 = IT_0007*IT_0008;
    const complex_t IT_0010 = std::pow(m_d, 2);
    const complex_t IT_0011 = 2*mX;
    const complex_t IT_0012 = std::pow(IT_0011, 2);
    const complex_t IT_0013 = std::pow(m_dtR1, 2);
    const complex_t IT_0014 = std::pow((-2)*s_24 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0010 + IT_0012 + -IT_0013, -1);
    const complex_t IT_0015 = (complex_t{0, 1})*IT_0009*IT_0014;
    const complex_t IT_0016 = lpp_201*U_dtR_12;
    const complex_t IT_0017 = lpp_202*U_dtR_22;
    const complex_t IT_0018 = IT_0016 + IT_0017;
    const complex_t IT_0019 = lpp_210*U_dtR_12;
    const complex_t IT_0020 = lpp_220*U_dtR_22;
    const complex_t IT_0021 = -IT_0019 + -IT_0020;
    const complex_t IT_0022 = IT_0018 + IT_0021;
    const complex_t IT_0023 = gw*gwdR*std::conj(U_dtR_02);
    const complex_t IT_0024 = IT_0022*IT_0023;
    const complex_t IT_0025 = std::pow(m_dtR3, 2);
    const complex_t IT_0026 = std::pow((-2)*s_24 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0010 + IT_0012 + -IT_0025, -1);
    const complex_t IT_0027 = (complex_t{0, 1})*IT_0024*IT_0026;
    const complex_t IT_0028 = IT_0015 + IT_0027;
    const complex_t IT_0029 = lpp_201*U_dtR_11;
    const complex_t IT_0030 = lpp_202*U_dtR_21;
    const complex_t IT_0031 = IT_0029 + IT_0030;
    const complex_t IT_0032 = lpp_210*U_dtR_11;
    const complex_t IT_0033 = lpp_220*U_dtR_21;
    const complex_t IT_0034 = -IT_0032 + -IT_0033;
    const complex_t IT_0035 = IT_0031 + IT_0034;
    const complex_t IT_0036 = gw*gwdR*std::conj(U_dtR_01);
    const complex_t IT_0037 = IT_0035*IT_0036;
    const complex_t IT_0038 = std::pow(m_dtR2, 2);
    const complex_t IT_0039 = std::pow((-2)*s_24 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0010 + IT_0012 + -IT_0038, -1);
    const complex_t IT_0040 = (complex_t{0, 1})*IT_0037*IT_0039;
    const complex_t IT_0041 = -IT_0040;
    const complex_t IT_0042 = std::conj(IT_0028) + -std::conj(IT_0041);
    const complex_t IT_0043 = std::pow((-2)*s_25 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0010 + IT_0012 + -IT_0013, -1);
    const complex_t IT_0044 = (complex_t{0, 1})*IT_0009*IT_0043;
    const complex_t IT_0045 = std::pow((-2)*s_25 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0010 + IT_0012 + -IT_0038, -1);
    const complex_t IT_0046 = (complex_t{0, 1})*IT_0037*IT_0045;
    const complex_t IT_0047 = std::pow((-2)*s_25 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0010 + IT_0012 + -IT_0025, -1);
    const complex_t IT_0048 = (complex_t{0, 1})*IT_0024*IT_0047;
    const complex_t IT_0049 = -IT_0044 + -IT_0046 + -IT_0048;
    const complex_t IT_0050 = s_25*s_34;
    const complex_t IT_0051 = s_23*s_45;
    const complex_t IT_0052 = -IT_0051;
    const complex_t IT_0053 = IT_0000 + IT_0050 + IT_0052;
    const complex_t IT_0054 = 0.5*IT_0053;
    const complex_t IT_0055 = (-0.5)*IT_0053;
    return IT_0000*(IT_0028 + -IT_0041)*IT_0042 + std::conj(IT_0049)*(IT_0041
      *IT_0054 + IT_0028*IT_0055) + IT_0049*(std::conj(IT_0049)*IT_0050 +
       std::conj(IT_0041)*IT_0054 + std::conj(IT_0028)*IT_0055);
}
} // End of namespace brparity
