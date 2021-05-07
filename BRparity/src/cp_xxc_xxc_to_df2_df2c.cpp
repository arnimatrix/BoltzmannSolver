#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_xxc_xxc_to_df2_df2c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_XXc_to_df2_df2c(
        param_t const &param
        )
{
    clearcache();
    auto const &gw = param.gw;
    auto const &mX = param.mX;
    auto const &m_s = param.m_s;
    auto const &GbLt = param.GbLt;
    auto const &GbRt = param.GbRt;
    auto const &GdLt = param.GdLt;
    auto const &GdRt = param.GdRt;
    auto const &GsLt = param.GsLt;
    auto const &GsRt = param.GsRt;
    auto const &gwdL = param.gwdL;
    auto const &gwdR = param.gwdR;
    auto const &s_23 = param.s_23;
    auto const &s_24 = param.s_24;
    auto const &s_25 = param.s_25;
    auto const &s_34 = param.s_34;
    auto const &s_35 = param.s_35;
    auto const &s_45 = param.s_45;
    auto const &m_dtL1 = param.m_dtL1;
    auto const &m_dtL2 = param.m_dtL2;
    auto const &m_dtL3 = param.m_dtL3;
    auto const &m_dtR1 = param.m_dtR1;
    auto const &m_dtR2 = param.m_dtR2;
    auto const &m_dtR3 = param.m_dtR3;
    auto const &U_dtL_10 = param.U_dtL_10;
    auto const &U_dtL_11 = param.U_dtL_11;
    auto const &U_dtL_12 = param.U_dtL_12;
    auto const &U_dtR_10 = param.U_dtR_10;
    auto const &U_dtR_11 = param.U_dtR_11;
    auto const &U_dtR_12 = param.U_dtR_12;
    const complex_t IT_0022 = gw*gwdL*U_dtL_10;
    const complex_t IT_0023 = gw*gwdL*std::conj(U_dtL_10);
    const complex_t IT_0024 = IT_0022*IT_0023;
    const complex_t IT_0003 = std::pow(m_s, 2);
    const complex_t IT_0004 = 2*mX;
    const complex_t IT_0005 = std::pow(IT_0004, 2);
    const complex_t IT_0025 = std::pow(m_dtL1, 2);
    const complex_t IT_0026 = std::pow((-2)*s_25 + (complex_t{0, 1})*GdLt
      *m_dtL1 + IT_0003 + IT_0005 + -IT_0025, -1);
    const complex_t IT_0027 = (complex_t{0, 1})*IT_0024*IT_0026;
    const complex_t IT_0028 = gw*gwdL*U_dtL_11;
    const complex_t IT_0029 = gw*gwdL*std::conj(U_dtL_11);
    const complex_t IT_0030 = IT_0028*IT_0029;
    const complex_t IT_0031 = std::pow(m_dtL2, 2);
    const complex_t IT_0032 = std::pow((-2)*s_25 + (complex_t{0, 1})*GsLt
      *m_dtL2 + IT_0003 + IT_0005 + -IT_0031, -1);
    const complex_t IT_0033 = (complex_t{0, 1})*IT_0030*IT_0032;
    const complex_t IT_0034 = gw*gwdL*U_dtL_12;
    const complex_t IT_0035 = gw*gwdL*std::conj(U_dtL_12);
    const complex_t IT_0036 = IT_0034*IT_0035;
    const complex_t IT_0037 = std::pow(m_dtL3, 2);
    const complex_t IT_0038 = std::pow((-2)*s_25 + (complex_t{0, 1})*GbLt
      *m_dtL3 + IT_0003 + IT_0005 + -IT_0037, -1);
    const complex_t IT_0039 = (complex_t{0, 1})*IT_0036*IT_0038;
    const complex_t IT_0040 = -IT_0027 + -IT_0033 + -IT_0039;
    const complex_t IT_0000 = gw*gwdR*U_dtR_10;
    const complex_t IT_0001 = gw*gwdR*std::conj(U_dtR_10);
    const complex_t IT_0002 = IT_0000*IT_0001;
    const complex_t IT_0006 = std::pow(m_dtR1, 2);
    const complex_t IT_0007 = std::pow((-2)*s_25 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0003 + IT_0005 + -IT_0006, -1);
    const complex_t IT_0008 = (complex_t{0, 1})*IT_0002*IT_0007;
    const complex_t IT_0009 = gw*gwdR*U_dtR_11;
    const complex_t IT_0010 = gw*gwdR*std::conj(U_dtR_11);
    const complex_t IT_0011 = IT_0009*IT_0010;
    const complex_t IT_0012 = std::pow(m_dtR2, 2);
    const complex_t IT_0013 = std::pow((-2)*s_25 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0003 + IT_0005 + -IT_0012, -1);
    const complex_t IT_0014 = (complex_t{0, 1})*IT_0011*IT_0013;
    const complex_t IT_0015 = gw*gwdR*U_dtR_12;
    const complex_t IT_0016 = gw*gwdR*std::conj(U_dtR_12);
    const complex_t IT_0017 = IT_0015*IT_0016;
    const complex_t IT_0018 = std::pow(m_dtR3, 2);
    const complex_t IT_0019 = std::pow((-2)*s_25 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0003 + IT_0005 + -IT_0018, -1);
    const complex_t IT_0020 = (complex_t{0, 1})*IT_0017*IT_0019;
    const complex_t IT_0021 = -IT_0008 + -IT_0014 + -IT_0020;
    const complex_t IT_0041 = IT_0003*IT_0005;
    const complex_t IT_0042 = 3*IT_0041;
    const complex_t IT_0043 = s_25*s_34;
    const complex_t IT_0044 = 3*IT_0043;
    const complex_t IT_0052 = s_45*IT_0005;
    const complex_t IT_0053 = 1.5*IT_0052;
    const complex_t IT_0054 = std::pow((-2)*s_24 + (complex_t{0, 1})*GdLt
      *m_dtL1 + IT_0003 + IT_0005 + -IT_0025, -1);
    const complex_t IT_0055 = (complex_t{0, 1})*IT_0024*IT_0054;
    const complex_t IT_0056 = std::pow((-2)*s_24 + (complex_t{0, 1})*GsLt
      *m_dtL2 + IT_0003 + IT_0005 + -IT_0031, -1);
    const complex_t IT_0057 = (complex_t{0, 1})*IT_0030*IT_0056;
    const complex_t IT_0058 = std::pow((-2)*s_24 + (complex_t{0, 1})*GbLt
      *m_dtL3 + IT_0003 + IT_0005 + -IT_0037, -1);
    const complex_t IT_0059 = (complex_t{0, 1})*IT_0036*IT_0058;
    const complex_t IT_0060 = IT_0055 + IT_0057 + IT_0059;
    const complex_t IT_0045 = std::pow((-2)*s_24 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0003 + IT_0005 + -IT_0006, -1);
    const complex_t IT_0046 = (complex_t{0, 1})*IT_0002*IT_0045;
    const complex_t IT_0047 = std::pow((-2)*s_24 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0003 + IT_0005 + -IT_0012, -1);
    const complex_t IT_0048 = (complex_t{0, 1})*IT_0011*IT_0047;
    const complex_t IT_0049 = std::pow((-2)*s_24 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0003 + IT_0005 + -IT_0018, -1);
    const complex_t IT_0050 = (complex_t{0, 1})*IT_0017*IT_0049;
    const complex_t IT_0051 = IT_0046 + IT_0048 + IT_0050;
    const complex_t IT_0061 = s_23*IT_0003;
    const complex_t IT_0062 = 1.5*IT_0061;
    const complex_t IT_0063 = s_24*s_35;
    const complex_t IT_0064 = 0.333333333333333*std::conj(IT_0021);
    const complex_t IT_0065 = IT_0040*(std::conj(IT_0021)*IT_0042 + std::conj
      (IT_0040)*IT_0044 + IT_0053*std::conj(IT_0060) + std::conj(IT_0051)
      *IT_0062) + IT_0021*(std::conj(IT_0040)*IT_0042 + std::conj(IT_0021)
      *IT_0044 + std::conj(IT_0051)*IT_0053 + std::conj(IT_0060)*IT_0062) + 3
      *IT_0051*(0.333333333333333*IT_0042*std::conj(IT_0060) + 0.333333333333333
      *std::conj(IT_0040)*IT_0062 + std::conj(IT_0051)*IT_0063 + IT_0053*IT_0064
      ) + 3*IT_0060*(0.333333333333333*IT_0042*std::conj(IT_0051) +
       0.333333333333333*std::conj(IT_0040)*IT_0053 + std::conj(IT_0060)*IT_0063
       + IT_0062*IT_0064);
    return IT_0065;
}
} // End of namespace brparity
