#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_xxc_xxc_to_uf1_uf1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XXc_XXc_to_uf1_uf1c(
        param_t const &param
        )
{
    clearcache();
    auto const &gw = param.gw;
    auto const &mX = param.mX;
    auto const &m_u = param.m_u;
    auto const &GcLt = param.GcLt;
    auto const &GcRt = param.GcRt;
    auto const &GtLt = param.GtLt;
    auto const &GtRt = param.GtRt;
    auto const &GuLt = param.GuLt;
    auto const &GuRt = param.GuRt;
    auto const &gwuL = param.gwuL;
    auto const &gwuR = param.gwuR;
    auto const &s_23 = param.s_23;
    auto const &s_24 = param.s_24;
    auto const &s_25 = param.s_25;
    auto const &s_34 = param.s_34;
    auto const &s_35 = param.s_35;
    auto const &s_45 = param.s_45;
    auto const &m_utL1 = param.m_utL1;
    auto const &m_utL2 = param.m_utL2;
    auto const &m_utL3 = param.m_utL3;
    auto const &m_utR1 = param.m_utR1;
    auto const &m_utR2 = param.m_utR2;
    auto const &m_utR3 = param.m_utR3;
    auto const &U_utL_00 = param.U_utL_00;
    auto const &U_utL_01 = param.U_utL_01;
    auto const &U_utL_02 = param.U_utL_02;
    auto const &U_utR_00 = param.U_utR_00;
    auto const &U_utR_01 = param.U_utR_01;
    auto const &U_utR_02 = param.U_utR_02;
    const complex_t IT_0022 = gw*gwuL*U_utL_00;
    const complex_t IT_0023 = gw*gwuL*std::conj(U_utL_00);
    const complex_t IT_0024 = IT_0022*IT_0023;
    const complex_t IT_0003 = std::pow(m_u, 2);
    const complex_t IT_0004 = 2*mX;
    const complex_t IT_0005 = std::pow(IT_0004, 2);
    const complex_t IT_0025 = std::pow(m_utL1, 2);
    const complex_t IT_0026 = std::pow((-2)*s_25 + (complex_t{0, 1})*GuLt
      *m_utL1 + IT_0003 + IT_0005 + -IT_0025, -1);
    const complex_t IT_0027 = (complex_t{0, 1})*IT_0024*IT_0026;
    const complex_t IT_0028 = gw*gwuL*U_utL_01;
    const complex_t IT_0029 = gw*gwuL*std::conj(U_utL_01);
    const complex_t IT_0030 = IT_0028*IT_0029;
    const complex_t IT_0031 = std::pow(m_utL2, 2);
    const complex_t IT_0032 = std::pow((-2)*s_25 + (complex_t{0, 1})*GcLt
      *m_utL2 + IT_0003 + IT_0005 + -IT_0031, -1);
    const complex_t IT_0033 = (complex_t{0, 1})*IT_0030*IT_0032;
    const complex_t IT_0034 = gw*gwuL*U_utL_02;
    const complex_t IT_0035 = gw*gwuL*std::conj(U_utL_02);
    const complex_t IT_0036 = IT_0034*IT_0035;
    const complex_t IT_0037 = std::pow(m_utL3, 2);
    const complex_t IT_0038 = std::pow((-2)*s_25 + (complex_t{0, 1})*GtLt
      *m_utL3 + IT_0003 + IT_0005 + -IT_0037, -1);
    const complex_t IT_0039 = (complex_t{0, 1})*IT_0036*IT_0038;
    const complex_t IT_0040 = -IT_0027 + -IT_0033 + -IT_0039;
    const complex_t IT_0000 = gw*gwuR*U_utR_00;
    const complex_t IT_0001 = gw*gwuR*std::conj(U_utR_00);
    const complex_t IT_0002 = IT_0000*IT_0001;
    const complex_t IT_0006 = std::pow(m_utR1, 2);
    const complex_t IT_0007 = std::pow((-2)*s_25 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0003 + IT_0005 + -IT_0006, -1);
    const complex_t IT_0008 = (complex_t{0, 1})*IT_0002*IT_0007;
    const complex_t IT_0009 = gw*gwuR*U_utR_01;
    const complex_t IT_0010 = gw*gwuR*std::conj(U_utR_01);
    const complex_t IT_0011 = IT_0009*IT_0010;
    const complex_t IT_0012 = std::pow(m_utR2, 2);
    const complex_t IT_0013 = std::pow((-2)*s_25 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0003 + IT_0005 + -IT_0012, -1);
    const complex_t IT_0014 = (complex_t{0, 1})*IT_0011*IT_0013;
    const complex_t IT_0015 = gw*gwuR*U_utR_02;
    const complex_t IT_0016 = gw*gwuR*std::conj(U_utR_02);
    const complex_t IT_0017 = IT_0015*IT_0016;
    const complex_t IT_0018 = std::pow(m_utR3, 2);
    const complex_t IT_0019 = std::pow((-2)*s_25 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0003 + IT_0005 + -IT_0018, -1);
    const complex_t IT_0020 = (complex_t{0, 1})*IT_0017*IT_0019;
    const complex_t IT_0021 = -IT_0008 + -IT_0014 + -IT_0020;
    const complex_t IT_0041 = IT_0003*IT_0005;
    const complex_t IT_0042 = 3*IT_0041;
    const complex_t IT_0043 = s_25*s_34;
    const complex_t IT_0044 = 3*IT_0043;
    const complex_t IT_0052 = s_45*IT_0005;
    const complex_t IT_0053 = 1.5*IT_0052;
    const complex_t IT_0054 = std::pow((-2)*s_24 + (complex_t{0, 1})*GuLt
      *m_utL1 + IT_0003 + IT_0005 + -IT_0025, -1);
    const complex_t IT_0055 = (complex_t{0, 1})*IT_0024*IT_0054;
    const complex_t IT_0056 = std::pow((-2)*s_24 + (complex_t{0, 1})*GcLt
      *m_utL2 + IT_0003 + IT_0005 + -IT_0031, -1);
    const complex_t IT_0057 = (complex_t{0, 1})*IT_0030*IT_0056;
    const complex_t IT_0058 = std::pow((-2)*s_24 + (complex_t{0, 1})*GtLt
      *m_utL3 + IT_0003 + IT_0005 + -IT_0037, -1);
    const complex_t IT_0059 = (complex_t{0, 1})*IT_0036*IT_0058;
    const complex_t IT_0060 = IT_0055 + IT_0057 + IT_0059;
    const complex_t IT_0045 = std::pow((-2)*s_24 + (complex_t{0, 1})*GuRt
      *m_utR1 + IT_0003 + IT_0005 + -IT_0006, -1);
    const complex_t IT_0046 = (complex_t{0, 1})*IT_0002*IT_0045;
    const complex_t IT_0047 = std::pow((-2)*s_24 + (complex_t{0, 1})*GcRt
      *m_utR2 + IT_0003 + IT_0005 + -IT_0012, -1);
    const complex_t IT_0048 = (complex_t{0, 1})*IT_0011*IT_0047;
    const complex_t IT_0049 = std::pow((-2)*s_24 + (complex_t{0, 1})*GtRt
      *m_utR3 + IT_0003 + IT_0005 + -IT_0018, -1);
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
