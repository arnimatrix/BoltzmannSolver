#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_xx_df3_to_uf1c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_df3_to_uf1c_df3c(
        param_t const &param
        )
{
    clearcache();
    auto const &gw = param.gw;
    auto const &mX = param.mX;
    auto const &m_b = param.m_b;
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
    auto const &lpp_002 = param.lpp_002;
    auto const &lpp_012 = param.lpp_012;
    auto const &lpp_020 = param.lpp_020;
    auto const &lpp_021 = param.lpp_021;
    auto const &U_dtR_00 = param.U_dtR_00;
    auto const &U_dtR_01 = param.U_dtR_01;
    auto const &U_dtR_02 = param.U_dtR_02;
    auto const &U_dtR_10 = param.U_dtR_10;
    auto const &U_dtR_11 = param.U_dtR_11;
    auto const &U_dtR_12 = param.U_dtR_12;
    auto const &U_dtR_20 = param.U_dtR_20;
    auto const &U_dtR_21 = param.U_dtR_21;
    auto const &U_dtR_22 = param.U_dtR_22;
    const complex_t IT_0000 = lpp_002*U_dtR_00;
    const complex_t IT_0001 = lpp_012*U_dtR_10;
    const complex_t IT_0002 = IT_0000 + IT_0001;
    const complex_t IT_0003 = lpp_020*U_dtR_00;
    const complex_t IT_0004 = lpp_021*U_dtR_10;
    const complex_t IT_0005 = -IT_0003 + -IT_0004;
    const complex_t IT_0006 = IT_0002 + IT_0005;
    const complex_t IT_0007 = gw*gwdR*std::conj(U_dtR_20);
    const complex_t IT_0008 = IT_0006*IT_0007;
    const complex_t IT_0009 = std::pow(m_b, 2);
    const complex_t IT_0010 = 2*mX;
    const complex_t IT_0011 = std::pow(IT_0010, 2);
    const complex_t IT_0012 = std::pow(m_dtR1, 2);
    const complex_t IT_0041 = std::pow(2*s_23 + (complex_t{0, 1})*GdRt*m_dtR1 
      + IT_0009 + IT_0011 + -IT_0012, -1);
    const complex_t IT_0042 = (complex_t{0, 1})*IT_0008*IT_0041;
    const complex_t IT_0015 = lpp_002*U_dtR_01;
    const complex_t IT_0016 = lpp_012*U_dtR_11;
    const complex_t IT_0017 = IT_0015 + IT_0016;
    const complex_t IT_0018 = lpp_020*U_dtR_01;
    const complex_t IT_0019 = lpp_021*U_dtR_11;
    const complex_t IT_0020 = -IT_0018 + -IT_0019;
    const complex_t IT_0021 = IT_0017 + IT_0020;
    const complex_t IT_0022 = gw*gwdR*std::conj(U_dtR_21);
    const complex_t IT_0023 = IT_0021*IT_0022;
    const complex_t IT_0024 = std::pow(m_dtR2, 2);
    const complex_t IT_0043 = std::pow(2*s_23 + (complex_t{0, 1})*GsRt*m_dtR2 
      + IT_0009 + IT_0011 + -IT_0024, -1);
    const complex_t IT_0044 = (complex_t{0, 1})*IT_0023*IT_0043;
    const complex_t IT_0027 = lpp_002*U_dtR_02;
    const complex_t IT_0028 = lpp_012*U_dtR_12;
    const complex_t IT_0029 = IT_0027 + IT_0028;
    const complex_t IT_0030 = lpp_020*U_dtR_02;
    const complex_t IT_0031 = lpp_021*U_dtR_12;
    const complex_t IT_0032 = -IT_0030 + -IT_0031;
    const complex_t IT_0033 = IT_0029 + IT_0032;
    const complex_t IT_0034 = gw*gwdR*std::conj(U_dtR_22);
    const complex_t IT_0035 = IT_0033*IT_0034;
    const complex_t IT_0036 = std::pow(m_dtR3, 2);
    const complex_t IT_0045 = std::pow(2*s_23 + (complex_t{0, 1})*GbRt*m_dtR3 
      + IT_0009 + IT_0011 + -IT_0036, -1);
    const complex_t IT_0046 = (complex_t{0, 1})*IT_0035*IT_0045;
    const complex_t IT_0047 = -IT_0042 + -IT_0044 + -IT_0046;
    const complex_t IT_0048 = s_23*s_45;
    const complex_t IT_0013 = std::pow((-2)*s_25 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0009 + IT_0011 + -IT_0012, -1);
    const complex_t IT_0014 = (complex_t{0, 1})*IT_0008*IT_0013;
    const complex_t IT_0025 = std::pow((-2)*s_25 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0009 + IT_0011 + -IT_0024, -1);
    const complex_t IT_0026 = (complex_t{0, 1})*IT_0023*IT_0025;
    const complex_t IT_0037 = std::pow((-2)*s_25 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0009 + IT_0011 + -IT_0036, -1);
    const complex_t IT_0038 = (complex_t{0, 1})*IT_0035*IT_0037;
    const complex_t IT_0039 = IT_0014 + IT_0026 + IT_0038;
    const complex_t IT_0040 = s_25*s_34;
    const complex_t IT_0049 = s_24*s_35;
    const complex_t IT_0050 = IT_0040 + IT_0048 + -IT_0049;
    const complex_t IT_0051 = -IT_0050;
    return 2*IT_0047*(std::conj(IT_0047)*IT_0048 + 0.5*std::conj(IT_0039)
      *IT_0051) + 2*IT_0039*(std::conj(IT_0039)*IT_0040 + 0.5*std::conj(IT_0047)
      *IT_0051);
}
} // End of namespace brparity
