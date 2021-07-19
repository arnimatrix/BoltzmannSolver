#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_XX_XX_to_df2_df2c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_XX_XX_to_df2_df2c(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t m_s = param.m_s;
    const real_t GbLt = param.GbLt;
    const real_t GbRt = param.GbRt;
    const real_t GdLt = param.GdLt;
    const real_t GdRt = param.GdRt;
    const real_t GsLt = param.GsLt;
    const real_t GsRt = param.GsRt;
    const real_t gwdL = param.gwdL;
    const real_t gwdR = param.gwdR;
    const real_t m_XX = param.m_XX;
    const real_t s_12 = param.s_12;
    const real_t s_13 = param.s_13;
    const real_t s_14 = param.s_14;
    const real_t s_23 = param.s_23;
    const real_t s_24 = param.s_24;
    const real_t s_34 = param.s_34;
    const real_t m_dtL1 = param.m_dtL1;
    const real_t m_dtL2 = param.m_dtL2;
    const real_t m_dtL3 = param.m_dtL3;
    const real_t m_dtR1 = param.m_dtR1;
    const real_t m_dtR2 = param.m_dtR2;
    const real_t m_dtR3 = param.m_dtR3;
    const real_t reg_prop = param.reg_prop;
    const complex_t U_dtL_10 = param.U_dtL_10;
    const complex_t U_dtL_11 = param.U_dtL_11;
    const complex_t U_dtL_12 = param.U_dtL_12;
    const complex_t U_dtR_10 = param.U_dtR_10;
    const complex_t U_dtR_11 = param.U_dtR_11;
    const complex_t U_dtR_12 = param.U_dtR_12;
    const complex_t IT_0031 = (complex_t{0, -1})*gw*gwdL*U_dtL_10;
    const complex_t IT_0032 = (complex_t{0, -1})*gw*gwdL*std::conj(U_dtL_10);
    const complex_t IT_0033 = IT_0031*IT_0032;
    const complex_t IT_0003 = std::pow(m_s, 2);
    const complex_t IT_0004 = std::pow(m_XX, 2);
    const complex_t IT_0034 = std::pow(m_dtL1, 2);
    const complex_t IT_0052 = std::pow((-2)*s_13 + (complex_t{0, 1})*GdLt
      *m_dtL1 + IT_0003 + IT_0004 + -IT_0034 + reg_prop, -1);
    const complex_t IT_0053 = (complex_t{0, 1})*IT_0033*IT_0052;
    const complex_t IT_0037 = (complex_t{0, -1})*gw*gwdL*U_dtL_11;
    const complex_t IT_0038 = (complex_t{0, -1})*gw*gwdL*std::conj(U_dtL_11);
    const complex_t IT_0039 = IT_0037*IT_0038;
    const complex_t IT_0040 = std::pow(m_dtL2, 2);
    const complex_t IT_0054 = std::pow((-2)*s_13 + (complex_t{0, 1})*GsLt
      *m_dtL2 + IT_0003 + IT_0004 + -IT_0040 + reg_prop, -1);
    const complex_t IT_0055 = (complex_t{0, 1})*IT_0039*IT_0054;
    const complex_t IT_0043 = (complex_t{0, -1})*gw*gwdL*U_dtL_12;
    const complex_t IT_0044 = (complex_t{0, -1})*gw*gwdL*std::conj(U_dtL_12);
    const complex_t IT_0045 = IT_0043*IT_0044;
    const complex_t IT_0046 = std::pow(m_dtL3, 2);
    const complex_t IT_0056 = std::pow((-2)*s_13 + (complex_t{0, 1})*GbLt
      *m_dtL3 + IT_0003 + IT_0004 + -IT_0046 + reg_prop, -1);
    const complex_t IT_0057 = (complex_t{0, 1})*IT_0045*IT_0056;
    const complex_t IT_0058 = 0.5*IT_0053 + 0.5*IT_0055 + 0.5*IT_0057;
    const complex_t IT_0029 = s_34*IT_0004;
    const complex_t IT_0030 = 6*IT_0029;
    const complex_t IT_0035 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdLt
      *m_dtL1 + IT_0003 + IT_0004 + -IT_0034 + reg_prop, -1);
    const complex_t IT_0036 = (complex_t{0, 1})*IT_0033*IT_0035;
    const complex_t IT_0041 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsLt
      *m_dtL2 + IT_0003 + IT_0004 + -IT_0040 + reg_prop, -1);
    const complex_t IT_0042 = (complex_t{0, 1})*IT_0039*IT_0041;
    const complex_t IT_0047 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbLt
      *m_dtL3 + IT_0003 + IT_0004 + -IT_0046 + reg_prop, -1);
    const complex_t IT_0048 = (complex_t{0, 1})*IT_0045*IT_0047;
    const complex_t IT_0049 = (-0.5)*IT_0036 + (-0.5)*IT_0042 + (-0.5)*IT_0048;
    const complex_t IT_0000 = (complex_t{0, -1})*gw*gwdR*U_dtR_10;
    const complex_t IT_0001 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_10);
    const complex_t IT_0002 = IT_0000*IT_0001;
    const complex_t IT_0005 = std::pow(m_dtR1, 2);
    const complex_t IT_0022 = std::pow((-2)*s_23 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0003 + IT_0004 + -IT_0005 + reg_prop, -1);
    const complex_t IT_0023 = (complex_t{0, 1})*IT_0002*IT_0022;
    const complex_t IT_0008 = (complex_t{0, -1})*gw*gwdR*U_dtR_11;
    const complex_t IT_0009 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_11);
    const complex_t IT_0010 = IT_0008*IT_0009;
    const complex_t IT_0011 = std::pow(m_dtR2, 2);
    const complex_t IT_0024 = std::pow((-2)*s_23 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0003 + IT_0004 + -IT_0011 + reg_prop, -1);
    const complex_t IT_0025 = (complex_t{0, 1})*IT_0010*IT_0024;
    const complex_t IT_0014 = (complex_t{0, -1})*gw*gwdR*U_dtR_12;
    const complex_t IT_0015 = (complex_t{0, -1})*gw*gwdR*std::conj(U_dtR_12);
    const complex_t IT_0016 = IT_0014*IT_0015;
    const complex_t IT_0017 = std::pow(m_dtR3, 2);
    const complex_t IT_0026 = std::pow((-2)*s_23 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0003 + IT_0004 + -IT_0017 + reg_prop, -1);
    const complex_t IT_0027 = (complex_t{0, 1})*IT_0016*IT_0026;
    const complex_t IT_0028 = (-0.5)*IT_0023 + (-0.5)*IT_0025 + (-0.5)*IT_0027;
    const complex_t IT_0050 = s_12*IT_0003;
    const complex_t IT_0051 = 6*IT_0050;
    const complex_t IT_0021 = s_13*s_24;
    const complex_t IT_0006 = std::pow((-2)*s_13 + (complex_t{0, 1})*GdRt
      *m_dtR1 + IT_0003 + IT_0004 + -IT_0005 + reg_prop, -1);
    const complex_t IT_0007 = (complex_t{0, 1})*IT_0002*IT_0006;
    const complex_t IT_0012 = std::pow((-2)*s_13 + (complex_t{0, 1})*GsRt
      *m_dtR2 + IT_0003 + IT_0004 + -IT_0011 + reg_prop, -1);
    const complex_t IT_0013 = (complex_t{0, 1})*IT_0010*IT_0012;
    const complex_t IT_0018 = std::pow((-2)*s_13 + (complex_t{0, 1})*GbRt
      *m_dtR3 + IT_0003 + IT_0004 + -IT_0017 + reg_prop, -1);
    const complex_t IT_0019 = (complex_t{0, 1})*IT_0016*IT_0018;
    const complex_t IT_0020 = 0.5*IT_0007 + 0.5*IT_0013 + 0.5*IT_0019;
    const complex_t IT_0059 = IT_0003*IT_0004;
    const complex_t IT_0060 = 12*IT_0059;
    const complex_t IT_0061 = s_14*s_23;
    const complex_t IT_0062 = 12*IT_0061;
    const complex_t IT_0063 = IT_0058*(IT_0030*std::conj(IT_0049) + std::conj
      (IT_0028)*IT_0051 + 12*IT_0021*std::conj(IT_0058) + std::conj(IT_0020)
      *IT_0060) + IT_0020*(12*std::conj(IT_0020)*IT_0021 + std::conj(IT_0028)
      *IT_0030 + std::conj(IT_0049)*IT_0051 + std::conj(IT_0058)*IT_0060) +
       IT_0028*(std::conj(IT_0020)*IT_0030 + IT_0051*std::conj(IT_0058) +
       std::conj(IT_0049)*IT_0060 + std::conj(IT_0028)*IT_0062) + IT_0049*
      (std::conj(IT_0020)*IT_0051 + IT_0030*std::conj(IT_0058) + std::conj
      (IT_0028)*IT_0060 + std::conj(IT_0049)*IT_0062);
    return IT_0063;
}
} // End of namespace brparity
