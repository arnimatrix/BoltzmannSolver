#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_utr2_to_df2c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_utR2_to_df2c_df3c(
        param_t const &param
        )
{
    clearcache();
    auto const &s_34 = param.s_34;
    auto const &lpp_012 = param.lpp_012;
    auto const &lpp_021 = param.lpp_021;
    auto const &lpp_112 = param.lpp_112;
    auto const &lpp_121 = param.lpp_121;
    auto const &lpp_212 = param.lpp_212;
    auto const &lpp_221 = param.lpp_221;
    auto const &U_utR_01 = param.U_utR_01;
    auto const &U_utR_11 = param.U_utR_11;
    auto const &U_utR_21 = param.U_utR_21;
    const complex_t IT_0000 = lpp_021*U_utR_01;
    const complex_t IT_0001 = lpp_121*U_utR_11;
    const complex_t IT_0002 = lpp_221*U_utR_21;
    const complex_t IT_0003 = IT_0000 + IT_0001 + IT_0002;
    const complex_t IT_0004 = lpp_012*U_utR_01;
    const complex_t IT_0005 = lpp_112*U_utR_11;
    const complex_t IT_0006 = lpp_212*U_utR_21;
    const complex_t IT_0007 = IT_0004 + IT_0005 + IT_0006;
    const complex_t IT_0008 = -IT_0007;
    const complex_t IT_0009 = IT_0003 + IT_0008;
    const complex_t IT_0010 = (complex_t{0, 1})*IT_0009;
    return 4*s_34*IT_0010*std::conj(IT_0010);
}
} // End of namespace brparity
