#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_utr3_to_df3c_df1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_utR3_to_df3c_df1c(
        param_t const &param
        )
{
    clearcache();
    auto const &s_34 = param.s_34;
    auto const &lpp_002 = param.lpp_002;
    auto const &lpp_020 = param.lpp_020;
    auto const &lpp_102 = param.lpp_102;
    auto const &lpp_120 = param.lpp_120;
    auto const &lpp_202 = param.lpp_202;
    auto const &lpp_220 = param.lpp_220;
    auto const &U_utR_02 = param.U_utR_02;
    auto const &U_utR_12 = param.U_utR_12;
    auto const &U_utR_22 = param.U_utR_22;
    const complex_t IT_0000 = lpp_020*U_utR_02;
    const complex_t IT_0001 = lpp_120*U_utR_12;
    const complex_t IT_0002 = lpp_220*U_utR_22;
    const complex_t IT_0003 = IT_0000 + IT_0001 + IT_0002;
    const complex_t IT_0004 = lpp_002*U_utR_02;
    const complex_t IT_0005 = lpp_102*U_utR_12;
    const complex_t IT_0006 = lpp_202*U_utR_22;
    const complex_t IT_0007 = IT_0004 + IT_0005 + IT_0006;
    const complex_t IT_0008 = -IT_0007;
    const complex_t IT_0009 = IT_0003 + IT_0008;
    const complex_t IT_0010 = (complex_t{0, 1})*IT_0009;
    return 4*s_34*IT_0010*std::conj(IT_0010);
}
} // End of namespace brparity
