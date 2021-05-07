#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_utr1_to_df2c_df1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_utR1_to_df2c_df1c(
        param_t const &param
        )
{
    clearcache();
    auto const &s_34 = param.s_34;
    auto const &lpp_001 = param.lpp_001;
    auto const &lpp_010 = param.lpp_010;
    auto const &lpp_101 = param.lpp_101;
    auto const &lpp_110 = param.lpp_110;
    auto const &lpp_201 = param.lpp_201;
    auto const &lpp_210 = param.lpp_210;
    auto const &U_utR_00 = param.U_utR_00;
    auto const &U_utR_10 = param.U_utR_10;
    auto const &U_utR_20 = param.U_utR_20;
    const complex_t IT_0000 = lpp_010*U_utR_00;
    const complex_t IT_0001 = lpp_110*U_utR_10;
    const complex_t IT_0002 = lpp_210*U_utR_20;
    const complex_t IT_0003 = IT_0000 + IT_0001 + IT_0002;
    const complex_t IT_0004 = lpp_001*U_utR_00;
    const complex_t IT_0005 = lpp_101*U_utR_10;
    const complex_t IT_0006 = lpp_201*U_utR_20;
    const complex_t IT_0007 = IT_0004 + IT_0005 + IT_0006;
    const complex_t IT_0008 = -IT_0007;
    const complex_t IT_0009 = IT_0003 + IT_0008;
    const complex_t IT_0010 = (complex_t{0, 1})*IT_0009;
    return 4*s_34*IT_0010*std::conj(IT_0010);
}
} // End of namespace brparity
