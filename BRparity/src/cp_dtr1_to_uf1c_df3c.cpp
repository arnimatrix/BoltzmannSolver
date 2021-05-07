#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_dtr1_to_uf1c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR1_to_uf1c_df3c(
        param_t const &param
        )
{
    clearcache();
    auto const &s_34 = param.s_34;
    auto const &lpp_002 = param.lpp_002;
    auto const &lpp_012 = param.lpp_012;
    auto const &lpp_020 = param.lpp_020;
    auto const &lpp_021 = param.lpp_021;
    auto const &U_dtR_00 = param.U_dtR_00;
    auto const &U_dtR_10 = param.U_dtR_10;
    const complex_t IT_0000 = lpp_002*U_dtR_00;
    const complex_t IT_0001 = lpp_012*U_dtR_10;
    const complex_t IT_0002 = IT_0000 + IT_0001;
    const complex_t IT_0003 = lpp_020*U_dtR_00;
    const complex_t IT_0004 = lpp_021*U_dtR_10;
    const complex_t IT_0005 = -IT_0003 + -IT_0004;
    const complex_t IT_0006 = IT_0002 + IT_0005;
    const complex_t IT_0007 = (complex_t{0, 1})*IT_0006;
    return 4*s_34*IT_0007*std::conj(IT_0007);
}
} // End of namespace brparity
