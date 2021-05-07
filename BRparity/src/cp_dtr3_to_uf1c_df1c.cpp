#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_dtr3_to_uf1c_df1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR3_to_uf1c_df1c(
        param_t const &param
        )
{
    clearcache();
    auto const &s_34 = param.s_34;
    auto const &lpp_001 = param.lpp_001;
    auto const &lpp_002 = param.lpp_002;
    auto const &lpp_010 = param.lpp_010;
    auto const &lpp_020 = param.lpp_020;
    auto const &U_dtR_12 = param.U_dtR_12;
    auto const &U_dtR_22 = param.U_dtR_22;
    const complex_t IT_0000 = lpp_001*U_dtR_12;
    const complex_t IT_0001 = lpp_002*U_dtR_22;
    const complex_t IT_0002 = IT_0000 + IT_0001;
    const complex_t IT_0003 = lpp_010*U_dtR_12;
    const complex_t IT_0004 = lpp_020*U_dtR_22;
    const complex_t IT_0005 = -IT_0003 + -IT_0004;
    const complex_t IT_0006 = IT_0002 + IT_0005;
    const complex_t IT_0007 = (complex_t{0, 1})*IT_0006;
    const complex_t IT_0008 = -IT_0007;
    return 4*s_34*IT_0008*std::conj(IT_0008);
}
} // End of namespace brparity
