#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_dtr2_to_uf2c_df2c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR2_to_uf2c_df2c(
        param_t const &param
        )
{
    clearcache();
    auto const &s_34 = param.s_34;
    auto const &lpp_101 = param.lpp_101;
    auto const &lpp_110 = param.lpp_110;
    auto const &lpp_112 = param.lpp_112;
    auto const &lpp_121 = param.lpp_121;
    auto const &U_dtR_01 = param.U_dtR_01;
    auto const &U_dtR_21 = param.U_dtR_21;
    const complex_t IT_0000 = lpp_101*U_dtR_01;
    const complex_t IT_0001 = lpp_121*U_dtR_21;
    const complex_t IT_0002 = IT_0000 + IT_0001;
    const complex_t IT_0003 = lpp_110*U_dtR_01;
    const complex_t IT_0004 = lpp_112*U_dtR_21;
    const complex_t IT_0005 = -IT_0003 + -IT_0004;
    const complex_t IT_0006 = IT_0002 + IT_0005;
    const complex_t IT_0007 = (complex_t{0, 1})*IT_0006;
    const complex_t IT_0008 = -IT_0007;
    return 4*s_34*IT_0008*std::conj(IT_0008);
}
} // End of namespace brparity
