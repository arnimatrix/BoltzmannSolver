#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_dtr1_to_uf3c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR1_to_uf3c_df3c(
        param_t const &param
        )
{
    clearcache();
    auto const &s_34 = param.s_34;
    auto const &lpp_202 = param.lpp_202;
    auto const &lpp_212 = param.lpp_212;
    auto const &lpp_220 = param.lpp_220;
    auto const &lpp_221 = param.lpp_221;
    auto const &U_dtR_00 = param.U_dtR_00;
    auto const &U_dtR_10 = param.U_dtR_10;
    const complex_t IT_0000 = lpp_202*U_dtR_00;
    const complex_t IT_0001 = lpp_212*U_dtR_10;
    const complex_t IT_0002 = IT_0000 + IT_0001;
    const complex_t IT_0003 = lpp_220*U_dtR_00;
    const complex_t IT_0004 = lpp_221*U_dtR_10;
    const complex_t IT_0005 = -IT_0003 + -IT_0004;
    const complex_t IT_0006 = IT_0002 + IT_0005;
    const complex_t IT_0007 = (complex_t{0, 1})*IT_0006;
    const complex_t IT_0008 = -IT_0007;
    return 4*s_34*IT_0008*std::conj(IT_0008);
}
} // End of namespace brparity
