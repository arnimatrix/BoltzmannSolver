#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_dtR1_to_uf2c_df2c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR1_to_uf2c_df2c(
        param_t const &param
        )
{
    clearcache();
    const real_t s_23 = param.s_23;
    const complex_t lpp_101 = param.lpp_101;
    const complex_t lpp_110 = param.lpp_110;
    const complex_t lpp_112 = param.lpp_112;
    const complex_t lpp_121 = param.lpp_121;
    const complex_t U_dtR_00 = param.U_dtR_00;
    const complex_t U_dtR_20 = param.U_dtR_20;
    const complex_t IT_0000 = lpp_101*U_dtR_00;
    const complex_t IT_0001 = lpp_110*U_dtR_00;
    const complex_t IT_0002 = lpp_112*U_dtR_20;
    const complex_t IT_0003 = lpp_121*U_dtR_20;
    const complex_t IT_0004 = (complex_t{0, 1})*(IT_0000 + -IT_0001 + -IT_0002
       + IT_0003);
    const complex_t IT_0005 = -IT_0004;
    return 4*s_23*IT_0005*std::conj(IT_0005);
}
} // End of namespace brparity
