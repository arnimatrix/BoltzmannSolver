#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_dtR3_to_uf2c_df1c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR3_to_uf2c_df1c(
        param_t const &param
        )
{
    clearcache();
    const real_t s_23 = param.s_23;
    const complex_t lpp_101 = param.lpp_101;
    const complex_t lpp_102 = param.lpp_102;
    const complex_t lpp_110 = param.lpp_110;
    const complex_t lpp_120 = param.lpp_120;
    const complex_t U_dtR_12 = param.U_dtR_12;
    const complex_t U_dtR_22 = param.U_dtR_22;
    const complex_t IT_0000 = lpp_101*U_dtR_12;
    const complex_t IT_0001 = lpp_102*U_dtR_22;
    const complex_t IT_0002 = lpp_110*U_dtR_12;
    const complex_t IT_0003 = lpp_120*U_dtR_22;
    const complex_t IT_0004 = (complex_t{0, -1})*(IT_0000 + IT_0001 + -IT_0002
       + -IT_0003);
    return 4*s_23*IT_0004*std::conj(IT_0004);
}
} // End of namespace brparity
