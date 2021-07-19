#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_dtR3_to_uf1c_df1c.h"
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
    const real_t s_23 = param.s_23;
    const complex_t lpp_001 = param.lpp_001;
    const complex_t lpp_002 = param.lpp_002;
    const complex_t lpp_010 = param.lpp_010;
    const complex_t lpp_020 = param.lpp_020;
    const complex_t U_dtR_12 = param.U_dtR_12;
    const complex_t U_dtR_22 = param.U_dtR_22;
    const complex_t IT_0000 = lpp_001*U_dtR_12;
    const complex_t IT_0001 = lpp_002*U_dtR_22;
    const complex_t IT_0002 = lpp_010*U_dtR_12;
    const complex_t IT_0003 = lpp_020*U_dtR_22;
    const complex_t IT_0004 = (complex_t{0, -1})*(IT_0000 + IT_0001 + -IT_0002
       + -IT_0003);
    return 4*s_23*IT_0004*std::conj(IT_0004);
}
} // End of namespace brparity
