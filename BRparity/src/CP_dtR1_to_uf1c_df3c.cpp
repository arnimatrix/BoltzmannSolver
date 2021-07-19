#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_dtR1_to_uf1c_df3c.h"
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
    const real_t s_23 = param.s_23;
    const complex_t lpp_002 = param.lpp_002;
    const complex_t lpp_012 = param.lpp_012;
    const complex_t lpp_020 = param.lpp_020;
    const complex_t lpp_021 = param.lpp_021;
    const complex_t U_dtR_00 = param.U_dtR_00;
    const complex_t U_dtR_10 = param.U_dtR_10;
    const complex_t IT_0000 = lpp_002*U_dtR_00;
    const complex_t IT_0001 = lpp_012*U_dtR_10;
    const complex_t IT_0002 = lpp_020*U_dtR_00;
    const complex_t IT_0003 = lpp_021*U_dtR_10;
    const complex_t IT_0004 = (complex_t{0, 1})*(IT_0000 + IT_0001 + -IT_0002 
      + -IT_0003);
    return 4*s_23*IT_0004*std::conj(IT_0004);
}
} // End of namespace brparity
