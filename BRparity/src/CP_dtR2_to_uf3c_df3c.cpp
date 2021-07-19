#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_dtR2_to_uf3c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR2_to_uf3c_df3c(
        param_t const &param
        )
{
    clearcache();
    const real_t s_23 = param.s_23;
    const complex_t lpp_202 = param.lpp_202;
    const complex_t lpp_212 = param.lpp_212;
    const complex_t lpp_220 = param.lpp_220;
    const complex_t lpp_221 = param.lpp_221;
    const complex_t U_dtR_01 = param.U_dtR_01;
    const complex_t U_dtR_11 = param.U_dtR_11;
    const complex_t IT_0000 = lpp_202*U_dtR_01;
    const complex_t IT_0001 = lpp_212*U_dtR_11;
    const complex_t IT_0002 = lpp_220*U_dtR_01;
    const complex_t IT_0003 = lpp_221*U_dtR_11;
    const complex_t IT_0004 = (complex_t{0, 1})*(IT_0000 + IT_0001 + -IT_0002 
      + -IT_0003);
    const complex_t IT_0005 = -IT_0004;
    return 4*s_23*IT_0005*std::conj(IT_0005);
}
} // End of namespace brparity
