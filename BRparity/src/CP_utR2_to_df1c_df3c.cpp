#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_utR2_to_df1c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_utR2_to_df1c_df3c(
        param_t const &param
        )
{
    clearcache();
    const real_t s_23 = param.s_23;
    const complex_t lpp_002 = param.lpp_002;
    const complex_t lpp_020 = param.lpp_020;
    const complex_t lpp_102 = param.lpp_102;
    const complex_t lpp_120 = param.lpp_120;
    const complex_t lpp_202 = param.lpp_202;
    const complex_t lpp_220 = param.lpp_220;
    const complex_t U_utR_01 = param.U_utR_01;
    const complex_t U_utR_11 = param.U_utR_11;
    const complex_t U_utR_21 = param.U_utR_21;
    const complex_t IT_0000 = lpp_002*U_utR_01;
    const complex_t IT_0001 = lpp_102*U_utR_11;
    const complex_t IT_0002 = lpp_202*U_utR_21;
    const complex_t IT_0003 = (complex_t{0, 1})*(IT_0000 + IT_0001 + IT_0002);
    const complex_t IT_0004 = -IT_0003;
    const complex_t IT_0005 = lpp_020*U_utR_01;
    const complex_t IT_0006 = lpp_120*U_utR_11;
    const complex_t IT_0007 = lpp_220*U_utR_21;
    const complex_t IT_0008 = (complex_t{0, 1})*(IT_0005 + IT_0006 + IT_0007);
    const complex_t IT_0009 = IT_0004 + IT_0008;
    const complex_t IT_0010 = -IT_0009;
    return 4*s_23*IT_0010*std::conj(IT_0010);
}
} // End of namespace brparity
