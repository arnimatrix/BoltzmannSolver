#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_utR3_to_df2c_df3c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_utR3_to_df2c_df3c(
        param_t const &param
        )
{
    clearcache();
    const real_t s_23 = param.s_23;
    const complex_t lpp_012 = param.lpp_012;
    const complex_t lpp_021 = param.lpp_021;
    const complex_t lpp_112 = param.lpp_112;
    const complex_t lpp_121 = param.lpp_121;
    const complex_t lpp_212 = param.lpp_212;
    const complex_t lpp_221 = param.lpp_221;
    const complex_t U_utR_02 = param.U_utR_02;
    const complex_t U_utR_12 = param.U_utR_12;
    const complex_t U_utR_22 = param.U_utR_22;
    const complex_t IT_0000 = lpp_012*U_utR_02;
    const complex_t IT_0001 = lpp_112*U_utR_12;
    const complex_t IT_0002 = lpp_212*U_utR_22;
    const complex_t IT_0003 = (complex_t{0, 1})*(IT_0000 + IT_0001 + IT_0002);
    const complex_t IT_0004 = -IT_0003;
    const complex_t IT_0005 = lpp_021*U_utR_02;
    const complex_t IT_0006 = lpp_121*U_utR_12;
    const complex_t IT_0007 = lpp_221*U_utR_22;
    const complex_t IT_0008 = (complex_t{0, 1})*(IT_0005 + IT_0006 + IT_0007);
    const complex_t IT_0009 = IT_0004 + IT_0008;
    return 4*s_23*IT_0009*std::conj(IT_0009);
}
} // End of namespace brparity
