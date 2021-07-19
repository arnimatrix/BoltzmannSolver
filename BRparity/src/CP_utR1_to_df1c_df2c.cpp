#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_utR1_to_df1c_df2c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_utR1_to_df1c_df2c(
        param_t const &param
        )
{
    clearcache();
    const real_t s_23 = param.s_23;
    const complex_t lpp_001 = param.lpp_001;
    const complex_t lpp_010 = param.lpp_010;
    const complex_t lpp_101 = param.lpp_101;
    const complex_t lpp_110 = param.lpp_110;
    const complex_t lpp_201 = param.lpp_201;
    const complex_t lpp_210 = param.lpp_210;
    const complex_t U_utR_00 = param.U_utR_00;
    const complex_t U_utR_10 = param.U_utR_10;
    const complex_t U_utR_20 = param.U_utR_20;
    const complex_t IT_0000 = lpp_001*U_utR_00;
    const complex_t IT_0001 = lpp_101*U_utR_10;
    const complex_t IT_0002 = lpp_201*U_utR_20;
    const complex_t IT_0003 = (complex_t{0, 1})*(IT_0000 + IT_0001 + IT_0002);
    const complex_t IT_0004 = -IT_0003;
    const complex_t IT_0005 = lpp_010*U_utR_00;
    const complex_t IT_0006 = lpp_110*U_utR_10;
    const complex_t IT_0007 = lpp_210*U_utR_20;
    const complex_t IT_0008 = (complex_t{0, 1})*(IT_0005 + IT_0006 + IT_0007);
    const complex_t IT_0009 = IT_0004 + IT_0008;
    const complex_t IT_0010 = -IT_0009;
    return 4*s_23*IT_0010*std::conj(IT_0010);
}
} // End of namespace brparity
