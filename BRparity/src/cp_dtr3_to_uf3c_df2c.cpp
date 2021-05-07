#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_dtr3_to_uf3c_df2c.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR3_to_uf3c_df2c(
        param_t const &param
        )
{
    clearcache();
    auto const &s_34 = param.s_34;
    auto const &lpp_201 = param.lpp_201;
    auto const &lpp_210 = param.lpp_210;
    auto const &lpp_212 = param.lpp_212;
    auto const &lpp_221 = param.lpp_221;
    auto const &U_dtR_02 = param.U_dtR_02;
    auto const &U_dtR_22 = param.U_dtR_22;
    const complex_t IT_0000 = lpp_201*U_dtR_02;
    const complex_t IT_0001 = lpp_221*U_dtR_22;
    const complex_t IT_0002 = IT_0000 + IT_0001;
    const complex_t IT_0003 = lpp_210*U_dtR_02;
    const complex_t IT_0004 = lpp_212*U_dtR_22;
    const complex_t IT_0005 = -IT_0003 + -IT_0004;
    const complex_t IT_0006 = IT_0002 + IT_0005;
    const complex_t IT_0007 = (complex_t{0, 1})*IT_0006;
    return 4*s_34*IT_0007*std::conj(IT_0007);
}
} // End of namespace brparity
