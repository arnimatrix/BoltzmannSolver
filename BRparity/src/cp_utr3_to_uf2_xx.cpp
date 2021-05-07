#include "clooptools.h"
#include "marty/looptools_init.h"
#include "cp_utr3_to_uf2_xx.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_utR3_to_uf2_XX(
        param_t const &param
        )
{
    clearcache();
    auto const &gw = param.gw;
    auto const &gwuR = param.gwuR;
    auto const &s_34 = param.s_34;
    auto const &U_utR_12 = param.U_utR_12;
    const complex_t IT_0000 = gw*gwuR*U_utR_12;
    const complex_t IT_0001 = (complex_t{0, 1})*IT_0000;
    const complex_t IT_0002 = -IT_0001;
    return 2*s_34*IT_0002*std::conj(IT_0002);
}
} // End of namespace brparity