#include "clooptools.h"
#include "marty/looptools_init.h"
#include "CP_dtR1_to_df3_XX.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t CP_dtR1_to_df3_XX(
        param_t const &param
        )
{
    clearcache();
    const real_t gw = param.gw;
    const real_t gwdR = param.gwdR;
    const real_t s_23 = param.s_23;
    const complex_t U_dtR_20 = param.U_dtR_20;
    const complex_t IT_0000 = (complex_t{0, -1})*gw*gwdR*U_dtR_20;
    return 2*s_23*IT_0000*std::conj(IT_0000);
}
} // End of namespace brparity