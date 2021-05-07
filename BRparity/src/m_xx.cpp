#include "clooptools.h"
#include "marty/looptools_init.h"
#include "m_xx.h"
#include "common.h"

#include "params.h"
#include "group_g.h"

#include "global.h"
#include "libcomplexop.h"

namespace brparity {

complex_t m_XX(
        param_t const &param
        )
{
    clearcache();
    auto const &mX = param.mX;
    const complex_t IT_0000 = 2*mX;
    return IT_0000;
}
} // End of namespace brparity
