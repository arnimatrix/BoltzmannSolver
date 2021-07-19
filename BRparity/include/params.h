#ifndef CSL_LIB_PARAM_H_INCLUDED
#define CSL_LIB_PARAM_H_INCLUDED

#include <array>
#include "common.h"
#include "libcomplexop.h"
#include "csl/initSanitizer.h"

namespace brparity {

struct param_t {

    ///////////////////////////////////////
    // Elementary parameters to be defined 
    ///////////////////////////////////////

    csl::InitSanitizer<real_t> gw { "gw" };
    csl::InitSanitizer<real_t> mX { "mX" };
    csl::InitSanitizer<real_t> m_b { "m_b" };
    csl::InitSanitizer<real_t> m_c { "m_c" };
    csl::InitSanitizer<real_t> m_d { "m_d" };
    csl::InitSanitizer<real_t> m_s { "m_s" };
    csl::InitSanitizer<real_t> m_t { "m_t" };
    csl::InitSanitizer<real_t> m_u { "m_u" };
    csl::InitSanitizer<real_t> GbLt { "GbLt" };
    csl::InitSanitizer<real_t> GbRt { "GbRt" };
    csl::InitSanitizer<real_t> GcLt { "GcLt" };
    csl::InitSanitizer<real_t> GcRt { "GcRt" };
    csl::InitSanitizer<real_t> GdLt { "GdLt" };
    csl::InitSanitizer<real_t> GdRt { "GdRt" };
    csl::InitSanitizer<real_t> GsLt { "GsLt" };
    csl::InitSanitizer<real_t> GsRt { "GsRt" };
    csl::InitSanitizer<real_t> GtLt { "GtLt" };
    csl::InitSanitizer<real_t> GtRt { "GtRt" };
    csl::InitSanitizer<real_t> GuLt { "GuLt" };
    csl::InitSanitizer<real_t> GuRt { "GuRt" };
    csl::InitSanitizer<real_t> gwdL { "gwdL" };
    csl::InitSanitizer<real_t> gwdR { "gwdR" };
    csl::InitSanitizer<real_t> gwuL { "gwuL" };
    csl::InitSanitizer<real_t> gwuR { "gwuR" };
    csl::InitSanitizer<real_t> s_12 { "s_12" };
    csl::InitSanitizer<real_t> s_13 { "s_13" };
    csl::InitSanitizer<real_t> s_14 { "s_14" };
    csl::InitSanitizer<real_t> s_23 { "s_23" };
    csl::InitSanitizer<real_t> s_24 { "s_24" };
    csl::InitSanitizer<real_t> s_34 { "s_34" };
    csl::InitSanitizer<real_t> reg_prop { "reg_prop" };
    csl::InitSanitizer<complex_t> mdL_00 { "mdL_00" };
    csl::InitSanitizer<complex_t> mdL_01 { "mdL_01" };
    csl::InitSanitizer<complex_t> mdL_02 { "mdL_02" };
    csl::InitSanitizer<complex_t> mdL_10 { "mdL_10" };
    csl::InitSanitizer<complex_t> mdL_11 { "mdL_11" };
    csl::InitSanitizer<complex_t> mdL_12 { "mdL_12" };
    csl::InitSanitizer<complex_t> mdL_20 { "mdL_20" };
    csl::InitSanitizer<complex_t> mdL_21 { "mdL_21" };
    csl::InitSanitizer<complex_t> mdL_22 { "mdL_22" };
    csl::InitSanitizer<complex_t> mdR_00 { "mdR_00" };
    csl::InitSanitizer<complex_t> mdR_01 { "mdR_01" };
    csl::InitSanitizer<complex_t> mdR_02 { "mdR_02" };
    csl::InitSanitizer<complex_t> mdR_10 { "mdR_10" };
    csl::InitSanitizer<complex_t> mdR_11 { "mdR_11" };
    csl::InitSanitizer<complex_t> mdR_12 { "mdR_12" };
    csl::InitSanitizer<complex_t> mdR_20 { "mdR_20" };
    csl::InitSanitizer<complex_t> mdR_21 { "mdR_21" };
    csl::InitSanitizer<complex_t> mdR_22 { "mdR_22" };
    csl::InitSanitizer<complex_t> muL_00 { "muL_00" };
    csl::InitSanitizer<complex_t> muL_01 { "muL_01" };
    csl::InitSanitizer<complex_t> muL_02 { "muL_02" };
    csl::InitSanitizer<complex_t> muL_10 { "muL_10" };
    csl::InitSanitizer<complex_t> muL_11 { "muL_11" };
    csl::InitSanitizer<complex_t> muL_12 { "muL_12" };
    csl::InitSanitizer<complex_t> muL_20 { "muL_20" };
    csl::InitSanitizer<complex_t> muL_21 { "muL_21" };
    csl::InitSanitizer<complex_t> muL_22 { "muL_22" };
    csl::InitSanitizer<complex_t> muR_00 { "muR_00" };
    csl::InitSanitizer<complex_t> muR_01 { "muR_01" };
    csl::InitSanitizer<complex_t> muR_02 { "muR_02" };
    csl::InitSanitizer<complex_t> muR_10 { "muR_10" };
    csl::InitSanitizer<complex_t> muR_11 { "muR_11" };
    csl::InitSanitizer<complex_t> muR_12 { "muR_12" };
    csl::InitSanitizer<complex_t> muR_20 { "muR_20" };
    csl::InitSanitizer<complex_t> muR_21 { "muR_21" };
    csl::InitSanitizer<complex_t> muR_22 { "muR_22" };
    csl::InitSanitizer<complex_t> lpp_001 { "lpp_001" };
    csl::InitSanitizer<complex_t> lpp_002 { "lpp_002" };
    csl::InitSanitizer<complex_t> lpp_010 { "lpp_010" };
    csl::InitSanitizer<complex_t> lpp_012 { "lpp_012" };
    csl::InitSanitizer<complex_t> lpp_020 { "lpp_020" };
    csl::InitSanitizer<complex_t> lpp_021 { "lpp_021" };
    csl::InitSanitizer<complex_t> lpp_101 { "lpp_101" };
    csl::InitSanitizer<complex_t> lpp_102 { "lpp_102" };
    csl::InitSanitizer<complex_t> lpp_110 { "lpp_110" };
    csl::InitSanitizer<complex_t> lpp_112 { "lpp_112" };
    csl::InitSanitizer<complex_t> lpp_120 { "lpp_120" };
    csl::InitSanitizer<complex_t> lpp_121 { "lpp_121" };
    csl::InitSanitizer<complex_t> lpp_201 { "lpp_201" };
    csl::InitSanitizer<complex_t> lpp_202 { "lpp_202" };
    csl::InitSanitizer<complex_t> lpp_210 { "lpp_210" };
    csl::InitSanitizer<complex_t> lpp_212 { "lpp_212" };
    csl::InitSanitizer<complex_t> lpp_220 { "lpp_220" };
    csl::InitSanitizer<complex_t> lpp_221 { "lpp_221" };


    ///////////////////////////////////////
    // Parameters functions of others  
    // through diagonalization or mass 
    // expressions, see updateSpectrum()  
    // in global.h or set them by hand  
    ///////////////////////////////////////

    csl::InitSanitizer<real_t> m_XX { "m_XX" };
    csl::InitSanitizer<real_t> m_dtL1 { "m_dtL1" };
    csl::InitSanitizer<real_t> m_dtL2 { "m_dtL2" };
    csl::InitSanitizer<real_t> m_dtL3 { "m_dtL3" };
    csl::InitSanitizer<real_t> m_dtR1 { "m_dtR1" };
    csl::InitSanitizer<real_t> m_dtR2 { "m_dtR2" };
    csl::InitSanitizer<real_t> m_dtR3 { "m_dtR3" };
    csl::InitSanitizer<real_t> m_utL1 { "m_utL1" };
    csl::InitSanitizer<real_t> m_utL2 { "m_utL2" };
    csl::InitSanitizer<real_t> m_utL3 { "m_utL3" };
    csl::InitSanitizer<real_t> m_utR1 { "m_utR1" };
    csl::InitSanitizer<real_t> m_utR2 { "m_utR2" };
    csl::InitSanitizer<real_t> m_utR3 { "m_utR3" };
    csl::InitSanitizer<complex_t> U_dtL_00 { "U_dtL_00" };
    csl::InitSanitizer<complex_t> U_dtL_01 { "U_dtL_01" };
    csl::InitSanitizer<complex_t> U_dtL_02 { "U_dtL_02" };
    csl::InitSanitizer<complex_t> U_dtL_10 { "U_dtL_10" };
    csl::InitSanitizer<complex_t> U_dtL_11 { "U_dtL_11" };
    csl::InitSanitizer<complex_t> U_dtL_12 { "U_dtL_12" };
    csl::InitSanitizer<complex_t> U_dtL_20 { "U_dtL_20" };
    csl::InitSanitizer<complex_t> U_dtL_21 { "U_dtL_21" };
    csl::InitSanitizer<complex_t> U_dtL_22 { "U_dtL_22" };
    csl::InitSanitizer<complex_t> U_dtR_00 { "U_dtR_00" };
    csl::InitSanitizer<complex_t> U_dtR_01 { "U_dtR_01" };
    csl::InitSanitizer<complex_t> U_dtR_02 { "U_dtR_02" };
    csl::InitSanitizer<complex_t> U_dtR_10 { "U_dtR_10" };
    csl::InitSanitizer<complex_t> U_dtR_11 { "U_dtR_11" };
    csl::InitSanitizer<complex_t> U_dtR_12 { "U_dtR_12" };
    csl::InitSanitizer<complex_t> U_dtR_20 { "U_dtR_20" };
    csl::InitSanitizer<complex_t> U_dtR_21 { "U_dtR_21" };
    csl::InitSanitizer<complex_t> U_dtR_22 { "U_dtR_22" };
    csl::InitSanitizer<complex_t> U_utL_00 { "U_utL_00" };
    csl::InitSanitizer<complex_t> U_utL_01 { "U_utL_01" };
    csl::InitSanitizer<complex_t> U_utL_02 { "U_utL_02" };
    csl::InitSanitizer<complex_t> U_utL_10 { "U_utL_10" };
    csl::InitSanitizer<complex_t> U_utL_11 { "U_utL_11" };
    csl::InitSanitizer<complex_t> U_utL_12 { "U_utL_12" };
    csl::InitSanitizer<complex_t> U_utL_20 { "U_utL_20" };
    csl::InitSanitizer<complex_t> U_utL_21 { "U_utL_21" };
    csl::InitSanitizer<complex_t> U_utL_22 { "U_utL_22" };
    csl::InitSanitizer<complex_t> U_utR_00 { "U_utR_00" };
    csl::InitSanitizer<complex_t> U_utR_01 { "U_utR_01" };
    csl::InitSanitizer<complex_t> U_utR_02 { "U_utR_02" };
    csl::InitSanitizer<complex_t> U_utR_10 { "U_utR_10" };
    csl::InitSanitizer<complex_t> U_utR_11 { "U_utR_11" };
    csl::InitSanitizer<complex_t> U_utR_12 { "U_utR_12" };
    csl::InitSanitizer<complex_t> U_utR_20 { "U_utR_20" };
    csl::InitSanitizer<complex_t> U_utR_21 { "U_utR_21" };
    csl::InitSanitizer<complex_t> U_utR_22 { "U_utR_22" };

    void reset()
    {
        using real_params = std::array<csl::InitSanitizer<real_t>*, 44>;
        using complex_params = std::array<csl::InitSanitizer<complex_t>*, 90>;

        for (auto &par : real_params{
                &gw, &mX, &m_b, &m_c, &m_d, 
                &m_s, &m_t, &m_u, &GbLt, &GbRt, &GcLt, 
                &GcRt, &GdLt, &GdRt, &GsLt, &GsRt, &GtLt, 
                &GtRt, &GuLt, &GuRt, &gwdL, &gwdR, &gwuL, 
                &gwuR, &s_12, &s_13, &s_14, &s_23, &s_24, 
                &s_34, &reg_prop, &m_XX, &m_dtL1, &m_dtL2, &m_dtL3, 
                &m_dtR1, &m_dtR2, &m_dtR3, &m_utL1, &m_utL2, &m_utL3, 
                &m_utR1, &m_utR2, &m_utR3, })
        {
            par->reset();
        }

        for (auto &par : complex_params{
                &mdL_00, &mdL_01, &mdL_02, &mdL_10, &mdL_11, 
                &mdL_12, &mdL_20, &mdL_21, &mdL_22, &mdR_00, &mdR_01, 
                &mdR_02, &mdR_10, &mdR_11, &mdR_12, &mdR_20, &mdR_21, 
                &mdR_22, &muL_00, &muL_01, &muL_02, &muL_10, &muL_11, 
                &muL_12, &muL_20, &muL_21, &muL_22, &muR_00, &muR_01, 
                &muR_02, &muR_10, &muR_11, &muR_12, &muR_20, &muR_21, 
                &muR_22, &lpp_001, &lpp_002, &lpp_010, &lpp_012, &lpp_020, 
                &lpp_021, &lpp_101, &lpp_102, &lpp_110, &lpp_112, &lpp_120, 
                &lpp_121, &lpp_201, &lpp_202, &lpp_210, &lpp_212, &lpp_220, 
                &lpp_221, &U_dtL_00, &U_dtL_01, &U_dtL_02, &U_dtL_10, &U_dtL_11, 
                &U_dtL_12, &U_dtL_20, &U_dtL_21, &U_dtL_22, &U_dtR_00, &U_dtR_01, 
                &U_dtR_02, &U_dtR_10, &U_dtR_11, &U_dtR_12, &U_dtR_20, &U_dtR_21, 
                &U_dtR_22, &U_utL_00, &U_utL_01, &U_utL_02, &U_utL_10, &U_utL_11, 
                &U_utL_12, &U_utL_20, &U_utL_21, &U_utL_22, &U_utR_00, &U_utR_01, 
                &U_utR_02, &U_utR_10, &U_utR_11, &U_utR_12, &U_utR_20, &U_utR_21, 
                &U_utR_22, })
        {
            par->reset();
        }
    }

    void print(std::ostream &out = std::cout)
    {
        using real_params = std::array<csl::InitSanitizer<real_t> const*, 44>;
        using complex_params = std::array<csl::InitSanitizer<complex_t> const*, 90>;

        out << "param_t struct:\n";
        out << "Real parameters\n";
        for (auto const &par : real_params{
                &gw, &mX, &m_b, &m_c, &m_d, 
                &m_s, &m_t, &m_u, &GbLt, &GbRt, &GcLt, 
                &GcRt, &GdLt, &GdRt, &GsLt, &GsRt, &GtLt, 
                &GtRt, &GuLt, &GuRt, &gwdL, &gwdR, &gwuL, 
                &gwuR, &s_12, &s_13, &s_14, &s_23, &s_24, 
                &s_34, &reg_prop, &m_XX, &m_dtL1, &m_dtL2, &m_dtL3, 
                &m_dtR1, &m_dtR2, &m_dtR3, &m_utL1, &m_utL2, &m_utL3, 
                &m_utR1, &m_utR2, &m_utR3, })
        {
            out << "  -> ";
            par->print(out);
        }

        out << "Complex parameters\n";
        for (auto const &par : complex_params{
                &mdL_00, &mdL_01, &mdL_02, &mdL_10, &mdL_11, 
                &mdL_12, &mdL_20, &mdL_21, &mdL_22, &mdR_00, &mdR_01, 
                &mdR_02, &mdR_10, &mdR_11, &mdR_12, &mdR_20, &mdR_21, 
                &mdR_22, &muL_00, &muL_01, &muL_02, &muL_10, &muL_11, 
                &muL_12, &muL_20, &muL_21, &muL_22, &muR_00, &muR_01, 
                &muR_02, &muR_10, &muR_11, &muR_12, &muR_20, &muR_21, 
                &muR_22, &lpp_001, &lpp_002, &lpp_010, &lpp_012, &lpp_020, 
                &lpp_021, &lpp_101, &lpp_102, &lpp_110, &lpp_112, &lpp_120, 
                &lpp_121, &lpp_201, &lpp_202, &lpp_210, &lpp_212, &lpp_220, 
                &lpp_221, &U_dtL_00, &U_dtL_01, &U_dtL_02, &U_dtL_10, &U_dtL_11, 
                &U_dtL_12, &U_dtL_20, &U_dtL_21, &U_dtL_22, &U_dtR_00, &U_dtR_01, 
                &U_dtR_02, &U_dtR_10, &U_dtR_11, &U_dtR_12, &U_dtR_20, &U_dtR_21, 
                &U_dtR_22, &U_utL_00, &U_utL_01, &U_utL_02, &U_utL_10, &U_utL_11, 
                &U_utL_12, &U_utL_20, &U_utL_21, &U_utL_22, &U_utR_00, &U_utR_01, 
                &U_utR_02, &U_utR_10, &U_utR_11, &U_utR_12, &U_utR_20, &U_utR_21, 
                &U_utR_22, })
        {
            out << "  -> ";
            par->print(out);
        }
        out << "\n";
    }
};


}

#endif
