#ifndef CSL_LIB_GLOBAL
#define CSL_LIB_GLOBAL
#include "params.h"
#include "common.h"

namespace brparity {

struct SpectrumInput;
struct SpectrumOutput;

SpectrumOutput updateDiagonalization(SpectrumInput const&);

void updateDiagonalization(param_t &params);

////////////////////////////////////////////////////
// Here are the parameters to set before calling    
// updateDiagonalization()                          
////////////////////////////////////////////////////
struct SpectrumInput {
    complex_t mdL_00;
    complex_t mdL_01;
    complex_t mdL_02;
    complex_t mdL_10;
    complex_t mdL_11;
    complex_t mdL_12;
    complex_t mdL_20;
    complex_t mdL_21;
    complex_t mdL_22;
    complex_t mdR_00;
    complex_t mdR_01;
    complex_t mdR_02;
    complex_t mdR_10;
    complex_t mdR_11;
    complex_t mdR_12;
    complex_t mdR_20;
    complex_t mdR_21;
    complex_t mdR_22;
    complex_t muL_00;
    complex_t muL_01;
    complex_t muL_02;
    complex_t muL_10;
    complex_t muL_11;
    complex_t muL_12;
    complex_t muL_20;
    complex_t muL_21;
    complex_t muL_22;
    complex_t muR_00;
    complex_t muR_01;
    complex_t muR_02;
    complex_t muR_10;
    complex_t muR_11;
    complex_t muR_12;
    complex_t muR_20;
    complex_t muR_21;
    complex_t muR_22;
};

////////////////////////////////////////////////////
// Here are the masses and mixings                 
// result of the diagonalization                   
////////////////////////////////////////////////////
struct SpectrumOutput {
    real_t m_dtL1;
    real_t m_dtL2;
    real_t m_dtL3;
    real_t m_dtR1;
    real_t m_dtR2;
    real_t m_dtR3;
    real_t m_utL1;
    real_t m_utL2;
    real_t m_utL3;
    real_t m_utR1;
    real_t m_utR2;
    real_t m_utR3;

    complex_t U_dtL_00;
    complex_t U_dtL_01;
    complex_t U_dtL_02;
    complex_t U_dtL_10;
    complex_t U_dtL_11;
    complex_t U_dtL_12;
    complex_t U_dtL_20;
    complex_t U_dtL_21;
    complex_t U_dtL_22;
    complex_t U_dtR_00;
    complex_t U_dtR_01;
    complex_t U_dtR_02;
    complex_t U_dtR_10;
    complex_t U_dtR_11;
    complex_t U_dtR_12;
    complex_t U_dtR_20;
    complex_t U_dtR_21;
    complex_t U_dtR_22;
    complex_t U_utL_00;
    complex_t U_utL_01;
    complex_t U_utL_02;
    complex_t U_utL_10;
    complex_t U_utL_11;
    complex_t U_utL_12;
    complex_t U_utL_20;
    complex_t U_utL_21;
    complex_t U_utL_22;
    complex_t U_utR_00;
    complex_t U_utR_01;
    complex_t U_utR_02;
    complex_t U_utR_10;
    complex_t U_utR_11;
    complex_t U_utR_12;
    complex_t U_utR_20;
    complex_t U_utR_21;
    complex_t U_utR_22;
};

////////////////////////////////////////////////////
// Here is a generic function to read results      
// of the diagonalization in a corresponding struct
////////////////////////////////////////////////////

template<class Type>
void readDiagonalizationInputs(
        SpectrumInput &diagData,
        Type    const &input
        )
{
    diagData.mdL_00 = input.mdL_00;
    diagData.mdL_01 = input.mdL_01;
    diagData.mdL_02 = input.mdL_02;
    diagData.mdL_10 = input.mdL_10;
    diagData.mdL_11 = input.mdL_11;
    diagData.mdL_12 = input.mdL_12;
    diagData.mdL_20 = input.mdL_20;
    diagData.mdL_21 = input.mdL_21;
    diagData.mdL_22 = input.mdL_22;
    diagData.mdR_00 = input.mdR_00;
    diagData.mdR_01 = input.mdR_01;
    diagData.mdR_02 = input.mdR_02;
    diagData.mdR_10 = input.mdR_10;
    diagData.mdR_11 = input.mdR_11;
    diagData.mdR_12 = input.mdR_12;
    diagData.mdR_20 = input.mdR_20;
    diagData.mdR_21 = input.mdR_21;
    diagData.mdR_22 = input.mdR_22;
    diagData.muL_00 = input.muL_00;
    diagData.muL_01 = input.muL_01;
    diagData.muL_02 = input.muL_02;
    diagData.muL_10 = input.muL_10;
    diagData.muL_11 = input.muL_11;
    diagData.muL_12 = input.muL_12;
    diagData.muL_20 = input.muL_20;
    diagData.muL_21 = input.muL_21;
    diagData.muL_22 = input.muL_22;
    diagData.muR_00 = input.muR_00;
    diagData.muR_01 = input.muR_01;
    diagData.muR_02 = input.muR_02;
    diagData.muR_10 = input.muR_10;
    diagData.muR_11 = input.muR_11;
    diagData.muR_12 = input.muR_12;
    diagData.muR_20 = input.muR_20;
    diagData.muR_21 = input.muR_21;
    diagData.muR_22 = input.muR_22;
}

template<class Type>
void readDiagonalizationOutputs(
        SpectrumOutput const &diagData,
        Type                 &output
        )
{
    output.m_dtL1 = diagData.m_dtL1;
    output.m_dtL2 = diagData.m_dtL2;
    output.m_dtL3 = diagData.m_dtL3;
    output.m_dtR1 = diagData.m_dtR1;
    output.m_dtR2 = diagData.m_dtR2;
    output.m_dtR3 = diagData.m_dtR3;
    output.m_utL1 = diagData.m_utL1;
    output.m_utL2 = diagData.m_utL2;
    output.m_utL3 = diagData.m_utL3;
    output.m_utR1 = diagData.m_utR1;
    output.m_utR2 = diagData.m_utR2;
    output.m_utR3 = diagData.m_utR3;
    output.U_dtL_00 = diagData.U_dtL_00;
    output.U_dtL_01 = diagData.U_dtL_01;
    output.U_dtL_02 = diagData.U_dtL_02;
    output.U_dtL_10 = diagData.U_dtL_10;
    output.U_dtL_11 = diagData.U_dtL_11;
    output.U_dtL_12 = diagData.U_dtL_12;
    output.U_dtL_20 = diagData.U_dtL_20;
    output.U_dtL_21 = diagData.U_dtL_21;
    output.U_dtL_22 = diagData.U_dtL_22;
    output.U_dtR_00 = diagData.U_dtR_00;
    output.U_dtR_01 = diagData.U_dtR_01;
    output.U_dtR_02 = diagData.U_dtR_02;
    output.U_dtR_10 = diagData.U_dtR_10;
    output.U_dtR_11 = diagData.U_dtR_11;
    output.U_dtR_12 = diagData.U_dtR_12;
    output.U_dtR_20 = diagData.U_dtR_20;
    output.U_dtR_21 = diagData.U_dtR_21;
    output.U_dtR_22 = diagData.U_dtR_22;
    output.U_utL_00 = diagData.U_utL_00;
    output.U_utL_01 = diagData.U_utL_01;
    output.U_utL_02 = diagData.U_utL_02;
    output.U_utL_10 = diagData.U_utL_10;
    output.U_utL_11 = diagData.U_utL_11;
    output.U_utL_12 = diagData.U_utL_12;
    output.U_utL_20 = diagData.U_utL_20;
    output.U_utL_21 = diagData.U_utL_21;
    output.U_utL_22 = diagData.U_utL_22;
    output.U_utR_00 = diagData.U_utR_00;
    output.U_utR_01 = diagData.U_utR_01;
    output.U_utR_02 = diagData.U_utR_02;
    output.U_utR_10 = diagData.U_utR_10;
    output.U_utR_11 = diagData.U_utR_11;
    output.U_utR_12 = diagData.U_utR_12;
    output.U_utR_20 = diagData.U_utR_20;
    output.U_utR_21 = diagData.U_utR_21;
    output.U_utR_22 = diagData.U_utR_22;
}

} // End of namespace brparity

#endif
