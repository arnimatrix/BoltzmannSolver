#include "global.h"
#include "libdiagonalization.h"
#include "brparity.h"
#include "libcomplexop.h"

namespace brparity {


void updateDiagonalization(param_t &params)
{
    SpectrumInput  inputs;
    readDiagonalizationInputs (inputs,  params);
    SpectrumOutput outputs = updateDiagonalization(inputs);
    readDiagonalizationOutputs(outputs, params);
}

SpectrumOutput updateDiagonalization(SpectrumInput const &inputs)
{
    auto const &mdL_00 = inputs.mdL_00;
    auto const &mdL_01 = inputs.mdL_01;
    auto const &mdL_02 = inputs.mdL_02;
    auto const &mdL_10 = inputs.mdL_10;
    auto const &mdL_11 = inputs.mdL_11;
    auto const &mdL_12 = inputs.mdL_12;
    auto const &mdL_20 = inputs.mdL_20;
    auto const &mdL_21 = inputs.mdL_21;
    auto const &mdL_22 = inputs.mdL_22;
    auto const &mdR_00 = inputs.mdR_00;
    auto const &mdR_01 = inputs.mdR_01;
    auto const &mdR_02 = inputs.mdR_02;
    auto const &mdR_10 = inputs.mdR_10;
    auto const &mdR_11 = inputs.mdR_11;
    auto const &mdR_12 = inputs.mdR_12;
    auto const &mdR_20 = inputs.mdR_20;
    auto const &mdR_21 = inputs.mdR_21;
    auto const &mdR_22 = inputs.mdR_22;
    auto const &muL_00 = inputs.muL_00;
    auto const &muL_01 = inputs.muL_01;
    auto const &muL_02 = inputs.muL_02;
    auto const &muL_10 = inputs.muL_10;
    auto const &muL_11 = inputs.muL_11;
    auto const &muL_12 = inputs.muL_12;
    auto const &muL_20 = inputs.muL_20;
    auto const &muL_21 = inputs.muL_21;
    auto const &muL_22 = inputs.muL_22;
    auto const &muR_00 = inputs.muR_00;
    auto const &muR_01 = inputs.muR_01;
    auto const &muR_02 = inputs.muR_02;
    auto const &muR_10 = inputs.muR_10;
    auto const &muR_11 = inputs.muR_11;
    auto const &muR_12 = inputs.muR_12;
    auto const &muR_20 = inputs.muR_20;
    auto const &muR_21 = inputs.muR_21;
    auto const &muR_22 = inputs.muR_22;

    SpectrumOutput outputs;

    Diagonalizer::applyDiagonalization(
        {
            -muL_00,
            -muL_01 + -muL_10,
            -muL_02 + -muL_20,
            0,
            -muL_11,
            -muL_12 + -muL_21,
            0,
            0,
            -muL_22,
        },
        {&outputs.U_utL_00, &outputs.U_utL_01, &outputs.U_utL_02, &outputs.U_utL_10, &outputs.U_utL_11, &outputs.U_utL_12, &outputs.U_utL_20, &outputs.U_utL_21, &outputs.U_utL_22, },
        {&outputs.m_utL1, &outputs.m_utL2, &outputs.m_utL3, }
        );

outputs.m_utL1 = std::sqrt(outputs.m_utL1);
outputs.m_utL2 = std::sqrt(outputs.m_utL2);
outputs.m_utL3 = std::sqrt(outputs.m_utL3);
    Diagonalizer::applyDiagonalization(
        {
            -muR_00,
            -muR_01 + -muR_10,
            -muR_02 + -muR_20,
            0,
            -muR_11,
            -muR_12 + -muR_21,
            0,
            0,
            -muR_22,
        },
        {&outputs.U_utR_00, &outputs.U_utR_01, &outputs.U_utR_02, &outputs.U_utR_10, &outputs.U_utR_11, &outputs.U_utR_12, &outputs.U_utR_20, &outputs.U_utR_21, &outputs.U_utR_22, },
        {&outputs.m_utR1, &outputs.m_utR2, &outputs.m_utR3, }
        );

outputs.m_utR1 = std::sqrt(outputs.m_utR1);
outputs.m_utR2 = std::sqrt(outputs.m_utR2);
outputs.m_utR3 = std::sqrt(outputs.m_utR3);
    Diagonalizer::applyDiagonalization(
        {
            -mdL_00,
            -mdL_01 + -mdL_10,
            -mdL_02 + -mdL_20,
            0,
            -mdL_11,
            -mdL_12 + -mdL_21,
            0,
            0,
            -mdL_22,
        },
        {&outputs.U_dtL_00, &outputs.U_dtL_01, &outputs.U_dtL_02, &outputs.U_dtL_10, &outputs.U_dtL_11, &outputs.U_dtL_12, &outputs.U_dtL_20, &outputs.U_dtL_21, &outputs.U_dtL_22, },
        {&outputs.m_dtL1, &outputs.m_dtL2, &outputs.m_dtL3, }
        );

outputs.m_dtL1 = std::sqrt(outputs.m_dtL1);
outputs.m_dtL2 = std::sqrt(outputs.m_dtL2);
outputs.m_dtL3 = std::sqrt(outputs.m_dtL3);
    Diagonalizer::applyDiagonalization(
        {
            -mdR_00,
            -mdR_01 + -mdR_10,
            -mdR_02 + -mdR_20,
            0,
            -mdR_11,
            -mdR_12 + -mdR_21,
            0,
            0,
            -mdR_22,
        },
        {&outputs.U_dtR_00, &outputs.U_dtR_01, &outputs.U_dtR_02, &outputs.U_dtR_10, &outputs.U_dtR_11, &outputs.U_dtR_12, &outputs.U_dtR_20, &outputs.U_dtR_21, &outputs.U_dtR_22, },
        {&outputs.m_dtR1, &outputs.m_dtR2, &outputs.m_dtR3, }
        );

outputs.m_dtR1 = std::sqrt(outputs.m_dtR1);
outputs.m_dtR2 = std::sqrt(outputs.m_dtR2);
outputs.m_dtR3 = std::sqrt(outputs.m_dtR3);
    Diagonalizer::applyDiagonalization(
        {
            -muL_00,
            -muL_01 + -muL_10,
            -muL_02 + -muL_20,
            0,
            -muL_11,
            -muL_12 + -muL_21,
            0,
            0,
            -muL_22,
        },
        {&outputs.U_utL_00, &outputs.U_utL_01, &outputs.U_utL_02, &outputs.U_utL_10, &outputs.U_utL_11, &outputs.U_utL_12, &outputs.U_utL_20, &outputs.U_utL_21, &outputs.U_utL_22, },
        {&outputs.m_utL1, &outputs.m_utL2, &outputs.m_utL3, }
        );

outputs.m_utL1 = std::sqrt(outputs.m_utL1);
outputs.m_utL2 = std::sqrt(outputs.m_utL2);
outputs.m_utL3 = std::sqrt(outputs.m_utL3);
    Diagonalizer::applyDiagonalization(
        {
            -muR_00,
            -muR_01 + -muR_10,
            -muR_02 + -muR_20,
            0,
            -muR_11,
            -muR_12 + -muR_21,
            0,
            0,
            -muR_22,
        },
        {&outputs.U_utR_00, &outputs.U_utR_01, &outputs.U_utR_02, &outputs.U_utR_10, &outputs.U_utR_11, &outputs.U_utR_12, &outputs.U_utR_20, &outputs.U_utR_21, &outputs.U_utR_22, },
        {&outputs.m_utR1, &outputs.m_utR2, &outputs.m_utR3, }
        );

outputs.m_utR1 = std::sqrt(outputs.m_utR1);
outputs.m_utR2 = std::sqrt(outputs.m_utR2);
outputs.m_utR3 = std::sqrt(outputs.m_utR3);
    Diagonalizer::applyDiagonalization(
        {
            -mdL_00,
            -mdL_01 + -mdL_10,
            -mdL_02 + -mdL_20,
            0,
            -mdL_11,
            -mdL_12 + -mdL_21,
            0,
            0,
            -mdL_22,
        },
        {&outputs.U_dtL_00, &outputs.U_dtL_01, &outputs.U_dtL_02, &outputs.U_dtL_10, &outputs.U_dtL_11, &outputs.U_dtL_12, &outputs.U_dtL_20, &outputs.U_dtL_21, &outputs.U_dtL_22, },
        {&outputs.m_dtL1, &outputs.m_dtL2, &outputs.m_dtL3, }
        );

outputs.m_dtL1 = std::sqrt(outputs.m_dtL1);
outputs.m_dtL2 = std::sqrt(outputs.m_dtL2);
outputs.m_dtL3 = std::sqrt(outputs.m_dtL3);
    Diagonalizer::applyDiagonalization(
        {
            -mdR_00,
            -mdR_01 + -mdR_10,
            -mdR_02 + -mdR_20,
            0,
            -mdR_11,
            -mdR_12 + -mdR_21,
            0,
            0,
            -mdR_22,
        },
        {&outputs.U_dtR_00, &outputs.U_dtR_01, &outputs.U_dtR_02, &outputs.U_dtR_10, &outputs.U_dtR_11, &outputs.U_dtR_12, &outputs.U_dtR_20, &outputs.U_dtR_21, &outputs.U_dtR_22, },
        {&outputs.m_dtR1, &outputs.m_dtR2, &outputs.m_dtR3, }
        );

outputs.m_dtR1 = std::sqrt(outputs.m_dtR1);
outputs.m_dtR2 = std::sqrt(outputs.m_dtR2);
outputs.m_dtR3 = std::sqrt(outputs.m_dtR3);
    return outputs;
}

} // End of namespace brparity

