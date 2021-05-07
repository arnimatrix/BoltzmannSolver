#include <fstream>
#include <iostream>
#include <istream>
#include <marty.h>
#include <stdlib.h>
#include <string>

using namespace std;
using namespace csl;
using namespace mty;

const std::string path_to_generated_library = ".";

Expr cc(Expr const& expr)
{
    return GetComplexConjugate(expr);
}

void calculateAndGenerateLib(std::string const&            initProcessName,
                             mty::Model&                   model,
                             std::vector<Insertion> const& insertions,
                             mty::Library&                 lib,
                             std::ofstream&                file,
                             bool                          show = false,
                             bool                          loop = false

)
{
    const auto ampl = (loop)
                        ? model.computeAmplitude(Order::OneLoop, insertions)
                        : model.computeAmplitude(Order::TreeLevel, insertions);

    if (show) {
        Display(ampl);
        if (!ampl.empty())
            Show(ampl);
    }
    const auto sampl            = model.computeSquaredAmplitude(ampl);
    const auto finalProcessName = FindProcessName(initProcessName, insertions);
    if (sampl != CSL_0) {
        lib.addFunction(finalProcessName, sampl);
        file << finalProcessName << std::endl;
    }
}

int main()
{

    // Model building

    Model toyModel;
    toyModel.addGaugedGroup(group::Type::SU, "C", 3);
    toyModel.addGaugedGroup(group::Type::U1, "Q");
    toyModel.addFlavorGroup("SM_flavor", 3, true);
    toyModel.init();

    toyModel.renameParticle("A_C", "G");
    toyModel.renameParticle("A_Q", "Q");

    Particle uLs = scalarboson_s("uLt", toyModel);
    uLs->setGroupRep("C", { 1, 0 });
    uLs->setGroupRep("Q", { 2, 3 });
    uLs->setFundamentalFlavorRep("SM_flavor");

    Particle dLs = scalarboson_s("dLt", toyModel);
    dLs->setGroupRep("C", { 1, 0 });
    dLs->setGroupRep("Q", { -1, 3 });
    dLs->setFundamentalFlavorRep("SM_flavor");

    Particle uRs = scalarboson_s("uRt", toyModel);
    uRs->setGroupRep("C", { 1, 0 });
    uRs->setGroupRep("Q", { 2, 3 });
    uRs->setFundamentalFlavorRep("SM_flavor");

    Particle dRs = scalarboson_s("dRt", toyModel);
    dRs->setGroupRep("C", { 1, 0 });
    dRs->setGroupRep("Q", { -1, 3 });
    dRs->setFundamentalFlavorRep("SM_flavor");

    Particle uL = weylfermion_s("u_L", toyModel, Chirality::Left);
    uL->setGroupRep("C", { 1, 0 });
    uL->setGroupRep("Q", { 2, 3 });
    uL->setFundamentalFlavorRep("SM_flavor");

    Particle dL = weylfermion_s("d_L", toyModel, Chirality::Left);
    dL->setGroupRep("C", { 1, 0 });
    dL->setGroupRep("Q", { -1, 3 });
    dL->setFundamentalFlavorRep("SM_flavor");

    Particle uR = weylfermion_s("u_R", toyModel, Chirality::Right);
    uR->setGroupRep("C", { 1, 0 });
    uR->setGroupRep("Q", { 2, 3 });
    uR->setFundamentalFlavorRep("SM_flavor");

    Particle dR = weylfermion_s("d_R", toyModel, Chirality::Right);
    dR->setGroupRep("C", { 1, 0 });
    dR->setGroupRep("Q", { -1, 3 });
    dR->setFundamentalFlavorRep("SM_flavor");

    Particle XX = diracfermion_s("XX", toyModel);
    XX->setSelfConjugate(true);
    XX->setGroupRep("Q", { 0 });

    toyModel.addParticle(uLs);
    toyModel.addParticle(uRs);
    toyModel.addParticle(dLs);
    toyModel.addParticle(dRs);
    toyModel.addParticle(uL);
    toyModel.addParticle(dL);
    toyModel.addParticle(uR);
    toyModel.addParticle(dR);
    toyModel.addParticle(XX);

    // SM flavor space
    Index I = toyModel.generateIndex("SM_flavor", "u_L");
    Index J = toyModel.generateIndex("SM_flavor", "u_L");
    Index K = toyModel.generateIndex("SM_flavor", "u_L");

    // SU(2)L space
    Index i = toyModel.generateIndex("C", "u_L");
    Index j = toyModel.generateIndex("C", "u_L");

    // SU(3)C space
    Index  a    = toyModel.generateIndex("C", "u_R");
    Index  b    = toyModel.generateIndex("C", "u_R");
    Index  c    = toyModel.generateIndex("C", "u_R");
    Tensor feps = a.getSpace()->getEpsilon();

    auto ii = toyModel.generateIndices(3, "C", "u_L");

    auto   flavorSpace = toyModel.getVectorSpace("SM_flavor", "u_L");
    Tensor Lpp = csl::Tensor("lpp", { flavorSpace, flavorSpace, flavorSpace });
    Lpp->setComplexProperty(csl::ComplexProperty::Complex);

    Tensor MuL = csl::Tensor("muL", { flavorSpace, flavorSpace });
    MuL->setComplexProperty(csl::ComplexProperty::Complex);

    Tensor MuR = csl::Tensor("muR", { flavorSpace, flavorSpace });
    MuR->setComplexProperty(csl::ComplexProperty::Complex);

    Tensor MdL = csl::Tensor("mdL", { flavorSpace, flavorSpace });
    MdL->setComplexProperty(csl::ComplexProperty::Complex);

    Tensor MdR = csl::Tensor("mdR", { flavorSpace, flavorSpace });
    MdR->setComplexProperty(csl::ComplexProperty::Complex);

    Expr   mX  = constant_s("mX");
    Tensor guL = csl::Tensor("guL", { flavorSpace, flavorSpace });
    Tensor guR = csl::Tensor("guR", { flavorSpace, flavorSpace });
    Tensor gdL = csl::Tensor("gdL", { flavorSpace, flavorSpace });
    Tensor gdR = csl::Tensor("gdR", { flavorSpace, flavorSpace });
    toyModel.addTensorCoupling(guL);
    toyModel.addTensorCoupling(guR);
    toyModel.addTensorCoupling(gdL);
    toyModel.addTensorCoupling(gdR);

    auto   al = DiracIndices(2);
    Tensor C  = DiracCMatrix();

    toyModel.addLagrangianTerm(
      -Lpp({ I, J, K }) * uRs({ I, a }) * feps({ a, b, c }) *
        dR({ J, b, al[0] }) * C({ al[0], al[1] }) * dR({ K, c, al[1] }),
      true // Add also the complex conjugate of this term
    );

    toyModel.addLagrangianTerm(
      -Lpp({ I, J, K }) * uR({ I, a, al[0] }) * feps({ a, b, c }) *
        C({ al[0], al[1] }) * dR({ J, b, al[1] }) * dRs({ K, c }),
      true // Add also the complex conjugate of this term
    );

    toyModel.addLagrangianTerm(
      -Lpp({ I, J, K }) * uR({ I, a, al[0] }) * C({ al[0], al[1] }) *
        feps({ a, b, c }) * dRs({ J, b }) * dR({ K, c, al[1] }),
      true // Add also the complex conjugate of this term
    );

    toyModel.addLagrangianTerm(
      -guL({ I, J }) * uLs({ I, i }) * cc(uL({ J, i, al[0] })) * XX(al[0]),
      true // Add also the complex conjugate of this term
    );

    toyModel.addLagrangianTerm(
      -guR({ I, J }) * uRs({ I, a }) * cc(uR({ J, a, al[0] })) * XX(al[0]),
      true // Add also the complex conjugate of this term
    );

    toyModel.addLagrangianTerm(
      -gdL({ I, J }) * dLs({ I, i }) * cc(dL({ J, i, al[0] })) * XX(al[0]),
      true // Add also the complex conjugate of this term
    );

    toyModel.addLagrangianTerm(
      -gdR({ I, J }) * dRs({ I, a }) * cc(dR({ J, a, al[0] })) * XX(al[0]),
      true // Add also the complex conjugate of this term
    );

    toyModel.addLagrangianTerm(-mX * cc(XX(al[0])) * XX(al[0]), true);

    toyModel.addLagrangianTerm(cc(uLs({ I, i })) * MuL({ I, J }) *
                               uLs({ J, i }));
    toyModel.addLagrangianTerm(cc(uRs({ I, a })) * MuR({ I, J }) *
                               uRs({ J, a }));
    toyModel.addLagrangianTerm(cc(dLs({ I, i })) * MdL({ I, J }) *
                               dLs({ J, i }));
    toyModel.addLagrangianTerm(cc(dRs({ I, a })) * MdR({ I, J }) *
                               dRs({ J, a }));

    toyModel.diagonalizeYukawa(
      "guL",                // diagonalizing Y_e
      { "gw", "gw", "gw" }, // Masses on the diagonal
      constant_s("gwuL")    // Global factor f in Y --> f*M
    );

    toyModel.diagonalizeYukawa(
      "gdL",                // diagonalizing Y_e
      { "gw", "gw", "gw" }, // Masses on the diagonal
      constant_s("gwdL")    // Global factor f in Y --> f*M
    );

    toyModel.diagonalizeYukawa(
      "guR",                // diagonalizing Y_e
      { "gw", "gw", "gw" }, // Masses on the diagonal
      constant_s("gwuR")    // Global factor f in Y --> f*M
    );

    toyModel.diagonalizeYukawa(
      "gdR",                // diagonalizing Y_e
      { "gw", "gw", "gw" }, // Masses on the diagonal
      constant_s("gwdR")    // Global factor f in Y --> f*M
    );

    csl::Expr m_u = constant_s("m_u");
    csl::Expr m_c = constant_s("m_c");
    csl::Expr m_t = constant_s("m_t");
    csl::Expr m_d = constant_s("m_d");
    csl::Expr m_s = constant_s("m_s");
    csl::Expr m_b = constant_s("m_b");

    cout << toyModel << endl;
    cout << endl << "=============================" << endl;
    // Model breaking

    toyModel.breakFlavorSymmetry("SM_flavor");
    toyModel.renameParticle("u_L_1", "uf1_L ; u_L");
    toyModel.renameParticle("u_L_2", "uf2_L ; c_L");
    toyModel.renameParticle("u_L_3", "uf3_L ; t_L");
    toyModel.renameParticle("u_R_1", "uf1_R ; u_R");
    toyModel.renameParticle("u_R_2", "uf2_R ; c_R");
    toyModel.renameParticle("u_R_3", "uf3_R ; t_R");
    toyModel.renameParticle("d_L_1", "df1_L ; d_L");
    toyModel.renameParticle("d_L_2", "df2_L ; s_L");
    toyModel.renameParticle("d_L_3", "df3_L ; b_L");
    toyModel.renameParticle("d_R_1", "df1_R ; d_R");
    toyModel.renameParticle("d_R_2", "df2_R ; s_R");
    toyModel.renameParticle("d_R_3", "df3_R ; b_R");

    toyModel.addFermionicMass("uf1_L", "uf1_R", m_u);
    toyModel.addFermionicMass("df1_L", "df1_R", m_d);
    toyModel.addFermionicMass("uf2_L", "uf2_R", m_c);
    toyModel.addFermionicMass("df2_L", "df2_R", m_s);
    toyModel.addFermionicMass("uf3_L", "uf3_R", m_t);
    toyModel.addFermionicMass("df3_L", "df3_R", m_b);

    cout << toyModel << endl;

    toyModel.renameParticle("uLt_1", "utL1 ; \\tilde{u}_L");
    toyModel.renameParticle("uLt_2", "utL2 ; \\tilde{c}_L");
    toyModel.renameParticle("uLt_3", "utL3 ; \\tilde{t}_L");
    toyModel.renameParticle("uRt_1", "utR1 ; \\tilde{u}_R");
    toyModel.renameParticle("uRt_2", "utR2 ; \\tilde{c}_R");
    toyModel.renameParticle("uRt_3", "utR3 ; \\tilde{t}_R");
    toyModel.renameParticle("dLt_1", "dtL1 ; \\tilde{d}_L");
    toyModel.renameParticle("dLt_2", "dtL2 ; \\tilde{s}_L");
    toyModel.renameParticle("dLt_3", "dtL3 ; \\tilde{b}_L");
    toyModel.renameParticle("dRt_1", "dtR1 ; \\tilde{u}_R");
    toyModel.renameParticle("dRt_2", "dtR2 ; \\tilde{c}_R");
    toyModel.renameParticle("dRt_3", "dtR3 ; \\tilde{t}_R");

    toyModel.rotateFields({ "utL1", "utL2", "utL3" }, true);
    Expr GuLt = constant_s("GuLt");
    toyModel.getParticle("utL1")->setWidth(GuLt);
    Expr GcLt = constant_s("GcLt");
    toyModel.getParticle("utL2")->setWidth(GcLt);
    Expr GtLt = constant_s("GtLt");
    toyModel.getParticle("utL3")->setWidth(GtLt);

    toyModel.rotateFields({ "utR1", "utR2", "utR3" }, true);
    Expr GuRt = constant_s("GuRt");
    toyModel.getParticle("utR1")->setWidth(GuRt);
    Expr GcRt = constant_s("GcRt");
    toyModel.getParticle("utR2")->setWidth(GcRt);
    Expr GtRt = constant_s("GtRt");
    toyModel.getParticle("utR3")->setWidth(GtRt);

    toyModel.rotateFields({ "dtL1", "dtL2", "dtL3" }, true);
    Expr GdLt = constant_s("GdLt");
    toyModel.getParticle("dtL1")->setWidth(GdLt);
    Expr GsLt = constant_s("GsLt");
    toyModel.getParticle("dtL2")->setWidth(GsLt);
    Expr GbLt = constant_s("GbLt");
    toyModel.getParticle("dtL3")->setWidth(GbLt);

    toyModel.rotateFields({ "dtR1", "dtR2", "dtR3" }, true);
    Expr GdRt = constant_s("GdRt");
    toyModel.getParticle("dtR1")->setWidth(GdRt);
    Expr GsRt = constant_s("GsRt");
    toyModel.getParticle("dtR2")->setWidth(GsRt);
    Expr GbRt = constant_s("GbRt");
    toyModel.getParticle("dtR3")->setWidth(GbRt);

    toyModel.refresh();
    std::cout << toyModel << std::endl;

    auto rules = toyModel.getFeynmanRules();
    Display(rules); // Displays expressions in terminal
    // Show(rules); // Shows diagrams in the application

    mty::Library myLib("BRparity", path_to_generated_library);

    std::vector<Particle> usR = toyModel.getPhysicalParticles([&](Particle p) {
        return p->getSpinDimension() == 1 // scalar
               && p->getName()[0] == 'u' &&
               p->getName()[2] == 'R'; // name starts with 'u'
    });
    std::vector<Particle> dsR = toyModel.getPhysicalParticles([&](Particle p) {
        return p->getSpinDimension() == 1 // scalar
               && p->getName()[0] == 'd' &&
               p->getName()[2] == 'R'; // name starts with 'd'
    });

    std::vector<Particle> uq = toyModel.getPhysicalParticles([&](Particle p) {
        return p->getSpinDimension() == 2 // fermion
               && p->getName()[0] == 'u'; // name starts with 'u'
    });
    std::vector<Particle> dq = toyModel.getPhysicalParticles([&](Particle p) {
        return p->getSpinDimension() == 2 // fermion
               && p->getName()[0] == 'd'; // name starts with 'd'
    });
    /**/
    ofstream DecayB("BRparity/script/DecayB.txt");
    ofstream DecayAB("BRparity/script/DecayAB.txt");
    ofstream ScattB("BRparity/script/ScattB.txt");
    ofstream ScattAB("BRparity/script/ScattAB.txt");
    ofstream Scatt("BRparity/script/Scatt.txt");
    ofstream Decay("BRparity/script/Decay.txt");

    for (size_t i = 0; i != usR.size(); ++i) {
        Particle usRi = usR[i];
        for (size_t j = 0; j != dq.size(); ++j) {
            Particle dfRj = dq[j];
            for (size_t k = 0; k != dq.size(); ++k) {
                Particle dfRk      = dq[k];
                auto     insertion = { Incoming(usRi),
                                   Outgoing(AntiPart(dfRj)),
                                   Outgoing(AntiPart(dfRk)) };
                calculateAndGenerateLib(
                  "CP", toyModel, insertion, myLib, Decay);
            }
        }
        for (size_t j = 0; j != uq.size(); ++j) {
            Particle ufRj  = uq[j];
            auto insertion = { Incoming(usRi), Outgoing(ufRj), Outgoing("XX") };
            auto insertion1 = { Incoming(usRi),
                                Outgoing(ufRj),
                                Outgoing(AntiPart("XX")) };
            calculateAndGenerateLib("CP", toyModel, insertion, myLib, Decay);
            calculateAndGenerateLib("CP", toyModel, insertion1, myLib, Decay);
        }
    }

    for (size_t i = 0; i != dsR.size(); ++i) {
        Particle dsRi = dsR[i];
        for (size_t j = 0; j != uq.size(); ++j) {
            Particle ufRj = uq[j];
            for (size_t k = 0; k != dq.size(); ++k) {
                Particle dfRk      = dq[k];
                auto     insertion = { Incoming(dsRi),
                                   Outgoing(AntiPart(ufRj)),
                                   Outgoing(AntiPart(dfRk)) };
                calculateAndGenerateLib(
                  "CP", toyModel, insertion, myLib, Decay);
            }
        }
        for (size_t j = 0; j != dq.size(); ++j) {
            Particle dfRj  = dq[j];
            auto insertion = { Incoming(dsRi), Outgoing(dfRj), Outgoing("XX") };
            auto insertion1 = { Incoming(dsRi),
                                Outgoing(dfRj),
                                Outgoing(AntiPart("XX")) };
            calculateAndGenerateLib("CP", toyModel, insertion, myLib, Decay);
            calculateAndGenerateLib("CP", toyModel, insertion1, myLib, Decay);
        }
    }

    for (size_t i = 0; i != uq.size(); ++i) {
        Particle ufRi = uq[i];
        for (size_t j = 0; j != dq.size(); ++j) {
            Particle dfRj = dq[j];
            for (size_t k = 0; k != dq.size(); ++k) {
                Particle dfRk       = dq[k];
                auto     insertion  = { Incoming(AntiPart("XX")),
                                   Outgoing(ufRi),
                                   Outgoing(dfRj),
                                   Outgoing(dfRk) };
                auto     insertion1 = { Incoming("XX"),
                                    Outgoing(ufRi),
                                    Outgoing(dfRj),
                                    Outgoing(dfRk) };
                calculateAndGenerateLib(
                  "CP", toyModel, insertion, myLib, DecayB);
                calculateAndGenerateLib(
                  "CP", toyModel, insertion1, myLib, DecayB);
                calculateAndGenerateLib(
                  "CP", toyModel, AntiPart(insertion), myLib, DecayAB);
                calculateAndGenerateLib(
                  "CP", toyModel, AntiPart(insertion1), myLib, DecayAB);
            }
        }
    }

    for (size_t i = 0; i != uq.size(); ++i) {
        Particle uLi        = uq[i];
        auto     insertion  = { Incoming("XX"),
                           Incoming(AntiPart("XX")),
                           Outgoing(uLi),
                           Outgoing(AntiPart(uLi)) };
        auto     insertion1 = { Incoming("XX"),
                            Incoming("XX"),
                            Outgoing(uLi),
                            Outgoing(AntiPart(uLi)) };
        auto     insertion2 = { Incoming(AntiPart("XX")),
                            Incoming(AntiPart("XX")),
                            Outgoing(uLi),
                            Outgoing(AntiPart(uLi)) };
        calculateAndGenerateLib("CP", toyModel, insertion, myLib, Scatt);
        calculateAndGenerateLib("CP", toyModel, insertion1, myLib, Scatt);
        calculateAndGenerateLib("CP", toyModel, insertion2, myLib, Scatt);
    }

    for (size_t i = 0; i != dq.size(); ++i) {
        Particle dLi        = dq[i];
        auto     insertion  = { Incoming("XX"),
                           Incoming(AntiPart("XX")),
                           Outgoing(dLi),
                           Outgoing(AntiPart(dLi)) };
        auto     insertion1 = { Incoming("XX"),
                            Incoming("XX"),
                            Outgoing(dLi),
                            Outgoing(AntiPart(dLi)) };
        auto     insertion2 = { Incoming(AntiPart("XX")),
                            Incoming(AntiPart("XX")),
                            Outgoing(dLi),
                            Outgoing(AntiPart(dLi)) };
        calculateAndGenerateLib("CP", toyModel, insertion, myLib, Scatt);
        calculateAndGenerateLib("CP", toyModel, insertion1, myLib, Scatt);
        calculateAndGenerateLib("CP", toyModel, insertion2, myLib, Scatt);
    }

    for (size_t i = 0; i != uq.size(); ++i) {
        Particle ufRi = uq[i];
        for (size_t j = 0; j != dq.size(); ++j) {
            Particle dfRj = dq[j];
            for (size_t k = 0; k != dq.size(); ++k) {
                Particle dfRk       = dq[k];
                auto     insertion  = { Incoming("XX"),
                                   Incoming(AntiPart(ufRi)),
                                   Outgoing(dfRj),
                                   Outgoing(dfRk) };
                auto     insertion1 = { Incoming(AntiPart("XX")),
                                    Incoming(AntiPart(ufRi)),
                                    Outgoing(dfRj),
                                    Outgoing(dfRk) };
                calculateAndGenerateLib(
                  "CP", toyModel, insertion, myLib, ScattB);
                calculateAndGenerateLib(
                  "CP", toyModel, insertion1, myLib, ScattB);
                calculateAndGenerateLib(
                  "CP", toyModel, AntiPart(insertion), myLib, ScattAB);
                calculateAndGenerateLib(
                  "CP", toyModel, AntiPart(insertion1), myLib, ScattAB);
            }
        }
    }

    for (size_t i = 0; i != uq.size(); ++i) {
        Particle ufRi = uq[i];
        for (size_t j = 0; j != dq.size(); ++j) {
            Particle dfRj = dq[j];
            for (size_t k = 0; k != dq.size(); ++k) {
                Particle dfRk       = dq[k];
                auto     insertion  = { Incoming("XX"),
                                   Incoming(AntiPart(dfRj)),
                                   Outgoing(ufRi),
                                   Outgoing(dfRk) };
                auto     insertion1 = { Incoming(AntiPart("XX")),
                                    Incoming(AntiPart(dfRj)),
                                    Outgoing(ufRi),
                                    Outgoing(dfRk) };
                calculateAndGenerateLib(
                  "CP", toyModel, insertion, myLib, DecayB);
                calculateAndGenerateLib(
                  "CP", toyModel, insertion1, myLib, DecayB);
                calculateAndGenerateLib(
                  "CP", toyModel, AntiPart(insertion), myLib, DecayAB);
                calculateAndGenerateLib(
                  "CP", toyModel, AntiPart(insertion1), myLib, DecayAB);
            }
        }
    }

    Decay.close();
    DecayB.close();
    DecayAB.close();
    Scatt.close();
    ScattB.close();
    ScattAB.close();

    myLib.generateSpectrum(toyModel);
    // myLib.importLHAModule(".");
    // myLib.addIPath("/usr/local/Cellar/boost/1.74.0/include");
    // myLib.addIPath("/usr/local/Cellar/eigen/3.3.8_1/include");
    // myLib.addLPath("/usr/local/Cellar/boost/1.74.0/lib");
    myLib.addLPath("/usr/local/lib");
    myLib.addIPath("/usr/local/include");

    myLib.generateSpectrum(toyModel);

    // Cleaning existing library
    myLib.cleanExistingSources();

    // Exporting the additional cpp and h files to the library
    [[maybe_unused]] int res = system("./export_to_lib.sh");

    myLib.build(4);

    return 0;
}
