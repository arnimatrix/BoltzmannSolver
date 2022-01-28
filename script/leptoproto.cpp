#include <fstream>
#include <iostream>
#include <istream>
#include <marty.h>
#include <stdlib.h>
#include <string>

#include "processdata.h"
#include "quantumnumberdata.h"
#include "libjsondata.h"

using namespace std;
using namespace csl;
using namespace mty;

// To change by the user !
std::string path_to_generated_library = ".";

Expr cc(Expr const &expr) { return GetComplexConjugate(expr); }

void calculateAndGenerateLib(
                             std::string            const &initProcessName,
                             mty::Model                        &model,
                             std::vector<Insertion> const &insertions,
                             mty::Library                      &lib,
                             QuantumNumberData      const &qData,
                             ProcessData                  &pData,
                             bool                          show = false,
                             bool                          asym = false
                             
                             )
{
    
    if(asym)
    {
        FeynOptions options;
        options.addFilters(
                [&](FeynmanDiagram const &diag) {
                    return !diag.isInLoop("lL_1")
                        && !diag.isInLoop("lL_2")
                    && !diag.isInLoop("lL_3")
                    && !diag.isInLoop("W")
                    && !diag.isInLoop("B");
                }
            );
        
        //mty::option::decomposeInOperators = true;
        //mty::option::decomposeInLocalOperator = true;
        
        
        const auto ampl1  = model.computeAmplitude(Order::OneLoop, insertions, options);
        
        const auto ampl0  = model.computeAmplitude(Order::TreeLevel, insertions);
        
       if (show) {
            Display(ampl0);
            if (!ampl0.empty())
                Show(ampl0);
            Display(ampl1);
            if (!ampl1.empty())
                Show(ampl1);
        }
        const auto sampl1 = model.computeSquaredAmplitude(ampl0,ampl1);
        
        const auto finalProcessName1 = FindProcessName(initProcessName + "L", insertions);
        if (sampl1 != CSL_0) {
            pData.addProcess(finalProcessName1, insertions, qData);
            lib.addFunction(finalProcessName1, sampl1);
        }
        
        const auto sampl0 = model.computeSquaredAmplitude(ampl0);
        
        const auto finalProcessName0 = FindProcessName(initProcessName + "T", insertions);
        if (sampl0 != CSL_0) {
            pData.addProcess(finalProcessName0, insertions, qData);
            lib.addFunction(finalProcessName0, sampl0);
        }
        
    }
    else
    {
        const auto ampl  = model.computeAmplitude(Order::TreeLevel, insertions);
        
        if (show) {
            Display(ampl);
            if (!ampl.empty())
                Show(ampl);
        }
        
        const auto sampl = model.computeSquaredAmplitude(ampl);
        const auto finalProcessName = FindProcessName(initProcessName, insertions);
        if (sampl != CSL_0) {
            pData.addProcess(finalProcessName, insertions, qData);
            lib.addFunction(finalProcessName, sampl);
            //file << finalProcessName << std::endl;
        }
        
    }
    
    
}

int main() {
    
    // Model building
    
    // Dirac space
    Index a = DiracIndex();
    Index b = DiracIndex();
    
    
    
    Model toyModel;
    toyModel.addGaugedGroup(group::Type::SU, "L", 2);
    toyModel.addGaugedGroup(group::Type::U1, "Y");
    toyModel.addFlavorGroup("SM_flavor", 3,true);
    toyModel.init();
    
    toyModel.renameParticle("A_L", "W");
    toyModel.renameParticle("A_Y", "B");
    
    
    Particle et = scalarboson_s("et", toyModel);
    et->setGroupRep("L", {1});
    et->setGroupRep("Y", {-1, 2});
    toyModel.addParticle(et);
    
    /*
     
     Particle H = scalarboson_s("H", toyModel);
     H->setGroupRep("L", {1});
     H->setGroupRep("Y", {-1, 2});
     toyModel.addParticle(H);
     
    Particle Li = weylfermion_s("L_L", toyModel, Chirality::Left);
    Li->setGroupRep("L", {1});
    Li->setGroupRep("Y", {-1, 2});
    Li->setFundamentalFlavorRep("SM_flavor");
    toyModel.addParticle(Li);
    
    Particle Ei = weylfermion_s("E_R", toyModel, Chirality::Right);
    Ei->setGroupRep("Y", {-1, 1});
    Ei->setFundamentalFlavorRep("SM_flavor");
    toyModel.addParticle(Ei);
    */
    
    Particle Li = diracfermion_s("L_L", toyModel);
    Li->setGroupRep("L", {1});
    Li->setGroupRep("Y", {-1, 2});
    Li->setFundamentalFlavorRep("SM_flavor");
    toyModel.addParticle(Li);
    
    Particle Ni = weylfermion_s("N", toyModel, Chirality::Right);
    Ni->setGroupRep("Y", {0});
    Ni->setFundamentalFlavorRep("SM_flavor");
    toyModel.addParticle(Ni);
    
    Particle SS = scalarboson_s("SS", toyModel);
    SS->setGroupRep("Y", {0});
    SS->setSelfConjugate(true);
    toyModel.addParticle(SS);
    
    
    // SM flavor space
    Index I = toyModel.generateIndex("SM_flavor", "L_L");
    Index J = toyModel.generateIndex("SM_flavor", "L_L");
    
    // SU(2)L space
    Index i = toyModel.generateIndex("L", "L_L");
    Index j = toyModel.generateIndex("L", "L_L");
    Tensor eps = i.getSpace()->getEpsilon();
    
    //Tensor pauli = toyModel.getGenerator("L", "L_L");
    Tensor pauli = GetGenerator(toyModel,"L", "L_L");
    auto ii = toyModel.generateIndices(3, "L", "L_L");
    
    auto flavorSpace = toyModel.getVectorSpace("SM_flavor", "L_L");
    Tensor Ye = csl::Tensor("Y_l", {flavorSpace, flavorSpace});
    Ye->setComplexProperty(csl::ComplexProperty::Complex);
    
    Tensor Yet = csl::Tensor("Y_et", {flavorSpace, flavorSpace});
    Yet->setComplexProperty(csl::ComplexProperty::Complex);
    
    Tensor Ynu = csl::Tensor("Y_nu", {flavorSpace, flavorSpace});
    Ynu->setComplexProperty(csl::ComplexProperty::Complex);
    
    Tensor YS = csl::Tensor("Y_S", {flavorSpace, flavorSpace});
    YS->setComplexProperty(csl::ComplexProperty::Complex);
    
    Tensor MN = csl::Tensor("M_N", {flavorSpace, flavorSpace});
    MN->setComplexProperty(csl::ComplexProperty::Real);
    
    Expr muet2 = constant_s("muet2");
    Expr muS2 = constant_s("muS2");
    Expr muSet = constant_s("muSet");
    muSet->setComplexProperty(csl::ComplexProperty::Complex);
    Expr let = constant_s("lamet");
    Expr lS = constant_s("lS");
    Expr lSet = constant_s("lSet");
    //Expr vh = constant_s("vh");
    //Expr vS = constant_s("vS");
    
    auto al = DiracIndices(3);
    Tensor C = dirac4.C_matrix;
    Tensor PL = dirac4.P_L;
    Tensor PR = dirac4.P_R;
    
    
    toyModel.addLagrangianTerm(-Ni({I, al[0]})*C({al[0],al[1]})*PL({al[1],al[2]})*Yet({I, J})*cc(et(i))*Li({J, i, al[2]}),
                               true // Add also the complex conjugate of this term
                               );
    
    toyModel.addLagrangianTerm(-Ni({I, al[0]})*C({al[0],al[1]})*MN({I, J})*Ni({J, al[1]}),true // Add also the complex conjugate of this term
                               );
    
    toyModel.addLagrangianTerm(-SS*Ni({I, al[0]})*C({al[0],al[1]})*YS({I, J})*Ni({J, al[1]}));
    
    toyModel.addLagrangianTerm(muet2*cc(et(i))*et(i));
    
    toyModel.addLagrangianTerm(muS2*SS*SS);
    
    toyModel.addLagrangianTerm(-lS*SS*SS*SS*SS);
    
    Expr prod1 = cc(et(i))*et(i);
    toyModel.addLagrangianTerm(let*prod1*RenamedIndices(prod1));
    
    toyModel.addLagrangianTerm(lSet*SS*SS*cc(et(i))*et(i));
    //toyModel.addLagrangianTerm(lSH*cc(SS)*SS*cc(H(i))*H(i));
    
    
    toyModel.addLagrangianTerm(muSet*SS*cc(et(i))*et(i));
    
    cout << toyModel << endl;
    
    // Model breaking
    
    //Expr m_le = constant_s("m_le");
    //Expr m_lmu = constant_s("m_lmu");
    //Expr m_ltau = constant_s("m_ltau");
    
    Expr m_N_1 = constant_s("m_N_1");
    Expr m_N_2 = constant_s("m_N_2");
    Expr m_N_3 = constant_s("m_N_3");
    
    //csl::Tensor M_l = csl::tensor_s("M_l",{flavorSpace, flavorSpace},csl::matrix_s({{m_le, 0, 0},{0, m_lmu, 0},{0, 0, m_ltau}}));
    csl::Tensor M_N = csl::tensor_s("M_N",{flavorSpace, flavorSpace},csl::matrix_s({{m_N_1, 0, 0},{0,  m_N_2, 0},{0, 0,   m_N_3}}));
    
    Index f_i = GetIndex(flavorSpace);
    Index f_j = GetIndex(flavorSpace);
    //toyModel.replace(Ye,sqrt_s(2)/(vh) * M_l({f_i, f_j}));
    toyModel.replace(MN,M_N({f_i, f_j}));
    //toyModel.replace(YS,sqrt_s(2)*M_N({f_i, f_j})/vS);
    
    
    toyModel.breakFlavorSymmetry("SM_flavor");
    toyModel.renameParticle("L_L_1","lL_1; e_L");
    toyModel.renameParticle("L_L_2","lL_2 ; \\mu_L");
    toyModel.renameParticle("L_L_3","lL_3 ; \\tau_L");
    /*
    
    toyModel.renameParticle("E_R_1","le_R ; e_R");
    toyModel.renameParticle("E_R_2","lmu_R ; \\mu_R");
    toyModel.renameParticle("E_R_3","ltau_R ; \\tau_R");
     toyModel.breakGaugeSymmetry("Y");
     
     toyModel.breakGaugeSymmetry("L");
     
     toyModel.renameParticle("L_L_1", "Nu_L");
     toyModel.renameParticle("L_L_2", "E_L");
     
     toyModel.renameParticle("H_1", "H0 ; H_0");
     toyModel.renameParticle("H_2", "Hm ; \\H^-");
     
     toyModel.renameParticle("et_1", "et0 ; \\eta_0");
     toyModel.renameParticle("et_2", "etm ; \\eta^-");
     
     toyModel.breakFlavorSymmetry("SM_flavor");
     toyModel.renameParticle("E_L_1","le_L");
     toyModel.renameParticle("E_L_2","lmu_L ; \\mu_L");
     toyModel.renameParticle("E_L_3","ltau_L ; \\tau_L");
     
     toyModel.renameParticle("E_R_1","le_R");
     toyModel.renameParticle("E_R_2","lmu_R ; \\mu_R");
     toyModel.renameParticle("E_R_3","ltau_R ; \\tau_R");
     
     // Higgs mechanism simulated
     Particle W_1 = toyModel.getParticle("W_1");
     Particle W_2 = toyModel.getParticle("W_2");
     Particle F_W_1 = W_1->getFieldStrength();
     Particle F_W_2 = W_2->getFieldStrength();
     Particle W = W_1->generateSimilar("W");
     W->setSelfConjugate(false);
     Index mu = MinkowskiIndex();
     Index nu = MinkowskiIndex();
     // W_1 goes to (W+ + W-) / sqrt(2)
     toyModel.replace(W_1, (W(mu) + GetComplexConjugate(W(mu))) / sqrt_s(2));
     toyModel.replace(F_W_1, (W({mu, nu}) + GetComplexConjugate(W({mu, nu}))) / sqrt_s(2));
     // W_2 goes to i*(W+ - W-) / sqrt(2)
     toyModel.replace(W_2, CSL_I * (W({mu}) - GetComplexConjugate(W({mu}))) / sqrt_s(2));
     toyModel.replace(F_W_2, CSL_I * (W({mu, nu}) - GetComplexConjugate(W({mu, nu}))) / sqrt_s(2));
     
     Particle S0 = toyModel.getParticle("SS");
     Particle SR = S0->generateSimilar("S_R");
     SR->setSelfConjugate(true);
     
     Particle SI = S0->generateSimilar("S_I");
     SI->setSelfConjugate(true);
     toyModel.replace(S0, (vS + SR + CSL_I*SI) / sqrt_s(2));
     
     //toyModel.replace(S0, (vS + SR) / sqrt_s(2));
     
     Particle H0 = toyModel.getParticle("H0");
     Particle hh = H0->generateSimilar("hh");
     hh->setSelfConjugate(true);
     Particle Ah = H0->generateSimilar("Ah");
     Ah->setSelfConjugate(true);
     toyModel.replace(H0, (vh + hh + CSL_I*Ah) / sqrt_s(2));
     
     Particle et0 = toyModel.getParticle("et0");
     
     Particle etR = et0->generateSimilar("eta_R ; \\eta^0_R");
     etR->setSelfConjugate(true);
     Particle etI = et0->generateSimilar("eta_I ; \\eta^0_I");
     etI->setSelfConjugate(true);
     
     toyModel.replace(et0, (etR + CSL_I*etI)/ sqrt_s(2));
     
     //toyModel.rotateFields({"eta_R","eta_I"},true);
     
     toyModel.rotateFields({"hh","S_R"},true);
     //Expr GSR = constant_s("GSR");
     //toyModel.getParticle("S_R")->setWidth(GSR);
     
     toyModel.rotateFields({"Ah","S_I"},true);
     //Expr GSI = constant_s("GSI");
     //toyModel.getParticle("S_I")->setWidth(GSI);
     Rename(toyModel,"Ah","G0");
     Rename(toyModel,"Hm","Gm");
     
     toyModel.getParticle("Nu_L_1")->setSelfConjugate(true);
     toyModel.getParticle("Nu_L_2")->setSelfConjugate(true);
     toyModel.getParticle("Nu_L_3")->setSelfConjugate(true);
     toyModel.promoteToMajorana("Nu_L_1");
     toyModel.promoteToMajorana("Nu_L_2");
     toyModel.promoteToMajorana("Nu_L_3");
     
     toyModel.renameParticle("Nu_L_1", "nu_e ; \\nu_{e L}");
     toyModel.renameParticle("Nu_L_2", "nu_mu ; \\nu _{\\mu L}");
     toyModel.renameParticle("Nu_L_3", "nu_tau ; \\nu _{\\tau L}");
     toyModel.rotateFields({"nu_e", "nu_mu", "nu_tau"}, true);
     */
    
    
    toyModel.getParticle("N_1")->setSelfConjugate(true);
    toyModel.getParticle("N_2")->setSelfConjugate(true);
    toyModel.getParticle("N_3")->setSelfConjugate(true);
    toyModel.promoteToMajorana("N_1");
    toyModel.promoteToMajorana("N_2");
    toyModel.promoteToMajorana("N_3");
    //toyModel.rotateFields({"N_1","N_2","N_3"},true);
    
    /*
     Expr GN1 = constant_s("GN1");
     toyModel.getParticle("N_1")->setWidth(GN1);
     Expr GN2 = constant_s("GN2");
     toyModel.getParticle("N_2")->setWidth(GN2);
     Expr GN3 = constant_s("GN3");
     toyModel.getParticle("N_3")->setWidth(GN3);
     */
    
    Expr thetaW = constant_s("thetaW");
    Expr cW = cos_s(thetaW);
    Expr sW = sin_s(thetaW);
    // We give the rotation matrix explicitly here, between double curly braces
    //toyModel.rotateFields({"W_3", "B"},{"Z"  , "A"},{{cW , sW},{-sW, cW}});
    
    //toyModel.diagonalizeSymbolically("B");
    //toyModel.renameParticle("B","A");
    //toyModel.renameParticle("W_3","Z");
    //Particle AA = toyModel.getParticle("A");
    //Particle Z = toyModel.getParticle("Z");
    
    //toyModel.promoteToGoldstone("Gm", "W");
    //toyModel.promoteToGoldstone("G0", "Z");
    
    Expr e = constant_s("e");
    Expr gY = toyModel.getScalarCoupling("g_Y");
    Expr gL = toyModel.getScalarCoupling("g_L");
    toyModel.replace(gY, e / cW);
    toyModel.replace(gL, e / sW);
    
    //Expr M_W = constant_s("M_W");
    //Expr M_Z = constant_s("M_Z");
    
    toyModel.refresh();
    std::cout << toyModel << std::endl;
    
    ///////////////////////////////////////////////////
    // Quantum number definitions, to complete
    ///////////////////////////////////////////////////
    
    QuantumNumberData qData(toyModel);
    ProcessData       pData(toyModel);
    qData.addQuantumNumber(
                           "L",
                           {"N_1", "N_2", "N_3", "et","SS","lL_1","lL_2","lL_3","W","B"},
                           { 0   ,  0    , 0    , 0  , 0  , 1    , 1    , 1    , 0 ,0  }
                           );
    qData.addQuantumNumber(
                           "Q",
                           {"N_1", "N_2", "N_3", "et","SS","lL_1","lL_2","lL_3","W","B"},
                           { 0   ,  0    , 0    , -1  , 0  , -1    , -1    , -1 , 0 ,0  }
                           );
    
    ///////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////
    
    
    auto rules = toyModel.getFeynmanRules();
    Display(rules); // Displays expressions in terminal
    //Show(rules); // Shows diagrams in the application
    //--------------------------
    
    mty::Library myLib("LeptoProto", path_to_generated_library);
    mty::option::decomposeInOperators = true;
    mty::option::decomposeInLocalOperator = false;
    
    //calculateAndGenerateLib("CP", toyModel, {Incoming("etm"),Incoming("etm"), Outgoing("e"), Outgoing("e")}, myLib,false,false);
    //calculateAndGenerateLib("CP", toyModel, {Incoming("etm"),AntiPart(Incoming("etm")), AntiPart(Outgoing("W")), Outgoing("W")}, myLib,true,false);
    
    
    auto XX = toyModel.getPhysicalParticles([&](Particle p) {
        const auto name = p->getName();
        return (name.find("N_") != std::string::npos || name == "et");
    });
    
    auto leptons = toyModel.getPhysicalParticles([&](Particle p) {
        const auto name = p->getName();
        return (name.find("lL_") != std::string::npos);
    });
    /*
     auto leptons = toyModel.getPhysicalParticles([&](Particle p) {
     const auto name = p->getName();
     return (name == "le" || name == "lmu" || name == "ltau");
     });
     
     auto neutrinos = toyModel.getPhysicalParticles([&](Particle p) {
     const auto name = p->getName();
     return (name.find("nu_") != std::string::npos);
     });
     */
    auto vectors = toyModel.getPhysicalParticles([&](Particle p) {
        return p->getSpinDimension() == 3
        && (p->getGaugedGroup() == toyModel.getGaugedGroup("L") || p->getGaugedGroup() == toyModel.getGaugedGroup("Y"));
    });
    
    auto smbath = toyModel.getPhysicalParticles([&](Particle p) {
        const auto name = p->getName();
        return ((name.find("lL_") != std::string::npos || name == "SS")||(p->getSpinDimension() == 3
        && (p->getGaugedGroup() == toyModel.getGaugedGroup("L") || p->getGaugedGroup() == toyModel.getGaugedGroup("Y"))));
    });
    /*
     for( const auto &dd : XX) {
     for( const auto &vv : vectors) {
     auto insertions = {
     Incoming(dd),
     AntiPart(Incoming(dd)),
     Outgoing(vv),
     Outgoing(AntiPart(vv))
     };
     calculateAndGenerateLib("Lepto", toyModel, insertions, myLib,qData,pData);
     }
     }
     */
    
    mty::option::excludeTadpoles = true;
    mty::option::excludeMassCorrections = true;
    mty::option::excludeTriangles = false;
    mty::option::excludeBoxes = true;
    
    for( const auto &dd : XX) {
        for( const auto &dd1 : XX) {
            for( const auto &sb : smbath) {
                for( const auto &sb1 : smbath) {
                    auto ins0 = {
                        Incoming(dd),
                        Incoming(dd1),
                        Outgoing(sb),
                        Outgoing(sb1)
                    };
                    calculateAndGenerateLib("Tree", toyModel, ins0, myLib,qData,pData);
                    
                    //calculateAndGenerateLib("Asym", toyModel,ins0, myLib,qData,pData,false,true);
                    //calculateAndGenerateLib("Asym", toyModel,AntiPart(ins0), myLib,qData,pData,false,true);
                    /**/
                    auto ins1 = {
                        Incoming(dd),
                        AntiPart(Incoming(sb)),
                        AntiPart(Outgoing(dd1)),
                        Outgoing(sb1)
                    };
                    calculateAndGenerateLib("Tree", toyModel, ins1, myLib,qData,pData);
                    
                    //calculateAndGenerateLib("Asym", toyModel,ins1, myLib,qData,pData,false,true);
                    //calculateAndGenerateLib("Asym", toyModel,AntiPart(ins1), myLib,qData,pData,false,true);
                    
                    auto ins2 = {
                        Incoming(dd),
                        AntiPart(Incoming(dd1)),
                        Outgoing(sb),
                        AntiPart(Outgoing(sb1))
                    };
                    
                    calculateAndGenerateLib("Tree", toyModel, ins2, myLib,qData,pData);
                    
                    calculateAndGenerateLib("Asym", toyModel,ins2, myLib,qData,pData,false,true);
                    calculateAndGenerateLib("Asym", toyModel,AntiPart(ins2), myLib,qData,pData,false,true);
                    
                }
                auto ins3 = {
                    Incoming(dd),
                    Outgoing(sb),
                    AntiPart(Outgoing(dd1))
                };
                calculateAndGenerateLib("Tree", toyModel, ins3, myLib,qData,pData);
                
                calculateAndGenerateLib("Asym", toyModel,ins3, myLib,qData,pData,false,true);
                calculateAndGenerateLib("Asym", toyModel,AntiPart(ins3), myLib,qData,pData,false,true);
                
            }
            /*
            for( const auto &vv : vectors) {
                auto ins4 = {
                    Incoming(dd),
                    AntiPart(Incoming(dd1)),
                    Outgoing(vv),
                    AntiPart(Outgoing(vv))
                };
                calculateAndGenerateLib("Tree", toyModel, ins4, myLib,qData,pData);
                
            }*/
        }
        
        /*
        for( const auto &sb : smbath) {
            for( const auto &sb1 : smbath) {
                auto ins4 = {
                    Incoming(dd),
                    Outgoing(sb),
                    AntiPart(Outgoing(sb1))
                };
                calculateAndGenerateLib("Tree", toyModel, ins4, myLib,qData,pData);
                
                //calculateAndGenerateLib("Asym", toyModel,ins4, myLib,qData,pData,false,true);
                //calculateAndGenerateLib("Asym", toyModel,AntiPart(ins4), myLib,qData,pData,false,true);
                
                auto ins5 = {
                    Incoming(dd),
                    Outgoing(sb),
                    Outgoing(sb1)
                };
                calculateAndGenerateLib("Tree", toyModel, ins5, myLib,qData,pData);
                
                //calculateAndGenerateLib("Asym", toyModel,ins5, myLib,qData,pData,false,true);
                //calculateAndGenerateLib("Asym", toyModel,AntiPart(ins5), myLib,qData,pData,false,true);
            }
        }
         */
    }
    /**/
    
    //============================================
    
    
    auto insertion = {Incoming("N_3"), Outgoing("lL_3"), AntiPart(Outgoing("et"))};
    
    //calculateAndGenerateLib("Asym", toyModel,insertion, myLib,qData,pData,true,true);
    
    //calculateAndGenerateLib("Asym", toyModel,AntiPart(insertion), myLib,qData,pData,true,true);
    
    //=====================================================================================
    /*
     auto insertions =  {Incoming("etm"),Incoming("etm"),Outgoing("le"),Outgoing("le")};
     
     calculateAndGenerateLib("Lepto", toyModel, insertions, myLib,qData,pData);
     
     calculateAndGenerateLib("Lepto", toyModel, AntiPart(insertions) , myLib,qData,pData);
     
     
     for (const auto &ll : leptons) {
     for (const auto &ll1 : leptons) {
     auto insertions =  {Incoming("etm"),Incoming("etm"),Outgoing(ll),Outgoing(ll1)};
     
     calculateAndGenerateLib("Lepto", toyModel, insertions, myLib,qData,pData);
     
     calculateAndGenerateLib("Lepto", toyModel, AntiPart(insertions) , myLib,qData,pData);
     }
     }
     
     for (const auto &ll : leptons) {
     
     auto ins0 = {Incoming("N_3"),Outgoing(ll),AntiPart(Outgoing("etm"))};
     
     calculateAndGenerateLib("Lepto", toyModel, ins0, myLib,qData,pData);
     
     calculateAndGenerateLib("Lepto", toyModel, AntiPart(ins0) , myLib,qData,pData);
     
     for( const auto &vv : vectors) {
     
     auto ins1 =  {Incoming("etm"),AntiPart(Incoming("N_3")),Outgoing(ll), Outgoing(vv)};
     
     calculateAndGenerateLib("Lepto", toyModel, ins1, myLib,qData,pData);
     
     calculateAndGenerateLib("Lepto", toyModel, AntiPart(ins1) , myLib,qData,pData);
     }
     
     }
     
     for (const auto &ll : neutrinos) {
     auto insertions = {
     Incoming("N_3"),
     Outgoing(ll),
     AntiPart(Outgoing("eta_R"))
     };
     calculateAndGenerateLib("Lepto", toyModel, insertions, myLib,qData,pData);
     
     calculateAndGenerateLib("Lepto", toyModel, AntiPart(insertions) , myLib,qData,pData);
     
     insertions = {
     Incoming("N_3"),
     Outgoing(ll),
     AntiPart(Outgoing("eta_I"))
     };
     calculateAndGenerateLib("Lepto", toyModel, insertions, myLib,qData,pData);
     
     calculateAndGenerateLib("Lepto", toyModel, AntiPart(insertions) , myLib,qData,pData);
     }
     */
    /*
     calculateAndGenerateLib("Lepto", toyModel, {Incoming("N_1"), Outgoing("le"), AntiPart(Outgoing("etm"))}, myLib,qData,pData);
     calculateAndGenerateLib("Lepto", toyModel, {Incoming("N_1"), Outgoing("lmu"), AntiPart(Outgoing("etm"))}, myLib,qData,pData);
     calculateAndGenerateLib("Lepto", toyModel, {Incoming("N_1"), Outgoing("ltau"), AntiPart(Outgoing("etm"))},myLib,qData,pData);
     
     calculateAndGenerateLib("Lepto", toyModel, {Incoming("N_1"), Outgoing("nu_e"), AntiPart(Outgoing("eta_R"))}, myLib,qData,pData);
     calculateAndGenerateLib("Lepto", toyModel, {Incoming("N_1"), Outgoing("nu_mu"), AntiPart(Outgoing("eta_R"))}, myLib,qData,pData);
     calculateAndGenerateLib("Lepto", toyModel, {Incoming("N_1"), Outgoing("nu_tau"), AntiPart(Outgoing("eta_R"))}, myLib,qData,pData);
     
     calculateAndGenerateLib("Lepto", toyModel, {Incoming("N_1"), Outgoing("nu_e"), AntiPart(Outgoing("eta_I"))}, myLib,qData,pData);
     calculateAndGenerateLib("Lepto", toyModel, {Incoming("N_1"), Outgoing("nu_mu"), AntiPart(Outgoing("eta_I"))}, myLib,qData,pData);
     calculateAndGenerateLib("Lepto", toyModel, {Incoming("N_1"), Outgoing("nu_tau"), AntiPart(Outgoing("eta_I"))}, myLib,qData,pData);
     */
    /*
     
     
     */
    
    
    ///////////////////////////////////////////////////
    // Saving JSON data to test.json
    ///////////////////////////////////////////////////
    
    saveParticleData("test.json", qData, pData);
    
    //myLib.applyDiagonalizationData(toyModel);
    myLib.generateSpectrum(toyModel);
    
    myLib.cleanExistingSources();
    
    myLib.setGccCompiler();
    myLib.addIPath("/usr/local/Cellar/boost/1.75.0_1/include");
    myLib.addIPath("/usr/local/Cellar/eigen/3.3.8_1/include");
    myLib.addLPath("/usr/local/Cellar/boost/1.75.0_1/lib");
    myLib.addLPath("/usr/local/lib");
    myLib.addIPath("/usr/local/include");
    
    myLib.build(4);
    
    myLib.print();
    
    
    
    // Exporting the additional cpp and h files to the library
    
    [[maybe_unused]] int res = system("./export_to_lib.sh");
    
    std::cout << " Hi!! " << std::endl;
    
    
    
    return 0;
}
