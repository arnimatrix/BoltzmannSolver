
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include "BEv1.h"
#include "brparity.h"

#include <ctime>

using namespace Eigen;
using Eigen::MatrixXd;
using namespace boost::math;

using namespace brparity;

using namespace std;

int main(int argc, char const**)
{

    param_t pp;

    double ee     = 1. / 137.;
    double sinthW = std::sqrt(0.23121);

    double mu, md, ms, mc, mb, mt;
    mu = 2.16E-03;
    md = 4.67E-3;
    ms = 9.3E-2;
    mc = 1.27;
    mb = 4.18;
    mt = 172.76;

    pp.m_u  = mu;
    pp.m_d  = md;
    pp.m_s  = ms;
    pp.m_c  = mc;
    pp.m_b  = mb;
    pp.m_t  = mt;
    pp.gw   = ee * sinthW;
    pp.gwuL = 1. / 6.;
    pp.gwdL = 1. / 6.;
    pp.gwuR = 2. / 3.;
    pp.gwdR = -1. / 3.;
    pp.mX   = 1.E+02;

    Eigen::VectorXcd th, mm;
    th    = Eigen::VectorXcd::Zero(3);
    mm    = Eigen::VectorXcd::Zero(3);
    th(0) = M_PI / 4. + 0.01 * 1i;
    mm(0) = -pow(4.5E+03, 2.);
    mm(1) = -pow(6.2E+03, 2.);
    mm(2) = -pow(7.6E+03, 2.);

    cout << endl << "mass =" << mm << endl;

    Eigen::MatrixXcd mds;

    mds = Mass(th, mm);

    pp.mdL_00 = mds(0, 0);
    pp.mdL_01 = mds(0, 1);
    pp.mdL_02 = mds(0, 2);
    pp.mdL_10 = mds(1, 0);
    pp.mdL_11 = mds(1, 1);
    pp.mdL_12 = mds(1, 2);
    pp.mdL_20 = mds(2, 0);
    pp.mdL_21 = mds(2, 1);
    pp.mdL_22 = mds(2, 2);

    pp.mdR_00 = mds(0, 0);
    pp.mdR_01 = mds(0, 1);
    pp.mdR_02 = mds(0, 2);
    pp.mdR_10 = mds(1, 0);
    pp.mdR_11 = mds(1, 1);
    pp.mdR_12 = mds(1, 2);
    pp.mdR_20 = mds(2, 0);
    pp.mdR_21 = mds(2, 1);
    pp.mdR_22 = mds(2, 2);

    mm(0) = -pow(3.3E+03, 2.);
    mm(1) = -pow(3.4E+03, 2.);
    mm(2) = -pow(4.5E+03, 2.);

    Eigen::MatrixXcd mus;

    mus = Mass(th, mm);

    pp.muL_00 = mus(0, 0);
    pp.muL_01 = mus(0, 1);
    pp.muL_02 = mus(0, 2);
    pp.muL_10 = mus(1, 0);
    pp.muL_11 = mus(1, 1);
    pp.muL_12 = mus(1, 2);
    pp.muL_20 = mus(2, 0);
    pp.muL_21 = mus(2, 1);
    pp.muL_22 = mus(2, 2);

    pp.muR_00 = mus(0, 0);
    pp.muR_01 = mus(0, 1);
    pp.muR_02 = mus(0, 2);
    pp.muR_10 = mus(1, 0);
    pp.muR_11 = mus(1, 1);
    pp.muR_12 = mus(1, 2);
    pp.muR_20 = mus(2, 0);
    pp.muR_21 = mus(2, 1);
    pp.muR_22 = mus(2, 2);

    /*
    pp.lpp_000 = 0.10 + 0.5*1i;
    pp.lpp_001 = 0.12 + 0.9*1i;
    pp.lpp_002 = 0.13 + 0.3*1i;
    pp.lpp_100 = 0.32 + 0.7*1i;
    pp.lpp_200 = 0.33 + 0.1*1i;
    pp.lpp_010 = 0.22 + 0.7*1i;
    pp.lpp_020 = 0.23 + 0.2*1i;
    */
    cout << endl
         << "The mass matrix Down squark is given as :" << endl
         << mds << endl;

    cout << endl
         << "The mass matrix UP squark is given as :" << endl
         << mus << endl;

    cout << endl << "The mass matrix is given as :" << endl << mds << endl;

    updateDiagonalization(pp);

    vector<pair<std::string, double>> sqDecay;
    nthDecay(sqDecay, pp);

    pp.GbLt = 0.;
    pp.GbRt = squarkD(sqDecay, "dtR3");
    pp.GcLt = 0.;
    pp.GcRt = squarkD(sqDecay, "dtR2");
    pp.GdLt = 0.;
    pp.GdRt = squarkD(sqDecay, "dtR1");
    pp.GsLt = 0.;
    pp.GsRt = squarkD(sqDecay, "utR2");
    pp.GtLt = 0.;
    pp.GtRt = squarkD(sqDecay, "utR3");
    pp.GuLt = 0.;
    pp.GuRt = squarkD(sqDecay, "utR1");

    cout << "mdtL1 = " << pp.m_dtL1 << endl;
    cout << "mdtL2 = " << pp.m_dtL2 << endl;
    cout << "mdtL3 = " << pp.m_dtL3 << endl;
    cout << "mdtR1 = " << pp.m_dtR1 << endl;
    cout << "mdtR2 = " << pp.m_dtR2 << endl;
    cout << "mdtR3 = " << pp.m_dtR3 << endl;

    cout << "mutL1 = " << pp.m_utL1 << endl;
    cout << "mutL2 = " << pp.m_utL2 << endl;
    cout << "mutL3 = " << pp.m_utL3 << endl;
    cout << "mutR1 = " << pp.m_utR1 << endl;
    cout << "mutR2 = " << pp.m_utR2 << endl;
    cout << "mutR3 = " << pp.m_utR3 << endl;

    cout << "mXX : " << prtlm("XX", pp) << endl;
    cout << "mXX = " << pp.mX << endl;
    std::array<std::string, 2> pname1 = { "CP_XX_to_uf1_df2_df1",
                                          "CP_XXc_to_uf1_df2_df1" };
    std::array<std::string, 2> pname2 = { "CP_XX_to_uf1c_df2c_df1c",
                                          "CP_XXc_to_uf1c_df2c_df1c" };

    std::array<std::string, 2> pname3 = { "CP_XX_XXc_to_uf1_uf1c",
                                          "CP_XX_XXc_to_uf2_uf2c" };

    //================ Function =========================

    complex_t decayB = 0., decayAntiB = 0.;
    double    gammaB = 0., gammaAntiB = 0., ann = 0.;
    int       pid;
    cout << "mXX = " << prtlm("XX", pp) << endl;
    for (size_t ll = 0; ll != pname1.size(); ++ll) {
        findAmp(pname1[ll], pid);
        decayB += Decay(pid, pp);
        gammaB += gamma(pid, prtlm("XX", pp), pp);
    }
    cout << "DecayB = " << decayB << endl;
    cout << "gammaB = " << gammaB << endl;

    for (size_t ll = 0; ll != pname2.size(); ++ll) {
        findAmp(pname2[ll], pid);
        decayAntiB += Decay(pid, pp);
        gammaAntiB += gamma(pid, prtlm("XX", pp), pp);
    }
    cout << "DecayAntiB = " << decayAntiB << endl;
    cout << "gammaAntiB = " << gammaAntiB << endl;

    cout << "Asym Decay : " << abs(decayB - decayAntiB) << endl;

    cout << "Asym Decay : " << abs(gammaB - gammaAntiB) << endl;
    /**/
    for (int ll = 0; ll != pname3.size(); ++ll) {
        // ann += gamma(pname3[ll],prtlm("XX",pp),pp);
    }

    vector<std::string> pn;

    pn.push_back("XX");
    pn.push_back("XXc");
    pn.push_back("uf1");
    pn.push_back("uf1c");

    cout << endl << "Pr name = " << Prname(pn, 2, 2) << endl;
    findAmp(Prname(pn, 2, 2), pid);
    cout << endl << "Gamma = " << gamma(pid, prtlm("XX", pp), pp) << endl;

    for (auto& pX : pXSM) {
        cout << "Name : " << pX.first << " Charge : " << pX.second << endl;
    }

    initialize();
    // To be coded in the functh

    //==== Input

    //============

    /**/
    double xi  = 1.0E-01;
    double xst = log10(xi);

    double xf  = m_XX(pp).real() / 0.001;
    double xfl = log10(xf);

    size_t           nn = pXSM.size() + 1;
    Eigen::VectorXcd yn(nn);

    double T = m_XX(pp).real() / xi;

    for (size_t pi = 0; pi != pXSM.size(); ++pi) {
        double gg = 1.;
        for (auto& gi : gptl) {
            gg *= (pXSM[pi].first.find(gi.first) != std::string::npos)
                    ? gi.second
                    : 1.;
        }
        yn(pi) = Yeq(gg, prtlm(pXSM[pi].first, pp), T);

        cout << "Yeq[" << pi << "]" << yn(pi) << endl;
    }
    yn(nn - 1) = 0.;

    ofstream ofile("result.txt");

    diff(ofile, xfl, xst, yn, pp);

    ofile.close();

    //============================================================

    return 0;
}
