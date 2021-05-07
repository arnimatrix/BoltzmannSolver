#include <iostream>
#include <fstream>
#include "brparity.h"
#include "BEv1.h"

using Eigen::MatrixXd;
using namespace boost::math;
using namespace brparity;
using namespace std;

plist XXSMSM[1][1];
plist XSMXSM[1][1];

plist BXXSMSM[1][1];
plist AntiBXXSMSM[1][1];

plist BXSMXSM[1][1];
plist AntiBXSMXSM[1][1];


plist XSMSMSM[1];
plist BXSMSMSM[1];
plist AntiBXSMSMSM[1];

double prtlm(std::string str,param_t const &pp)
{
    double mass = 0.;
    
    mass += (str.find("uf1")!=std::string::npos) ? pp.m_u : 0.;
    mass += (str.find("uf2")!=std::string::npos) ? pp.m_c : 0.;
    mass += (str.find("uf3")!=std::string::npos) ? pp.m_t : 0.;
    mass += (str.find("df1")!=std::string::npos) ? pp.m_d : 0.;
    mass += (str.find("df2")!=std::string::npos) ? pp.m_s : 0.;
    mass += (str.find("df3")!=std::string::npos) ? pp.m_b : 0.;
    mass += (str.find("utR1")!=std::string::npos) ? pp.m_utR1 : 0.;
    mass += (str.find("utR2")!=std::string::npos) ? pp.m_utR2 : 0.;
    mass += (str.find("utR3")!=std::string::npos) ? pp.m_utR3 : 0.;
    mass += (str.find("dtR1")!=std::string::npos) ? pp.m_dtR1 : 0.;
    mass += (str.find("dtR2")!=std::string::npos) ? pp.m_dtR2 : 0.;
    mass += (str.find("dtR3")!=std::string::npos) ? pp.m_dtR3 : 0.;
    mass += (str.find("XX")!=std::string::npos) ? m_XX(pp).real() : 0.;
    
    //cout << "Mass :" << str << " => " << mass << endl;
    return mass;
    
}

void facB(int &pid,pair <double,double> &fac,double const &T)
{
    std::string pname = f_G[pid].name;
    vector<pair <size_t,double>> pp;
    
    double sp = geffT(T)*2*pow(M_PI,2.)/45.;
    
    double Ne = 0.,Nu = 0.,Nd = 0.;
    
    for(auto &mm : ml)
        Ne += (T>mm) ? 1. : 0.;
    
    for(auto &mm : muq)
        Nu += (T>mm) ? 1. : 0.;
    
    for(auto &mm : mdq)
        Nd += (T>mm) ? 1. : 0.;
    
    double mu_u = 0.,mu_d = 0.;
    
    
    mu_u += (T>100.) ? -10./28. : 6.*(1./3. + 0.5/Nd + 0.5/Ne)/(1. + 3.*Nu/Ne + Nu/Nd + 2.*Nu);
    
    
    mu_d += (T>100.) ? 38./28.: 6.*(1. - 1./3.*Nu*mu_u)/(2.*Nd);
    
    
    size_t found = pname.find("uf");
    while (found!=std::string::npos)
    {
        pp.push_back(make_pair(found,sp*mu_u));
        found = pname.find("uf",found+1);
    }
    
    found = pname.find("df");
    while (found!=std::string::npos)
    {
        pp.push_back(make_pair(found,sp*mu_d));
        found = pname.find("df",found+1);
    }
    
    found = pname.find("XX");
    while (found!=std::string::npos)
    {
        pp.push_back(make_pair(found,0.));
        found = pname.find("XX",found+1);
    }
    found = pname.find("to");
    
    double in = 0.,out = 0.;
    for(auto & pv : pp)
    {
        (pv.first<found) ? in+=pv.second : out+=pv.second;
    }
    
    fac = (pid < ((int) f_G.size())) ? make_pair(in,out) : make_pair(0.,0.);
    
}



complex_t Complex(int const &x,double const &y)
{
    return x + y*1i;
}
std::string Prname(vector<std::string> const &ptln,int const &inp,int const &outp)
{
    std::string res = "CP";
    
    if(inp+outp ==((int) ptln.size()))
    {
        for (int ii = 0; ii<inp; ii++)
        {
            res += sp + ptln[ii] ;
        }
        res += sp + to;
        for (int ii = inp; ii<inp+outp; ii++)
        {
            res += sp + ptln[ii];
        }
    }
    return res;
}


double geffT(double x)
{
    int size = temp.size();
    
    int i = 0;                                                                  // find left end of interval for interpolation
    if ( x >= temp[size - 2] )                                                 // special case: beyond right end
    {
        i = size - 2;
    }
    else
    {
        while ( x > temp[i+1] ) i++;
    }
    double xL = temp[i], yL = geff[i], xR = geff[i+1], yR = geff[i+1];      // points on either side (unless beyond ends)
    if ( x < xL ) yR = yL;
    if ( x > xR ) yL = yR;
    
    
    double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient
    
    return yL + dydx * ( x - xL );                                              // linear interpolation
}

double Yeq(double gi,double mi,double T)
{
    double yy = 0.;
    if(mi>0.)
    {
        yy = 45.*gi/geffT(T)*pow(mi/(2.*M_PI*T),2.)*boost::math::cyl_bessel_k(2, mi/T);
    }
    if(mi==0.)
    {
        yy = 45./2.*gi/geffT(T)*pow(M_PI,-4.); //cout << "h1!!" << endl;
    }
    return yy;
}

double squarkD(vector<pair <std::string,double>> sg,std::string str)
{
    double decay = 0.;
    
    for(auto &data : sg){
        if(data.first==str) {
            decay = data.second;
        }
    }
    return decay;
}



Eigen::MatrixXcd Mass(Eigen::VectorXcd &th,Eigen::VectorXcd &mm)
{
    Eigen::MatrixXcd res;
    Eigen::MatrixXcd O12(3,3);
    Eigen::MatrixXcd O13(3,3);
    Eigen::MatrixXcd O23(3,3);
    Eigen::MatrixXcd MM(3,3);
    
    cout << endl << th << endl;
    cout << endl << mm << endl;
    Eigen::MatrixXcd temp(3,3);
    
    O12 << cos(th(0)),-sin(th(0)),0.,
    sin(th(0)),cos(th(0)),0.,
    0.,0.,1.;
    O13 << cos(th(1)),0.,-sin(th(1)),
    0.,1.,0.,
    sin(th(1)),0.,cos(th(1));
    O23 << 1.,0.,0.,
    0.,cos(th(2)),-sin(th(2)),
    0.,sin(th(2)),cos(th(2));
    
    cout << endl << O12 << endl;
    cout << endl << O13 << endl;
    cout << endl << O23 << endl;
    MM << mm(0), 0., 0.,0.,mm(1),0.,0.,0.,mm(2);
    
    temp = O12 * O13 * O23;
    
    //res = temp.inverse() * MM * temp;
    
    res = (temp.conjugate()).transpose() * MM * temp;
    
    return res;
}



Eigen::MatrixXcd expM(Eigen::MatrixXd M,int &dim)
{
    Eigen::MatrixXcd res,ev;
    //Eigen::MatrixXcd evt;
    
    Eigen::EigenSolver<MatrixXd> es(M);
    
    //cout << "M : " << endl << M << endl << endl;
    ev = Eigen::MatrixXcd::Zero(dim,dim);
    //evt = Eigen::MatrixXcd::Zero(dim,dim);
    //cout << "ev : " << endl << ev << endl << endl;
    int i;
    for(i=0; i<dim; i++)
    {
        ev(i,i) = exp(es.eigenvalues()[i]);
        //evt(i,i) = es.eigenvalues()[i];
    }
    //cout << "ev : " << endl << ev << endl << endl;
    Eigen::MatrixXcd V = es.eigenvectors();
    
    //cout << "evt : " << endl << V*evt*V.inverse() << endl << endl;
    //cout << "M : " << endl << M << endl << endl;
    res = V*ev*V.inverse();
    return res;
}

complex_t L(complex_t xx,complex_t yy,complex_t zz)
{
    return xx*xx + yy*yy + zz*zz - 2.*xx*yy - 2.*xx*zz - 2.*yy*zz;
}

double jac(complex_t sqss,Eigen::VectorXcd &mm,vector<double> const &rr,int const &intp,int const &outp)
{
    double res;
    if((intp==1)&&(outp==3)){
        double mD = mm(0).real();
        double m1 = mm(1).real();
        double m2 = mm(2).real();
        double m3 = mm(3).real();
        double s23 = pow(m2 + m3 + (-m1 - m2 - m3 + mD)*rr[0],2.);
        
        res = 0.5/mD*1./(4.*M_PI)*sqrt(L(mD*mD,m1*m1,s23).real())/(8.*M_PI*mD*mD)*sqrt(L(s23,m2*m2,m3*m3)).real()/(8.*M_PI*s23);
    }
    if((intp==1)&&(outp==2)){
        double mD = mm(0).real();
        double m1 = mm(1).real();
        double m2 = mm(2).real();
        
        res = 1./(16.*M_PI*pow(mD,3.))*sqrt(L(mD*mD,m1*m1,m2*m2)).real();
    }
    if((intp==2)&&(outp==2)){
        complex_t ss = pow(sqss,2.);
        double ma = mm(0).real();
        double mb = mm(1).real();
        double m2 = mm(2).real();
        double m3 = mm(3).real();
        res = 1./(64.*M_PI*M_PI*ss.real())*sqrt(L(ss,m2*m2,m3*m3)).real()/sqrt(L(ss,ma*ma,mb*mb)).real();
    }
    return res;
}


Eigen::MatrixXcd pmom(
        complex_t sqss,
        Eigen::VectorXcd &mm,
        vector<double> const &cthz,
        vector<double> const &rr,
        int const &inp, 
        int const &outp)
{
    Eigen::MatrixXcd res(inp+outp,inp+outp);
    
    std::array<Eigen::Vector3cd,4> pmom;
    
    for (int ii=0; ii!=4; ++ii){pmom[ii] = Eigen::VectorXcd::Zero(3);}
    
    if((inp==1)&&(outp==3))
    {
        double mD = mm(0).real();
        double m1 = mm(1).real();
        double m2 = mm(2).real();
        double m3 = mm(3).real();
        double cth1 = cthz[0];
        double cth2 = cthz[1];
        
        double s23 = pow(m2 + m3 + (-m1 - m2 - m3 + mD)*rr[0],2.);
        
        res = Eigen::MatrixXcd::Zero(4,4);
        
        pmom[0] << mD, 0.,0.;
        pmom[1] << (sqrt(pow(mD,2))*(1 + pow(m1,2)/pow(mD,2) - s23/pow(mD,2)))/2.,
        (Complex(0,-0.5)*sqrt(1 - pow(cth1,2))*sqrt(L(pow(mD,2),pow(m1,2),s23)))/sqrt(pow(mD,2)),
        (Complex(0,-0.5)*cth1*sqrt(L(pow(mD,2),pow(m1,2),s23)))/sqrt(pow(mD,2));
        pmom[2] << (-((pow(m1,2) - pow(mD,2) - s23)*(pow(m2,2) - pow(m3,2) + s23)) +
                    (cth1*cth2 + sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2)))*sqrt(L(pow(mD,2),pow(m1,2),s23))*sqrt(L(s23,pow(m2,2),pow(m3,2))))/
        (4.*sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(sqrt(1 - pow(cth1,2))*(pow(m2,2) - pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) +
                                                   (pow(cth1,2)*sqrt(1 - pow(cth2,2))*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) +
                                                    sqrt(1 - pow(cth2,2))*(-pow(m1,2) + pow(mD,2) + s23) +
                                                    cth1*sqrt(1 - pow(cth1,2))*cth2*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23))*sqrt(L(s23,pow(m2,2),pow(m3,2)))
                                                   ))/(sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(cth1*(pow(m2,2) - pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) +
                                                                                              (2*cth2*sqrt(pow(mD,2))*sqrt(s23) + pow(cth1,2)*cth2*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23) +
                                                                                               cth1*sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2))*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23))*
                                                                                              sqrt(L(s23,pow(m2,2),pow(m3,2)))))/(sqrt(pow(mD,2))*s23);
        pmom[3] << (-((pow(m2,2) - pow(m3,2) - s23)*(-pow(m1,2) + pow(mD,2) + s23)) -
                    (cth1*cth2 + sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2)))*sqrt(L(pow(mD,2),pow(m1,2),s23))*sqrt(L(s23,pow(m2,2),pow(m3,2))))/
        (4.*sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(sqrt(1 - pow(cth1,2))*(-pow(m2,2) + pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) -
                                                   (pow(cth1,2)*sqrt(1 - pow(cth2,2))*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) +
                                                    sqrt(1 - pow(cth2,2))*(-pow(m1,2) + pow(mD,2) + s23) +
                                                    cth1*sqrt(1 - pow(cth1,2))*cth2*(-pow(m1,2) + pow(mD,2) - 2*sqrt(pow(mD,2))*sqrt(s23) + s23))*sqrt(L(s23,pow(m2,2),pow(m3,2)))
                                                   ))/(sqrt(pow(mD,2))*s23),(Complex(0,0.25)*(cth1*(-pow(m2,2) + pow(m3,2) + s23)*sqrt(L(pow(mD,2),pow(m1,2),s23)) +
                                                                                              (pow(cth1,2)*cth2*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) +
                                                                                               cth1*sqrt(1 - pow(cth1,2))*sqrt(1 - pow(cth2,2))*(pow(m1,2) - pow(mD,2) + 2*sqrt(pow(mD,2))*sqrt(s23) - s23) -
                                                                                               2*cth2*sqrt(pow(mD,2))*sqrt(s23))*sqrt(L(s23,pow(m2,2),pow(m3,2)))))/(sqrt(pow(mD,2))*s23);
        
        res << pmom[0].transpose() * pmom[0], pmom[0].transpose() * pmom[1],      pmom[0].transpose() * pmom[2], pmom[0].transpose() * pmom[3],
        pmom[1].transpose() * pmom[0], pmom[1].transpose() * pmom[1], pmom[1].transpose() * pmom[2], pmom[1].transpose() * pmom[3],
        pmom[2].transpose() * pmom[0], pmom[2].transpose() * pmom[1], pmom[2].transpose() * pmom[2], pmom[2].transpose() * pmom[3],
        pmom[3].transpose() * pmom[0], pmom[3].transpose() * pmom[1], pmom[3].transpose() * pmom[2], pmom[3].transpose() * pmom[3];
    }
    if((inp==1)&&(outp==2)) {
        double mD = mm(0).real();
        double m1 = mm(1).real();
        double m2 = mm(2).real();
        
        pmom[0] << mD, 0.,0.;
        pmom[1] << (pow(m1,2) - pow(m2,2) + pow(mD,2))/(2.*sqrt(pow(mD,2))),0,
        -0.5*(sqrt(L(pow(mD,2),pow(m1,2),pow(m2,2))))/sqrt(pow(mD,2))*1i;
        pmom[2] << (pow(m1,2) - pow(m2,2) + pow(mD,2))/(2.*sqrt(pow(mD,2))),0,
        0.5*(sqrt(L(pow(mD,2),pow(m1,2),pow(m2,2))))/sqrt(pow(mD,2))*1i;
        
        res = Eigen::MatrixXcd::Zero(3,3);
        
        res << pmom[0].transpose() * pmom[0], pmom[0].transpose() * pmom[1],      pmom[0].transpose() * pmom[2],
        pmom[1].transpose() * pmom[0], pmom[1].transpose() * pmom[1], pmom[1].transpose() * pmom[2],
        pmom[2].transpose() * pmom[0], pmom[2].transpose() * pmom[1], pmom[2].transpose() * pmom[2];
    }
    if((inp==2)&&(outp==2)) {
        double ss = pow(sqss.real(),2.);
        double m1 = mm(0).real();
        double m2 = mm(1).real();
        double m3 = mm(2).real();
        double m4 = mm(3).real();
        double cth = cthz[0];
        
        pmom[0] << ((1 + pow(m1,2.)/ss - pow(m2,2.)/ss)*sqrt(ss))/2.,0,-sqrt(L(ss,pow(m1,2.),pow(m2,2.)))/(2.*sqrt(ss))*1i;
        pmom[1] << ((1 - pow(m1,2.)/ss + pow(m2,2.)/ss)*sqrt(ss))/2.,0,sqrt(L(ss,pow(m1,2),pow(m2,2.)))/(2.*sqrt(ss))*1i;
        pmom[2] << ((1 + pow(m3,2)/ss - pow(m4,2)/ss)*sqrt(ss))/2.,-0.5*(sqrt(1 - pow(cth,2))*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i,
        -0.5*(cth*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i;
        pmom[3] << ((1 - pow(m3,2.)/ss + pow(m4,2.)/ss)*sqrt(ss))/2.,0.5*(sqrt(1 - pow(cth,2.))*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i,
        0.5*(cth*sqrt(L(ss,pow(m3,2.),pow(m4,2.))))/sqrt(ss)*1i;
        
        res = Eigen::MatrixXcd::Zero(4,4);
        
        res << pmom[0].transpose() * pmom[0], pmom[0].transpose() * pmom[1],      pmom[0].transpose() * pmom[2], pmom[0].transpose() * pmom[3],
        pmom[1].transpose() * pmom[0], pmom[1].transpose() * pmom[1], pmom[1].transpose() * pmom[2], pmom[1].transpose() * pmom[3],
        pmom[2].transpose() * pmom[0], pmom[2].transpose() * pmom[1], pmom[2].transpose() * pmom[2], pmom[2].transpose() * pmom[3],
        pmom[3].transpose() * pmom[0], pmom[3].transpose() * pmom[1], pmom[3].transpose() * pmom[2], pmom[3].transpose() * pmom[3];
    }
    return res;
}

void findAmp(std::string pname,int &id)
{
    id = f_G.size() + 1;
    for(int ll=0; ll!=f_G.size();++ll)
    {
        if (f_G[ll].name == pname) {
            id = ll;
        }
    }
    
}

double ptlSM(std::string str)
{
    double res = 0.;
    for (auto &chptl : chSM) {
        res += (chptl.first == str) ? chptl.second : 0.;
    }
    return res;
}

double ptlB(std::string str)
{
    double res = 0.;
    for (auto &chptl : chB) {
        res += (chptl.first == str) ? chptl.second : 0.;
    }
    return res;
}


void nthDecay(vector<pair <std::string,double>> &sqDecay,param_t &pp)
{
    int pid;
    double ch1,ch2,ch3,charge;
    
    std::string f1,f2;
    for(size_t qq = 0;qq!=pnthSM.size();++qq)
    {
        double decay = 0.;
        ch1 = ptlSM(pnthSM[qq].first);
        std::string squark = pnthSM[qq].first;
        for(size_t qi = 0;qi!=pthSM.size();++qi){
            for(size_t qj = 0;qj!=pthSM.size();++qj){
                
                ch2 = ptlSM(pthSM[qi].first);
                ch3 = ptlSM(pthSM[qj].first);
                
                f1 = pthSM[qi].first;
                f2 = pthSM[qj].first;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn1;
                    pn1.push_back(squark);
                    pn1.push_back(f1);
                    pn1.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn1,1,2) << endl;
                    findAmp(Prname(pn1,1,2),pid);
                    decay += (pid<((int) f_G.size())) ? Decay(pid,pp).real() : 0.;
                    
                }
                
                ch2 = ptlSM(pthSM[qi].second);
                ch3 = ptlSM(pthSM[qj].first);
                
                f1 = pthSM[qi].second;
                f2 = pthSM[qj].first;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn1;
                    pn1.push_back(squark);
                    pn1.push_back(f1);
                    pn1.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn1,1,2) << endl;
                    findAmp(Prname(pn1,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                    
                }
                
                ch2 = ptlSM(pthSM[qi].first);
                ch3 = ptlSM(pthSM[qj].second);
                
                f1 = pthSM[qi].first;
                f2 = pthSM[qj].second;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn1;
                    pn1.push_back(squark);
                    pn1.push_back(f1);
                    pn1.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn1,1,2) << endl;
                    findAmp(Prname(pn1,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;}
                
                ch2 = ptlSM(pthSM[qi].second);
                ch3 = ptlSM(pthSM[qj].second);
                
                f1 = pthSM[qi].second;
                f2 = pthSM[qj].second;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn1;
                    pn1.push_back(squark);
                    pn1.push_back(f1);
                    pn1.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn1,1,2) << endl;
                    findAmp(Prname(pn1,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;}
            }
            /**/
            for(size_t qj = 0;qj!=pXSM.size();++qj){
                
                ch2 = ptlSM(pthSM[qi].first);
                ch3 = ptlSM(pXSM[qj].first);
                
                f1 = pthSM[qi].first;
                f2 = pXSM[qj].first;
                
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn2;
                    pn2.push_back(squark);
                    pn2.push_back(f1);
                    pn2.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn2,1,2) << endl;
                    findAmp(Prname(pn2,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                }
                
                ch2 = ptlSM(pthSM[qi].second);
                ch3 = ptlSM(pXSM[qj].first);
                
                f1 = pthSM[qi].second;
                f2 = pXSM[qj].first;
                charge = ch1 - ch2 - ch3;
                
                if (charge == 0) {
                    vector<std::string> pn2;
                    pn2.push_back(squark);
                    pn2.push_back(f1);
                    pn2.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn2,1,2) << endl;
                    findAmp(Prname(pn2,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                }
                
                ch2 = ptlSM(pthSM[qi].first);
                ch3 = ptlSM(pXSM[qj].second);
                
                f1 = pthSM[qi].first;
                f2 = pXSM[qj].second;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn2;
                    pn2.push_back(squark);
                    pn2.push_back(f1);
                    pn2.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn2,1,2) << endl;
                    findAmp(Prname(pn2,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                }
                
                ch2 = ptlSM(pthSM[qi].second);
                ch3 = ptlSM(pXSM[qj].second);
                
                f1 = pthSM[qi].second;
                f2 = pXSM[qj].second;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn2;
                    pn2.push_back(squark);
                    pn2.push_back(f1);
                    pn2.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn2,1,2) << endl;
                    findAmp(Prname(pn2,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                    
                }
            }
        }
        
        squark = pnthSM[qq].second;
        for(size_t qi = 0;qi!=pthSM.size();++qi){
            for(size_t qj = 0;qj!=pthSM.size();++qj){
                
                ch2 = ptlSM(pthSM[qi].first);
                ch3 = ptlSM(pthSM[qj].first);
                
                f1 = pthSM[qi].first;
                f2 = pthSM[qj].first;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn1;
                    pn1.push_back(squark);
                    pn1.push_back(f1);
                    pn1.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn1,1,2) << endl;
                    findAmp(Prname(pn1,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                    
                }
                
                ch2 = ptlSM(pthSM[qi].second);
                ch3 = ptlSM(pthSM[qj].first);
                
                f1 = pthSM[qi].second;
                f2 = pthSM[qj].first;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn1;
                    pn1.push_back(squark);
                    pn1.push_back(f1);
                    pn1.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn1,1,2) << endl;
                    findAmp(Prname(pn1,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                    
                }
                
                ch2 = ptlSM(pthSM[qi].first);
                ch3 = ptlSM(pthSM[qj].second);
                
                f1 = pthSM[qi].first;
                f2 = pthSM[qj].second;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn1;
                    pn1.push_back(squark);
                    pn1.push_back(f1);
                    pn1.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn1,1,2) << endl;
                    findAmp(Prname(pn1,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                    
                }
                
                ch2 = ptlSM(pthSM[qi].second);
                ch3 = ptlSM(pthSM[qj].second);
                
                f1 = pthSM[qi].second;
                f2 = pthSM[qj].second;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn1;
                    pn1.push_back(squark);
                    pn1.push_back(f1);
                    pn1.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn1,1,2) << endl;
                    findAmp(Prname(pn1,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                    
                }
                
            }
            /**/
            for(size_t qj = 0;qj!=pXSM.size();++qj){
                
                ch2 = ptlSM(pthSM[qi].first);
                ch3 = ptlSM(pXSM[qj].first);
                
                f1 = pthSM[qi].first;
                f2 = pXSM[qj].first;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn2;
                    pn2.push_back(squark);
                    pn2.push_back(f1);
                    pn2.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn2,1,2) << endl;
                    findAmp(Prname(pn2,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;}
                
                ch2 = ptlSM(pthSM[qi].second);
                ch3 = ptlSM(pXSM[qj].first);
                
                f1 = pthSM[qi].second;
                f2 = pXSM[qj].first;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn2;
                    pn2.push_back(squark);
                    pn2.push_back(f1);
                    pn2.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn2,1,2) << endl;
                    findAmp(Prname(pn2,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;}
                
                ch2 = ptlSM(pthSM[qi].first);
                ch3 = ptlSM(pXSM[qj].second);
                
                f1 = pthSM[qi].first;
                f2 = pXSM[qj].second;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn2;
                    pn2.push_back(squark);
                    pn2.push_back(f1);
                    pn2.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn2,1,2) << endl;
                    findAmp(Prname(pn2,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;}
                
                ch2 = ptlSM(pthSM[qi].second);
                ch3 = ptlSM(pXSM[qj].second);
                
                f1 = pthSM[qi].second;
                f2 = pXSM[qj].second;
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn2;
                    pn2.push_back(squark);
                    pn2.push_back(f1);
                    pn2.push_back(f2);
                    cout << endl << "Pr name = "<< Prname(pn2,1,2) << endl;
                    findAmp(Prname(pn2,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;}
            }
        }
        sqDecay.push_back(std::make_pair(squark,decay));
        
        cout << "Decay of sq[ "<< "]"<< qq <<" = "<< squark << " = " << decay << endl;
    }
    
    
}

complex_t Decay(int &pid,param_t &pp)
{
    complex_t res = 0.;
    
    std::string pname = f_G[pid].name;
    
    vector< pair <int,double> > vect;
    vect.reserve(ptl.size());
    size_t found;
    
    for(int pt=0;pt!=ptl.size();++pt){
        found = pname.find(ptl[pt]);
        while (found!=std::string::npos) {
            vect.push_back( make_pair(found,prtlm(ptl[pt],pp)) );
            found = pname.find(ptl[pt],found+1);
        }
    }
    
    sort(vect.begin(), vect.end());
    int in = 0,out = 0;
    
    for(size_t tt = 0;tt!=vect.size();++tt){
        (vect[tt].first<pname.find("to")) ? in++ : out++;
    }
    
    //cout << "The process is " << in << "-> " << out << endl;
    
    complex_t jacD = 0.,decay = 0.;
    Eigen::MatrixXcd pmat(in+out,in+out);
    Eigen::VectorXcd mass(in+out);
    
    for(int tt = 0;tt!=vect.size();++tt){
        mass(tt) = vect[tt].second;
    }
    
    //cout << "The Masses are: "<< endl << mass  << endl;
    if(in==1){
        
        if(out==3){
            decay = 0.;
            for(int ii = 0;ii<13;ii++){
                for(int jj = 0;jj<13;jj++){
                    for(int kk = 0;kk<13;kk++){
                        
                        vector<double> cth,rr;
                        cth.push_back(xs[jj]);
                        cth.push_back(xs[kk]);
                        rr.push_back(xs1[ii]);
                        pmat = pmom(mass(0),mass,cth,rr,in,out);
                        pp.s_12 = pmat(0,1).real();
                        pp.s_13 = pmat(0,2).real();
                        pp.s_14 = pmat(0,3).real();
                        pp.s_23 = pmat(1,2).real();
                        pp.s_24 = pmat(1,3).real();
                        pp.s_34 = pmat(2,3).real();
                        
                        jacD = jac(mass(0),mass,rr,in,out);
                        /*
                         for(int ll=0; ll!=f_G.size();++ll){
                         if(f_G[ll].name == pname){
                         decay += f_G[ll](pp)*jacD*as[jj]*as[kk]*as1[ii];
                         }
                         }*/
                        decay += f_G[pid](pp)*jacD*as[jj]*as[kk]*as1[ii];
                    }
                }
            }
        }
        if(out==2){
            decay = 0.;
            vector<double> cth,rr;
            cth.push_back(1.);
            rr.push_back(1.);
            pmat = pmom(mass(0),mass,cth,rr,in,out);
            pp.s_12 = pmat(0,1).real();
            pp.s_13 = pmat(0,2).real();
            pp.s_23 = pmat(1,2).real();
            
            jacD = jac(mass(0),mass,rr,in,out);
            /*
             for(int kk=0; kk!=f_G.size();++kk){
             if(f_G[kk].name == pname){
             decay += f_G[kk](pp)*jacD;}
             }*/
            decay += f_G[pid](pp)*jacD;
            
        }
    }
    res = decay;
    return res;
}

double gamma(int &pid,double const &T,param_t &data)
{
    double result = 0;
    
    double gg = 1.;
    vector< pair <int,double> > vect;
    vector< pair <int,double> > gvect;
    size_t found;
    
    std::string pname = f_G[pid].name;
    
    for(int pt = 0;pt!=ptl.size();++pt){
        found = pname.find(ptl[pt]);
        while (found!=std::string::npos) {
            vect.push_back( make_pair(found,prtlm(ptl[pt],data)) );
            found = pname.find(ptl[pt],found+1);
        }
    }
    
    sort(vect.begin(), vect.end());
    int in = 0,out = 0;
    
    for(int tt = 0;tt!=vect.size();++tt){
        (vect[tt].first<pname.find("to")) ? in++ : out++;
        
    }
    
    for(auto &gi : gptl){
        found = pname.find(gi.first);
        while (found<pname.find("to")) {
            gg *= gi.second;
            found = pname.find(gi.first,found+1);
        }
    }
    
    Eigen::MatrixXcd pmat(in+out,in+out);
    Eigen::VectorXcd mass(in+out);
    
    
    vector<double> mmax;
    for(int tt = 0;tt!=vect.size();++tt){
        mass(tt) = vect[tt].second;
        mmax.push_back(vect[tt].second);
    }
    
    if((in==2)&&(out==2)){
        
        double s,pin,pout;
        
        double ma,mb,m1,m2,ms,Tsc;
        
        ms = *max_element(mmax.begin(), mmax.end());
        Tsc = T/ms;
        ma = abs(mass(0)/ms);
        mb = abs(mass(1)/ms);
        m1 = abs(mass(2)/ms);
        m2 = abs(mass(3)/ms);
        int ii,jj;
        double jac = 0.,mi=0.,mj=0.;
        complex_t amp = 0.;
        //cout << " ma = " << ma << " mb = " << mb << " m1 = " << m1 << " m2 = " << m2 << endl << endl;
        if((pow(ma+mb, 2.))>=(pow(m1+m2, 2.)))
        {
            mi=ma;
            mj=mb;
            
        }
        if((pow(ma+mb, 2.))<(pow(m1+m2, 2.)))
        {
            mi=m1;
            mj=m2;
        }
        for (ii=0; ii<12; ii++)
        {
            for (jj=0; jj<13; jj++)
            {
                s = pow(mi-mj,2.) + 2.*mi*mj*(2. + xki[ii]);//xs1[ii](1.- xs1[ii]);
                pin = 0.5*sqrt(L(s,ma*ma,mb*mb)/s).real();
                pout = 0.5*sqrt(L(s,m1*m1,m2*m2)/s).real();
                
                vector<double> cth,rr;
                cth.push_back(xs[jj]);
                rr.push_back(1.);
                pmat = pmom(sqrt(s)*ms,mass,cth,rr,in,out);
                
                data.s_12 = pmat(0,1).real();
                data.s_13 = pmat(0,2).real();
                data.s_14 = pmat(0,3).real();
                data.s_23 = pmat(1,2).real();
                data.s_24 = pmat(1,3).real();
                data.s_34 = pmat(2,3).real();
                
                jac = gg*pow(ms,4.)*2.*mi*mj*aki[ii]*as[jj];//*(1./(pow(xs1[ii] - 1.,2.)));
                amp = 0.;
                /*
                 for(int kk=0; kk!=f_G.size();++kk){
                 if(f_G[kk].name == pname){
                 amp += f_G[kk](data)*jac;}
                 }*/
                amp += f_G[pid](data)*jac;
                result += exp(xki[ii])*Tsc/(256.*pow(M_PI,5.))*pin*pout*amp.real()*1./sqrt(s)*jac*boost::math::cyl_bessel_k(1, sqrt(s)/Tsc);
                
            }
        }
    }
    
    if(in==1)
    {
        complex_t decay = Decay(pid,data);
        
        double mD = mass(0).real();
        result = gg*mD*mD/(2.*pow(M_PI,2.))*T*boost::math::cyl_bessel_k(1, abs(mD)/T)*decay.real();
    }
    //cout << "\n  Res = " << result << endl;
    return result;
}

void initialize()
{
    int pid;
    double charge,chargeB;
    double ch1,ch2,ch3,ch4;
    double chB1,chB2,chB3,chB4;
    vector<std::string> pX;
    
    for(size_t ii = 0;ii!=pXSM.size();++ii){
        for(size_t jj = 0;jj!=pXSM.size();++jj){
            
            /**/
            XXSMSM[ii][jj].prid.clear();
            BXXSMSM[ii][jj].prid.clear();
            AntiBXXSMSM[ii][jj].prid.clear();
            
            XSMXSM[ii][jj].prid.clear();
            BXSMXSM[ii][jj].prid.clear();
            AntiBXSMXSM[ii][jj].prid.clear();
            
            for(size_t kk = 0;kk!=pthSM.size();++kk){
                for(size_t ll = 0;ll!=pthSM.size();++ll){
                    
                    //cout << "XX["<< kk <<"]["<<ll<<"]" << endl;
                    
                    //cout << f_G[pid].name << "1 : " << pid << endl;
                    //================== FFFF
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pXSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pXSM[jj].first);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    
                    //cout << f_G[pid].name << ">> : " << Prname(pX,2,2) << endl;
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[kk].first);
                    ch3 = ptlSM(pXSM[jj].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[kk].first);
                    chB3 = ptlB(pXSM[jj].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== FSFF
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pXSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pXSM[jj].second);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[kk].second);
                    ch3 = ptlSM(pXSM[jj].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[kk].second);
                    chB3 = ptlB(pXSM[jj].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    
                    //================== SFFF
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pXSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pXSM[jj].first);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[kk].first);
                    ch3 = ptlSM(pXSM[jj].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[kk].first);
                    chB3 = ptlB(pXSM[jj].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== SSFF
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pXSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pXSM[jj].second);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[kk].second);
                    ch3 = ptlSM(pXSM[jj].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[kk].second);
                    chB3 = ptlB(pXSM[jj].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    
                    //================== FFFS
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pXSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pXSM[jj].first);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[kk].first);
                    ch3 = ptlSM(pXSM[jj].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[kk].first);
                    chB3 = ptlB(pXSM[jj].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== FSFS
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pXSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pXSM[jj].second);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[kk].second);
                    ch3 = ptlSM(pXSM[jj].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[kk].second);
                    chB3 = ptlB(pXSM[jj].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    
                    //================== SFFS
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pXSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pXSM[jj].first);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[kk].first);
                    ch3 = ptlSM(pXSM[jj].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[kk].first);
                    chB3 = ptlB(pXSM[jj].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== SSFS
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pXSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pXSM[jj].second);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[kk].second);
                    ch3 = ptlSM(pXSM[jj].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[kk].second);
                    chB3 = ptlB(pXSM[jj].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[ll].second);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== FFSF
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pXSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pXSM[jj].first);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[kk].first);
                    ch3 = ptlSM(pXSM[jj].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[kk].first);
                    chB3 = ptlB(pXSM[jj].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== FSSF
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pXSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pXSM[jj].second);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[kk].second);
                    ch3 = ptlSM(pXSM[jj].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[kk].second);
                    chB3 = ptlB(pXSM[jj].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(charge == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    
                    //================== SFSF
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pXSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pXSM[jj].first);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[kk].first);
                    ch3 = ptlSM(pXSM[jj].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[kk].first);
                    chB3 = ptlB(pXSM[jj].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== SSSF
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pXSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pXSM[jj].second);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[kk].second);
                    ch3 = ptlSM(pXSM[jj].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[kk].second);
                    chB3 = ptlB(pXSM[jj].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[ll].first);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== FFSS
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pXSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pXSM[jj].first);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[kk].first);
                    ch3 = ptlSM(pXSM[jj].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[kk].first);
                    chB3 = ptlB(pXSM[jj].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== FSSS
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pXSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pXSM[jj].second);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[kk].second);
                    ch3 = ptlSM(pXSM[jj].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[kk].second);
                    chB3 = ptlB(pXSM[jj].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    
                    //================== SFSS
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pXSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pXSM[jj].first);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[kk].first);
                    ch3 = ptlSM(pXSM[jj].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[kk].first);
                    chB3 = ptlB(pXSM[jj].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== SSSS
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pXSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pXSM[jj].second);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXXSMSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[kk].second);
                    ch3 = ptlSM(pXSM[jj].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[kk].second);
                    chB3 = ptlB(pXSM[jj].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pXSM[jj].second);
                    pX.push_back(pthSM[ll].second);
                    
                    findAmp(Prname(pX,2,2),pid);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMXSM[ii][jj].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    
                }
            }
            
        }
    }
    
    
    for(size_t ii = 0;ii!=pXSM.size();++ii){
        XSMSMSM[ii].prid.clear();
        BXSMSMSM[ii].prid.clear();
        AntiBXSMSMSM[ii].prid.clear();
        
        for(size_t jj = 0;jj!=pthSM.size();++jj){
            for(size_t kk = 0;kk!=pthSM.size();++kk){
                for(size_t ll = 0;ll!=pthSM.size();++ll){
                    
                    
                    //================== FFFF
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[jj].first);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[jj].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].first);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== FSFF
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[jj].second);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[jj].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].first);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    pX.clear();
                    
                    //================== SFFF
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[jj].first);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[jj].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].first);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== SSFF
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[jj].second);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[jj].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].first);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== FFFS
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[jj].first);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[jj].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].second);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    //================== FSFS
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[jj].second);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[jj].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].second);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== SFFS
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[jj].first);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[jj].first);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].second);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== SSFS
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].first);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[jj].second);
                    chB3 = ptlB(pthSM[kk].first);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[jj].second);
                    pX.push_back(pthSM[kk].first);
                    pX.push_back(pthSM[ll].second);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== FFSF
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[jj].first);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pXSM[jj].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].first);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== FSSF
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[jj].second);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[jj].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].first);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== SFSF
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[jj].first);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[jj].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].first);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    //================== SSSF
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].first);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[jj].second);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].first);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[jj].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].first);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== FFSS
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[jj].first);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[jj].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].second);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                    //================== FSSS
                    ch1 = ptlSM(pXSM[ii].first);
                    ch2 = ptlSM(pthSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].first);
                    chB2 = ptlB(pthSM[jj].second);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].first);
                    pX.push_back(pthSM[jj].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].second);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    
                    //================== SFSS
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[jj].first);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[jj].first);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[jj].first);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].second);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    pX.clear();
                    
                    //================== SSSS
                    ch1 = ptlSM(pXSM[ii].second);
                    ch2 = ptlSM(pthSM[jj].second);
                    ch3 = ptlSM(pthSM[kk].second);
                    ch4 = ptlSM(pthSM[ll].second);
                    
                    chB1 = ptlB(pXSM[ii].second);
                    chB2 = ptlB(pthSM[jj].second);
                    chB3 = ptlB(pthSM[kk].second);
                    chB4 = ptlB(pthSM[ll].second);
                    
                    pX.push_back(pXSM[ii].second);
                    pX.push_back(pthSM[jj].second);
                    pX.push_back(pthSM[kk].second);
                    pX.push_back(pthSM[ll].second);
                    
                    charge = ch1 + ch2 - ch3 - ch4;
                    chargeB = chB1 + chB2 - chB3 - chB4;
                    findAmp(Prname(pX,2,2),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    charge = ch1 - ch2 - ch3 - ch4;
                    chargeB = chB1 - chB2 - chB3 - chB4;
                    findAmp(Prname(pX,1,3),pid);
                    if((charge == 0)&&(pid< ((int) f_G.size())))
                    {
                        XSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    
                    if((charge == 0)&&(chargeB == 1)&&(pid< ((int) f_G.size())))
                    {
                        BXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    if((charge == 0)&&(chargeB == -1)&&(pid< ((int) f_G.size())))
                    {
                        AntiBXSMSMSM[ii].prid.push_back(pid);
                        //cout << f_G[pid].name << " : " << pid << endl;
                    }
                    
                    pX.clear();
                    
                }
            }
            
        }
    }
    
    
}


void functh(param_t &pp,double hh,double zl,Eigen::VectorXcd &yn)
{
    double xx = pow(10., zl);
    
    double T = m_XX(pp).real()/xx;
    double const &mX = m_XX(pp).real();
    
    size_t nn = pXSM.size() + 1;
    Eigen::VectorXcd tempf(nn),bth(nn),yth(nn),ypth(nn),yeq(nn);
    
    Eigen::MatrixXcd tempfp(nn,nn),athinv(nn,nn),exath(nn,nn);
    Eigen::MatrixXd ath(nn,nn);
    
    Eigen::MatrixXcd Gam2(nn,nn);
    Eigen::VectorXcd Gam1(nn);
    
    tempf = Eigen::VectorXcd::Zero(nn);
    tempfp = Eigen::MatrixXcd::Zero(nn,nn);
    
    yth = yn;
    
    for(size_t pi = 0;pi!=pXSM.size();++pi){
        double gg = 1.;
        for(auto &gi : gptl){
            gg *= (gi.first==pXSM[pi].first) ? gi.second : 1.;
        }
        yeq(pi) = Yeq(gg,prtlm(pXSM[pi].first,pp),T);
        cout << "\tΔ["<< pi << "] = " << yth(pi)/yeq(pi) - 1.;
    }
    cout << endl;
    yeq(nn-1) = Yeq(2.,0.,T);
    
    double Hb = sqrt(4.*pow(M_PI,3.)*geffT(T)/45.)*pow(T,2.)/mpl;
    
    double s = geffT(T)*2*pow(M_PI,2.)/45.*pow(T,3.);
    
    complex<double> Xsec = 0., Dec = 0.,pr = 0.;
    
    Gam1 = Eigen::VectorXcd::Zero(nn);
    Gam2 = Eigen::MatrixXcd::Zero(nn,nn);
    
    for(size_t ii = 0;ii!=pXSM.size();++ii){
        for(size_t jj = 0;jj!=pXSM.size();++jj){
            Xsec = 0.;
            auto &pXX = XXSMSM[ii][jj].prid;
            for(auto &pid : pXX)
            {
                pr = (pid < ((int) f_G.size())) ? gamma(pid,T,pp)/(s*Hb*xx) : 0.;
                Xsec += (pr.real()>0.) ? pr : 0.;
                //cout << endl << "Pr name1 = "<< f_G[pid].name << ": " << pr << endl;
            }
            
            Gam2(ii,jj) += Xsec;
            if (ii<jj) {
                Gam2(jj,ii) += Xsec;
            }
            tempf(ii) += -log(10)*xx*((yth(ii)*yth(jj)/(yeq(ii)*yeq(jj)) - 1.)*Gam2(ii,jj));
        }
        
    }
    
    
    for(size_t ii = 0;ii != pXSM.size();++ii){
        for(size_t kk = 0;kk!=pXSM.size();++kk){
            for(size_t jj = 0;jj!=pXSM.size();++jj){
                double del_ki = (kk==ii) ? 1. : 0;
                double del_kj = (kk==jj) ? 1. : 0;
                tempfp(ii,kk) += -log(10)*xx*((del_ki*yth(jj)/(yeq(ii)*yeq(jj)) + del_kj*yth(ii)/(yeq(ii)*yeq(jj)))*Gam2(ii,jj));
            }
        }
    }
    
    Xsec = 0.;
    Gam2 = Eigen::MatrixXcd::Zero(nn,nn);
    
    for(size_t ii = 0;ii!=pXSM.size();++ii){
        for(size_t jj = 0;jj!=pXSM.size();++jj){
            Xsec = 0.;
            if(ii<=jj){
                
                auto &pXX = XSMXSM[ii][jj].prid;
                for(auto &pid : pXX)
                {
                    pr = (pid < ((int) f_G.size())) ? gamma(pid,T,pp)/(s*Hb*xx) : 0.;
                    Xsec += (pr.real()>0.) ? pr : 0.;
                    //cout << endl << "Pr name1 = "<< f_G[pid].name << ": " << pr << endl;
                }
                
                Gam2(ii,jj) += Xsec;
                
                Gam2(jj,ii) += (ii<jj) ? Xsec : 0.;
                
                tempf(ii) += -log(10)*xx*((yth(ii)/yeq(ii) - yth(jj)/yeq(jj))*Gam2(ii,jj));
            }
        }
    }
    
    for(size_t ii = 0;ii != pXSM.size();++ii){
        for(size_t kk = 0;kk!=pXSM.size();++kk){
            for(size_t jj = 0;jj!=pXSM.size();++jj){
                double del_ki = (kk==ii) ? 1. : 0;
                double del_kj = (kk==jj) ? 1. : 0;
                tempfp(ii,kk) += -log(10)*xx*((del_ki/yeq(ii) - del_kj/yeq(jj))*Gam2(ii,jj));
            }
        }
    }
    //cout << "Gam2 :" << endl << Gam2 << endl << endl;
    Gam1 = Eigen::VectorXcd::Zero(nn);
    
    for(size_t ii = 0;ii!=pXSM.size();++ii){
        Xsec = 0.;
        
        auto &pXX = XSMSMSM[ii].prid;
        for(auto &pid : pXX)
        {
            pr = (pid < ((int) f_G.size())) ? gamma(pid,T,pp)/(s*Hb*xx) : 0.;
            Xsec += (pr.real()>0.) ? pr : 0.;
            //cout << endl << "Pr name1 = "<< f_G[pid].name << ": " << pr << endl;
        }
        
        Gam1(ii) += Xsec;
        tempf(ii) += -log(10)*xx*((yth(ii)/yeq(ii) - 1.)*Gam1(ii));
        
        for(size_t kk = 0;kk!=pXSM.size();++kk){
            double del_ki = (kk==ii) ? 1. : 0;
            tempfp(ii,kk) += -log(10)*xx*((del_ki/yeq(ii))*Gam1(ii));
        }
    }
    if(nn > pXSM.size()){
        //cout << "Gam1 :" << endl << Gam1.transpose() << endl << endl;
        Gam1 = Eigen::VectorXcd::Zero(nn);
        
        //cout << "Gam1 :" << endl << Gam1.transpose() << endl << endl;
        /**/
        complex_t XsecB = 0.,XsecAB = 0.;
        complex_t DecB = 0.,DecAB = 0.;
        
        Eigen::MatrixXcd GCPOdd1(nn,nn),GCPOdd2(nn,nn),GCPOdd3(nn,nn);
        Eigen::MatrixXcd GCPEven1a(nn,nn),GCPEven2a(nn,nn);
        Eigen::MatrixXcd GCPEven1b(nn,nn),GCPEven2b(nn,nn);
        
        GCPOdd1 = Eigen::MatrixXcd::Zero(nn,nn);
        GCPOdd2 = Eigen::MatrixXcd::Zero(nn,nn);
        GCPEven1a = Eigen::MatrixXcd::Zero(nn,nn);
        GCPEven2a = Eigen::MatrixXcd::Zero(nn,nn);
        GCPEven1b = Eigen::MatrixXcd::Zero(nn,nn);
        GCPEven2b = Eigen::MatrixXcd::Zero(nn,nn);
        
        complex_t CPOdd;
        pair <complex_t,complex_t> CPEven;
        
        for(size_t ii = 0;ii!=pXSM.size();++ii){
            CPOdd = 0.;
            CPEven = make_pair(0.,0.);
            pair <double,double> fB = make_pair(0.,0.);
            pair <double,double> fAB = make_pair(0.,0.);
            
            auto &pXB = BXSMSMSM[ii].prid;
            for(auto &pid : pXB)
            {
                pr = (pid < ((int) f_G.size())) ? gamma(pid,T,pp)/(s*Hb*xx) : 0.;
                XsecB += (pr.real()>0.) ? pr : 0.;
                //cout << endl << "Pr name1 = "<< f_G[pid].name << ": " << pr << endl;
                facB(pid,fB,T);
            }
            
            auto &pXAB = AntiBXSMSMSM[ii].prid;
            for(auto &pid : pXAB)
            {
                pr = (pid < ((int) f_G.size())) ? gamma(pid,T,pp)/(s*Hb*xx) : 0.;
                XsecAB += (pr.real()>0.) ? pr : 0.;
                facB(pid,fAB,T);
                
            }
            
            CPOdd += XsecB-XsecAB;
            CPEven.first += fB.first*XsecB + fAB.first*XsecAB;
            CPEven.second += fB.second*XsecB + fAB.second*XsecAB;
            
            cout << " fB.first : " << fB.first << " fAB.first : " << fAB.first << endl;
            cout << " fB.second : " << fB.second << " fAB.second : " << fAB.second << endl;
            cout << "XsecB : " << XsecB;
            cout << "\t XsecAB : " << XsecAB << endl;
            cout << "GCPEvena : " << CPEven.first;
            cout << "\t GCPEvenb : " << CPEven.second << endl;
            
            
            tempf(nn-1) += log(10)*xx*( (yth(ii)/yeq(ii) - 1.)*CPOdd - yth(nn-1)/yeq(nn-1)*(CPEven.first*yth(ii)/yeq(ii) + CPEven.second));
            tempfp(nn-1,nn-1) += -log(10)*xx/yeq(nn-1)*((CPEven.first*yth(ii)/yeq(ii) + CPEven.second));
            tempfp(nn-1,ii) += log(10)*xx*(CPOdd/yeq(ii) - yth(nn-1)/yeq(nn-1)*(CPEven.first/yeq(ii)));
        }
        
        
        for(size_t ii= 0;ii!=pXSM.size();++ii){
            for(size_t jj= 0;jj!=pXSM.size();++jj){
                if(ii<=jj){
                    XsecB = 0.;
                    XsecAB = 0.;
                    
                    pair <double,double> fB = make_pair(0.,0.);
                    pair <double,double> fAB = make_pair(0.,0.);
                    
                    auto &pXXB = BXXSMSM[ii][jj].prid;
                    for(auto &pid : pXXB){
                        pr = (pid < ((int) f_G.size())) ? gamma(pid,T,pp)/(s*Hb*xx) : 0.;
                        XsecB += (pr.real()>0.) ? pr : 0.;
                        //cout << endl << "Pr name1 = "<< f_G[pid].name << ": " << pr << endl;
                        facB(pid,fB,T);
                    }
                    
                    auto &pXXAB = AntiBXXSMSM[ii][jj].prid;
                    for(auto &pid : pXXAB){
                        pr = (pid < ((int) f_G.size())) ? gamma(pid,T,pp)/(s*Hb*xx) : 0.;
                        XsecAB += (pr.real()>0.) ? pr : 0.;
                        //cout << endl << "Pr name1 = "<< f_G[pid].name << ": " << pr << endl;
                        facB(pid,fAB,T);
                    }
                    
                    GCPOdd1(ii,jj) += XsecB-XsecAB;
                    GCPEven1a(ii,jj) += fB.first*XsecB+fAB.first*XsecAB;
                    GCPEven1b(ii,jj) += fB.second*XsecB+fAB.second*XsecAB;
                    
                    if (ii<jj) {
                        GCPOdd1(jj,ii) += XsecB-XsecAB;
                        GCPEven1a(jj,ii) += fB.first*XsecB+fAB.first*XsecAB;
                        GCPEven1b(jj,ii) += fB.second*XsecB+fAB.second*XsecAB;
                    }
                    
                    XsecB = 0.;
                    XsecAB = 0.;
                    fB = make_pair(0.,0.);
                    fAB = make_pair(0.,0.);
                    
                    auto &pXSMB = BXSMXSM[ii][jj].prid;
                    for(auto &pid : pXSMB){
                        pr = (pid< ((int) f_G.size())) ? gamma(pid,T,pp)/(s*Hb*xx) : 0.;
                        XsecB += (pr.real()>0.) ? pr : 0.;
                        //cout << endl << "Pr name1 = "<< f_G[pid].name << ": " << pr << endl;
                        facB(pid,fB,T);
                    }
                    
                    auto &pXSMAB = AntiBXSMXSM[ii][jj].prid;
                    
                    for(auto &pid : pXSMAB){
                        pr = (pid < ((int) f_G.size())) ? gamma(pid,T,pp)/(s*Hb*xx) : 0.;
                        XsecAB += (pr.real()>0.) ? pr : 0.;
                        //cout << endl << "Pr name1 = "<< f_G[pid].name << ": " << pr << endl;
                        facB(pid,fAB,T);}
                    
                    
                    GCPOdd2(ii,jj) += XsecB-XsecAB;
                    GCPEven2a(ii,jj) += fB.first*XsecB+fAB.first*XsecAB;
                    GCPEven2b(ii,jj) += fB.second*XsecB+fAB.second*XsecAB;
                    
                    if(ii<jj){
                        GCPOdd2(jj,ii) += XsecB-XsecAB;
                        GCPEven2a(jj,ii) += fB.first*XsecB+fAB.first*XsecAB;
                        GCPEven2b(jj,ii) += fB.second*XsecB+fAB.second*XsecAB;
                    }
                   
                }
                
                
                tempf(nn-1) += log(10)*xx*( (yth(ii)*yth(jj)/(yeq(ii)*yeq(jj)) - 1.)*GCPOdd1(ii,jj) - yth(nn-1)/yeq(nn-1)*(GCPEven1a(ii,jj)*yth(ii)*yth(jj)/(yeq(ii)*yeq(jj)) + GCPEven1b(ii,jj)));
                
                tempfp(nn-1,nn-1) += -log(10)*xx/yeq(nn-1)*((GCPEven1a(ii,jj)*yth(ii)*yth(jj)/(yeq(ii)*yeq(jj)) + GCPEven1b(ii,jj)));
                
                tempf(nn-1) += log(10)*xx*( (yth(ii)/yeq(ii) - yth(jj)/yeq(jj))*GCPOdd2(ii,jj) - yth(nn-1)/yeq(nn-1)*(GCPEven2a(ii,jj)*yth(ii)/yeq(ii) + GCPEven2b(ii,jj)*yth(jj)/yeq(jj)));
                
                tempfp(nn-1,nn-1) += -log(10)*xx/yeq(nn-1)*(GCPEven2a(ii,jj)*yth(ii)/yeq(ii) + GCPEven2b(ii,jj)*yth(jj)/yeq(jj));
                
            }
        }
        for(size_t kk = 0; kk!=pXSM.size(); ++kk){
            for(size_t ii= 0;ii!=pXSM.size();++ii){
                for(size_t jj= 0;jj!=pXSM.size();++jj){
                    
                    double del_ki = (kk==ii) ? 1. : 0.;
                    double del_kj = (kk==jj) ? 1. : 0.;
                    tempfp(nn-1,kk) += log(10)*xx*((del_ki*yth(jj)/(yeq(ii)*yeq(jj)) + yth(ii)*del_kj/(yeq(ii)*yeq(jj)))*GCPOdd1(ii,jj) - yth(nn-1)/yeq(nn-1)*(GCPEven1a(ii,jj)*del_ki*yth(jj)/(yeq(ii)*yeq(jj)) + GCPEven1a(ii,jj)*yth(ii)*del_kj/(yeq(ii)*yeq(jj))));
                    
                    tempfp(nn-1,kk) += log(10)*xx*( (del_ki/yeq(ii) - del_kj/yeq(jj))*GCPOdd2(ii,jj) - yth(nn-1)/yeq(nn-1)*(GCPEven2a(ii,jj)*del_ki/yeq(ii) + GCPEven2b(ii,jj)*del_kj/yeq(jj)));
                    
                }
            }
            
        }
    }
    ath = tempfp.real();
    bth = tempf - tempfp*yth;
    
    int dim = pXSM.size() + 1;
    
    if((ath.determinant()!=0)&&(!isnan(ath.determinant()))&&(!isinf(ath.determinant())))
    {
        ypth = expM(ath*hh, dim)*yth + tempfp.inverse()*(expM(ath*hh, dim) -  MatrixXd::Identity(dim, dim))*bth;
    }
    else
    {
        ypth = yth;
    }
    
    yn = ypth;
}


void diff(std::ofstream &file,double zmax,double zin,Eigen::VectorXcd &yn,param_t &p)
{
    int i,iter = MAXSTEPS,nn;
    
    double h = (zmax - zin)/((double) MAXSTEPS);
    double z = zin;
    //int tot = ((int) pXSM.size()) + 1;
    
    double T = m_XX(p).real()/pow(10., z);
    
    size_t tot = pXSM.size() + 1;
    
    cout << "Tot =" << tot << endl;
    Eigen::VectorXcd yeq(tot);
    
    for(size_t pi = 0;pi!=pXSM.size();++pi){
        double gg = 1.;
        for(auto &gi : gptl){
            gg *= (gi.first==pXSM[pi].first) ? gi.second : 1.;
        }
        yeq(pi) = Yeq(gg,prtlm(pXSM[pi].first,p),T);
    }
    yeq(tot-1) = Yeq(2.,0.,T);
    
    file << T << " ";
    cout<< T << " ";
    for (size_t yi = 0; yi!=tot; ++yi) {
        file << yn(yi).real() << " " << yeq(yi).real() << " ";
        cout<< sbyrho*yn(yi).real() << " ";
    }
    file << endl;
    cout << endl;
    
    for(i=0;i<iter;i++)
    {
        T = m_XX(p).real()/pow(10., z);
        
        for(size_t pi = 0;pi!=pXSM.size();++pi){
            double gg = 1.;
            for(auto &gi : gptl){
                gg *= (gi.first==pXSM[pi].first) ? gi.second : 1.;
            }
            yeq(pi) = Yeq(gg,prtlm(pXSM[pi].first,p),T);
        }
        yeq(tot-1) = Yeq(2.,0.,T);
        
        functh(p,h,z,yn);  //printf("\n");
        z = z + h;
        
        file << T << " ";
        cout<< T << " ";
        for (size_t yi = 0; yi!=tot; ++yi) {
            file << yn(yi).real() << " " << yeq(yi).real() << " ";
            cout<< sbyrho*yn(yi).real() << " ";
        }
        file << endl;
        cout << endl;
        /*
         cout << pow(10., z[i]) << "\t" << m_XX(p)/pow(10., z[i]) <<"\t"<<sqrt(p.muS2.real())*sbyrho*Yeq(1., sqrt(p.muS2.real()), T)<<"\t"<<sqrt(p.muS2.real())*sbyrho*(yn[0][i]);
         cout <<"\t"<<p.m_X.real()*sbyrho*Yeq(2., p.m_X.real(), T) << "\t" << (p.m_X.real())*sbyrho*(yn[1][i]) << " " << 2./3.*sbyrho*yn[2][i] << endl;
         cout << "Ωh^2 =  " << p.m_X.real()*sbyrho*(yn[1][i]) + sqrt(p.muS2.real())*sbyrho*(yn[0][i])<< endl;
         
         cout << pow(10., z[i]) << "\t" << p.m_X.real()/pow(10., z[i]) <<"\t"<<sqrt(p.muS2.real())*sbyrho*Yeq(1., sqrt(p.muS2.real()), T)<<"\t"<<sqrt(p.muS2.real())*sbyrho*(yn[0][i])*Yeq(1., sqrt(p.muS2.real()),T);
         cout <<"\t"<<p.m_X.real()*sbyrho*Yeq(2., p.m_X.real(), T) << "\t" << (p.m_X.real())*sbyrho*(yn[1][i])*Yeq(2., p.m_X.real(), T) << " " << 2./3.*sbyrho*yn[2][i] << endl;
         cout << "Ωh^2 =  " << p.m_X.real()*sbyrho*(yn[1][i])*Yeq(2., p.m_X.real(), T) + sqrt(p.muS2.real())*sbyrho*(yn[0][i])*Yeq(1., sqrt(p.muS2.real()),T) << endl;*/
        
    }
}
