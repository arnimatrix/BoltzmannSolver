#include "BE.h"


#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>

using namespace std;

using Eigen::MatrixXd;
using namespace Eigen;


using namespace mty::lib;




namespace leptoproto{


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

//=====================================================================

complex_t Complex(int const &x,double const &y)
{
    return x + y*1i;
}
/*
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
*/

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

double Massp(std::string name,param_t &pp)
{
    double mi;
    std::string ptl_mass = "m_";
    
    bool cm0 = false;
    for(const auto &m0 : Massless){
        if((name==m0)||(name==m0+"c"))
        {
            cm0 = true;
            
        }
    }
    if(cm0)
    {
        mi = 0.;
    }
    else{
        auto pt_mass = pp.realParams[ptl_mass + name];
        mi = *pt_mass;
        
    }
    
    return mi;
}

double Yeq(std::string name,param_t &pp,double T)
{
    double yy = 0.;
    double mi = Massp(name,pp);
    double gi = 1.;
    for(const auto &ga : gptl){
        gi *= (name==ga.first) ? ga.second : 1.;
        //std::cout << "name = " << ga.first << std::endl;
    }
    if(mi>0.)
    {
        yy = 45.*gi/geffT(T)*pow(mi/(2.*M_PI*T),2.)*cyl_bessel_k(2, mi/T);
    }
    if(mi==0.)
    {
        yy = 45./2.*gi/geffT(T)*pow(M_PI,-4.); //cout << "h1!!" << endl;
    }
    return yy;
}
/*

double nthermalD(vector<pair <std::string,double>> sg,std::string str)
{
    double decay = 0.;
    
    for(auto &data : sg){
        if(data.first==str) {
            decay = data.second;
        }
    }
    //cout << "Decay nth : =" << decay << endl;
    return decay;
}
*/




Eigen::MatrixXcd expM(Eigen::MatrixXd M,int &dim)
{
    Eigen::MatrixXcd res,ev;
    //Eigen::MatrixXcd evt;
    
    Eigen::EigenSolver<MatrixXd> es(M);
    
    //cout << "es : " << es.eigenvalues().transpose() << endl << endl;
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
        
        res = ((mD - m1 - m2 - m3) > 0.) ? 0.5/mD*1./(4.*M_PI)*sqrt(L(mD*mD,m1*m1,s23).real())/(8.*M_PI*mD*mD)*sqrt(L(s23,m2*m2,m3*m3)).real()/(8.*M_PI*s23) : 0.;
    }
    if((intp==1)&&(outp==2)){
        double mD = mm(0).real();
        double m1 = mm(1).real();
        double m2 = mm(2).real();
        
        res = ((mD - m1 - m2) > 0.) ? 1./(16.*M_PI*pow(mD,3.))*sqrt(L(mD*mD,m1*m1,m2*m2)).real() : 0.;
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


Eigen::MatrixXcd pmom(complex_t sqss,Eigen::VectorXcd &mm,vector<double> const &cthz,vector<double> const &rr,int const &inp, int const &outp)
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

/*
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
        ch1 = ptlSM(pnthSM[qq]);
        std::string squark = pnthSM[qq];
        for(size_t qi = 0;qi!=pthSM.size();++qi){
            for(size_t qj = 0;qj!=pthSM.size();++qj){
                
                ch2 = ptlSM(pthSM[qi]);
                ch3 = ptlSM(pthSM[qj]);
                
                f1 = pthSM[qi];
                f2 = pthSM[qj];
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn1;
                    pn1.push_back(squark);
                    pn1.push_back(f1);
                    pn1.push_back(f2);
                    //cout << endl << "Pr name = "<< Prname(pn1,1,2) << endl;
                    findAmp(Prname(pn1,1,2),pid);
                    decay += (pid<((int) f_G.size())) ? Decay(pid,pp).real() : 0.;
                    
                }
            }
            
            for(size_t qj = 0;qj!=pXSM.size();++qj){
                
                ch2 = ptlSM(pthSM[qi]);
                ch3 = ptlSM(pXSM[qj]);
                
                f1 = pthSM[qi];
                f2 = pXSM[qj];
                
                charge = ch1 - ch2 - ch3;
                if (charge == 0) {
                    vector<std::string> pn2;
                    pn2.push_back(squark);
                    pn2.push_back(f1);
                    pn2.push_back(f2);
                    //cout << endl << "Pr name = "<< Prname(pn2,1,2) << endl;
                    findAmp(Prname(pn2,1,2),pid);
                    decay += (pid<((int)f_G.size())) ? Decay(pid,pp).real() : 0.;
                }
                
            }
        }
        sqDecay.push_back(std::make_pair(squark,decay));
        //cout << "Decay of sq[ "<< "]"<< qq <<" = "<< squark << " = " << decay << endl;
    }
    
    
}
 */

double epsilon1(std::string ll,std::string ptl,param_t &pp)
{
    std::vector<Process> procL;
    ParticleData pdata;
    pdata.loadFile("script/test.json");
    procL = pdata.getProcessesFromQNumberViolation(ll);
    std::string ptl_mass = "m_";
    complex<double> ampLp = 0.;
    complex<double> ampLap = 0.;
    
    complex<double> tree = 0.; // Only for testing
    for(const auto &prc : procL){
        int in_proc = prc.inParticles.size();
        int out_proc = prc.outParticles.size();
        double DeltaL=prc.qNumbers.begin()->second;
        if((in_proc<2)&&(out_proc>1)){
            if(prc.inParticles[0].name.find(ptl)!=std::string::npos){
                
                complex_t decay = 0.;
                Eigen::MatrixXcd pmat(1+out_proc,1+out_proc);
                Eigen::VectorXcd mass(1+out_proc);
                
                vector<double> mmax;
                size_t tt_max = prc.inParticles.size() + prc.outParticles.size();
                size_t tt = 0;
                bool cm0 = false;
                for(const auto &prname : prc.inParticles){
                    for(const auto &m0 : Massless){
                        if((prname.name.find(m0)!=std::string::npos))
                        {
                            cm0 = true;
                        }
                    }
                    if(cm0)
                    {
                        //std::cout << "Mass of Photon "<< prname.name << " = 0" << std::endl;
                        mass(tt) = 0.;
                    }
                    else{
                        auto pt_mass = pp.realParams[ptl_mass + prname.name];
                        //std::cout << "Mass of "<< prname.name << " = " << *pt_mass << std::endl;
                        mmax.push_back(*pt_mass);
                        mass(tt) = *pt_mass;
                        
                    }
                    ++tt;
                }
                
                for(const auto &prname : prc.outParticles){
                    for(const auto &m0 : Massless){
                        if((prname.name.find(m0)!=std::string::npos))
                        {
                            cm0 = true;
                        }
                    }
                    if(cm0)
                    {
                        //std::cout << "Mass of Photon "<< prname.name << " = 0" << std::endl;
                        mass(tt) = 0.;
                    }
                    else{
                        auto pt_mass = pp.realParams[ptl_mass + prname.name];
                        //std::cout << "Mass of "<< prname.name << " = " << *pt_mass << std::endl;
                        mmax.push_back(*pt_mass);
                        mass(tt) = *pt_mass;
                        
                    }
                    ++tt;
                }
                
                vector<double> cth,rr;
                cth.push_back(1.);
                rr.push_back(1.);
                
                pmat = pmom(mass(0),mass,cth,rr,1,out_proc);
                pp.s_12 = pmat(0,1).real();
                pp.s_13 = pmat(0,2).real();
                pp.s_23 = pmat(1,2).real();
                int pid;
                if(prc.name.find("AsymL")!=std::string::npos){
                    findAmp(prc.name,pid);
                    if(DeltaL==1){
                        //std::cout << " Process = " << prc.name << std::endl;
                        //std::cout << " Qnumber  = " << prc.qNumbers.begin()->second << std::endl;
                        ampLp += f_G[pid](pp);
                        //std::cout << "Ampl Loop P = " << ampLp << std::endl;
                    }
                }
                if(prc.name.find("AsymT")!=std::string::npos){
                    findAmp(prc.name,pid);
                    if(DeltaL==1){
                        //std::cout << " Process = " << prc.name << std::endl;
                        //std::cout << " Qnumber  = " << prc.qNumbers.begin()->second << std::endl;
                        tree += f_G[pid](pp);
                        //std::cout << "Ampl Tree P = " << tree << std::endl;
                    }
                }
                if(prc.name.find("AsymL")!=std::string::npos){
                    findAmp(prc.name,pid);
                    if(DeltaL==-1){
                        //std::cout << " Process = " << prc.name << std::endl;
                        //std::cout << " Qnumber  = " << prc.qNumbers.begin()->second << std::endl;
                        ampLap += f_G[pid](pp);
                        //std::cout << "Ampl Loop AP = " << ampLap << std::endl;
                    }
                }
                if(prc.name.find("AsymT")!=std::string::npos){
                    findAmp(prc.name,pid);
                    if(DeltaL==-1){
                        //std::cout << " Process = " << prc.name << std::endl;
                        //std::cout << " Qnumber  = " << prc.qNumbers.begin()->second << std::endl;
                        tree += f_G[pid](pp);
                        //std::cout << "Ampl Tree P = " << tree << std::endl;
                    }
                }
            }
        }
    }
    
    return 2.*((ampLp - ampLap)).real()/(tree.real());
}


complex_t Decay(int &pid,param_t &pp,Process prc,int loop)
{
    complex_t res = 0.;
    
    std::string pname = f_G[pid].name;
    std::string ptl_mass = "m_";
    //std::cout << "Name =" << pname << std::endl;
    size_t found;
    double pmass;
    int out = prc.outParticles.size();
    
    /*
    for(size_t tt = 0;tt!=vect_ptl.size();++tt){
        (vect_ptl[tt].first<pname.find("to")) ? in++ : out++;
        std::cout << ">" << vect_ptl[tt].second << std::endl;
    }
    */
    //std::cout << "The process is 1->" << out << std::endl;
    
    complex_t jacD = 0.,decay = 0.;
    Eigen::MatrixXcd pmat(1+out,1+out);
    Eigen::VectorXcd mass(1+out);
    
    
    vector<double> mmax;
    size_t tt_max = prc.inParticles.size() + prc.outParticles.size();
    size_t tt = 0;
    bool cm0 = false;
    for(const auto &prname : prc.inParticles){
        for(const auto &m0 : Massless){
            if((prname.name.find(m0)!=std::string::npos))
            {
                cm0 = true;
            }
        }
        if(cm0)
        {
            //std::cout << "Mass of Photon "<< prname.name << " = 0" << std::endl;
            mass(tt) = 0.;
        }
        else{
            auto pt_mass = pp.realParams[ptl_mass + prname.name];
            //std::cout << "Mass of "<< prname.name << " = " << *pt_mass << std::endl;
            mmax.push_back(*pt_mass);
            mass(tt) = *pt_mass;
            
        }
        ++tt;
    }
    
    for(const auto &prname : prc.outParticles){
        for(const auto &m0 : Massless){
            if((prname.name.find(m0)!=std::string::npos))
            {
                cm0 = true;
            }
        }
        if(cm0)
        {
            //std::cout << "Mass of Photon "<< prname.name << " = 0" << std::endl;
            mass(tt) = 0.;
        }
        else{
            auto pt_mass = pp.realParams[ptl_mass + prname.name];
            //std::cout << "Mass of "<< prname.name << " = " << *pt_mass << std::endl;
            mmax.push_back(*pt_mass);
            mass(tt) = *pt_mass;
            
        }
        ++tt;
    }
    
    
    //cout << "The Masses are: "<< endl << mass  << endl;
    if(out==3){
        decay = 0.;
        for(int ii = 0;ii<13;ii++){
            for(int jj = 0;jj<13;jj++){
                for(int kk = 0;kk<13;kk++){
                    
                    vector<double> cth,rr;
                    cth.push_back(xs[jj]);
                    cth.push_back(xs[kk]);
                    rr.push_back(xs1[ii]);
                    pmat = pmom(mass(0),mass,cth,rr,1,out);
                    pp.s_12 = pmat(0,1).real();
                    pp.s_13 = pmat(0,2).real();
                    pp.s_14 = pmat(0,3).real();
                    pp.s_23 = pmat(1,2).real();
                    pp.s_24 = pmat(1,3).real();
                    pp.s_34 = pmat(2,3).real();
                    
                    jacD = jac(mass(0),mass,rr,1,out);
                    /*
                     for(int ll=0; ll!=f_G.size();++ll){
                     if(f_G[ll].name == pname){
                     decay += f_G[ll](pp)*jacD*as[jj]*as[kk]*as1[ii];
                     }
                     }*/
                    
                    if(loop==0){
                        if(f_G[pid].name.find("Tree")!=std::string::npos){
                            decay += f_G[pid](pp)*jacD*as[jj]*as[kk]*as1[ii];
                        }
                    }
                    if(loop==1){
                        if(f_G[pid].name.find("AsymL")!=std::string::npos){
                            decay += f_G[pid](pp)*jacD*as[jj]*as[kk]*as1[ii];
                        }
                    }
                    //cout << "AMp : " << f_G[pid].name << " : " << f_G[pid](pp) << endl;
                }
            }
        }
    }
    if(out==2){
        decay = 0.;
        vector<double> cth,rr;
        cth.push_back(1.);
        rr.push_back(1.);
        pmat = pmom(mass(0),mass,cth,rr,1,out);
        pp.s_12 = pmat(0,1).real();
        pp.s_13 = pmat(0,2).real();
        pp.s_23 = pmat(1,2).real();
        
        jacD = jac(mass(0),mass,rr,1,out);
        
        //std::cout << "jac = " << jacD << std::endl;
        //std::cout << "Proc = " << f_G[pid].name << std::endl;
        //pp.print();
        /*
         for(int kk=0; kk!=f_G.size();++kk){
         if(f_G[kk].name == pname){
         decay += f_G[kk](pp)*jacD;}
         }*/
        
        if(loop==0){
            if(f_G[pid].name.find("Tree")!=std::string::npos){
                decay += f_G[pid](pp)*jacD;
            }
        }
        if(loop==1){
            if(f_G[pid].name.find("AsymL")!=std::string::npos){
                decay += f_G[pid](pp)*jacD;
            }
        }
        
        
        
    }
    res = decay;
    return res;
}

double gamma(int &pid,double const &T,param_t &data,Process prc,int loop)
{
    double result = 0;
    
    double gg = 1.;
    vector< pair <int,double> > gvect;
    
    std::string pname = f_G[pid].name;
    std::string ptl_mass = "m_";
    
    //std::cout << "Name =" << pname << std::endl;
    double pmass;
        
    int in = prc.inParticles.size(),out = prc.outParticles.size();
    /*
    for(int tt = 0;tt!=vect.size();++tt){
        (vect[tt].first<pname.find("to")) ? in++ : out++;
        
    }
    */
    size_t found;
    
    for(const auto &inpp : prc.inParticles)
    {
        for(auto &gi : gptl){
            gg *= (inpp.name == gi.first) ? gi.second : 1.;
        }
    }
    
    
    Eigen::MatrixXcd pmat(in+out,in+out);
    Eigen::VectorXcd mass(in+out);
    
    vector<double> mmax;
    size_t tt_max = prc.inParticles.size() + prc.outParticles.size();
    size_t tt = 0;
    bool cm0 = false;
    for(const auto &prname : prc.inParticles){
        for(const auto &m0 : Massless){
            if((prname.name==m0)||(prname.name==m0+"c"))
            {
                cm0 = true;
            }
        }
        if(cm0)
        {
            //std::cout << "Mass of Photon "<< prname.name << " = 0" << std::endl;
            mass(tt) = 0.;
        }
        else{
            auto pt_mass = data.realParams[ptl_mass + prname.name];
            //std::cout << "Mass of "<< prname.name << " = " << *pt_mass << std::endl;
            mmax.push_back(*pt_mass);
            mass(tt) = *pt_mass;
            
        }
        ++tt;
    }
    
    for(const auto &prname : prc.outParticles){
        for(const auto &m0 : Massless){
            if((prname.name==m0)||(prname.name==m0+"c"))
            {
                cm0 = true;
            }
        }
        if(cm0)
        {
            //std::cout << "Mass of Photon "<< prname.name << " = 0" << std::endl;
            mass(tt) = 0.;
        }
        else{
            auto pt_mass = data.realParams[ptl_mass + prname.name];
            //std::cout << "Mass of "<< prname.name << " = " << *pt_mass << std::endl;
            mmax.push_back(*pt_mass);
            mass(tt) = *pt_mass;
            
        }
        ++tt;
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
                if(loop==0){
                    if(f_G[pid].name.find("Tree")!=std::string::npos){
                        amp += f_G[pid](data)*jac;
                    }
                }
                if(loop==1){
                    if(f_G[pid].name.find("AsymL")!=std::string::npos){
                        amp += f_G[pid](data)*jac;
                    }
                }
                result += exp(xki[ii])*Tsc/(256.*pow(M_PI,5.))*pin*pout*amp.real()*1./sqrt(s)*cyl_bessel_k(1, sqrt(s)/Tsc);
                cth.clear();
                rr.clear();
            }
        }
    }
    
    if(in==1)
    {
        complex_t decay = Decay(pid,data,prc,loop);
        
        double mD = mass(0).real();
        result = gg*mD*mD/(2.*pow(M_PI,2.))*T*cyl_bessel_k(1, abs(mD)/T)*decay.real();
    }
    //cout << "\n  Res = " << result << endl;
    return result;
}

void BESolver(std::ofstream &ofile,param_t &pp)
{
    std::vector<Process> procL,procQ;
    
    ParticleData pdata;
    
    pdata.loadFile("script/test.json");
    
    int N_proc;
    int out_proc;
    std::cout << "Process for L != 0" << std::endl;
    procL = pdata.getProcessesFromQNumberViolation("L");
    procQ = pdata.getProcessesFromQNumberConservation("Q");
    int id = 0;
    int DeltaL1 = 1,DeltaL2 = 2;
    double gammaL1 = 0.,gammaL2 = 0.;
    
    vector<double> mass_scale;
    
    for(auto &name : Xptl){
        mass_scale.push_back(Massp(name,pp));
    }
    
    double Mscale = *max_element(mass_scale.begin(), mass_scale.end());;
    size_t nn = Xptl.size() + 1;
    
    Eigen::VectorXcd yn(nn),eps1(Xptl.size());
    
    Eigen::MatrixXcd eps2a(Xptl.size(),Xptl.size()),eps2b(Xptl.size(),Xptl.size());
    
    Eigen::VectorXcd tempf(nn),bth(nn),yth(nn),ypth(nn),yeq(nn);
    
    Eigen::MatrixXcd tempfp(nn,nn),athinv(nn,nn),exath(nn,nn);
    Eigen::MatrixXd ath(nn,nn);
    
    Eigen::MatrixXcd Gam2(nn,nn);
    Eigen::VectorXcd Gam1(nn);
    
    
    yn = Eigen::VectorXcd::Zero(nn);
    /*
    size_t yi = 0;
    for(const auto &mname : Xptl){
        yn(yi) = Yeq(mname,pp,pp.m_N_3*10.);
        ++yi;
    }*/
    //yn << 1.,1.,1.,1.,0.;
    yn << 0.,0.,0.,0.,0.;
    
    
    
    eps1 << epsilon1("L","N_3",pp),epsilon1("L","N_2",pp),epsilon1("L","N_1",pp),0.;
    
    std::cout << "eps =  " << eps1(0) << " " << eps1(1) << " " << eps1(2) << std::endl;
    eps2a << 0., 0.,0.,0.,
             0., 0.,0.,0.,
             0., 0.,0.,0.,
             0., 0.,0.,0.;
    
    eps2b << 0., 0.,0.,0.,
             0., 0.,0.,0.,
             0., 0.,0.,0.,
             0., 0.,0.,0.;

    double steps = 1000.;
    double sph = 28./79.;
    double hh = 3./steps;
    double zi = 0.;
    double TT;
    double tolerance = 1.E-05;
    do{
        double zz = pow(10., -1. + hh*zi);
        TT = Mscale/zz;
        
        /*
        std::cout << TT << " ";
        for(const auto &mname : Xptl){
            std::cout << Yeq(mname,pp,TT) << " ";
        }
        std::cout << std::endl;
        */
        complex<double> Xsec = 0., pr = 0.;
        
        gammaL1 = 0;
        gammaL2 = 0.;
        double Hb = sqrt(4.*pow(M_PI,3.)*geffT(TT)/45.)*pow(TT,2.)/mpl;
        
        double s = geffT(TT)*2*pow(M_PI,2.)/45.*pow(TT,3.);
        
        double sT = geffT(TT)*2*pow(M_PI,2.)/45.;
        
        double yeql = 4./pow(M_PI,2.)*pow(sT,-1.);
        /**/
        tempf = Eigen::VectorXcd::Zero(nn);
        tempfp = Eigen::MatrixXcd::Zero(nn,nn);
        
        Gam1 = Eigen::VectorXcd::Zero(nn);
        Gam2 = Eigen::MatrixXcd::Zero(nn,nn);
        
        yth = yn;
        
        size_t si = 0;
        
        auto k1byk2 = [&](const double& zx) { return zx/(1.5 + zx); };
        
        for(const auto &mname : Xptl){
            yeq(si) = Yeq(mname,pp,TT);
            //tempf(si) += (yth(si)+1.)*cyl_bessel_k(1, Massp(mname,pp)/TT)/cyl_bessel_k(2, Massp(mname,pp)/TT);
            //tempfp(si,si) += cyl_bessel_k(1, Massp(mname,pp)/TT)/cyl_bessel_k(2, Massp(mname,pp)/TT);
            tempf(si) += (yth(si)+1.)*k1byk2(Massp(mname,pp)/TT);
            tempfp(si,si) += k1byk2(Massp(mname,pp)/TT);
            ++si;
        }
        
        //std::cout << "tempfp > " << tempf.transpose() << std::endl;
        //std::cout << "tempfp > " << std::endl << tempfp << std::endl;
        std::cout << std::endl;
        
        yeq(nn-1) = yeql;
        
        for(size_t ii = 0;ii!=Xptl.size();++ii){
            for(size_t jj = 0;jj!=Xptl.size();++jj){
                for(const auto &prname : procQ){
                    N_proc = prname.inParticles.size();
                    if(N_proc>1){
                        Xsec = 0.;
                        if((prname.inParticles[0].name.find(Xptl[ii]) != std::string::npos)&&(prname.inParticles[1].name.find(Xptl[jj] ) != std::string::npos))
                        {
                            findAmp(prname.name,id);
                            pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                            Xsec += (pr.real()>tolerance*yeq(ii).real()) ? pr : 0.;
                        }
                        
                        tempf(ii) += -log(10)*zz/yeq(ii)*((yth(ii)*yth(jj) + yth(ii) + yth(jj))*Xsec);
                        
                        for(size_t kk = 0;kk!=Xptl.size();++kk){
                            
                            double del_ki = (kk==ii) ? 1. : 0;
                            double del_kj = (kk==jj) ? 1. : 0;
                            tempfp(ii,kk) += -log(10)*zz/yeq(ii)*(del_ki*yth(jj) + del_ki + del_kj + del_kj*yth(ii))*Xsec;
                        }
                        
                    }
                //}
                //for(const auto &prname : procL){
                    for(const auto &pin : prname.inParticles){
                        for(const auto &pout : prname.outParticles){
                            if((pin.name.find(Xptl[ii]) != std::string::npos)&&(pout.name.find(Xptl[jj] ) != std::string::npos)){
                                Xsec = 0.;
                                findAmp(prname.name,id);
                                pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                Xsec += (pr.real()>tolerance*yeq(ii).real()) ? pr : 0.;
                                
                                tempf(ii) += -log(10)*zz/yeq(ii)*((yth(ii) - yth(jj))*Xsec);
                                
                                for(size_t kk = 0;kk!=Xptl.size();++kk){
                                    double del_ki = (kk==ii) ? 1. : 0;
                                    double del_kj = (kk==jj) ? 1. : 0;
                                    tempfp(ii,kk) += -log(10)*zz/yeq(ii)*((del_ki - del_kj)*Xsec);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        
        if(nn > Xptl.size()){
            //Gam1 = Eigen::VectorXcd::Zero(nn);
            complex<double> XsecL = 0.;
            
            for(const auto &prname : procL){
                double DeltaL = 1.;
                //====================================
                
                DeltaL=prname.qNumbers.begin()->second;
                N_proc = prname.inParticles.size();
                out_proc = prname.outParticles.size();
                if((N_proc<2)&&(out_proc<3)){
                    for(size_t ii = 0;ii!=Xptl.size();++ii){
                        for(size_t jj = 0;jj!=Xptl.size();++jj){
                            for(const auto &outptl : prname.outParticles){
                                if((prname.inParticles[0].name.find(Xptl[ii])!=std::string::npos)&&(outptl.name.find(Xptl[jj])!=std::string::npos)){
                                    findAmp(prname.name,id);
                                    Xsec = 0.;
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                    Xsec += (pr.real()>tolerance*yeql) ? pr : 0.;
                                    
                                    /**/
                                    XsecL = 0.;
                                    if(DeltaL>0){
                                        pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                        XsecL += (pr.real()>tolerance*yeql) ? 2.*pr : 0.;
                                    }
                                    if(DeltaL<0){
                                        pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                        XsecL -= (pr.real()>tolerance*yeql) ? 2.*pr : 0.;
                                    }
                                    
                                    //tempf(nn-1) += log(10)*zz*((yth(ii) - yth(jj))*eps1(ii)*Xsec - DeltaL*yth(nn-1)/yeql*Xsec*(yth(jj) + 1.) );
                                    tempf(nn-1) += log(10)*zz*((yth(ii) - yth(jj))*XsecL - DeltaL*yth(nn-1)/yeql*Xsec*(yth(jj) + 1.) );
                                    tempfp(nn-1,nn-1) += -log(10)*zz*DeltaL/yeql*Xsec*(yth(jj) + 1.);
                                    //tempfp(nn-1,ii) += log(10)*zz*eps1(ii)*Xsec;
                                    tempfp(nn-1,ii) += log(10)*zz*XsecL;
                                    //tempfp(nn-1,jj) += -log(10)*zz*(eps1(ii)*Xsec + yth(nn-1)*DeltaL/yeql*Xsec);
                                    tempfp(nn-1,jj) += -log(10)*zz*(XsecL + yth(nn-1)*DeltaL/yeql*Xsec);
                                }
                            }
                        }
                    }
                }
                /**/
                if((N_proc>1)&&(out_proc<3)){
                    bool in_lepto = false;
                    bool out_lepto = false;
                    for(const auto &inptl : prname.inParticles)
                        if(inptl.name.find("lL_") != std::string::npos) in_lepto = true;
                    for(const auto &outptl : prname.outParticles)
                        if(outptl.name.find("lL_") != std::string::npos) out_lepto = true;
                    
                    complex<double> ytemp = 0.,yX = 0.;
                    double deli = 0.;
                    double delj = 0.;
                    for(size_t ii = 0;ii!=Xptl.size();++ii){
                        for(size_t jj = 0;jj!=Xptl.size();++jj){
                            for(const auto &inptl : prname.inParticles){
                                for(const auto &outptl : prname.outParticles){
                                    if((inptl.name.find(Xptl[ii]) != std::string::npos)&&(outptl.name.find(Xptl[jj]) != std::string::npos)){
                                        findAmp(prname.name,id);
                                        Xsec = 0.;
                                        pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                        Xsec += (pr.real()>tolerance*yeql) ? pr : 0.;
                                        
                                        XsecL = 0.;
                                        /*
                                        if(DeltaL>0){
                                            pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                            XsecL += (pr.real()>tolerance) ? 2.*pr : 0.;
                                        }
                                        if(DeltaL<0){
                                            pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                            XsecL -= (pr.real()>tolerance) ? 2.*pr : 0.;
                                        }*/
                                         
                                        ytemp = (in_lepto) ? yth(ii) + 1. : 0.;
                                        yX += ytemp;
                                        ytemp = (out_lepto) ? yth(jj) + 1. : 0.;
                                        yX += ytemp;
                                        
                                        deli = (in_lepto) ? 1. : 0.;
                                        delj = (out_lepto) ? 1. : 0.;
                                        
                                        //std::cout << "XsecL = " << XsecL << std::endl;
                                        
                                        //tempf(nn-1) += log(10)*zz*( (yth(ii) - yth(jj))*eps2a(ii,jj)*Xsec - yth(nn-1)*DeltaL/yeql*Xsec*yX);
                                        tempf(nn-1) += log(10)*zz*( (yth(ii) - yth(jj))*XsecL - yth(nn-1)*DeltaL/yeql*Xsec*yX);
                                        tempfp(nn-1,nn-1) += -log(10)*zz*DeltaL/yeql*Xsec*yX;
                                        
                                        //tempfp(nn-1,ii) += log(10)*zz*(eps2a(ii,jj)*Xsec - yth(nn-1)*DeltaL/yeql*Xsec*deli);
                                        //tempfp(nn-1,jj) += -log(10)*zz*(eps2a(ii,jj)*Xsec + yth(nn-1)*DeltaL/yeql*Xsec*delj);
                                        tempfp(nn-1,ii) += log(10)*zz*(XsecL - yth(nn-1)*DeltaL/yeql*Xsec*deli);
                                        tempfp(nn-1,jj) += -log(10)*zz*(XsecL + yth(nn-1)*DeltaL/yeql*Xsec*delj);
                                        
                                    }
                                }
                            }
                        }
                        
                    }
                    
                    for(size_t ii = 0;ii!=Xptl.size();++ii){
                        for(size_t jj = 0;jj!=Xptl.size();++jj){
                            for(const auto &inptl1 : prname.inParticles){
                                for(const auto &inptl2 : prname.inParticles){
                                    if((inptl1.name.find(Xptl[ii]) != std::string::npos )&&(inptl2.name.find(Xptl[jj])!= std::string::npos)){
                                        findAmp(prname.name,id);
                                        Xsec = 0.;
                                        pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                        Xsec += (pr.real()>tolerance*yeql) ? pr : 0.;
                                        
                                        XsecL = 0.;
                                        /*
                                        if(DeltaL>0){
                                         pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                         XsecL += (pr.real()>tolerance) ? 2.*pr : 0.;
                                        }
                                        if(DeltaL<0){
                                         pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                         XsecL -= (pr.real()>tolerance) ? 2.*pr : 0.;
                                        }
                                         */
                                        //tempf(nn-1) += log(10)*zz*( (yth(ii)*yth(jj) + yth(ii) + yth(jj) )*eps2b(ii,jj)*Xsec - yth(nn-1)*DeltaL/yeql*Xsec);
                                        //tempfp(nn-1,ii) += log(10)*zz*(yth(jj) + 1.)*eps2b(ii,jj)*Xsec;
                                        //tempfp(nn-1,jj) += log(10)*zz*(yth(ii) + 1.)*eps2b(ii,jj)*Xsec;
                                        tempf(nn-1) += log(10)*zz*( (yth(ii)*yth(jj) + yth(ii) + yth(jj) )*XsecL - yth(nn-1)*DeltaL/yeql*Xsec);
                                        tempfp(nn-1,jj) += log(10)*zz*(yth(ii) + 1.)*XsecL;
                                        tempfp(nn-1,ii) += log(10)*zz*(yth(jj) + 1.)*XsecL;
                                        tempfp(nn-1,nn-1) += -log(10)*zz*DeltaL/yeql*Xsec;
                                    }
                                }
                            }
                        }
                    }
                }
                
                //======================================
            }
        }
        
        /*
        std::cout << "tempfp > " << tempf.transpose() << std::endl;
        std::cout << "tempfp > " << std::endl << tempfp << std::endl;
        std::cout << std::endl;
        std::cout << "tempfp^(-1) > " << std::endl << tempfp.inverse() << std::endl;
        std::cout << std::endl;
        std::cout << "ath^(-1)" << std::endl << ath.determinant() << std::endl;
         */
        ath = tempfp.real();
        bth = tempf - tempfp*yth;
        int dim = Xptl.size() + 1;
        //std::cout << "Exp Matrix >> " << std::endl << expM(ath*hh, dim) << std::endl << std::endl;
        
        if((ath.determinant()!=0)&&(!isnan(ath.determinant()))&&(!isinf(ath.determinant())))
        {
            ypth = expM(ath*hh, dim)*yth + tempfp.inverse()*(expM(ath*hh, dim) -  MatrixXd::Identity(dim, dim))*bth;
        }
        else
        {
            ypth = yth;
        }
          
        
        /*
        std::cout << "===============================================" << std::endl;
        std::cout << "ΔL = " << DeltaL1 <<std::endl;
        std::cout << "===============================================" << std::endl;
        */
        gammaL1 = 0.;
        gammaL2 = 0.;
        for(const auto &prname : procL)
        {
            
            if(DeltaL1==prname.qNumbers.begin()->second){
                findAmp(prname.name,id);
                if(!(prname.name.find("OneLoop")!=std::string::npos)){
                    double gamma_temp = gamma(id,TT,pp,prname,0);
                    gammaL1 += gamma_temp;
                }
                //gammaL1 += gamma(id,TT,pp,prname);
                //std::cout << "fG = " << id;
                //std::cout << "name = " << prname.name << " γ("<< TT << ") = " << gamma_temp << std::endl;
            }
            if(DeltaL2==prname.qNumbers.begin()->second){
                findAmp(prname.name,id);
                if(!(prname.name.find("OneLoop")!=std::string::npos)){
                    double gamma_temp = gamma(id,TT,pp,prname,0);
                    gammaL2 += gamma_temp;
                    //std::cout << "fG = " << id;
                    //std::cout << " name = " << prname.name << std::endl;
                    
                    
                }
            }
        }
        
        
        ofile << TT << "\t"<< zz << "\t" << DeltaL1*gammaL1/(s*Hb*zz*yeql) << "\t" << DeltaL2*gammaL2/(s*Hb*zz*yeql)
        << "\t" << yn(nn-1).real() << std::endl;
        std::cout << TT << "\t"<< zz << "\t" << DeltaL1*gammaL1/(s*Hb*zz*yeql) << "\t" << DeltaL2*gammaL2/(s*Hb*zz*yeql) << "\t" << sph*yn(nn-1).real()/(8.6E-11) << std::endl;
        
        
        
        yn = ypth;
        
        std::cout << yn.transpose() << std::endl;
      
        
        ++zi;
    }while((zi<=steps));
    
}


double washout(double zz,double &Mscale,std::vector<Process> &procL,param_t &pp)
{
    double TT = Mscale/zz;
    //std::cout << "TT = " << TT << std::endl;
    double Hb = sqrt(4.*pow(M_PI,3.)*geffT(TT)/45.)*pow(TT,2.)/mpl;
    double sT = geffT(TT)*2*pow(M_PI,2.)/45.;
    double yeql = 4./pow(M_PI,2.)*pow(sT,-1.);
    double s = geffT(TT)*2*pow(M_PI,2.)/45.*pow(TT,3.);
    double gammaL1 = 0.;
    double gammaL2 = 0.;
    int DeltaL1 = 1,DeltaL2 = 2;
    int id = 0;
    for(const auto &prname : procL)
    {
        
        if(DeltaL1==prname.qNumbers.begin()->second){
            findAmp(prname.name,id);
            if((prname.name.find("Tree")!=std::string::npos)){
                double gamma_temp = gamma(id,TT,pp,prname,0);
                gammaL1 += gamma_temp;
            }
        }
        if(DeltaL2==prname.qNumbers.begin()->second){
            findAmp(prname.name,id);
            if((prname.name.find("Tree")!=std::string::npos)){
                double gamma_temp = gamma(id,TT,pp,prname,0);
                gammaL2 += gamma_temp;
            }
        }
    }
    return DeltaL1*gammaL1/(s*Hb*zz*yeql) + DeltaL2*gammaL2/(s*Hb*zz*yeql);
    
}

double bisection(double a, double b,double z_cut,param_t &pp,double &Mscale,asym_param ap)
{
    double fa,fb,fc;
    
    std::vector<Process> procL,procQ;
    
    //ParticleData pdata;
    
    //pdata.loadFile("script/test.json");
    
    procL = ap.procL;//pdata.getProcessesFromQNumberViolation("L");
    
    
    fa = washout(pow(10.,a),Mscale,procL,pp) - z_cut;
    
    fb = washout(pow(10.,b),Mscale,procL,pp) - z_cut;
    
    if (fa * fb >= 0) {
        //cout << "You have not assumed right a and b\n";
        return 0;
    }
    double cc = a;
    while ((b-a) >= EP) {
        // Find middle point
        cc = (a+b)/2;
        // Check if middle point is root
        fc = washout(pow(10.,cc),Mscale,procL,pp) - z_cut;
        if (fc == 0.0)
            break;
        // Decide the side to repeat the steps
        else if (fc*fa < 0)
            b = cc;
        else
            a = cc;
    }
    //cout << "The value of root is : " << cc << endl;
    
    return cc;
}


double freezeout(param_t &pp,asym_param ap)
{
    std::vector<Process> procL,procQ;
    
    //ParticleData pdata;
    
    //pdata.loadFile("script/test.json");
    
    procL = ap.procL;//pdata.getProcessesFromQNumberViolation("L");
    procQ = ap.procQ;//pdata.getProcessesFromQNumberConservation("Q");
    
    //double yB_ratio = 0.;
    
    vector<double> mass_scale;
    
    for(auto &name : Xptl){
        mass_scale.push_back(Massp(name,pp));
    }
    
    double Mscale = *max_element(mass_scale.begin(), mass_scale.end());;
    
    //double z_result = bisection(-1., 2.,pp,Mscale);
    
    //std::cout << "====== \t z_reasult =  " << z_result << std::endl;
    double steps = 10.;
    double sph = 28./79.;
    double hh = 3./steps;
    
    double TT;
    int id = 0;
    double zi = 0.;
    double zz_in = pow(10.,bisection(-1., 2.,10.,pp,Mscale,ap));//zi;
    double zz_out = pow(10.,bisection(-1., 2.,0.1,pp,Mscale,ap));//zi;
     
    double tolerance = 1.E-05;
    
    int N_proc;
    int out_proc;
    size_t nn = Xptl.size() + 1;
    
    Eigen::VectorXcd yn(nn),eps1(Xptl.size());
    
    Eigen::MatrixXcd eps2a(Xptl.size(),Xptl.size()),eps2b(Xptl.size(),Xptl.size());
    
    Eigen::VectorXcd tempf(nn),bth(nn),yth(nn),ypth(nn),yeq(nn);
    
    Eigen::MatrixXcd tempfp(nn,nn),athinv(nn,nn),exath(nn,nn);
    Eigen::MatrixXd ath(nn,nn);
    
    Eigen::MatrixXcd Gam2(nn,nn);
    Eigen::VectorXcd Gam1(nn);
    
    
    yn = Eigen::VectorXcd::Zero(nn);
    yn << 0.,0.,0.,0.,0.;
    
    hh = (log10(zz_out) - log10(zz_in))/steps;
    
    zi = 0.;
    do{
        double zz = pow(10., log10(zz_in) + hh*zi);
        TT = Mscale/zz;
        
        complex<double> Xsec = 0., pr = 0.;
        
        double Hb = sqrt(4.*pow(M_PI,3.)*geffT(TT)/45.)*pow(TT,2.)/mpl;
        double sT = geffT(TT)*2*pow(M_PI,2.)/45.;
        double yeql = 4./pow(M_PI,2.)*pow(sT,-1.);
        double s = geffT(TT)*2*pow(M_PI,2.)/45.*pow(TT,3.);
        
        tempf = Eigen::VectorXcd::Zero(nn);
        tempfp = Eigen::MatrixXcd::Zero(nn,nn);
        
        Gam1 = Eigen::VectorXcd::Zero(nn);
        Gam2 = Eigen::MatrixXcd::Zero(nn,nn);
        
        yth = yn;
        
        size_t si = 0;
        
        auto k1byk2 = [&](const double& zx) { return zx/(1.5 + zx); };
        
        for(const auto &mname : Xptl){
            yeq(si) = Yeq(mname,pp,TT);
            tempf(si) += (yth(si)+1.)*k1byk2(Massp(mname,pp)/TT);
            tempfp(si,si) += k1byk2(Massp(mname,pp)/TT);
            ++si;
        }
        
        
        yeq(nn-1) = yeql;
        
        for(size_t ii = 0;ii!=Xptl.size();++ii){
            for(size_t jj = 0;jj!=Xptl.size();++jj){
                for(const auto &prname : procQ){
                    N_proc = prname.inParticles.size();
                    if(N_proc>1){
                        Xsec = 0.;
                        if((prname.inParticles[0].name.find(Xptl[ii]) != std::string::npos)&&(prname.inParticles[1].name.find(Xptl[jj] ) != std::string::npos))
                        {
                            findAmp(prname.name,id);
                            pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                            Xsec += (pr.real()>tolerance*yeq(ii).real()) ? pr : 0.;
                        }
                        
                        tempf(ii) += -log(10)*zz/yeq(ii)*((yth(ii)*yth(jj) + yth(ii) + yth(jj))*Xsec);
                        
                        for(size_t kk = 0;kk!=Xptl.size();++kk){
                            double del_ki = (kk==ii) ? 1. : 0;
                            double del_kj = (kk==jj) ? 1. : 0;
                            tempfp(ii,kk) += -log(10)*zz/yeq(ii)*(del_ki*yth(jj) + del_ki + del_kj + del_kj*yth(ii))*Xsec;
                        }
                        
                    }
                //}
                //for(const auto &prname : procL){
                    for(const auto &pin : prname.inParticles){
                        for(const auto &pout : prname.outParticles){
                            if((pin.name.find(Xptl[ii]) != std::string::npos)&&(pout.name.find(Xptl[jj] ) != std::string::npos)){
                                Xsec = 0.;
                                findAmp(prname.name,id);
                                pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                Xsec += (pr.real()>tolerance*yeq(ii).real()) ? pr : 0.;
                                
                                tempf(ii) += -log(10)*zz/yeq(ii)*((yth(ii) - yth(jj))*Xsec);
                                
                                for(size_t kk = 0;kk!=Xptl.size();++kk){
                                    double del_ki = (kk==ii) ? 1. : 0;
                                    double del_kj = (kk==jj) ? 1. : 0;
                                    tempfp(ii,kk) += -log(10)*zz/yeq(ii)*((del_ki - del_kj)*Xsec);
                                }
                            }
                        }
                    }
                }
            }
        }
         
        
        
        if(nn > Xptl.size()){
            //Gam1 = Eigen::VectorXcd::Zero(nn);
            complex<double> XsecL = 0.;
            
            for(const auto &prname : procL){
                double DeltaL = 1.;
                //====================================
                
                DeltaL=prname.qNumbers.begin()->second;
                N_proc = prname.inParticles.size();
                out_proc = prname.outParticles.size();
                if((N_proc<2)&&(out_proc<3)){
                    for(size_t ii = 0;ii!=Xptl.size();++ii){
                        for(size_t jj = 0;jj!=Xptl.size();++jj){
                            for(const auto &outptl : prname.outParticles){
                                if((prname.inParticles[0].name.find(Xptl[ii])!=std::string::npos)&&(outptl.name.find(Xptl[jj])!=std::string::npos)){
                                    findAmp(prname.name,id);
                                    Xsec = 0.;
                                    pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                    Xsec += (pr.real()>tolerance*yeql) ? pr : 0.;
                                    
                                    /**/
                                    XsecL = 0.;
                                    if(DeltaL>0){
                                        pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                        XsecL += (pr.real()>tolerance*yeql) ? 2.*pr : 0.;
                                    }
                                    if(DeltaL<0){
                                        pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                        XsecL -= (pr.real()>tolerance*yeql) ? 2.*pr : 0.;
                                    }
                                    
                                    //tempf(nn-1) += log(10)*zz*((yth(ii) - yth(jj))*eps1(ii)*Xsec - DeltaL*yth(nn-1)/yeql*Xsec*(yth(jj) + 1.) );
                                    tempf(nn-1) += log(10)*zz*((yth(ii) - yth(jj))*XsecL - DeltaL*yth(nn-1)/yeql*Xsec*(yth(jj) + 1.) );
                                    tempfp(nn-1,nn-1) += -log(10)*zz*DeltaL/yeql*Xsec*(yth(jj) + 1.);
                                    //tempfp(nn-1,ii) += log(10)*zz*eps1(ii)*Xsec;
                                    tempfp(nn-1,ii) += log(10)*zz*XsecL;
                                    //tempfp(nn-1,jj) += -log(10)*zz*(eps1(ii)*Xsec + yth(nn-1)*DeltaL/yeql*Xsec);
                                    tempfp(nn-1,jj) += -log(10)*zz*(XsecL + yth(nn-1)*DeltaL/yeql*Xsec);
                                }
                            }
                        }
                    }
                }
                /**/
                if((N_proc>1)&&(out_proc<3)){
                    bool in_lepto = false;
                    bool out_lepto = false;
                    for(const auto &inptl : prname.inParticles)
                        if(inptl.name.find("lL_") != std::string::npos) in_lepto = true;
                    for(const auto &outptl : prname.outParticles)
                        if(outptl.name.find("lL_") != std::string::npos) out_lepto = true;
                    
                    complex<double> ytemp = 0.,yX = 0.;
                    double deli = 0.;
                    double delj = 0.;
                    for(size_t ii = 0;ii!=Xptl.size();++ii){
                        for(size_t jj = 0;jj!=Xptl.size();++jj){
                            for(const auto &inptl : prname.inParticles){
                                for(const auto &outptl : prname.outParticles){
                                    if((inptl.name.find(Xptl[ii]) != std::string::npos)&&(outptl.name.find(Xptl[jj]) != std::string::npos)){
                                        findAmp(prname.name,id);
                                        Xsec = 0.;
                                        pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                        Xsec += (pr.real()>tolerance*yeql) ? pr : 0.;
                                        
                                        XsecL = 0.;
                                         
                                        /*
                                        if(DeltaL>0){
                                            pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                            XsecL += (pr.real()>tolerance) ? 2.*pr : 0.;
                                        }
                                        if(DeltaL<0){
                                            pr = (id < ((int) f_G.size())) ? abs(gamma(id,TT,pp,prname,1)/(s*Hb*zz)) : 0.;
                                            XsecL -= (pr.real()>tolerance) ? 2.*pr : 0.;
                                        }*/
                                        ytemp = (in_lepto) ? yth(ii) + 1. : 0.;
                                        yX += ytemp;
                                        ytemp = (out_lepto) ? yth(jj) + 1. : 0.;
                                        yX += ytemp;
                                        
                                        deli = (in_lepto) ? 1. : 0.;
                                        delj = (out_lepto) ? 1. : 0.;
                                        
                                        tempf(nn-1) += log(10)*zz*( (yth(ii) - yth(jj))*XsecL - yth(nn-1)*DeltaL/yeql*Xsec*yX);
                                        tempfp(nn-1,nn-1) += -log(10)*zz*DeltaL/yeql*Xsec*yX;
                                        
                                        tempfp(nn-1,ii) += log(10)*zz*(XsecL - yth(nn-1)*DeltaL/yeql*Xsec*deli);
                                        tempfp(nn-1,jj) += -log(10)*zz*(XsecL + yth(nn-1)*DeltaL/yeql*Xsec*delj);
                                        
                                    }
                                }
                            }
                        }
                        
                    }
                    
                    for(size_t ii = 0;ii!=Xptl.size();++ii){
                        for(size_t jj = 0;jj!=Xptl.size();++jj){
                            for(const auto &inptl1 : prname.inParticles){
                                for(const auto &inptl2 : prname.inParticles){
                                    if((inptl1.name.find(Xptl[ii]) != std::string::npos )&&(inptl2.name.find(Xptl[jj])!= std::string::npos)){
                                        findAmp(prname.name,id);
                                        Xsec = 0.;
                                        pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,0)/(s*Hb*zz) : 0.;
                                        Xsec += (pr.real()>tolerance*yeql) ? pr : 0.;
                                        
                                        XsecL = 0.;
                                        /*
                                        if(DeltaL>0){
                                         pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                         XsecL += (pr.real()>tolerance) ? 2.*pr : 0.;
                                        }
                                        if(DeltaL<0){
                                         pr = (id < ((int) f_G.size())) ? gamma(id,TT,pp,prname,1)/(s*Hb*zz) : 0.;
                                         XsecL -= (pr.real()>tolerance) ? 2.*pr : 0.;
                                        }
                                         */
                                        tempf(nn-1) += log(10)*zz*( (yth(ii)*yth(jj) + yth(ii) + yth(jj) )*XsecL - yth(nn-1)*DeltaL/yeql*Xsec);
                                        tempfp(nn-1,jj) += log(10)*zz*(yth(ii) + 1.)*XsecL;
                                        tempfp(nn-1,ii) += log(10)*zz*(yth(jj) + 1.)*XsecL;
                                        tempfp(nn-1,nn-1) += -log(10)*zz*DeltaL/yeql*Xsec;
                                    }
                                }
                            }
                        }
                    }
                }
                
                //======================================
            }
        }
        
        ath = tempfp.real();
        bth = tempf - tempfp*yth;
        int dim = Xptl.size() + 1;
        
        if((ath.determinant()!=0)&&(!isnan(ath.determinant()))&&(!isinf(ath.determinant())))
        {
            ypth = expM(ath*hh, dim)*yth + tempfp.inverse()*(expM(ath*hh, dim) -  MatrixXd::Identity(dim, dim))*bth;
        }
        else
        {
            ypth = yth;
        }
        //ofile << TT << "\t"<< zz << "\t" << DeltaL1*gammaL1/(s*Hb*zz*yeql) << "\t" << DeltaL2*gammaL2/(s*Hb*zz*yeql) << "\t" << yn(nn-1).real() << std::endl;
        std::cout << TT << "\t"<< zz << "\t" << washout(zz,Mscale,procL,pp) << "\t" << sph*yn(nn-1).real()/(8.6E-11) << std::endl;
        
        yn = ypth;
        //std::cout << yn.transpose() << std::endl;
      
        ++zi;
    }while(zi<=steps);//while((TT>=100.));//
    
    return sph*yn(nn-1).real()/(8.6E-11);
    
}
}
