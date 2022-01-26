
#include "leptoproto.h"
#include "BE.h"

using namespace std;
using namespace csl;
using namespace mty;

using namespace leptoproto;

using namespace mty::lib;


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

/**/
void analysis(int i_start,int i_end,asym_param ap)
{
    double ii = ((double) i_start);
    
    std::string filename = "output-" + std::to_string(i_start) + std::to_string(i_end) + ".txt";
    std::ofstream output_file(filename);
    
    std::vector<Process> procL,procQ;
    
    //ParticleData pdata;
    
    //pdata.loadFile("script/test.json");
    
    procL = ap.procL;//pdata.getProcessesFromQNumberViolation("L");
    procQ = ap.procQ;//pdata.getProcessesFromQNumberConservation("Q");
    
    /*
    asym_param ap;
    
    ap.procL = procL;
    ap.procQ = procQ;
    */
    
    //srand(time(NULL));
    
    
    do
    {
        //double MN3 = ((double) ii)*1000.;
        //for(auto jj=1;jj!=6;++jj)
        {
            //int xx = rand()%20 + 1;
            //int yy = rand()%10 + 1;
            
            //double MN3 = ((double) xx)*1000.;
            double MN3 = ((double) ii)*100. + 1000.;
            
            double diff = 100.;//((double) jj)*100.;
            
            for(auto kk=1;kk!=20;++kk)
            {
                Eigen::VectorXd MN(3);
                
                MN << 100.,MN3 - 2.*diff,MN3;
                double meta = MN3 - diff;
                double muSet = ((double) kk)*0.5*MN3 + MN3;
                std::cout << MN3 << "\t" << muSet/MN3 << "\t" << diff << "\t" << std::endl;
                //double muSet = ((double) yy)*MN3;
                
                //std::cout << "K1 = " << cyl_bessel_k(1, 100.) << std::endl;
                double alph_e = 1./137.;
                double vev = 246.;
                double lam5 = 1.e+00;
                double lamet = 0.5;
                
                auto F = [&](double const& x){ return x/(x-1.)*log(x);};
                
                double th12 = M_PI*33.8/180.;
                double th13 = M_PI*8.61/180.;
                double th23 = M_PI*46.2/180.;
                double del = 281.2/180.*M_PI;
                
                Eigen::MatrixXcd mnu(3,3);
                mnu = Eigen::MatrixXd::Zero(3,3);
                
                double mnu1 = 1.0E-18*1.0E-09;
                double Dm21 = 7.5E-05*1.0E-18;
                double meff = 2.4*1.0E-3*1.0E-18;
                double mnu2 = std::sqrt(Dm21 + mnu1*mnu1);
                double mnu3N = std::sqrt(meff + Dm21 + mnu1*mnu1);
                //ouble mnu3I = std::sqrt(meff + mnu1*mnu1);
                
                // NH
                mnu << sqrt(mnu1), 0., 0.,
                0., sqrt(mnu2), 0.,
                0., 0., sqrt(mnu3N);
                //  IH
                // mnu << std::sqrt(mnu3I), 0., 0.,
                // 0., std::sqrt(mnu1), 0.,
                // 0., 0., std::sqrt(mnu2);
                
                //std::cout << "mnu >> " << std::endl << mnu<< std::endl;
                Eigen::MatrixXcd O12(3,3);
                Eigen::MatrixXcd O13(3,3);
                Eigen::MatrixXcd O23(3,3);
                Eigen::MatrixXcd MM(3,3);
                
                Eigen::MatrixXcd UPMNS(3,3);
                
                O12 << cos(th12),-sin(th12),0.,
                sin(th12),cos(th12),0.,
                0.,0.,1.;
                O13 << cos(th13),0.,-sin(th13)*std::exp(-del*1i),
                0.,1.,0.,
                sin(th13)*exp(del*1i),0.,cos(th13);
                O23 << 1.,0.,0.,
                0.,cos(th23),-sin(th23),
                0.,sin(th23),cos(th23);
                
                Eigen::MatrixXcd sqLam(3,3);
                
                double MetaR = pow(meta,2.) + 2.*lam5*vev*vev;
                double MetaI = pow(meta,2.) - 2.*lam5*vev*vev;
                
                complex<double> lam1,lam2,lam3;
                
                lam1 = MN(0)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(0),2.))
                                                 - F(MetaI/pow(MN(0),2.)));
                
                lam2 = MN(1)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(2),2.))
                                                 - F(MetaI/pow(MN(2),2.)));
                
                lam3 = MN(2)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(0),2.))
                                                 - F(MetaI/pow(MN(0),2.)));
                
                sqLam << sqrt(lam1), 0., 0.,
                0., sqrt(lam2), 0.,
                0., 0., sqrt(lam3);
                
                Eigen::MatrixXcd RO23(3,3);
                complex<double> Rth23 = M_PI/8. - M_PI/8.i;
                //complex<double> Rth23 = -M_PI/4.i;
                RO23 << 1.,0.,0.,
                0.,cos(Rth23),-sin(Rth23),
                0.,sin(Rth23),cos(Rth23);
                
                UPMNS = O23 * O13 * O12;
                Eigen::MatrixXcd Yuk;
                Yuk = (sqLam.inverse()*RO23* mnu * (UPMNS.conjugate()).transpose());
                
                param_t pp;
                
                pp.e = std::sqrt(4.*M_PI*alph_e);
                pp.m_N_1 = MN(0);
                pp.m_N_2 = MN(1);
                pp.m_N_3 = MN(2);
                pp.muet2 = -std::pow(meta,2.);
                
                pp.thetaW = std::asin(std::sqrt(0.23116));
                pp.reg_prop = 1.;
                pp.muS2 = -1.;
                pp.lSet = 0.1;
                pp.lamet = lamet;
                pp.Finite = 1.;
                pp.muSet = muSet;//12000.;
                pp.Y_et_00 = Yuk(0,0);
                pp.Y_et_01 = Yuk(0,1);
                pp.Y_et_10 = Yuk(1,0);//0.001*factor;
                pp.Y_et_02 = Yuk(0,2);
                pp.Y_et_20 = Yuk(2,0);//factor;
                pp.Y_et_11 = Yuk(1,1);//0.01i*factor;
                pp.Y_et_12 = Yuk(1,2);//0.01i*factor;
                pp.Y_et_21 = Yuk(2,1);//factor;
                pp.Y_et_22 = Yuk(2,2);//factor;
                
                pp.Y_S_00 = 0.;
                pp.Y_S_01 = 0.0;
                pp.Y_S_10 = 0.0;
                pp.Y_S_02 = 0.0;
                pp.Y_S_20 = 0.0;
                pp.Y_S_11 = 0.;
                pp.Y_S_12 = 0.5;//*std::exp(pi/4.i);
                pp.Y_S_21 = 0.5;//*std::exp(pi/4.i);
                pp.Y_S_22 = 0.0;
                
                updateSpectrum(pp);
                
                //pp.print();
                
                //double yB_ratio = freezeout(pp);
                double yB_ratio = freezeout(pp,ap);
                
                output_file << MN3 << "\t" << muSet/MN3 << "\t" << diff << "\t" << yB_ratio << std::endl;
                
            }
        }
        
        ++ii;
        
    }while(ii<=i_end);
    output_file.close();
}

int main(int argc, char* argv[]) {
    
    Eigen::VectorXd MN(3);
    
    MN << 100.,atof(argv[3]),atof(argv[1]);
    double meta = atof(argv[2]);
    double muSet = atof(argv[4]);
    
    std::cout << "K1 = " << cyl_bessel_k(1, 100.) << std::endl;
    double alph_e = 1./137.;
    double vev = 246.;
    double lam5 = 1.e+00;
    double lamet = 0.5;
    
    auto F = [&](double const& x){ return x/(x-1.)*log(x);};
    
    double th12 = M_PI*33.8/180.;
    double th13 = M_PI*8.61/180.;
    double th23 = M_PI*46.2/180.;
    double del = 281.2/180.*M_PI;
    
    Eigen::MatrixXcd mnu(3,3);
    mnu = Eigen::MatrixXd::Zero(3,3);
    
    double mnu1 = 1.0E-18*1.0E-09;
    double Dm21 = 7.5E-05*1.0E-18;
    double meff = 2.4*1.0E-3*1.0E-18;
    double mnu2 = std::sqrt(Dm21 + mnu1*mnu1);
    double mnu3N = std::sqrt(meff + Dm21 + mnu1*mnu1);
    //double mnu3I = std::sqrt(meff + mnu1*mnu1);
    
    /* NH*/
    mnu << sqrt(mnu1), 0., 0.,
              0., sqrt(mnu2), 0.,
          0., 0., sqrt(mnu3N);
    /*  IH
    mnu << std::sqrt(mnu3I), 0., 0.,
              0., std::sqrt(mnu1), 0.,
              0., 0., std::sqrt(mnu2);*/
    
    std::cout << "mnu >> " << std::endl << mnu<< std::endl;
    Eigen::MatrixXcd O12(3,3);
    Eigen::MatrixXcd O13(3,3);
    Eigen::MatrixXcd O23(3,3);
    Eigen::MatrixXcd MM(3,3);
    
    Eigen::MatrixXcd UPMNS(3,3);
    
    O12 << cos(th12),-sin(th12),0.,
    sin(th12),cos(th12),0.,
    0.,0.,1.;
    O13 << cos(th13),0.,-sin(th13)*std::exp(-del*1i),
    0.,1.,0.,
    sin(th13)*exp(del*1i),0.,cos(th13);
    O23 << 1.,0.,0.,
    0.,cos(th23),-sin(th23),
    0.,sin(th23),cos(th23);
    
    Eigen::MatrixXcd sqLam(3,3);
    
    double MetaR = pow(meta,2.) + 2.*lam5*vev*vev;
    double MetaI = pow(meta,2.) - 2.*lam5*vev*vev;
    
    complex<double> lam1,lam2,lam3;
    
    lam1 = MN(0)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(0),2.))
     - F(MetaI/pow(MN(0),2.)));
    
    lam2 = MN(1)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(2),2.))
     - F(MetaI/pow(MN(2),2.)));
    
    lam3 = MN(2)/(32.*pow(M_PI,2.))*(F(MetaR/pow(MN(0),2.))
     - F(MetaI/pow(MN(0),2.)));
    
    sqLam << sqrt(lam1), 0., 0.,
     0., sqrt(lam2), 0.,
     0., 0., sqrt(lam3);
    
    Eigen::MatrixXcd RO23(3,3);
    complex<double> Rth23 = M_PI/8. - M_PI/8.i;
    //complex<double> Rth23 = -M_PI/4.i;
    RO23 << 1.,0.,0.,
    0.,cos(Rth23),-sin(Rth23),
    0.,sin(Rth23),cos(Rth23);
    
    UPMNS = O23 * O13 * O12;
    Eigen::MatrixXcd Yuk;
    Yuk = (sqLam.inverse()*RO23* mnu * (UPMNS.conjugate()).transpose());
    
    std::cout << "sqLam >> " << std::endl << sqLam << std::endl;
    std::cout << "UPMNS >> " << std::endl << UPMNS << std::endl;
    
    std::cout << "sqLam >> " << std::endl << sqLam.inverse() << std::endl;
    std::cout << "Yukawa >> " << std::endl << Yuk << std::endl;
    
    param_t pp;
    
    pp.e = std::sqrt(4.*M_PI*alph_e);
    pp.m_N_1 = MN(0);
    pp.m_N_2 = MN(1);
    pp.m_N_3 = MN(2);
    pp.muet2 = -std::pow(meta,2.);
    
    pp.thetaW = std::asin(std::sqrt(0.23116));
    pp.reg_prop = 1.;
    pp.muS2 = -1.;
    pp.lSet = 0.1;
    pp.lamet = lamet;
    pp.Finite = 1.;
    pp.muSet = muSet;//12000.;
    pp.Y_et_00 = Yuk(0,0);
    pp.Y_et_01 = Yuk(0,1);
    pp.Y_et_10 = Yuk(1,0);//0.001*factor;
    pp.Y_et_02 = Yuk(0,2);
    pp.Y_et_20 = Yuk(2,0);//factor;
    pp.Y_et_11 = Yuk(1,1);//0.01i*factor;
    pp.Y_et_12 = Yuk(1,2);//0.01i*factor;
    pp.Y_et_21 = Yuk(2,1);//factor;
    pp.Y_et_22 = Yuk(2,2);//factor;
    
    pp.Y_S_00 = 0.;
    pp.Y_S_01 = 0.0;
    pp.Y_S_10 = 0.0;
    pp.Y_S_02 = 0.0;
    pp.Y_S_20 = 0.0;
    pp.Y_S_11 = 0.;
    pp.Y_S_12 = 0.5;//*std::exp(pi/4.i);
    pp.Y_S_21 = 0.5;//*std::exp(pi/4.i);
    pp.Y_S_22 = 0.0;
    
    updateSpectrum(pp);
    
    pp.print();
    
    std::cout << "eps = " << epsilon1("L","N_3",pp) << std::endl;
    
    std::vector<Process> procL,procQ;
    
    ParticleData pdata;
    
    pdata.loadFile("script/test.json");
    
    procL = pdata.getProcessesFromQNumberViolation("L");
    procQ = pdata.getProcessesFromQNumberConservation("Q");
    
    /**/
    asym_param ap;
    
    ap.procL = procL;
    ap.procQ = procQ;
    
    
    auto startTime = std::chrono::high_resolution_clock::now();
    double yB_ratio = freezeout(pp,ap);
    //double yB_ratio = freezeout(pp);
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime);
    cout << "Sec: " << duration.count()/1000000 << endl;
    std::cout << "yB =  " << yB_ratio << std::endl;
    
    
    
    /**/
    
    //std::string pname = argv[5];
    //std::ofstream ofile(pname);
    
    //BESolver(ofile,pp);
    //ofile.close();
     
    
    //analysis(0,10);
    /*
    thread t1(analysis,1,50,ap);
    thread t2(analysis,51,100,ap);
    thread t3(analysis,101,150,ap);
    thread t4(analysis,151,200,ap);
    thread t5(analysis,21,25);
    thread t6(analysis,50,59);
    thread t7(analysis,60,69);
    thread t8(analysis,70,79);
    thread t9(analysis,80,89);
    thread t10(analysis,90,99);
     */
     
    /*
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t7.join();
    t8.join();
    t9.join();
    t10.join();
               */
    
    //std::cout << "Process for Q = 0" << std::endl;
    //proc = pdata.getProcessesFromQNumberConservation("Q");
    
    //for(const auto &prname : proc)
    //std::cout << "name = " << prname.name << std::endl;
    
    return 0;
}
