// We compute the electric components -(1)tr(G0jG0j) = Nelec*g^2 to
// lowest order with SF boundary conditions.
// We use Martin Luscher's formulation based on an orbifolding symmetry.
// Gauge parameters for the flow (alpha) or the action (lambda) can be
// varied.

// initial code by Argia Rubeo with some corrections by S. Sint (July 2016)
// S. Sint (August 2016; added LW action, Symanzik flow, used symmetry of spatial momenta,
//  some simplifications)

// A. Rubeo May 2017;   further modifications
//                      moving definitions of the global variables in the library NelecDerCb.hpp
//                      adding a new input parameter c_b
//                      cuputing the analytical derivitative with respect to cb (de_E_el/de_cb)
//                      adding epsilon (eps) as input parameter
//                      check the analytical derivative with the symmetric numerical (der_num_sym)
//                      adding continuum value of N_el for c=0.2,0.3,0.4 (computed using continuum_energy_density/)
//                      computation of second dervitive
//                      continuum limit of the second derivative

//Note: in this program cbMat = Id4 + c_b * wilsonKernel


#include "Nelec.hpp"

void update_p (int n_0, int n_1, int n_2, int n_3);
void init_PmuPnu (arma::mat& P_muP_nu);
void init_cbMat (arma::mat& cbMat,arma::mat& P_muP_nu, double c_b);
void init_map(double);

int main(int argc, char *argv[])
{
//    cout << argc << argv[0] << endl;
    if (argc != 9 ) {
        cout << "Usage: Nelec eps ob fl ac c L/a T/a c_b \n" << endl;
        cout << "EXAMPLE: Nelec 0.01 imp z pl 0.3 8 8 0.02 \n" << endl;
        cout << "eps = parameter to compute numerical derivative \n" << endl;
        cout << "c= sqrt(8t)/L  \n" << endl;
        cout << "ob = observable, options: pl (plaq not symm),pls (plaq symm), cl (clover (symm by def)),plcl () \n" << endl;
        cout << "fl=flow, options: fl = w,z,s \n" << endl;
        cout << "ac=action, options: ac = pl,lw \n" << endl;
        cout << "L/a = spatial lattice size \n" << endl;
        cout << "T/a = temporal lattice size \n" << endl;
        cout << "c_b = imporvement coeff. in the initial condition for the flow \n" << endl;
        exit(1);
    }

    // input parameters
    eps = (double) atof(argv[1]);
    string ob = argv[2];   // options pl,cl,lw,plcl
    string fl = argv[3];   // options z or w
    string ac = argv[4];   // options pl or lw
    double c = (double) atof(argv[5]);
    L = (double) atof(argv[6]);
    T = (double) atof(argv[7]);
    c_b = (double) atof(argv[8]); // imporvement coeff. in the initial condition for the flow


    // flow time from c
    double t= c*c*L*L/8.0;

    // declare various matrices:

    // matrix for Wilson kernel
    arma::mat K_pl;

    // generic placeholder matrices
    arma::mat Kob;
    arma::mat Dprop;
    arma::mat Hker;
    arma::mat Dbar;
    arma::mat de_Dbar_de_cb;

    // projectors
    arma::mat TRAN;
    arma::mat LONG;

    //identity
    arma::mat Id4;
    // matrix K_Zeuthen
    arma::mat K_z;
    arma::mat Corrfact;
    // matrix LW kernel with alpha
    arma::mat K_lw;
     // matrix LW kernel with lambda
    arma::mat K_lw_lam;

    // matrix kernels for observables
    arma::mat K_cl_obs;
    arma::mat K_lw_obs;
    arma::mat K_pl_obs;
    arma::mat K_plcl_obs;

    //cbMat
    arma::mat cbMat(4,4);
    arma::mat P_muP_nu(4,4);
    arma::mat de_Dprop_de_cb(4,4);
    //numerical deriv
    arma::mat cbMat_p_eps(4,4);
    arma::mat Dbar_p_eps;
    arma::mat cbMat_m_eps(4,4);
    arma::mat Dbar_m_eps;


    Id4.eye(4,4);

    for (int x=1;x<T;x++){

        Nel[x]=0.0;
        secder_Nel[x]=0.0;
        firstder_Nel[x]=0.0;
        El[x]=0.0;
        Elpl[x]=0.0;
        Elpls[x]=0.0;
        Elcl[x]=0.0;
        Elplcl[x]=0.0;
        Elimp[x]=0.0;

        secder_Elpls[x]=0.0;
        secder_Elcl[x]=0.0;
        secder_Elpl[x]=0.0;
        secder_Elplcl[x]=0.0;
        secder_Elimp[x]=0.0;
        secder_El[x]=0.0;

        firstder_Elpls[x]=0.0;
        firstder_Elcl[x]=0.0;
        firstder_Elpl[x]=0.0;
        firstder_Elplcl[x]=0.0;
        firstder_Elimp[x]=0.0;
        firstder_El[x]=0.0;

        Elpl_p_eps[x]=0.0;
        Elcl_p_eps[x]=0.0;
        Elpl_m_eps[x]=0.0;
        Elcl_m_eps[x]=0.0;


        de_Nel[x]=0.0;
        de_El[x]=0.0;
        de_Elpl[x]=0.0;
        de_Elpls[x]=0.0;
        de_Elcl[x]=0.0;
        de_Elplcl[x]=0.0;
        de_Elimp[x]=0.0;


    }


    for(int n_0 =0; n_0< T; n_0++){
        for(int n_1 =0; n_1< L; n_1++){
            for(int n_2 =0; n_2< L; n_2++){
                for(int n_3 =0; n_3< L; n_3++){
                 //if (n_0+n_3!= 0) {
                    //updating p values for each n_0 n_1 n_2 n_3
                    update_p(n_0, n_1, n_2, n_3);

                     double degfact = 0.0;
                     degfact=1.;
                     //if ( (n_1 != n_2) && (n_1 != n_3) && (n_2 != n_3) ){ degfact = 6.;}
                     //else if ((n_1 == n_2) && (n_1 == n_3) && (n_2 == n_3)){ degfact = 1.;}
                     //else { degfact = 3.;}
                     //if (n_0 == 0) {degfact=degfact/2.;}

                     double Elpl_p=0.0;
                     double Elcl_p=0.0;

                     double Elpl_p_p_eps=0.0;
                     double Elcl_p_p_eps=0.0;
                     double Elpl_p_m_eps=0.0;
                     double Elcl_p_m_eps=0.0;

                     double de_Elpl_p=0.0;
                     double de_Elcl_p=0.0;


                     // Transverse and Longitudinal projectors
                     LONG
                     << p0_2nd << p0p1   << p0p2   << p0p3   << arma::endr
                     << p1p0   << p1_2nd << p1p2   << p1p3   << arma::endr
                     << p2p0   << p2p1   << p2_2nd << p2p3   << arma::endr
                     << p3p0   << p3p1   << p3p2   << p3_2nd << arma::endr;

                     LONG= (1.0/p_hat_2nd) * LONG;
                     TRAN = Id4 - LONG;

                     // correction factor for Zeuthen flow
                     Corrfact
                     << prefact0 << 0.0      << 0.0      << 0.0      << arma::endr
                     << 0.0      << prefact1 << 0.0      << 0.0      << arma::endr
                     << 0.0      << 0.0      << prefact2 << 0.0      << arma::endr
                     << 0.0      << 0.0      << 0.0      << prefact3 << arma::endr;

                     // LW kernel
                     K_lw
                     // 1st row
                     << p_hat_2nd - p0_2nd + (1./12.)*(p_hat_4th+ p_hat_2nd*p0_2nd - 2.0 * p0_4th)
                     << -p0p1 +(1./12.)*(-p0p1*(p0_2nd+p1_2nd))
                     << -p0p2 +(1./12.)*(-p0p2*(p0_2nd+p2_2nd))
                     << -p0p3 +(1./12.)*(-p0p3*(p0_2nd+p3_2nd))
                     << arma::endr
                     // 2nd row
                     << -p1p0 +(1./12.)*(-p1p0*(p1_2nd+p0_2nd))
                     << p_hat_2nd - p1_2nd + (1./12.)*(p_hat_4th+ p_hat_2nd*p1_2nd - 2.0 * p1_4th)
                     << -p1p2 +(1./12.)*(-p1p2*(p1_2nd+p2_2nd))
                     << -p1p3 +(1./12.)*(-p1p3*(p1_2nd+p3_2nd))
                     << arma::endr
                     // 3rd row
                     << -p2p0 +(1./12.)*(-p2p0*(p2_2nd+p0_2nd))
                     << -p2p1 +(1./12.)*(-p2p1*(p2_2nd+p1_2nd))
                     << p_hat_2nd - p2_2nd + (1./12.)*(p_hat_4th+ p_hat_2nd*p2_2nd - 2.0 * p2_4th)
                     << -p2p3 +(1./12.)*(-p2p3*(p2_2nd+p3_2nd))
                     << arma::endr
                     // 4th row
                     << -p3p0 +(1./12.)*(-p3p0*(p3_2nd+p0_2nd))
                     << -p3p1 +(1./12.)*(-p3p1*(p3_2nd+p1_2nd))
                     << -p3p2 +(1./12.)*(-p3p2*(p3_2nd+p2_2nd))
                     << p_hat_2nd - p3_2nd + (1./12.)*(p_hat_4th+ p_hat_2nd*p3_2nd - 2.0 * p3_4th)
                     << arma::endr;

                     // Zeuthen kernel
                     K_z = Corrfact * K_lw + alpha * p_hat_2nd * LONG;

                     K_lw_lam = K_lw + lambda*p_hat_2nd*LONG; // (to define LW propagator)
                     K_lw = K_lw + alpha*p_hat_2nd*LONG;  // (to define Symanzik flow)

                     init_cbMat (cbMat, P_muP_nu, c_b);
                     init_cbMat(cbMat_p_eps, P_muP_nu, c_b+eps);
                     init_cbMat(cbMat_m_eps, P_muP_nu, c_b-eps);

//                     cout << "cbMat_p_eps" << cbMat_p_eps << endl;
//                     cout << "cbMat" << cbMat << endl;



                     if (ac == "pl" ){Dprop = (1.0/ p_hat_2nd) * (TRAN + (1/lambda)* LONG);}
                     else if (ac == "lw" ){Dprop = (K_lw_lam.i());}   // inversion could be avoided using ML's formula (careful with gauge fixing!)
                     else{
                        cout << "Options for action are pl (plaq) or lw (Luscher-Weisz), not: " << ac << endl;
                        return 0;
                        }
                     if (fl == "w" ){Hker = exp(-t * p_hat_2nd) * TRAN + exp(-t* alpha* p_hat_2nd) * LONG;}
                     else if (fl == "z" ){Hker = expmat(-t * K_z);}
                     else if (fl == "s" ){Hker = expmat(-t * K_lw);}
                     else{ cout << "Options for flow are w (Wilson), z (Zeuthen) or s (Symanzik), not: " << fl << endl;
                        return 0;
                        }

                     //Dbar = Hker * Dprop * (Hker.t());

                     //the new Dbar contains the cb dependence
                     Dbar = Hker * (cbMat * Dprop * (cbMat.t()) ) * (Hker.t());

                     Dbar_p_eps = Hker * (cbMat_p_eps * Dprop * (cbMat_p_eps.t()) ) * (Hker.t());
                     Dbar_m_eps = Hker * (cbMat_m_eps * Dprop * (cbMat_m_eps.t()) ) * (Hker.t());

                     // PmuPnu matrix needed to introduce (de Dbar / de cb)
                     init_PmuPnu (P_muP_nu);

                     de_Dprop_de_cb = (((p_hat_2nd * Id4) - P_muP_nu) * Dprop * (cbMat.t())) + (cbMat * Dprop * ((p_hat_2nd * Id4) - P_muP_nu));

                     de_Dbar_de_cb = Hker * ( de_Dprop_de_cb ) * (Hker.t());
                     //de_Dbar_de_cb = (Hker * ((p_hat_2nd * Id4) - P_muP_nu) * Dprop * cbMat.t() * Hker.t()) + (Hker * cbMat * Dprop * ((p_hat_2nd * Id4) - P_muP_nu) * Hker.t());
                     // calculate observables for electric components

                     //                    for(int i =1; i< 4; i++){Elpl_p += p0_2nd*Dbar.at(i,i);}
                     Elpl_p += p0_2nd*(Dbar.at(1,1)+Dbar.at(2,2)+Dbar.at(3,3));
                     Elpl_p += p_hat_2nds*Dbar.at(0,0) - 2.*(p0p1*Dbar.at(0,1)+p0p2*Dbar.at(0,2)+p0p3*Dbar.at(0,3));
                     Elcl_p += (1. - p0_2nd/4.)* (1. - p1_2nd/4.)*(p0_2nd*Dbar.at(1,1) + p1_2nd*Dbar.at(0,0)- 2.*p0p1*Dbar.at(0,1));
                     Elcl_p += (1. - p0_2nd/4.)* (1. - p2_2nd/4.)*(p0_2nd*Dbar.at(2,2) + p2_2nd*Dbar.at(0,0)- 2.*p0p2*Dbar.at(0,2));
                     Elcl_p += (1. - p0_2nd/4.)* (1. - p3_2nd/4.)*(p0_2nd*Dbar.at(3,3) + p3_2nd*Dbar.at(0,0)- 2.*p0p3*Dbar.at(0,3));


                     Elpl_p_p_eps += p0_2nd*(Dbar_p_eps.at(1,1)+Dbar_p_eps.at(2,2)+Dbar_p_eps.at(3,3));
                     Elpl_p_p_eps += p_hat_2nds*Dbar_p_eps.at(0,0) - 2*(p0p1*Dbar_p_eps.at(0,1)+p0p2*Dbar_p_eps.at(0,2)+p0p3*Dbar_p_eps.at(0,3));
                     Elcl_p_p_eps += (1. - p0_2nd/4.)* (1. - p1_2nd/4.)*(p0_2nd*Dbar_p_eps.at(1,1) + p1_2nd*Dbar_p_eps.at(0,0)- 2.*p0p1*Dbar_p_eps.at(0,1));
                     Elcl_p_p_eps += (1. - p0_2nd/4.)* (1. - p2_2nd/4.)*(p0_2nd*Dbar_p_eps.at(2,2) + p2_2nd*Dbar_p_eps.at(0,0)- 2.*p0p2*Dbar_p_eps.at(0,2));
                     Elcl_p_p_eps += (1. - p0_2nd/4.)* (1. - p3_2nd/4.)*(p0_2nd*Dbar_p_eps.at(3,3) + p3_2nd*Dbar_p_eps.at(0,0)- 2.*p0p3*Dbar_p_eps.at(0,3));

                     Elpl_p_m_eps += p0_2nd*(Dbar_m_eps.at(1,1)+Dbar_m_eps.at(2,2)+Dbar_m_eps.at(3,3));
                     Elpl_p_m_eps += p_hat_2nds*Dbar_m_eps.at(0,0) - 2*(p0p1*Dbar_m_eps.at(0,1)+p0p2*Dbar_m_eps.at(0,2)+p0p3*Dbar_m_eps.at(0,3));
                     Elcl_p_m_eps += (1. - p0_2nd/4.)* (1. - p1_2nd/4.)*(p0_2nd*Dbar_m_eps.at(1,1) + p1_2nd*Dbar_m_eps.at(0,0)- 2.*p0p1*Dbar_m_eps.at(0,1));
                     Elcl_p_m_eps += (1. - p0_2nd/4.)* (1. - p2_2nd/4.)*(p0_2nd*Dbar_m_eps.at(2,2) + p2_2nd*Dbar_m_eps.at(0,0)- 2.*p0p2*Dbar_m_eps.at(0,2));
                     Elcl_p_m_eps += (1. - p0_2nd/4.)* (1. - p3_2nd/4.)*(p0_2nd*Dbar_m_eps.at(3,3) + p3_2nd*Dbar_m_eps.at(0,0)- 2.*p0p3*Dbar_m_eps.at(0,3));


                     de_Elpl_p += p0_2nd*(de_Dbar_de_cb.at(1,1)+de_Dbar_de_cb.at(2,2)+de_Dbar_de_cb.at(3,3));
                     de_Elpl_p += p_hat_2nds*de_Dbar_de_cb.at(0,0) - 2.*(p0p1*de_Dbar_de_cb.at(0,1)+p0p2*de_Dbar_de_cb.at(0,2)+p0p3*de_Dbar_de_cb.at(0,3));
                     de_Elcl_p += (1. - p0_2nd/4.)* (1. - p1_2nd/4.)*(p0_2nd*de_Dbar_de_cb.at(1,1) + p1_2nd*de_Dbar_de_cb.at(0,0)- 2.*p0p1*de_Dbar_de_cb.at(0,1));
                     de_Elcl_p += (1. - p0_2nd/4.)* (1. - p2_2nd/4.)*(p0_2nd*de_Dbar_de_cb.at(2,2) + p2_2nd*de_Dbar_de_cb.at(0,0)- 2.*p0p2*de_Dbar_de_cb.at(0,2));
                     de_Elcl_p += (1. - p0_2nd/4.)* (1. - p3_2nd/4.)*(p0_2nd*de_Dbar_de_cb.at(3,3) + p3_2nd*de_Dbar_de_cb.at(0,0)- 2.*p0p3*de_Dbar_de_cb.at(0,3));

                     // add time dependent terms for plaquette and clover compbinations
                     for (int x=1; x<T; x++){
                         Elpl[x] +=  degfact*sin(p_0*x+p_0/2.)*sin(p_0*x+p_0/2.)* Elpl_p;
                         Elcl[x] +=  degfact*sin(p_0*x)*sin(p_0*x)*Elcl_p;

                         //first deriv
                         firstder_Elpl[x] += degfact*Elpl_p* ((1./2.)*( sin(p_0*(x+1)+(p_0/2.)) * sin(p_0*(x+1)+(p_0/2.))\
                                                                    + sin(p_0*(x-1)+(p_0/2.)) * sin(p_0*(x-1)+(p_0/2.)) ));

                         // second derivative NOT IMPROVED
                         //secder_Elpl[x] += degfact*Elpl_p*((1./4.)*( sin(p_0*(x+2)+(p_0/2.)) * sin(p_0*(x+2)+(p_0/2.))\
                                                                    - 2 *sin(p_0*x+p_0/2.)*sin(p_0*x+p_0/2.)\
                                                                    + sin(p_0*(x-2)+(p_0/2.)) * sin(p_0*(x-2)+(p_0/2.)) ));
                         secder_Elpl[x] += degfact*Elpl_p*( - double(1./12.)*sin(p_0*(x+2)+(p_0/2.)) * sin(p_0*(x+2)+(p_0/2.)) \
                                                            + double(4./3.)* sin(p_0*(x+1)+(p_0/2.)) * sin(p_0*(x+1)+(p_0/2.)) \
                                                            - double(5./2.)*(sin(p_0*x+p_0/2. )*sin(p_0*x+p_0/2.)) \
                                                            + double(4./3.)*(sin(p_0*(x-1)+(p_0/2.))* sin(p_0*(x-1)+(p_0/2.))) \
                                                            - double(1./12.)*(sin(p_0*(x-2)+(p_0/2.))*sin(p_0*(x-2)+(p_0/2.))));
                        firstder_Elcl[x] += degfact*Elpl_p*((1./2.)*( sin(p_0*(x+1)) * sin(p_0*(x+1))\
                                                                    + sin(p_0*(x-1)) * sin(p_0*(x-1)) ));
                        //NOT IMPROVED
                        //secder_Elcl[x] += degfact*Elpl_p*((1./4.)*( sin(p_0*(x+2)) * sin(p_0*(x+2))\
                                                               - 2 *sin(p_0*x)*sin(p_0*x)\
                                                              + sin(p_0*(x-2)) * sin(p_0*(x-2)) ));
                         secder_Elcl[x] +=  degfact*Elcl_p*( - double(1./12.)*sin(p_0*(x+2)) * sin(p_0*(x+2)) \
                                                             + double(4./3.)* sin(p_0*(x+1)) * sin(p_0*(x+1)) \
                                                             - double(5./2.)*(sin(p_0*x )*sin(p_0*x)) \
                                                             + double(4./3.)*(sin(p_0*(x-1))* sin(p_0*(x-1))) \
                                                             - double(1./12.)*(sin(p_0*(x-2))*sin(p_0*(x-2))));

                         Elpl_p_eps[x] +=  degfact*sin(p_0*x+p_0/2.)*sin(p_0*x+p_0/2.)* Elpl_p_p_eps;
                         Elcl_p_eps[x] +=  degfact*sin(p_0*x)*sin(p_0*x)*Elcl_p_p_eps;

                         Elpl_m_eps[x] +=  degfact*sin(p_0*x+p_0/2.)*sin(p_0*x+p_0/2.)* Elpl_p_m_eps;
                         Elcl_m_eps[x] +=  degfact*sin(p_0*x)*sin(p_0*x)*Elcl_p_m_eps;


                         de_Elpl[x] +=  degfact*sin(p_0*x+p_0/2.)*sin(p_0*x+p_0/2.)* de_Elpl_p;
                         de_Elcl[x] +=  degfact*sin(p_0*x)*sin(p_0*x)*de_Elcl_p;

                     }
                        //}  // if n3 !=0 (momentum 0 excluded)
                    }
                }
            }
        }

    for (int x=1; x<T; x++){
        Elpls[x] = (Elpl[x]+Elpl[x-1])/2.;
        Elplcl[x]+= (4./3.)*Elpls[x]-(1./3.)*Elcl[x];

        firstder_Elpls[x] = (firstder_Elpl[x]+firstder_Elpl[x-1])/2.;
        firstder_Elplcl[x]+= (4./3.)*firstder_Elpls[x]-(1./3.)*firstder_Elcl[x];

        secder_Elpls[x] = (secder_Elpl[x]+secder_Elpl[x-1])/2.;
        secder_Elplcl[x]+= (4./3.)*secder_Elpls[x]-(1./3.)*secder_Elcl[x];

        Elpls_p_eps[x] = (Elpl_p_eps[x]+Elpl_p_eps[x-1])/2.;
        Elplcl_p_eps[x]+= (4./3.)*Elpls_p_eps[x]-(1./3.)*Elcl_p_eps[x];

        Elpls_m_eps[x] = (Elpl_m_eps[x]+Elpl_m_eps[x-1])/2.;
        Elplcl_m_eps[x]+= (4./3.)*Elpls_m_eps[x]-(1./3.)*Elcl_m_eps[x];


        de_Elpls[x] = (de_Elpl[x]+de_Elpl[x-1])/2.;
        de_Elplcl[x]+= (4./3.)*de_Elpls[x]-(1./3.)*de_Elcl[x];
    }
    for (int x=1; x<T; x++){
        Elimp[x] = (4./3.)*Elplcl[x] -(1./6.)*(Elplcl[x-1]+Elplcl[x+1]);
        //Elimp[x] = Elplcl[x] +(1./6.)*(-0.0160266241215); //sum range 20
        //Elimp[x] = Elplcl[x] +(1./6.)*(-0.0160266283632); //sum range 40
        secder_Elimp[x] = (4./3.)* secder_Elplcl[x] -(1./6.)*(secder_Elplcl[x-1]+secder_Elplcl[x+1]);
        firstder_Elimp[x] = (4./3.)* firstder_Elplcl[x] -(1./6.)*(firstder_Elplcl[x-1]+firstder_Elplcl[x+1]);

        de_Elimp[x] = (4./3.)*de_Elplcl[x] -(1./6.)*(de_Elplcl[x-1]+de_Elplcl[x+1]);

        Elimp_p_eps[x] = (4./3.)*Elplcl_p_eps[x] -(1./6.)*(Elplcl_p_eps[x-1]+Elplcl_p_eps[x+1]);
        Elimp_m_eps[x] = (4./3.)*Elplcl_m_eps[x] -(1./6.)*(Elplcl_m_eps[x-1]+Elplcl_m_eps[x+1]);

        }

    if (ob == "pl") {for (int x=1; x<T; x++){El[x] = Elpl[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "pls") {for (int x=1; x<T; x++){El[x] = Elpls[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "cl") {for (int x=1; x<T; x++){El[x] = Elcl[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "plcl"){for (int x=1; x<T; x++){El[x] = Elplcl[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "imp"){for (int x=1; x<T; x++){El[x] = Elimp[x] *(8.0/ ((double)T*L*L*L));}}

    //second derivative
    if (ob == "pl") {for (int x=1; x<T; x++){secder_El[x] = secder_Elpl[x] *(8.0/ ((double)T*L));}}
    else if (ob == "pls") {for (int x=1; x<T; x++){secder_El[x] = secder_Elpls[x] *(8.0/ ((double)T*L));}}
    else if (ob == "cl") {for (int x=1; x<T; x++){secder_El[x] = secder_Elcl[x] *(8.0/ ((double)T*L));}}
    else if (ob == "plcl"){for (int x=1; x<T; x++){secder_El[x] = secder_Elplcl[x] *(8.0/ ((double)T*L));}}
    else if (ob == "imp"){for (int x=1; x<T; x++){secder_El[x] = secder_Elimp[x] *(8.0/ ((double)T*L));}}

    //first derivative
    if (ob == "pl") {for (int x=1; x<T; x++){firstder_El[x] = firstder_Elpl[x] *(8.0/ ((double)T*L*L));}}
    else if (ob == "pls") {for (int x=1; x<T; x++){firstder_El[x] = firstder_Elpls[x] *(8.0/ ((double)T*L*L));}}
    else if (ob == "cl") {for (int x=1; x<T; x++){firstder_El[x] = firstder_Elcl[x] *(8.0/ ((double)T*L*L));}}
    else if (ob == "plcl"){for (int x=1; x<T; x++){firstder_El[x] = firstder_Elplcl[x] *(8.0/ ((double)T*L*L));}}
    else if (ob == "imp"){for (int x=1; x<T; x++){firstder_El[x] = firstder_Elimp[x] *(8.0/ ((double)T*L*L));}}

    if (ob == "pl") {for (int x=1; x<T; x++){El_p_eps[x] = Elpl_p_eps[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "pls") {for (int x=1; x<T; x++){El_p_eps[x] = Elpls_p_eps[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "cl") {for (int x=1; x<T; x++){El_p_eps[x] = Elcl_p_eps[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "plcl"){for (int x=1; x<T; x++){El_p_eps[x] = Elplcl_p_eps[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "imp"){for (int x=1; x<T; x++){El_p_eps[x] = Elimp_p_eps[x] *(8.0/ ((double)T*L*L*L));}}

    if (ob == "pl") {for (int x=1; x<T; x++){El_m_eps[x] = Elpl_m_eps[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "pls") {for (int x=1; x<T; x++){El_m_eps[x] = Elpls_m_eps[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "cl") {for (int x=1; x<T; x++){El_m_eps[x] = Elcl_m_eps[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "plcl"){for (int x=1; x<T; x++){El_m_eps[x] = Elplcl_m_eps[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "imp"){for (int x=1; x<T; x++){El_m_eps[x] = Elimp_m_eps[x] *(8.0/ ((double)T*L*L*L));}}


    if (ob == "pl") {for (int x=1; x<T; x++){de_El[x] = de_Elpl[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "pls") {for (int x=1; x<T; x++){de_El[x] = de_Elpls[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "cl") {for (int x=1; x<T; x++){de_El[x] = de_Elcl[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "plcl"){for (int x=1; x<T; x++){de_El[x] = de_Elplcl[x] *(8.0/ ((double)T*L*L*L));}}
    else if (ob == "imp"){for (int x=1; x<T; x++){de_El[x] = de_Elimp[x] *(8.0/ ((double)T*L*L*L));}}

    else{ cout << "Options for obs are pl (plaq), cl (clover), pls (sym. plaq) plcl (plaq-clover), imp (imprvd), not: " << ob << endl;
    return 0;
    }


    //continuum values from cont_energyDen.cpp (depending on c and x_0_T)
    for (int x=1; x<T; x++){

        if ((double)x/(double)T == 0.5) {

            if (c == 0.3){

                E_cont_el_T2 = 8.7406217269326777e-03;
                secder_cont_el_T2 = -0.016026628363208;

            }
            else if (c == 0.4){

                E_cont_el_T2 = 0.00724640965499;
                secder_cont_el_T2 = -0.0006010427791;

            }
            else if (c == 0.2){

                E_cont_el_T2 = 0.00934853264901;
                secder_cont_el_T2 = -0.00340266;
            }
            //c=0.5
            else{

                E_cont_el_T2 = 0.00510967817985;
                cout << "ERROR --computing second derivative only for c = 0.2, 0.3, 0.4 " << endl;

            }
            if (c != 0.3 && c != 0.4 && c != 0.5 && c != 0.25 && c != 0.2 && c != 0.35){
                cout << "computing only c=0.3, 0.4, 0.5, 0.25 ,0.2 , 0.35" << endl;
                exit(1);
            }
        }

        //x_0 = T/4
        if ((double)x/(double)T == 0.25) {

            if (c == 0.3){

                E_cont_el_T4 = 0.00855663273721;

            }
            else if (c == 0.4){

                E_cont_el_T4 = 0.00920113799832;

            }
            else if (c == 0.2){

                E_cont_el_T4 = 0.00900850006105;
            }
            //c=0.5
            else{

                E_cont_el_T4 = 0.0104747861051;

            }
            if (c != 0.3 && c != 0.4 && c != 0.5 && c != 0.2 ){
                cout << "computing only c=0.3, 0.4, 0.5 ,0.2 " << endl;
                exit(1);
            }
        }

    }


    //cout << "#L T  E_el_norm_to_cont        de_E_el/de_cb        ((E_lat/E_cont) -1)  cb                    E_lat                x_0/T" << endl;

    for (int x=1; x<T; x++){

        Nel[x] = t*t*El[x];
        secder_Nel[x] = t*t*secder_El[x];
        firstder_Nel[x] = t*t*firstder_El[x];
        //Nel[x] = 0.5 * (Nel[x] + Nel[x+1]);
        Nel_p_eps[x] = t*t*El_p_eps[x];
        Nel_m_eps[x] = t*t*El_m_eps[x];
        de_Nel[x] = t*t* de_El[x];

        double der_num = (Nel_p_eps[x] - Nel[x])*(1.0/eps);
        double der_num_sym = (Nel_p_eps[x] - Nel_m_eps[x]) * (1.0/(2.*eps));

//        if ((double)x/(double)T == 0.25) {

//            cout << setfill('0') << setw(2) << L
            //<< DELIMITER << setfill('0') << setw(2) << T
            //<< DELIMITER<< DELIMITER<< DELIMITER << fixed << setprecision(18) << Nel[x]
            //<< DELIMITER<< DELIMITER<< DELIMITER << fixed << setprecision(18) << Nel[x]/E_cont_el_T4
            //<< DELIMITER << de_Nel[x]
//            << DELIMITER << (Nel[x] / E_cont_el_T4)-1.0
//            << DELIMITER << (double)x/(double)T
//            << DELIMITER << c_b
//            << DELIMITER << setprecision(16) << scientific << Nel[x]
//            << DELIMITER << setprecision(16) << scientific << secder_Nel[x]
//            << endl;

//        }
        if ((double)x/(double)T == 0.5) {

           cout << setfill('0') << setw(2) << L
            //<< DELIMITER << (double)x/(double)T
            << DELIMITER << c_b
            << DELIMITER << setprecision(16) << scientific << Nel[x]
            //<< DELIMITER << setprecision(16) << scientific << ((Nel[x] / E_cont_el_T2)-1.0)
            << DELIMITER << setprecision(16) << scientific << secder_Nel[x]
            //<< DELIMITER << setprecision(16) << scientific << ((secder_Nel[x] / secder_cont_el_T2)-1.0)
            //<< DELIMITER << setfill('0') << setw(2) << T
            //<< DELIMITER << setprecision(16) << scientific << firstder_Nel[x]
            << endl;

        }


    }
    return 0;
}



// FUNCTIONS USED IN THE MAIN


void update_p (int n_0, int n_1, int n_2, int n_3)
{

    p_1= (2*M_PI/L) * n_1;
    p_2= (2*M_PI/L) * n_2;
    p_3= (2*M_PI/L) * n_3;
    p_0= (M_PI/T)   * (n_0+0.5);


    p_hat_0 = 2.0 * sin (0.5 * p_0);
    p_hat_1 = 2.0 * sin (0.5 * p_1);
    p_hat_2 = 2.0 * sin (0.5 * p_2);
    p_hat_3 = 2.0 * sin (0.5 * p_3);


    p0p1 = p_hat_0 * p_hat_1;
    p1p0 = p0p1;
    p0p2 = p_hat_0 * p_hat_2;
    p2p0 = p0p2;
    p0p3 = p_hat_0 * p_hat_3;
    p3p0 = p0p3;
    p1p2 = p_hat_1 * p_hat_2;
    p2p1 = p1p2;
    p1p3 = p_hat_1 * p_hat_3;
    p3p1 = p1p3;
    p2p3 = p_hat_2 * p_hat_3;
    p3p2 = p2p3;

    p0_2nd= p_hat_0 * p_hat_0;
    p0_4th= p0_2nd * p0_2nd;
    p1_2nd= p_hat_1 * p_hat_1;
    p1_4th= p1_2nd * p1_2nd;
    p2_2nd= p_hat_2 * p_hat_2;
    p2_4th= p2_2nd * p2_2nd;
    p3_2nd= p_hat_3 * p_hat_3;
    p3_4th= p3_2nd * p3_2nd;


    p_hat_2nd = p0_2nd + p1_2nd + p2_2nd + p3_2nd;
    p_hat_2nds = p1_2nd + p2_2nd + p3_2nd;

    p_hat_4th  = p0_4th + p1_4th + p2_4th + p3_4th;
    p_hat_4ths = p1_4th + p2_4th + p3_4th;

    prefact0 = 1.0 - p0_2nd/12.0;
    prefact1 = 1.0 - p1_2nd/12.0;
    prefact2 = 1.0 - p2_2nd/12.0;
    prefact3 = 1.0 - p3_2nd/12.0;

}

void init_PmuPnu (arma::mat& P_muP_nu){

    P_muP_nu
    << p_hat_0*p_hat_0 << p_hat_0*p_hat_1 << p_hat_0*p_hat_2 << p_hat_0* p_hat_3<< arma::endr
    << p_hat_1*p_hat_0 << p_hat_1*p_hat_1 << p_hat_1*p_hat_2 << p_hat_1* p_hat_3<< arma::endr
    << p_hat_2*p_hat_0 << p_hat_2*p_hat_1 << p_hat_2*p_hat_2 << p_hat_2* p_hat_3<< arma::endr
    << p_hat_3*p_hat_0 << p_hat_3*p_hat_1 << p_hat_3*p_hat_2 << p_hat_3* p_hat_3<< arma::endr;

}

void init_cbMat(arma::mat& cbMat, arma::mat& P_muP_nu, double c_b){

    arma::mat Id4;
    Id4.eye(4,4);

    arma::mat wilsonKernel;
    wilsonKernel.eye(4,4);

    wilsonKernel = p_hat_2nd * Id4 - P_muP_nu;

    //cbMat = Id4 + c_b * wilsonKernel;
    cbMat = expmat(c_b * wilsonKernel);

}
