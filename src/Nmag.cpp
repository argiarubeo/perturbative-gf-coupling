// We compute the magnetic components -(1/2)tr(GijGij) = Nmag*g^2 to
//  lowest order with SF boundary conditions.
// We use Martin Luscher's formulation based on an orbifolding symmetry.
// Gauge parameters for the flow (alpha) or the action (lambda) can be
// varied.

// initial code by Argia Rubeo with some corrections by S. Sint (July 2016)
// S. Sint (August 2016; kernels for observables, added LW action, Symanzik flow, used symmetry of spatial momenta,
//  some simplifications)

//A.Rubeo: further modifications,
//adding cb parameter and computing the derivative of the obs with respect to it
//tau shift
//second derivative improved

//Note:
//in this program cbMat = Id4 + c_b * wilsonKernel;
//plcl is the imprved combination for the magnetic component.

// compiling: g++ -O2 NmagNew.cpp -o NmagNew -larmadillo


#include <iomanip>
#include <cmath>

#include <iostream>
#include <armadillo>

using namespace std;

//global variables
char DELIMITER = ' ';

int L = 0;
int T = 0;
double c_b = 0.0;
double tau = 0.0;
double E_cont_mag = 0.0;
double secder_cont_mag = 0.0;
double E_cont_el = 0.0;
double x_0_T = 0.0;

//momenta
double p_1 = 0.0;
double p_2 = 0.0;
double p_3 = 0.0;
double p_0 = 0.0;

//momenta on the lattice (p hat)
double p_hat_0 = 0.0;
double p_hat_1 = 0.0;
double p_hat_2 = 0.0;
double p_hat_3 = 0.0;

// p hat squared
double p_hat_2nd = 0;
double p_hat_2nds = 0;

//definition to make shorter expression for zeuthen flow
double prefact0= 0;
double prefact1= 0;
double prefact2= 0;
double prefact3= 0;

// p hat to the fourth
double p_hat_4th = 0.0;
double p_hat_4ths = 0.0;

//product of momenta on the lattice

double p0_2nd = 0.0;
double p1_2nd = 0.0;
double p2_2nd = 0.0;
double p3_2nd = 0.0;
double p0_4th = 0.0;
double p1_4th = 0.0;
double p2_4th = 0.0;
double p3_4th = 0.0;

double p0p1 = 0.0;
double p1p0 = 0.0;
double p0p2 = 0.0;
double p2p0 = 0.0;
double p0p3 = 0.0;
double p3p0 = 0.0;
double p1p2 = 0.0;
double p2p1 = 0.0;
double p1p3 = 0.0;
double p3p1 = 0.0;
double p2p3 = 0.0;
double p3p2 = 0.0;

//gauge parameter for heat kernel (should be positive)
double alpha = 1.0;
//gauge parameter for action
double lambda = 1.0;

double Nmag[100];
double Emag[100];
double secder_Emag[100];
double secder_Nmag[100];
double firstder_Emag[100];
double firstder_Nmag[100];

void init_PmuPnu (arma::mat& P_muP_nu);
void init_cbMat (arma::mat& cbMat,arma::mat& P_muP_nu, double c_b);


void update_p (int n_0, int n_1, int n_2, int n_3)
{

    p_1= (2*M_PI/L) * n_1;
    p_2= (2*M_PI/L) * n_2;
    p_3= (2*M_PI/L) * n_3;
    p_0= (M_PI/T)   * n_0;


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

int main(int argc, char *argv[])
{
//    cout << argc << argv[0] << endl;
    if (argc != 9 ) {
        cout << "Usage: Nmag ob fl ac c L/a T/a cb tau" << endl;
        cout << "EXAMPLE: Nmag plcl z pl 0.3 16 16 0.1 0" << endl;
        cout << "ob = observable, options: pl,cl,plcl,lw " << endl;
        cout << "fl=flow, options: fl = w,z,s " << endl;
        cout << "ac=action, options: ac = pl,lw " << endl;
        cout << "c= sqrt(8t)/L  " << endl;
        cout << "L/a = spatial lattice size " << endl;
        cout << "T/a = temporal lattice size " << endl;
        cout << "cb = imporvement coeff. in the initial condition for the flow "<< endl;
        cout << "tau =(shift in flow time)" << endl;
        exit(1);
    }

    // input parameters
    string ob = argv[1];   // options pl,cl,lw,plcl
    string fl = argv[2];   // options z or w
    string ac = argv[3];   // options pl or lw
    double c = (double) atof(argv[4]);
    L = (double) atof(argv[5]);
    T = (double) atof(argv[6]);
    c_b = (double) atof(argv[7]); // imporvement coeff. in the initial condition for the flow
    tau = (double) atof(argv[8]); //tau shift

    // flow time from c
    double t= c*c*L*L/8.0;

    // declaration of various matrices:

    // matrix for Wilson kernel
    arma::mat K_pl;

    // generic placeholder matrices
    arma::mat Kob;
    arma::mat Dprop;
    arma::mat Hker;
    arma::mat Dbar;

    // projectors into transverse and longitudinal modes
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

    //matrices with cb dependence
    arma::mat D_wil_cb;
    arma::mat de_D_wil_bar_z_de_cb;
    arma::mat de_D_wil_bar_w_de_cb;
    arma::mat de_D_wil_bar_s_de_cb;
    arma::mat de_D_wil_cb;
    arma::mat D_lw_cb(4,4);
    arma::mat D_lw_cb_act(4,4);
    arma::mat de_D_lw_bar_s_de_cb;
    arma::mat de_D_lw_bar_z_de_cb;
    arma::mat de_D_lw_bar_w_de_cb;
    arma::mat de_D_lw_cb;

    arma::mat cbMat(4,4);

    arma::mat P_muP_nu(4,4);

    
    Id4.eye(4,4);


    for (int x=1;x<T;x++){
         Emag[x]=0.0;
         Nmag[x]=0.0;
         secder_Emag[x]=0.0;
         secder_Nmag[x]=0.0;
         firstder_Emag[x]=0.0;
         firstder_Nmag[x]=0.0;
    }


    for(int n_0 =1; n_0< T; n_0++){
        for(int n_1 =0; n_1< L; n_1++){
            for(int n_2 =n_1; n_2< L; n_2++){
                for(int n_3 =n_2; n_3< L; n_3++){

                    //updating p values for each n_0 n_1 n_2 n_3
                    update_p (n_0, n_1, n_2, n_3);

                    int degfact = 0;
                    if ( (n_1 != n_2) && (n_1 != n_3) && (n_2 != n_3) ){ degfact = 6;}
                    else if ((n_1 == n_2) && (n_1 == n_3) && (n_2 == n_3)){ degfact = 1;}
                    else { degfact = 3;}


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
                    << p_hat_2nd + (1./12.)*(p_hat_4th+ p_hat_2nd*p0_2nd - 2.0 * p0_4th)
                    << -p0p1 +(1./12.)*(-p0p1*(p0_2nd+p1_2nd))
                    << -p0p2 +(1./12.)*(-p0p2*(p0_2nd+p2_2nd))
                    << -p0p3 +(1./12.)*(-p0p3*(p0_2nd+p3_2nd))
                    << arma::endr
                    // 2nd row
                    << -p1p0 +(1./12.)*(-p1p0*(p1_2nd+p0_2nd))
                    << p_hat_2nd  + (1./12.)*(p_hat_4th+ p_hat_2nd*p1_2nd - 2.0 * p1_4th)
                    << -p1p2 +(1./12.)*(-p1p2*(p1_2nd+p2_2nd))
                    << -p1p3 +(1./12.)*(-p1p3*(p1_2nd+p3_2nd))
                    << arma::endr
                    // 3rd row
                    << -p2p0 +(1./12.)*(-p2p0*(p2_2nd+p0_2nd))
                    << -p2p1 +(1./12.)*(-p2p1*(p2_2nd+p1_2nd))
                    << p_hat_2nd  + (1./12.)*(p_hat_4th+ p_hat_2nd*p2_2nd - 2.0 * p2_4th)
                    << -p2p3 +(1./12.)*(-p2p3*(p2_2nd+p3_2nd))
                    << arma::endr
                    // 4th row
                    << -p3p0 +(1./12.)*(-p3p0*(p3_2nd+p0_2nd))
                    << -p3p1 +(1./12.)*(-p3p1*(p3_2nd+p1_2nd))
                    << -p3p2 +(1./12.)*(-p3p2*(p3_2nd+p2_2nd))
                    << p_hat_2nd  + (1./12.)*(p_hat_4th+ p_hat_2nd*p3_2nd - 2.0 * p3_4th)
                    << arma::endr;

                    // Zeuthen kernel
                    K_z = Corrfact * K_lw + alpha * p_hat_2nd * LONG;

                    K_lw_lam = K_lw + lambda*p_hat_2nd*LONG; // (to define LW propagator)
                    K_lw = K_lw + alpha*p_hat_2nd*LONG;  // (to define Symanzik flow)

                    // Kernels for observables
                    K_lw_obs  //(3x3) matrix for magnetic observable)
                    // 1st row
                    << p_hat_2nds - p1_2nd + (1./12.)*(p_hat_4ths+ p_hat_2nds*p1_2nd - 2.0 * p1_4th)
                    << -p1p2 +(1./12.)*(-p1p2*(p1_2nd+p2_2nd))
                    << -p1p3 +(1./12.)*(-p1p3*(p1_2nd+p3_2nd))
                    << arma::endr
                    // 2nd row
                    << -p2p1 +(1./12.)*(-p2p1*(p2_2nd+p1_2nd))
                    << p_hat_2nds - p2_2nd + (1./12.)*(p_hat_4ths+ p_hat_2nds*p2_2nd - 2.0 * p2_4th)
                    << -p2p3 +(1./12.)*(-p2p3*(p2_2nd+p3_2nd))
                    << arma::endr
                    // 3rd row
                    << -p3p1 +(1./12.)*(-p3p1*(p3_2nd+p1_2nd))
                    << -p3p2 +(1./12.)*(-p3p2*(p3_2nd+p2_2nd))
                    << p_hat_2nds - p3_2nd + (1./12.)*(p_hat_4ths+ p_hat_2nds*p3_2nd - 2.0 * p3_4th)
                    << arma::endr;

                    K_cl_obs   //(3x3) matrix for magnetic observable)
                    // 1st row
                    << (p_hat_2nds - p_hat_4ths/4.) * (1 - p1_2nd/4.) - p1_2nd*(1 - p1_2nd/4.)* (1 - p1_2nd/4.)
                    << -p1p2 *(1 - p1_2nd/4.)* (1 - p2_2nd/4.)
                    << -p1p3 *(1 - p1_2nd/4.)* (1 - p3_2nd/4.)
                    << arma::endr
                    // 2nd row
                    << -p2p1 *(1 - p2_2nd/4.)* (1 - p1_2nd/4.)
                    << (p_hat_2nds - p_hat_4ths/4.) * (1 - p2_2nd/4.) - p2_2nd*(1 - p2_2nd/4.)* (1 - p2_2nd/4.)
                    << -p2p3 *(1 - p2_2nd/4.)* (1 - p3_2nd/4.)
                    << arma::endr
                    // 3rd row
                    << -p3p1 *(1 - p3_2nd/4.)* (1 - p1_2nd/4.)
                    << -p3p2 *(1 - p3_2nd/4.)* (1 - p2_2nd/4.)
                    << (p_hat_2nds - p_hat_4ths/4.) * (1 - p3_2nd/4.) - p3_2nd*(1 - p3_2nd/4.)* (1 - p3_2nd/4.)
                    << arma::endr;


                    K_pl_obs  //(3x3) matrix for magnetic observable)
                    <<  p_hat_2nds - p1_2nd << -p1p2                 << -p1p3                << arma::endr
                    << -p2p1                <<  p_hat_2nds - p2_2nd  << -p2p3                << arma::endr
                    << -p3p1                << -p3p2                 <<  p_hat_2nds - p3_2nd << arma::endr;

                    K_plcl_obs = (4./3.)*K_pl_obs - (1./3.)*K_cl_obs;

                    init_PmuPnu (P_muP_nu);

                    //initialisation matrix c_b dependence
                    init_cbMat(cbMat,P_muP_nu, c_b);


                    if (ac == "pl" ){Dprop = (1.0/ p_hat_2nd) * (TRAN + (1/lambda)* LONG);}
                    else if (ac == "lw" ){Dprop = (K_lw_lam.i());}   // inversion could be avoided using ML's formula (careful with gauge fixing!)
                    else{
                        cout << "Options for action are pl (plaq) or lw (Luscher-Weisz), not: " << ac << endl;
                        return 0;
                        }
                    //tau shift in the flow : t --> (t + tau)
                    if (fl == "w" ){Hker = exp(- (t + tau)  * p_hat_2nd) * TRAN + exp(- (t + tau)* alpha* p_hat_2nd) * LONG;}
                    else if (fl == "z" ){Hker = expmat(-(t + tau) * K_z);}
                    else if (fl == "s" ){Hker = expmat(-(t + tau) * K_lw);}
                    else{ cout << "Options for flow are w (Wilson), z (Zeuthen) or s (Symanzik), not: " << fl << endl;
                        return 0;
                        }

                    //Dbar = Hker * Dprop * (Hker.t()); //(this line does not include the cb dependence)
                    Dbar = Hker * cbMat * Dprop * (Hker.t()) * (cbMat.t());

                    if (ob == "pl"){Kob = K_pl_obs;}
                    else if (ob == "cl"){Kob = K_cl_obs;}
                    else if (ob == "lw"){Kob = K_lw_obs;}
                    else if (ob == "plcl"){Kob = K_plcl_obs;}
                    else{ cout << "Options for observable are pl (plaquette), cl (clover), lw (Luscher-Weisz) plcl (plaq-clover), not: " << ob << endl;
                        return 0;
                        }

                    // contraction of the spatial indices to get magnetic term
                    // need to shift indices in Kob by 1 to account for C++ counting from 0
                    double Emag_p=0.0;

                    for(int i =1; i< 4; i++){
                        for(int j =1; j< 4; j++){
                            Emag_p += Kob.at(i-1,j-1)*Dbar.at(j,i);
                        }
                    }
                    for (int x=1; x<T; x++){
                    Emag[x] +=  degfact*sin(p_0*x )*sin(p_0*x)* Emag_p;
                    //first derivative
                    firstder_Emag[x] +=  degfact* Emag_p* (1./2.)* ( sin(p_0*(x+1) )*sin(p_0*(x+1)) - sin(p_0*(x-1) )*sin(p_0*(x-1)) );
                    //symm (standard) second derivative on the lattice
                    //secder_Emag[x] += degfact* Emag_p* ((sin (p_0* (x-1)) * sin (p_0* (x-1))) + (sin (p_0* (x+1)) * sin (p_0* (x+1))) - 2. *sin(p_0*x )*sin(p_0*x));

                    // improved second derivative
                    secder_Emag[x] +=  degfact* Emag_p*( double(-1./12.)*sin(p_0*(x+2)) * sin(p_0*(x+2)) \
                                                              + (4./3.)* sin(p_0*(x+1)) * sin(p_0*(x+1)) \
                                                              -(5./2.)*(sin(p_0*x )*sin(p_0*x)) \
                                                              + (4./3.)*(sin(p_0*(x-1))* sin(p_0*(x-1))) \
                                                              - (1./12.)*(sin(p_0*(x-2))*sin(p_0*(x-2))));

                    }
                }
            }
        }
    }
    for (int x=1; x<T; x++){
        Emag[x] = Emag[x] *(8.0/ ((double)T*L*L*L));
        Nmag[x] = t*t*Emag[x];
        secder_Emag[x] = secder_Emag[x] *(8.0/ ((double)T*L));
        secder_Nmag[x] = t*t*secder_Emag[x];
        firstder_Emag[x] = firstder_Emag[x] *(8.0/ ((double)T*L));
        firstder_Nmag[x] = t*t*firstder_Emag[x];

        x_0_T = (double)x/(double)T;

        if (x_0_T == 0.5) {
            //continuum values hardcoded from cont_energyDen.cpp (depending on the value of c)
            if (c == 0.3){

                E_cont_mag = 0.0085637412877;
                secder_cont_mag = -0.01709793688072;

            }
            else if (c == 0.4){

                E_cont_mag = 0.00671373472268;
                secder_cont_mag = -0.04326680245534;

            }
            else if (c == 0.2){

                E_cont_mag = 0.00931389344276;
                secder_cont_mag = -0.003402671;

            }
            //c=0.5
            else{

                E_cont_mag = 0.00413266156849;
                secder_cont_mag = -1.3193308065547363e-01;

            }
            if (c !=0.2 && c !=0.3 && c != 0.4 && c != 0.5 ){
                cout << "computing only c=0.2, 0.3, 0.4, 0.5" << endl;
                exit(1);
            }
        }
        if (x_0_T == 0.125) {
            //continuum values hardcoded from cont_energyDen.cpp (depending on the value of c)
            if (c == 0.3){

                E_cont_mag = 3.2947533904262327e-03;
                secder_cont_mag = -9.5319840480145614e-03;

            }
        }

        if (x_0_T == 0.25) {
            //continuum values hardcoded from cont_energyDen.cpp (depending on the value of c)
            if (c == 0.3){

                E_cont_mag = 0.00719900196433;
                secder_cont_mag = -1.9888715285744810e-01;

            }
            else if (c == 0.4){

                E_cont_mag = 0.00471420819893;
                secder_cont_mag = -1.0549645680401466e-01;

            }
            else if (c == 0.2){

                E_cont_mag = 0.00893688246541;
                secder_cont_mag = -1.3193308065547363e-01;

            }
            //c=0.5
            else{

                E_cont_mag = 0.00248499321041;
                cout << "WARNING: missing continuum value for sec derivative at c=0.5 x0/T=0.25 (cont_energyDen)" << endl;

            }
            if (c !=0.2 && c !=0.3 && c != 0.4 && c != 0.5 ){
                cout << "computing only c=0.2, 0.3, 0.4, 0.5" << endl;
                exit(1);
            }

        }


        if (x_0_T == 0.375) {

            if (c == 0.3){

                E_cont_mag = 8.3878564019138520e-03;
                secder_cont_mag = -5.4180065662296654e-02;

            }
        }


        //print
        if(x_0_T == 0.5){
            cout << setfill('0') << setw(2) << L
            << DELIMITER << x_0_T
            << DELIMITER << c_b
            << DELIMITER << fixed << setprecision(16) << scientific << Nmag[x]
            //<< DELIMITER << fixed << setprecision(16) << scientific << ((Nmag[x]/E_cont_mag))
            << DELIMITER << fixed << setprecision(16) << scientific << secder_Nmag[x]
            //<< DELIMITER << fixed << setprecision(16) << scientific << ((secder_Nmag[x]/secder_cont_mag))
            //<< DELIMITER << setfill('0') << setw(2) << T
            << endl;

        }

    }
    return 0;
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

    wilsonKernel = p_hat_2nd*Id4 - P_muP_nu;

    //cbMat = Id4 + c_b * wilsonKernel;
    cbMat = expmat(c_b * wilsonKernel);

}
