/********************************************************************
 ********************************************************************
 Alternating Direction Graph Matching
 Version: 0.1 (pre-release)
 Written by: D. Khuê Lê-Huu (khue.le@centralesupelec.fr)
 
 If you use any part of this code, please cite:
 
D. Khuê Lê-Huu and Nikos Paragios. Alternating Direction Graph Matching. IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2017.

BibTeX:

@inproceedings{lehuu2017adgm,
 title={Alternating Direction Graph Matching},
 author={L{\^e}-Huu, D. Khu{\^e} and Paragios, Nikos},
 booktitle = {Proceedings of the {IEEE} Conference on Computer Vision and Pattern Recognition ({CVPR})},
 year = {2017}
}


*********************
NOTE: 

- This is a preliminary re-implementation in C++ Eigen, 
I haven't tested it on the full benchmark yet and thus the performance 
is not guaranteed (in particular some caching part has not been implemented
yet, so it should be slower). Before using this software, please check if 
any updated version is available on my website www.khue.fr.
- In the current version, only third-order potentials are supported.

 ** This file is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied 
 ** warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 **
 ********************************************************************
 ********************************************************************/
 
#include <Eigen/Dense>

#include "outils.h"

using namespace std;
using namespace Eigen;

#define NONE 0
#define EQUALITY 1
#define INEQUALITY 2
#define ROW 1
#define COL 2
#define VARIANT_1 1
#define VARIANT_2 2


void update_step_3rd(int variant, double rho, int n1, int n2, int Nt3, const double* valH3,
                     const int* I1, const int* I2, const int* I3,
                     VectorXd &x1, VectorXd &x2, VectorXd &x3, VectorXd &y2, VectorXd &y3,
                     double &r, double &s)
{
    int A = n1*n2;
    if(variant == VARIANT_1){
        /// Step 1: update x1 = argmin 0.5||x||^2 - c1^T*x where
        /// Compute c1 = 0.5*(x2 + x3) - (y2 + y3 + Fx2x3)/(2.0*rho);
        VectorXd Fx2x3 = VectorXd::Zero(A);
        for(int ind = 0; ind < Nt3; ind++){
            if(x2[I2[ind]] > 0 && x3[I3[ind]] > 0)
                Fx2x3(I1[ind]) += valH3[ind]*x2[I2[ind]]*x3[I3[ind]];
        }

        VectorXd c1 = 0.5*(x2 + x3) - (y2 + y3 + Fx2x3)/(2.0*rho);
        Map<MatrixXd> C1(c1.data(), n1, n2);

        VectorXd x1_old = x1;
        MatrixXd X1T;
        if(n1 <= n2)
            X1T = lemma_EQUALITY(C1.transpose()); // Transpose because the lemma solves for each column
        else
            X1T = lemma_INEQUALITY(C1.transpose());

        Matrix<double,Dynamic,Dynamic,RowMajor> X1(X1T);
        x1 = VectorXd::Map(X1.data(), A);

        /// Step 2: update x2 = argmin 0.5||x||^2 - c2^T*x where
        /// c2 = x1 + (y2 - Fx1x3)/rho;
        VectorXd Fx1x3 = VectorXd::Zero(A);
        for(int ind = 0; ind < Nt3; ind++){
            if(x1[I1[ind]] > 0 && x3[I3[ind]] > 0)
                Fx1x3(I2[ind]) += valH3[ind]*x1[I1[ind]]*x3[I3[ind]];
        }
        VectorXd c2 = x1 + (y2 - Fx1x3)/rho;
        Map<MatrixXd> C2(c2.data(), n1, n2);

        VectorXd x2_old = x2;
        MatrixXd X2;
        if(n2 <= n1)
            X2 = lemma_EQUALITY(C2);
        else
            X2 = lemma_INEQUALITY(C2);

        x2 = VectorXd::Map(X2.data(), A);


        /// Step 3: update x3 = argmin 0.5||x||^2 - c2^T*x where
        /// c3 = x1 + (y3 - Fx1x2)/rho;
        VectorXd Fx1x2 = VectorXd::Zero(A);
        for(int ind = 0; ind < Nt3; ind++){
            if(x1[I1[ind]] > 0 && x2[I2[ind]] > 0)
                Fx1x2(I3[ind]) += valH3[ind]*x1[I1[ind]]*x2[I2[ind]];

        }
        VectorXd c3 = x1 + (y3 - Fx1x2)/rho;
        Map<MatrixXd> C3(c3.data(), n1, n2);

        VectorXd x3_old = x3;
        if(n1 <= n2){
            MatrixXd X3T = lemma_EQUALITY(C3.transpose());
            Matrix<double,Dynamic,Dynamic,RowMajor> X3(X3T);
            x3 = VectorXd::Map(X3.data(), A);
        }
        else{
            MatrixXd X3 = lemma_EQUALITY(C3);
            x3 = VectorXd::Map(X3.data(), A);
        }


        /// Step 4: update y
        y2 += rho*(x1 - x2);
        y3 += rho*(x1 - x3);


        /// Step 5: compute the residuals and update rho

        r = (x1 - x2).squaredNorm() + (x1 - x3).squaredNorm();
        s = (x1 - x1_old).squaredNorm() + (x2 - x2_old).squaredNorm() + (x3 - x3_old).squaredNorm();
    }
    else{
        /// Step 1: update x1 = argmin 0.5||x||^2 - c1^T*x where
        /// Compute c1 = x2 - (y2 + Fx2x3)/(rho);
        VectorXd Fx2x3 = VectorXd::Zero(A);
        for(int ind = 0; ind < Nt3; ind++){
            if(x2[I2[ind]] > 0 && x3[I3[ind]] > 0)
                Fx2x3(I1[ind]) += valH3[ind]*x2[I2[ind]]*x3[I3[ind]];
        }

        VectorXd c1 = x2 - (y2 + Fx2x3)/rho;
        Map<MatrixXd> C1(c1.data(), n1, n2);

        VectorXd x1_old = x1;
        MatrixXd X1T;
        if(n1 <= n2)
            X1T = lemma_EQUALITY(C1.transpose()); // Transpose because the lemma solves for each column
        else
            X1T = lemma_INEQUALITY(C1.transpose());

        Matrix<double,Dynamic,Dynamic,RowMajor> X1(X1T);
        x1 = VectorXd::Map(X1.data(), A);

        /// Step 2: update x2 = argmin 0.5||x||^2 - c2^T*x where
        /// Compute c2 = 0.5*(x1 + x3) + (y2 - y3 - Fx1x3)/(2*rho);
        VectorXd Fx1x3 = VectorXd::Zero(A);
        for(int ind = 0; ind < Nt3; ind++){
            if(x1[I1[ind]] > 0 && x3[I3[ind]] > 0)
                Fx1x3(I2[ind]) += valH3[ind]*x1[I1[ind]]*x3[I3[ind]];
        }
        VectorXd c2 = 0.5*(x1 + x3) + (y2 - y3 - Fx1x3)/(2.0*rho);
        Map<MatrixXd> C2(c2.data(), n1, n2);

        VectorXd x2_old = x2;
        MatrixXd X2;
        if(n2 <= n1)
            X2 = lemma_EQUALITY(C2);
        else
            X2 = lemma_INEQUALITY(C2);

        x2 = VectorXd::Map(X2.data(), A);


        /// Step 3: update x3 = argmin 0.5||x||^2 - c2^T*x where
        /// c3 = x2 + (y3 - Fx1x2)/rho;
        VectorXd Fx1x2 = VectorXd::Zero(A);
        for(int ind = 0; ind < Nt3; ind++){
            if(x1[I1[ind]] > 0 && x2[I2[ind]] > 0)
                Fx1x2(I3[ind]) += valH3[ind]*x1[I1[ind]]*x2[I2[ind]];

        }
        VectorXd c3 = x2 + (y3 - Fx1x2)/rho;
        Map<MatrixXd> C3(c3.data(), n1, n2);

        VectorXd x3_old = x3;
        if(n1 <= n2){
            MatrixXd X3T = lemma_EQUALITY(C3.transpose());
            Matrix<double,Dynamic,Dynamic,RowMajor> X3(X3T);
            x3 = VectorXd::Map(X3.data(), A);
        }
        else{
            MatrixXd X3 = lemma_EQUALITY(C3);
            x3 = VectorXd::Map(X3.data(), A);
        }

        /// Step 4: update y
        y2 += rho*(x1 - x2);
        y3 += rho*(x2 - x3);

        /// Step 5: compute the residuals and update rho
        r = (x1 - x2).squaredNorm() + (x2 - x3).squaredNorm();
        s = (x1 - x1_old).squaredNorm() + (x2 - x2_old).squaredNorm() + (x3 - x3_old).squaredNorm();
    }
}


int ADGM_3rdORDER(double* &x, double* &residuals, double &rho, const double* _x0, int n1, int n2,
                                const int* indH3, const double* _valH3, int Nt3, const double* rhos, int nrho,
                                int MAX_ITER, bool verb, bool restart, int iter1, int iter2, int variant)
/// Return the number of iterations
{
    assert(variant == VARIANT_1 || variant == VARIANT_2);
    //bool debug = false;
    int A = n1*n2;


    /// Normalize, and extract the diagonal
    double* valH3 = normalize(_valH3, Nt3);
    const int* I1 = &indH3[0];
    const int* I2 = &indH3[Nt3];
    const int* I3 = &indH3[2*Nt3];

    /// Initialization
    VectorXd x1 = VectorXd::Map(_x0, A) ;
    VectorXd x2 = VectorXd::Map(_x0, A) ;
    VectorXd x3 = VectorXd::Map(_x0, A) ;
    VectorXd y2 = VectorXd::Zero(A,1);
    VectorXd y3 = VectorXd::Zero(A,1);

    std::vector<double> res(MAX_ITER, 0.0);

    int iter;

    double res_best_so_far = DBL_MAX;

    int counter = 0;
    int i_rho = 0;
    double energy_best = DBL_MAX;
    VectorXd x_best;

    int iter1_cumulated = iter1;

    rho = rhos[i_rho];
    double r, s;

    for(iter=1; iter <= MAX_ITER; iter++)
    {

        update_step_3rd(variant, rho, n1, n2, Nt3, valH3,
                        I1, I2, I3, x1, x2, x3, y2, y3, r, s);
        double R = r + s;
        res[iter - 1] = R;

        if(verb && iter%1 == 0)
            cout<<iter<<")"<<" residual = "<<R<<endl;
        /// If convergence
        if(R <= 1e-6){
            VectorXd x_bin = discritize(x3, n1, n2);
            double energy = computeEnergy(x_bin, indH3, valH3, Nt3);
            if(verb){
                cout<<"Converged!!!!!! energy = "<<energy<<" (best was "<<energy_best<<")"<<endl;
            }
            if(energy < energy_best){
                energy_best = energy;
                x_best = x_bin;
            }
            x = x_best.data();
            residuals = new double[iter];
            std::copy(res.begin(), res.begin() + iter, residuals);
            return iter;
        }

        /// Only after iter1 iterations that we start to track the best residuals
        if(iter >= iter1_cumulated){
            // res_best_so_far has changed then restart the counter
            if(R <= 0.9*res_best_so_far){
                res_best_so_far = R;
                counter = 0;
            }else{
                counter++;
            }
            // If the best_so_far residual has not changed during iter2 iterations, then update rho
            if(counter > iter2){
                i_rho++;
                if(i_rho < nrho){
                    rho = rhos[i_rho];
                    counter = 0;
                    iter1_cumulated = iter + iter1;
                    // Keep track of the current x before updating rho
                    // Convert to discrete solution and compute the energy
                    VectorXd x_bin = discritize(x3, n1, n2);
                    double energy = computeEnergy(x_bin, indH3, valH3, Nt3);
                    if(verb){
                        cout<<"!!!!!! energy = "<<energy<<" (best was "<<energy_best<<"), Update rho = "<<rho<<endl;
                    }
                    if(energy < energy_best){
                        energy_best = energy;
                        x_best = x_bin;
                    }
                }

                if(restart){
                    if(verb)
                        cout<<"Restart"<<endl;
                    x1 = VectorXd::Map(_x0, A);
                    x2 = VectorXd::Map(_x0, A);
                    x3 = VectorXd::Map(_x0, A);
                }
            }
        }
    }

    VectorXd x_bin = discritize(x3, n1, n2);
    double energy = computeEnergy(x_bin, indH3, valH3, Nt3);
    if(verb){
        cout<<"Converged!!!!!! energy = "<<energy<<" (best was "<<energy_best<<")"<<endl;
    }
    if(energy < energy_best){
        energy_best = energy;
        x_best = x_bin;
    }


    x1 = VectorXd::Map(_x0, A);
    x2 = VectorXd::Map(_x0, A);
    x3 = VectorXd::Map(_x0, A);
    y2 = VectorXd::Zero(A,1);
    y3 = VectorXd::Zero(A,1);
    for(int iter2=1; iter2 <= MAX_ITER; iter2++)
    {

        update_step_3rd(variant, rho, n1, n2, Nt3, valH3,
                        I1, I2, I3, x1, x2, x3, y2, y3, r, s);
        double R = r + s;

        if(R <= 1e-6){
            VectorXd x_bin = discritize(x3, n1, n2);
            double energy = computeEnergy(x_bin, indH3, valH3, Nt3);
            if(verb){
                cout<<"Converged!!!!!! energy = "<<energy<<" (best was "<<energy_best<<")"<<endl;
            }
            if(energy < energy_best){
                energy_best = energy;
                x_best = x_bin;
            }
            break;
        }
    }
    x = x_best.data();



    residuals = new double[iter - 1];
    std::copy(res.begin(), res.begin() + iter - 1, residuals);
    return (iter - 1);
}



int ADGM_3rdORDER_SYMMETRIC(double* &x, double* &residuals, double &rho, const double* _x0, int n1, int n2,
                                const int* _indH3, const double* _valH3, int _Nt3, const double* rhos, int nrho,
                                int MAX_ITER, bool verb, bool restart, int iter1, int iter2, int variant)
/// Return the number of iterations
{
    assert(variant == VARIANT_1 || variant == VARIANT_2);
    //bool debug = true;
    int A = n1*n2;

    /// Taking into account the super-symmetry of the tensor _valH3
    assert(_Nt3 % 6 == 0);
    const int* _I1 = &_indH3[0];
    const int* _I2 = &_indH3[_Nt3];
    const int* _I3 = &_indH3[2*_Nt3];

    int Nt3 = _Nt3/6;
    double *valH3 = new double[Nt3];
    int *indH3 = new int[3*Nt3];

    int ind = 0;
    for(int _ind = 0; _ind < _Nt3; _ind++){
        if(_I1[_ind] < _I2[_ind] && _I2[_ind] < _I3[_ind]){
            valH3[ind]          = 6.0*_valH3[_ind];
            indH3[ind]          = _I1[_ind];
            indH3[ind + Nt3]    = _I2[_ind];
            indH3[ind + 2*Nt3]  = _I3[_ind];
            ind++;
        }
    }
    assert(ind == (Nt3 -1));

    const int* I1 = &indH3[0];
    const int* I2 = &indH3[Nt3];
    const int* I3 = &indH3[2*Nt3];

    /// Normalize, and extract the diagonal
    normalize(valH3, Nt3);

    /// Initialization
    VectorXd x1 = VectorXd::Map(_x0, A) ;
    VectorXd x2 = VectorXd::Map(_x0, A) ;
    VectorXd x3 = VectorXd::Map(_x0, A) ;
    VectorXd y2 = VectorXd::Zero(A,1);
    VectorXd y3 = VectorXd::Zero(A,1);

    std::vector<double> res(MAX_ITER, 0.0);

    int iter;

    double res_best_so_far = DBL_MAX;

    int counter = 0;
    int i_rho = 0;
    double energy_best = DBL_MAX;
    VectorXd x_best;

    int iter1_cumulated = iter1;

    rho = rhos[i_rho];
    double r, s;

    for(iter=1; iter <= MAX_ITER; iter++)
    {

        update_step_3rd(variant, rho, n1, n2, Nt3, valH3,
                        I1, I2, I3, x1, x2, x3, y2, y3, r, s);
        double R = r + s;
        res[iter - 1] = R;

        if(verb && iter%1 == 0)
            cout<<iter<<")"<<" residual = "<<R<<endl;
        /// If convergence
        if(R <= 1e-6){
            VectorXd x_bin = discritize(x3, n1, n2);
            double energy = computeEnergy(x_bin, indH3, valH3, Nt3);
            if(verb){
                cout<<"Converged!!!!!! energy = "<<energy<<" (best was "<<energy_best<<")"<<endl;
            }
            if(energy < energy_best){
                energy_best = energy;
                x_best = x_bin;
            }
            x = x_best.data();
            residuals = new double[iter];
            std::copy(res.begin(), res.begin() + iter, residuals);
            return iter;
        }

        /// Only after iter1 iterations that we start to track the best residuals
        if(iter >= iter1_cumulated){
            // res_best_so_far has changed then restart the counter
            if(R <= 0.9*res_best_so_far){
                res_best_so_far = R;
                counter = 0;
            }else{
                counter++;
            }
            // If the best_so_far residual has not changed during iter2 iterations, then update rho
            if(counter > iter2){
                i_rho++;
                if(i_rho < nrho){
                    rho = rhos[i_rho];
                    counter = 0;
                    iter1_cumulated = iter + iter1;
                    // Keep track of the current x before updating rho
                    // Convert to discrete solution and compute the energy
                    VectorXd x_bin = discritize(x3, n1, n2);
                    double energy = computeEnergy(x_bin, indH3, valH3, Nt3);
                    if(verb){
                        cout<<"!!!!!! energy = "<<energy<<" (best was "<<energy_best<<"), Update rho = "<<rho<<endl;
                    }
                    if(energy < energy_best){
                        energy_best = energy;
                        x_best = x_bin;
                    }
                }

                if(restart){
                    if(verb)
                        cout<<"Restart"<<endl;
                    x1 = VectorXd::Map(_x0, A);
                    x2 = VectorXd::Map(_x0, A);
                    x3 = VectorXd::Map(_x0, A);
                }
            }
        }
    }
    
    VectorXd x_bin = discritize(x3, n1, n2);
    double energy = computeEnergy(x_bin, indH3, valH3, Nt3);
    if(verb){
        cout<<"Converged!!!!!! energy = "<<energy<<" (best was "<<energy_best<<")"<<endl;
    }
    if(energy < energy_best){
        energy_best = energy;
        x_best = x_bin;
    }


    x1 = VectorXd::Map(_x0, A);
    x2 = VectorXd::Map(_x0, A);
    x3 = VectorXd::Map(_x0, A);
    y2 = VectorXd::Zero(A,1);
    y3 = VectorXd::Zero(A,1);
    for(int iter2=1; iter2 <= MAX_ITER; iter2++)
    {

        update_step_3rd(variant, rho, n1, n2, Nt3, valH3,
                        I1, I2, I3, x1, x2, x3, y2, y3, r, s);
        double R = r + s;

        if(R <= 1e-6){
            VectorXd x_bin = discritize(x3, n1, n2);
            double energy = computeEnergy(x_bin, indH3, valH3, Nt3);
            if(verb){
                cout<<"Converged!!!!!! energy = "<<energy<<" (best was "<<energy_best<<")"<<endl;
            }
            if(energy < energy_best){
                energy_best = energy;
                x_best = x_bin;
            }
            break;
        }
    }
    x = x_best.data();


    residuals = new double[iter - 1];
    std::copy(res.begin(), res.begin() + iter - 1, residuals);
    return (iter - 1);
}
