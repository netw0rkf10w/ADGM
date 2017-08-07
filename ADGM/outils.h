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
 
#ifndef OUTILS_H
#define OUTILS_H

#include <iostream>
#include <algorithm>
#include <assert.h>
#include <vector>
#include <Eigen/Dense>

#include "hungarian.h"
#include "condat_simplexproj.h"

using namespace std;
using namespace Eigen;


double* normalize(const double* v, int n)
{
    double* v1 = new double[n];

    double vmax = 0.0;
    for(size_t i = 0; i < n; i++){
        if(vmax < abs(v[i]))
            vmax = abs(v[i]);
    }
    for(size_t i = 0; i < n; i++){
        v1[i] = v[i]/vmax;
    }
    return v1;
}


void normalize(double* &v, int n)
{
    double vmax = 0.0;
    for(size_t i = 0; i < n; i++){
        if(vmax < abs(v[i]))
            vmax = abs(v[i]);
    }
    for(size_t i = 0; i < n; i++){
        v[i] = v[i]/vmax;
    }
}


double** vector_to_matrix(const VectorXd &v, int rows, int cols) {
    assert((int)v.size() == rows*cols);
    double** M;
    M = (double**)calloc(rows,sizeof(double*));
    for(int i = 0; i < rows; i++){
        M[i] = (double*)calloc(cols,sizeof(double));
        for(int j = 0; j < cols; j++)
            M[i][j] = v[i + rows*j];
    }
  return M;
}


VectorXd discritize(const VectorXd&X, int n1, int n2) {
    double** matrix = vector_to_matrix(X, n1, n2);
    VectorXd Y = VectorXd::Zero(X.size());
    hungarian_problem_t hunger;
    hungarian_init(&hunger, matrix, n1, n2, HUNGARIAN_MODE_MAXIMIZE_UTIL) ;
    hungarian_solve(&hunger);
    for(int i = 0; i < n1; i++){
        for(int j = 0; j < n2; j++)
            Y(i + n1*j) = (hunger.assignment[i][j] > 0)?1:0;
    }
    hungarian_free(&hunger);
    if (matrix != NULL) {
        for(int i = 0; i < n1; i++){
            if (matrix[i] != NULL)  free(matrix[i]);
        }
        free(matrix);
    }
    return Y;
}

double computeEnergy(const VectorXd &X, const int* indH3, const double* valH3, int Nt3)
{
    double energy = 0.0;
    for(int i = 0; i < Nt3; i++){
        int i1 = indH3[i];
        int i2 = indH3[i+Nt3];
        int i3 = indH3[i+2*Nt3];
        if(X(i1) > 0 && X(i2) > 0 && X(i3) > 0)
            energy += valH3[i]*X(i1)*X(i2)*X(i3);
        }
    return energy;
}


double computeEnergy_2nd(const VectorXd &x, const MatrixXd &M, VectorXd &d)
{
    double energy = 0.0;
    for(int i = 0; i < x.size(); i++){
        if(x(i) > 0)
            energy += x(i)*d(i);
    }
    for(int i = 0; i < x.size(); i++){
        for(int j = 0; j < x.size(); j++){
        if(j != i && x(i) > 0 && x(j) > 0)
            energy += M(i,j)*x(i)*x(j);
        }
    }
    return energy;
}


VectorXd sublemma_EQUALITY(const VectorXd &c)
/// Solve   min     1/2*||x||^2 - c^T*x
///         subject to      1^T*x = 1 and x >= 0.
/// Input c: n*1 vector
/// Output x: n*1 vector
{
    const double* y = c.data();
    unsigned int n = c.size();
    double *z = new double[n];
    simplexproj_Condat(y, z, n, 1.0);
    Map<VectorXd> x(z, n); // no need to free z
    return x;
}


VectorXd sublemma_INEQUALITY(const VectorXd &c)
/// Solve   min     1/2*||x||^2 - c^T*x
///         subject to      1^T*x <= 1 and x >= 0.
/// Input c: n*1 vector
/// Output x: n*1 vector
{
    VectorXd x = VectorXd::Zero(c.size(), 1);
    double s = 0.0;
    for(size_t i = 0; i < x.size(); i++)
    {
        if(c(i) > 0)
        {
            x(i)= c(i);
            s += c(i);
        }
    }
    if(s <= 1)
        return x;

//    VectorXd x = (c.array() > 0).select(c, 0);
//    double  s = x.sum();
//    if(s <= 1)
//        return x;
    return sublemma_EQUALITY(c);
}


MatrixXd lemma_EQUALITY(const MatrixXd &C)
/// Solve for each COLUMN   min     1/2*||x||^2 - c^T*x
///         subject to      1^T*x = 1 and x >= 0.
/// Input c: n*1 vector
/// Output x: n*1 vector
{
    MatrixXd X(C.rows(), C.cols());
    //#pragma omp parallel for
    for(int i = 0; i < X.cols(); i++){
        X.col(i) = sublemma_EQUALITY(C.col(i));
    }
    return X;
}


MatrixXd lemma_INEQUALITY(const MatrixXd &C)
/// Solve for each COLUMN   min     1/2*||x||^2 - c^T*x
///         subject to      1^T*x = 1 and x >= 0.
/// Input c: n*1 vector
/// Output x: n*1 vector
{
    // Thresholding C: if negative then set to 0
    //MatrixXd X = (C.array() > 0).select(C, 0);
    MatrixXd X = C.cwiseMax(0.0); // cheaper than (...).select(C,0)

    // Compute the sum of each column of P (i.e. the
    // sum of positive elements in C)
    MatrixXd s = X.colwise().sum();

    for(int j = 0; j < s.cols(); j++){
        if(s(j) > 1)
            X.col(j) = sublemma_EQUALITY(C.col(j));
    }
    return X;
}

#endif // OUTILS_H
