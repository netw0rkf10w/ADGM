/*************************
C++ MEX-file for the Alternating Direction Graph Matching (ADGM) algorithm
Written by D. Khue Le-Huu
khue.le@centralesupelec.fr
September 2016

Matlab call: [x, residuals, rho] = mexADGM_3rdORDER_SYMMETRIC(x0, indH1, valH1, indH2, valH2, indH3, valH3, rhos, MAX_ITER, verb, restart, iter1, iter2, variant)
**************************/

#include "mex.h"
#include "ADGM.h"

using namespace std;


///////////////////////////////////////////////////////////////////////////
// Check inputs
///////////////////////////////////////////////////////////////////////////
void checkInputs(int nrhs, const mxArray *prhs[])
{
    const int * tdims, * fdims;

    // Check number of inputs
    if (nrhs < 14)
    {
        mexErrMsgTxt("Incorrect number of inputs. Function expects at least 15 inputs.");
    }

    if (mxGetNumberOfDimensions(prhs[0])>2)
    {
        mexErrMsgTxt("Incorrect number of dimensions. First input must be a matrix.");
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
		
	// Check inputs to mex function
    checkInputs(nrhs, prhs);
	
	// read the input
	double* x0 = mxGetPr(prhs[0]);
	mwSize* indH1 = (mwSize*)mxGetPr(prhs[1]);
	double* valH1 = mxGetPr(prhs[2]);
	mwSize* indH2 = (mwSize*)mxGetPr(prhs[3]);
	double* valH2 = mxGetPr(prhs[4]);
	mwSize* indH3 = (mwSize*)mxGetPr(prhs[5]);
	double* valH3 = mxGetPr(prhs[6]);
	double* rhos = mxGetPr(prhs[7]);
	mwSize MAX_ITER = mxGetScalar(prhs[8]);
	bool verb = mxGetScalar(prhs[9]);
	bool restart = mxGetScalar(prhs[10]);
	mwSize iter1 = mxGetScalar(prhs[11]);
	mwSize iter2 = mxGetScalar(prhs[12]);
	mwSize variant = mxGetScalar(prhs[13]);
	if(verb)
		cout<<"Started C++ ADGM..."<<endl;	
	
	// get the dimensions of the input matrices
	mwSize n1 = mxGetM(prhs[0]);
	mwSize n2 = mxGetN(prhs[0]);
	mwSize Nt1 = mxGetM(prhs[1]);
	mwSize Nt2 = mxGetM(prhs[3]);
	mwSize Nt3 = mxGetM(prhs[5]);
	mwSize n_rho = max(mxGetN(prhs[7]), mxGetM(prhs[7]));
	
	// Create the output matrix and get a pointer to the real data in that matrix
	//plhs[0] = mxCreateDoubleMatrix(n1,n2,mxREAL);
	//double* x = mxGetPr(plhs[0]);
	double* x; 
	double* residuals;
	double rho; // the final rho will be returned as well
	int N = ADGM_3rdORDER_SYMMETRIC(x, residuals, rho, x0, n1, n2, indH3, valH3, Nt3,
                    rhos, n_rho, MAX_ITER, verb, restart, iter1, iter2, variant);
    
    // Put the data back into the output MATLAB array
    int A = n1*n2;
    mxArray * mx1 = mxCreateDoubleMatrix(A, 1, mxREAL);
    std::copy(x, x + A, mxGetPr(mx1));
    plhs[0] = mx1;
    mxArray * mx = mxCreateDoubleMatrix(N, 1, mxREAL);
    std::copy(residuals, residuals + N, mxGetPr(mx));
    plhs[1] = mx;
    plhs[2] = mxCreateDoubleScalar(rho);
    
}
