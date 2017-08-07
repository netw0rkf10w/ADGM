#include "mex.h"
#include "mexGetMatchingScore.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double* x = mxGetPr(prhs[0]);
	mwSize* indH = (mwSize*)mxGetPr(prhs[1]);
	double* valH = mxGetPr(prhs[2]);
	mwSize N = mxGetM(prhs[1]);
	
    double score = getMatchingScore(x, indH, valH, N);
    plhs[0] = mxCreateDoubleScalar(score);
}
