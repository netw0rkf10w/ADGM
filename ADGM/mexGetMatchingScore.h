/// D. Khuê Lê-Huu on October 19, 2016
double getMatchingScore(const double* X, const int* indH, const double* valH, int N)
{
    double score = 0.0;
    for(int i = 0; i < N; i++){
    	int i1 = indH[i];
    	int i2 = indH[i+N];
    	int i3 = indH[i+2*N];
    	if(X[i1] > 0 && X[i2] > 0 && X[i3] > 0)
        	score += valH[i]*X[i1]*X[i2]*X[i3];
        }
    return score;
}
