#ifndef CONDAT_H
#define CONDAT_H

/*
 #  File            : condat_simplexproj.c 
 #
 #  Version History : 1.0, Aug. 15, 2014 
 #
 #  Author          : Laurent Condat, PhD, CNRS research fellow in France.
 #
 #  Description     : This file contains an implementation in the C language
 #                    of algorithms described in the research paper:
 #	
 #                    L. Condat, "Fast Projection onto the Simplex and the
 #                    l1 Ball", preprint Hal-01056171, 2014.
 #
 #                    This implementation comes with no warranty: due to the
 #                    limited number of tests performed, there may remain
 #                    bugs. In case the functions would not do what they are
 #                    supposed to do, please email the author (contact info
 #                    to be found on the web).
 #
 #                    If you use this code or parts of it for any purpose,
 #                    the author asks you to cite the paper above or, in 
 #                    that event, its published version. Please email him if 
 #                    the proposed algorithms were useful for one of your 
 #                    projects, or for any comment or suggestion.
 #
 #  Usage rights    : Copyright Laurent Condat.
 #                    This file is distributed under the terms of the CeCILL
 #                    licence (compatible with the GNU GPL), which can be
 #                    found at the URL "http://www.cecill.info".
 #
 #  This software is governed by the CeCILL license under French law and
 #  abiding by the rules of distribution of free software. You can  use,
 #  modify and or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL :
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
*/


/* This code was compiled using
gcc -march=native -O2 condat_simplexproj.c -o main  -I/usr/local/include/ 
	-lm -lgsl  -L/usr/local/lib/
On my machine, gcc is actually a link to the compiler Apple LLVM version 5.1 
(clang-503.0.40) */


/* The following functions are implemented:
1) simplexproj_algo1
2) simplexproj_algo2
3) simplexproj_algo3_pivotrand
4) simplexproj_algo3_pivotmedian
5) simplexproj_Duchi
6) simplexproj_algo4
7) simplexproj_Condat (proposed algorithm)
All these functions take the same parameters. They project the vector y onto
the closest vector x of same length (parameter N in the paper) with x[n]>=0,
n=0..N-1, and sum_{n=0}^{N-1}x[n]=a. 
We can have x==y (projection done in place). If x!=y, the arrays x and y must
not overlap, as x is used for temporary calculations before y is accessed.
We must have length>=1 and a>0. 
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define datatype double /* type of the elements in y */


/* A simple implementation of quicksort to sort in decreasing order 
The elements equal to the pivot are shared equally bewteen the two 
subvectors */
static void quicksort(datatype *y, const int lo, const int hi) {
	int i=lo, j=hi;
	datatype temp, 
	pivot=y[lo+(int)(rand()/(((double)RAND_MAX+1.0)/(hi-lo+1)))];
	for (;;) {    
		while (y[i] > pivot) i++; 
		while (y[j] < pivot) j--;
		if (i>=j) break;
		temp=y[i];
		y[i++]=y[j];
		y[j--]=temp;
	}   
	if (i-1>lo) quicksort(y, lo, i-1);
	if (hi>j+1) quicksort(y, j+1, hi);
}


/* Algorithm 1 in the paper, using sorting */
static void simplexproj_algo1(const datatype* y, datatype* x,
const unsigned int length, const double a) {
	datatype*	aux = (x==y ? (datatype*)malloc(length*sizeof(datatype)) : x);
	int 	i=1;
	memcpy(aux,y,length*sizeof(datatype));
	quicksort(aux,0,length-1);
	double	tau=*aux-a;
	while ((i<length) && (aux[i]>tau)) {
		tau+=(aux[i]-tau)/(i+1);
		i++;
	}
	for (i=0; i<length; i++)
		x[i]=(y[i]>tau ? y[i]-tau : 0.0); 
	if (x==y) free(aux);
}


static void heap_sift(int root, int lastChild, datatype* x) {
    int child;
    datatype c;
    for (; (child = (root<<1) + 1) <= lastChild; root = child) {
		if (child < lastChild)
		    if (x[child] < x[child+1]) child++;	
		if (x[child] <= x[root]) break;
		c=x[root];
		x[root]=x[child];
		x[child]=c;
    }
}


static int heap_del_max(int numElems, datatype* x) {
    int lastChild = numElems - 1;
    datatype c = *x;
    *x = x[lastChild];
    x[lastChild--] = c;
    heap_sift(0, lastChild, x);  /* restore the heap property */
    return numElems - 1;
}


/* Algorithm 2 in the paper, using a heap */
static void simplexproj_algo2(const datatype* y, datatype* x,
const unsigned int length, const double a) {
	int 	j=length, i=length>>1;
   	double 	tau;
   	datatype*	aux = (x==y ? (datatype*)malloc(length*sizeof(datatype)) : x);
	memcpy(aux,y,length*sizeof(datatype));	
   	for (; i >= 0; i--) heap_sift(i, length-1, aux);
   	tau = *aux - a;          
	j = heap_del_max(j, aux);
	i = 1;
   	while ((i<length) && (*aux>tau)) {
		tau += (*aux-tau) / ++i;
		j = heap_del_max(j, aux);
	}
	for (i=0; i<length; i++) 
		x[i]=(y[i]>tau ? y[i]-tau : 0.0); 
	if (x==y) free(aux);
}


/* Algorithm 3 in the paper, using partitioning with respect to a pivot, 
chosen randomly */
static void simplexproj_algo3_pivotrand(const datatype* y, datatype* x,
const unsigned int length, const double a) {
	datatype*	auxlower = (x==y ? (datatype*)malloc(length*sizeof(datatype)) : x);
	datatype*	auxupper = (datatype*)malloc(length*sizeof(datatype));
	datatype*	aux;
	int 	auxlowerlength=0;
	int		auxupperlength=0;
	int		upperlength;
	int		auxlength;
	int		i=0;
	int 	nbpivoteq=0;
	datatype 	pivot=y[(int)(rand() / (((double)RAND_MAX+1.0)/length))];
	double	tauupper;
	double	tau=0.0;
	for (;i<length;i++) 
		if (y[i]<pivot) 
			auxlower[auxlowerlength++]=y[i];
		else if (y[i]==pivot) 
				nbpivoteq++;
		else {
			auxupper[auxupperlength++]=y[i];
			tau+=(y[i]-tau)/auxupperlength;
		}
	tau+=(pivot-tau)*((double)nbpivoteq/(double)(nbpivoteq+auxupperlength))
		-a/(nbpivoteq+auxupperlength);
	if (tau<pivot) {
		upperlength=auxupperlength+nbpivoteq;
		tauupper=tau;	
		auxlength=auxlowerlength;
		aux=auxlower;
	} else {		
		tauupper=0.0;
		upperlength=0;
		aux=auxupper;
		auxlength=auxupperlength;
	}
	while (auxlength>0) {
		pivot=aux[(int)(rand() / (((double)RAND_MAX+1.0)/auxlength))];
		tau=0.0;
		auxlowerlength=0;
		auxupperlength=0;
		nbpivoteq=0;
		for (i=0; i<auxlength; i++) 
			if (aux[i]<pivot) 
				auxlower[auxlowerlength++]=aux[i];
			else if (aux[i]==pivot) 
					nbpivoteq++;
			else {
				auxupper[auxupperlength++]=aux[i];
				tau+=(aux[i]-tau)/auxupperlength;
			}	
		if (upperlength==0) 
			tau+=(pivot-tau)*((double)nbpivoteq/(double)(nbpivoteq+
				auxupperlength))-a/(nbpivoteq+auxupperlength);
		else {
			tau+=(pivot-tau)*((double)nbpivoteq/(double)(nbpivoteq+
				auxupperlength));
			tau+=(tauupper-tau)*((double)upperlength/(double)(nbpivoteq+
				auxupperlength+upperlength));
		}	
		if (tau<pivot) {
			upperlength+=auxupperlength+nbpivoteq;
			tauupper=tau;	
			auxlength=auxlowerlength;
			aux=auxlower;
		} else {
			aux=auxupper;
			auxlength=auxupperlength;
		}
	}
	for (i=0; i<length; i++)
		x[i]=(y[i]>tau ? y[i]-tauupper : 0.0); 
	if (x==y) free(auxlower);
	free(auxupper);
}


/* C code written by L. Condat, translated from Fortran 77 code 
written by K. C. Kiwiel on 8 March 2006.
The algorithm is described in: 
K. C. Kiwiel, "On Floyd and Rivest's SELECT algorithm",
Theoretical Computer Science, vol. 347, pp. 214–238, 2005.
*/
/*
Selects the smallest k elements of the array x[0:r-1].
The input array is permuted so that the smallest k elements of
x are x(i), i = 0,...,k-1, (in arbitrary order) and x(k-1) is the
kth smallest element.
*/ 
void Kiwiel_select(datatype* x, int r, int k) {     
	#define 	nstack 10     
	/*These two arrays permit up to nstack levels of recursion.
	For standard parameters cs <= 1 and cutoff >= 600,
	nstack = 5 suffices for n up to 2**31-1 (maximum integer*4).*/
	int			stack1[nstack];
	int			stack2[nstack];
	#define		cutoff 600
	#define 	cs 0.5
	#define		csd 0.5
	int			i, j, m, s, sd, jstack=0, l=0;
	double		dm, swap, v, z;
	r--;  k--;
	/*entry to SELECT( x, length, l, r, k)
	SELECT will rearrange the values of the array segment x[l:r] so
	that x(k) (for some given k; l <= k <= r) will contain the
	(k-l+1)-th smallest value, l <= i <= k will imply x(i) <= x(k),
	and k <= i <= r will imply x(k) <= x(i).*/
	for (;;) {
		/*The test below prevents stack overflow.*/
		if ((r-l>cutoff) && (jstack<nstack)) {		
		/*Use SELECT recursively on a sample of length s to get an
		estimate for the (k-l+1)-th smallest element into x(k),
		biased slightly so that the (k-l+1)-th element is
		expected to lie in the smaller set after partitioning.*/
			s=cs*exp(2.0*(z=log(dm=m=r-l+1))/3.0)+0.5;
			sd=sqrt(z*s*(1.0-s/dm))*((i=k-l+1)<dm/2.0?-csd:csd)+0.5;
			if (i==m/2) sd=0;
			/*Push the current l and r on the stack.*/
			stack1[jstack]=l;
			stack2[jstack++]=r;
			/*Find new l and r for the next recursion.*/
			swap=k-i*(s/dm)+sd;
			if (swap>l) l=swap+0.5; 
			if (swap+s<r) r=swap+s+0.5; 
		} else {
			if (l>=r) {
				/*Exit if the stack is empty.*/
				if (jstack==0) return;
				/*Pop l and r from the stack.*/
				l=stack1[--jstack];
				r=stack2[jstack];			
			} 
			/*Partition x[l:r] about the pivot v := x(k).*/
			v=x[k];
			/*Initialize pointers for partitioning.*/
			i=l;
			j=r;
			x[k]=x[l];
			x[l]=v;
			if (v<x[r]) {
	            x[l]=x[r];
	            x[r]=v;
			}
			while (i < j) {
	            swap=x[j];
	            x[j--]=x[i];
	            x[i++]=swap;
				/*Scan up to find element >= v.*/
	            while (x[i]<v) i++;
				/*Scan down to find element <= v.*/
				while (x[j]>v) j--;
			}
			if (x[l]==v) {
				x[l]=x[j];
				x[j]=v;
			} else {
				swap=x[++j];
				x[j]=x[r];
				x[r]=swap;
			}
			/*Now adjust l, r so that they surround the subset
			containing the (k-l+1)-th smallest element.*/
			if (j<=k) l=j+1;
			if (k<=j) r=j-1;
		}
	}       
}


/* Algorithm 3 in the paper, using partitioning with respect to a pivot, 
chosen as the median of the list v, found using Kiwiel's implementation
of select. */
static void simplexproj_algo3_pivotmedian(const datatype* y, datatype* x,
const unsigned int length, const double a) {
	datatype*	aux0 = (x==y ? (datatype*)malloc(length*sizeof(datatype)) : x);
	datatype*	aux;
	datatype 	pivot;
	int		auxupperlength;
	int		upperlength;
	int		auxlength;
	int		i=0;
	int 	pospivot=length>>1;	
	datatype*  auxlower=aux0;
	datatype*	auxupper=aux0+pospivot+1;
	double	tauupper;
	double	tau;
	memcpy(aux0,y,length*sizeof(datatype));
	Kiwiel_select(aux0, length, pospivot+1);	
	tau=(pivot=aux0[pospivot])-a;
	auxupperlength=length-pospivot-1;
	for (i=0; i<auxupperlength; i++) 
		tau+=(auxupper[i]-tau)/(i+2);
	if (tau<pivot) {
		upperlength=auxupperlength+1;
		tauupper=tau;	
		auxlength=pospivot;
		aux=auxlower;
	} else {		
		tauupper=0.0;
		upperlength=0;
		aux=auxupper;
		auxlength=auxupperlength;
	}
	while (auxlength>0) {
		pospivot=auxlength>>1;
		Kiwiel_select(aux, auxlength, pospivot+1);	
		auxlower=aux;
		auxupper=aux+pospivot+1;
		if (upperlength==0)
			tau=(pivot=aux[pospivot])-a;
		else 
			tau=tauupper+((pivot=aux[pospivot])-tauupper)/(upperlength+1);
		auxupperlength=auxlength-pospivot-1;
		for (i=0; i<auxupperlength; i++) 
			tau+=(auxupper[i]-tau)/(i+upperlength+2);
		if (tau<pivot) {
			upperlength+=auxupperlength+1;
			tauupper=tau;	
			auxlength=pospivot;
			aux=auxlower;
		} else {
			aux=auxupper;
			auxlength=auxupperlength;
		}
	}
	for (i=0; i<length; i++)
		x[i]=(y[i]>tau ? y[i]-tauupper : 0.0); 
	if (x==y) free(aux0);
}


/* Algorithm using partitioning with respect to a pivot, chosen randomly, 
as given by Duchi et al. in "Efficient Projections onto the l1-Ball for 
Learning in High Dimensions" */
static void simplexproj_Duchi(const datatype* y, datatype* x,
const unsigned int length, const double a) {
	datatype*	auxlower = (x==y ? (datatype*)malloc(length*sizeof(datatype)) : x);
	datatype*	auxupper = (datatype*)malloc(length*sizeof(datatype));
	datatype*	aux=auxlower;
	datatype 	pivot;
	int 	auxlowerlength=0;
	int		auxupperlength=1;
	int		upperlength;
	int		auxlength;
	int		i=0;
	int 	pospivot=(int)(rand() / (((double)RAND_MAX+1.0)/length));	
	double	tauupper;
	double	tau=(pivot=y[pospivot])-a;
	while (i<pospivot) 
		if (y[i]<pivot) 
			auxlower[auxlowerlength++]=y[i++];
		else {
			auxupper[auxupperlength++]=y[i];
			tau+=(y[i++]-tau)/auxupperlength;
		}
	i++;
	while (i<length) 
		if (y[i]<pivot) 
			auxlower[auxlowerlength++]=y[i++];
		else {
			auxupper[auxupperlength++]=y[i];
			tau+=(y[i++]-tau)/auxupperlength;
		}
	if (tau<pivot) {
		upperlength=auxupperlength;
		tauupper=tau;	
		auxlength=auxlowerlength;
	} else {		
		tauupper=0.0;
		upperlength=0;
		aux=auxupper+1;
		auxlength=auxupperlength-1;
	}
	while (auxlength>0) {
		pospivot=(int)(rand() / (((double)RAND_MAX+1.0)/auxlength));
		if (upperlength==0)
			tau=(pivot=aux[pospivot])-a;
		else 
			tau=tauupper+((pivot=aux[pospivot])-tauupper)/(upperlength+1);
		i=0;
		auxlowerlength=0;
		auxupperlength=1;
		while (i<pospivot) 
			if (aux[i]<pivot) 
				auxlower[auxlowerlength++]=aux[i++];
			else {
				auxupper[auxupperlength++]=aux[i];
				tau+=(aux[i++]-tau)/(upperlength+auxupperlength);
			}
		i++;
		while (i<auxlength) 
			if (aux[i]<pivot) 
				auxlower[auxlowerlength++]=aux[i++];
			else {
				auxupper[auxupperlength++]=aux[i];
				tau+=(aux[i++]-tau)/(upperlength+auxupperlength);
			}
		if (tau<pivot) {
			upperlength+=auxupperlength;
			tauupper=tau;	
			auxlength=auxlowerlength;
			aux=auxlower;
		} else {
			aux=auxupper+1;
			auxlength=auxupperlength-1;
		}
	}
	for (i=0; i<length; i++)
		x[i]=(y[i]>tau ? y[i]-tauupper : 0.0); 
	if (x==y) free(auxlower);
	free(auxupper);
}


/* Algorithm 4 in the paper */
static void simplexproj_algo4(const datatype* y, datatype* x,
const unsigned int length, const double a) {
	datatype*	aux = (x==y ? (datatype*)malloc(length*sizeof(datatype)) : x);
	int		auxlength=1; 
	int		auxlengthold=length;
	double	tau=*y-a;
	double tauold;
	int 	i=1;
	for (; i<length; i++) 
		tau+=(y[i]-tau)/(i+1);
	i=0;
	while (y[i]<=tau) i++;
	tauold=tau;
	tau=(*aux=y[i++])-a;
	for (; i<length; i++)
		if (y[i]>tauold) { 
			tau+=((aux[auxlength]=y[i])-tau)/(auxlength+1);
			auxlength++;
		}
	while (auxlength<auxlengthold) {
		auxlengthold=auxlength;
		auxlength=1;	
		i=0;
		while (aux[i]<=tau) i++;
		tauold=tau;
		tau=(*aux=aux[i++])-a;
		for (; i<auxlengthold; i++)
			if (aux[i]>tauold) { 
				tau+=((aux[auxlength]=aux[i])-tau)/(auxlength+1);
				auxlength++;			
			}
	}
	for (i=0; i<length; i++)
		x[i]=(y[i]>tau ? y[i]-tau : 0.0); 
	if (x==y) free(aux);
}


/* Proposed algorithm */
static void simplexproj_Condat(const datatype* y, datatype* x,
const unsigned int length, const double a) {
	datatype*	aux = (x==y ? (datatype*)malloc(length*sizeof(datatype)) : x);
	datatype*  aux0=aux;
	int		auxlength=1; 
	int		auxlengthold=-1;	
	double	tau=(*aux=*y)-a;
	int 	i=1;
	for (; i<length; i++) 
		if (y[i]>tau) {
			if ((tau+=((aux[auxlength]=y[i])-tau)/(auxlength-auxlengthold))
			<=y[i]-a) {
				tau=y[i]-a;
				auxlengthold=auxlength-1;
			}
			auxlength++;
		} 
	if (auxlengthold>=0) {
		auxlength-=++auxlengthold;
		aux+=auxlengthold;
		while (--auxlengthold>=0) 
			if (aux0[auxlengthold]>tau) 
				tau+=((*(--aux)=aux0[auxlengthold])-tau)/(++auxlength);
	}
	do {
		auxlengthold=auxlength-1;
		for (i=auxlength=0; i<=auxlengthold; i++)
			if (aux[i]>tau) 
				aux[auxlength++]=aux[i];	
			else 
				tau+=(tau-aux[i])/(auxlengthold-i+auxlength);
	} while (auxlength<=auxlengthold);
	for (i=0; i<length; i++)
		x[i]=(y[i]>tau ? y[i]-tau : 0.0); 
	if (x==y) free(aux0);
} 


/* Code used in the paper to obtain the computation times. 
Here, this corresponds to Experiment 3 with N=10^6. */
//int main() {
//	datatype*	y;
//	datatype*	x;
//	double	val=0.0;
//	double 	aux;
//	double 	epsi;
//    int 	i,j;
//    int 	Nbrea;
//    unsigned int length;
//    clock_t	start, end;
//    const gsl_rng_type *T;
//  	gsl_rng *r;
//  	gsl_rng_env_setup();
//	gsl_rng_default_seed=32067;
//	T = gsl_rng_default;
//	r = gsl_rng_alloc(T);
//	double timemean=0.0, timevar=0.0, timedelta;
//	const double a = 1.0;
//	srand((unsigned int)pow(time(NULL)%100,3));
//    double*	 timetable;
    
//    length = 1000000;
//	Nbrea = 100;
//	y=(datatype*)malloc(length*sizeof(datatype));
//	x=(datatype*)malloc(length*sizeof(datatype));
//	timetable=(double*)malloc((Nbrea+1)*sizeof(double));
//	for (j=0; j<=Nbrea; j++) {
//		for (i=0; i<length; i++) {
//	    	y[i]=gsl_ran_gaussian(r, 0.001);
//	    }
//	    y[(int)(rand() / (((double)RAND_MAX+1.0)/length))]+=a;
//	    start=clock();
//	    simplexproj_Condat(y,x,length,a);
//	    end=clock();
//	    timetable[j]=(double)(end-start)/CLOCKS_PER_SEC;
//	}
//	/* we discard the first value, because the computation time is higher,
//	probably because of cache operations */
//	for (j=1; j<=Nbrea; j++) {
//		timedelta=timetable[j]-timemean;
//		timemean+=timedelta/j;
//		timevar+=timedelta*(timetable[j]-timemean);
//	}
//	timevar=sqrt(timevar/Nbrea);
//	printf("av. time: %e std dev: %e\n",timemean,timevar);
//	free(x);
//    free(y);
//    return 0;
//}





#endif // CONDAT_H
