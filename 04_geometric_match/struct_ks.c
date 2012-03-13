# include "struct.h"
/***************************************/
 void nrerror(char error_text[]){
     fprintf ( stderr, "%s\n", error_text);
     exit(1);
}
/***************************************/

void kstwo(double *data1, int n1, double * data2, int n2,
	   double *d, double *prob)
/* Given an array data1[1..n1], and an array data2[1..n2], this routine returns the K­S */
/* statistic d, and the significance level prob for the null hypothesis that the data sets are */
/* drawn from the same distribution. Small values of prob show that the cumulative distribution */
/* function of data1 is significantly different from that of data2. The arrays data1 and data2 */
/* are modified by being sorted into ascending order. */
{
    double probks(double alam);
    void sort(int n, double arr[]);
    int j1=1,j2=1;
    double d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;
    sort(n1,data1);
    sort(n2,data2);
    en1=n1;
    en2=n2;
    *d=0.0;
    while (j1 <= n1 && j2 <= n2) {    /* If we are not done... */
	if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=j1++/en1;   /*  Next step is in data1. */
	if (d2 <= d1) fn2=j2++/en2;   /*  Next step is in data2. */
	if ((dt=fabs(fn2-fn1)) > *d) *d=dt;
    }
    en=sqrt(en1*en2/(en1+en2));
    *prob=probks((en+0.12+0.11/en)*(*d));    /* Compute significance. */
}

#define EPS1 0.001
#define EPS2 1.0e-8
double probks(double alam)
/* Kolmogorov-Smirnov probability function. */
{
    int j;
    double a2,fac=2.0,sum=0.0,term,termbf=0.0;
    a2 = -2.0*alam*alam;
    for (j=1;j<=100;j++) {
	term=fac*exp(a2*j*j);
	sum += term;
	if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
	fac = -fac;   /*  Alternating signs in sum. */
	termbf=fabs(term);
    }
    return 1.0;    /* Get here only by failing to converge. */
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50
/* Here M is the size of subarrays sorted by straight insertion and NSTACK is the required auxiliary */
/* storage. */
void sort(int n, double arr[])
/* Sorts an array arr[1..n] into ascending numerical order using the Quicksort algorithm. n is */
/* input; arr is replaced on output by its sorted rearrangement. */
{
    int i,ir=n,j,k,l=1,*istack;
    int jstack=0;
    double  a,temp;
    istack= (int*) calloc (NSTACK, sizeof(int) );
    for (;;) {  /*   Insertion sort when subarray small enough. */
	if (ir-l < M) {
	    for (j=l+1;j<=ir;j++) {
		a=arr[j];
		for (i=j-1;i>=l;i--) {
		    if (arr[i] <= a) break;
		    arr[i+1]=arr[i];
		}
		arr[i+1]=a;
	    }
	    if (jstack == 0) break;
	    ir=istack[jstack--];   /*  Pop stack and begin a new round of parti- */
	    l=istack[jstack--];   /*  tioning. */
	} else {
	    k=(l+ir) >> 1;   /*  Choose median of left, center, and right el- */
	    SWAP(arr[k],arr[l+1])   /*  ements as partitioning element a. Also */
		if (arr[l] > arr[ir]) {   /*  rearrange so that a[l]  a[l+1]  a[ir]. */
		    SWAP(arr[l],arr[ir])
			}
	    if (arr[l+1] > arr[ir]) {
		SWAP(arr[l+1],arr[ir])
		    }
	    if (arr[l] > arr[l+1]) {
		SWAP(arr[l],arr[l+1])
		    }
	    i=l+1;   /*  Initialize pointers for partitioning. */
	    j=ir;
	    a=arr[l+1];    /* Partitioning element. */
	    for (;;) {  /*   Beginning of innermost loop. */
		do i++; while (arr[i] < a);    /* Scan up to find element > a. */
		do j--; while (arr[j] > a);   /*  Scan down to find element < a. */
		if (j < i) break;    /* Pointers crossed. Partitioning complete. */
		SWAP(arr[i],arr[j]);  /*   Exchange elements. */
	    }   /*  End of innermost loop. */
	    arr[l+1]=arr[j];    /* Insert partitioning element. */
	    arr[j]=a;
	    jstack += 2;
/* Push pointers to larger subarray on stack, process smaller subarray immediately. */
	    if (jstack > NSTACK) nrerror ("NSTACK too small in sort.");
	
	    if (ir-i+1 >= j-l) {
		istack[jstack]=ir;
		istack[jstack-1]=i;
		ir=j-1;
	    } else {
		istack[jstack]=j-1;
		istack[jstack-1]=l;
		l=i;
	    }
	}
    }
    free(istack);
}

