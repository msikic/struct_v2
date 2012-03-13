/* find transformation which takes x to y*/
# include "struct.h"

int  map2rotation (Protein *protein1, Protein *protein2, int *map_i2j,
		   double **x, double **y, double *q, double *T, double *rmsd){
    
    int mappedCoords2rotation (double **x, double **y, int map_size,
			       double q[4], double T[3], double *rmsd);

    int resctr1, resctr2, map_count = 0;
    int no_res_1 = protein1->length, no_res_2= protein2->length;
    int i, map_size;
    double ca1[3], ca2[3];
     
    
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	map_count ++;
    }
    map_size = map_count;

    /*  mapping is stored in terms of residues ids -
	turn this into coordinates of the corresponding Calphas */
    
    map_count = 0;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	
	find_Calpha ( protein1, resctr1, ca1 );
	find_Calpha ( protein2, resctr2, ca2 );

	if (map_count >= no_res_1+no_res_2 ) {
	    fprintf (stderr, "bleep\n"); exit (1);
	}
	
	for (i=0; i<3; i++) {
	    x[i][map_count] = ca1[i];
	    y[i][map_count] = ca2[i];
	}
	map_count ++;
    }

    mappedCoords2rotation (x, y, map_size, q, T, rmsd);
    
    return 0;
}

/*********************************************************/
/*********************************************************/

int mappedCoords2rotation (double **x, double **y, int map_size,
		   double q[4], double T[3], double *rmsd) {

    double x_mp[3], y_mp[3];
    int  ctr;
    int  no_vectors;
    int  i, j;
 
    double ATA     [4][4] = {{0.0}};
    double prev_ATA[4][4] = {{0.0}};
    double ATA_sum [4][4] = {{0.0}};
    double a[3] = {0.0}, b[3] = {0.0};
    double **R;
    
    int add_matrices  (double matrix1[4][4],double matrix2[4][4],
		       double result[4][4]);
    int construct_ATA (double ATA[4][4], double a[3], double  b[3]);

    /* note how we pass the matrix: pointer to the first element in the block */
    void dsyev_ (char * jobz, char *uplo,  int *n,
		  double *A, int * lda, double * w, double * work, int * lwork, int *info);

    if (!map_size) {
	*rmsd = -1;
	return 1;
    }

    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    memset ( &(q[0]), 0, 4*sizeof(double) );
    /* turn the matching atoms into vectors x and y - use only c-alphas*/
    no_vectors = map_size;
    
    /* check: */
    if (0) {
	printf (" Number of vectors read in: %d. \n", no_vectors);
	for ( ctr =0; ctr < no_vectors; ctr++ ) {
	    printf ("\t x%1d   %10.4lf  %10.4lf  %10.4lf   ",
		    ctr, x[0][ctr], x[1][ctr], x[2][ctr]);
	    printf ("\t y%1d   %10.4lf  %10.4lf  %10.4lf \n",
		    ctr, y[0][ctr], y[1][ctr], y[2][ctr]);
	}
	exit (1);
    }

    /* find the meanpoints: */
    for ( i =0; i < 3; i++ ) {
	x_mp[i] = 0.0;
	y_mp[i] = 0.0;
    }
    for ( ctr =0; ctr < no_vectors; ctr++ ) {
	for ( i =0; i < 3; i++ ) {
	    x_mp[i] += x[i][ctr];
	    y_mp[i] += y[i][ctr];
	}
    }
    for ( i =0; i < 3; i++ ) {
	x_mp[i] /= no_vectors;
	y_mp[i] /= no_vectors;
    }
    /* subtract them from x, y */
    for ( ctr =0; ctr < no_vectors; ctr++ ) {
	for ( i =0; i < 3; i++ ) {
	    x[i][ctr] -= x_mp[i];
	    y[i][ctr] -= y_mp[i];
	}
    }
    /* B = ATA_sum matrix to diagonalize in order to get the quaternion */
    for ( ctr =0; ctr < no_vectors; ctr++ ) {
   	for (i=0; i<3; i++ ) {
	    a[i] = y[i][ctr] + x[i][ctr];
	    b[i] = y[i][ctr] - x[i][ctr];
	}
 	construct_ATA (ATA, a, b);
	add_matrices (prev_ATA, ATA, ATA_sum);
	memcpy (prev_ATA[0], ATA_sum[0], 4*4*sizeof(double));
    }
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA_sum[i][j] /= no_vectors;
	}
    }
    /* diagonalize ATA_sum - the eigenvector corresponsing to the
       smallest lambda is the quaternion we are looking for; the
       eigenvalue is the rmsd*/
    /* use the nomenclature from dsyev*/
    char jobz= 'V'; /*Compute eigenvalues and eigenvectors.*/
    char uplo= 'U'; /* Upper triangle of A (the matrix we are diagonalizing) is stored; */
    int  n = 4;     /* order and the leading dimension of A */
    int  lda = 4;
    double ** A;
    int  info;
    int  lwork = 200;
    double w [4];
    double work[200];
    
    if ( !( A=dmatrix(4,4) ) ) exit (1);
    memcpy (A[0], ATA_sum[0], 4*4*sizeof(double));


   /* note how we pass the matrix: */
    dsyev_ ( &jobz, &uplo,  &n, A[0], &lda, w, work, &lwork, &info);
    if (  ! info) {
	*rmsd = sqrt (w[0]);
	for (i=0; i<4; i++ ) q[i] = A[0][i];
	if (0) {
	    /* w contains the eigenvalues */
	    printf ("\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", w[i]);
	    printf ("\nrmsd: %8.3lf \n", *rmsd);
	    printf ("quat:\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", q[i]);
	    printf ("\n");
	    /* printf (" opt lwork: %d\n", (int) work[0]); */
	}
    } else {
	fprintf (stderr, "Error in dsyev().\n");
	exit (1);
    }
    
    /* construct the rotation matrix R */
    quat_to_R (q,R);
    /* T = y_mp - R x_mp */
    for (i=0; i<3; i++ ) {
	T[i] = y_mp[i];
	for (j=0; j<3; j++ ) {
	    T[i] -= R[i][j]*x_mp[j];
	}
    }
   
    
    free_dmatrix(A);
    free_dmatrix(R);

 
    return 0;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double braket  (double ATA[4][4],double  q[4]){
    int i,j;
    double value = 0.0;

    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    value += q[i]*ATA[i][j]*q[j];
	}
    }
   
    return value;
}
 
/*******************************************************************************/
int construct_ATA (double ATA[4][4], double a[3], double  b[3]){

    int i,j,k;
    double A[4][4] = {{ 0.0, -b[0], -b[1], -b[2]},
		      {b[0],   0.0, -a[2],  a[1]},
		      {b[1],  a[2],   0.0, -a[0]},
		      {b[2], -a[1],  a[0],   0.0}};
    
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA[i][j] = 0.0;
	    for (k=0; k<4; k++ ) {
		ATA[i][j] += A[k][i]*A[k][j];
	    }
	}
    }
     
    return 0;
}

/*******************************************************************************/
int add_matrices  (double matrix1[4][4],double matrix2[4][4],
		   double result[4][4]){
    int i,j;

    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    result[i][j] = matrix1[i][j] + matrix2[i][j];
	}
    }
    return 0;

}

    
