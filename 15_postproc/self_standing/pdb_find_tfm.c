/* find transformation which takes one pdb into another - the matching order
 is assumed to be given by the input pdb's*/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include "pdb.h"
# include "utils.h"

# define N  3
# define NO_VECTORS  21
# define BUFSIZE    150

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int layer;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 40

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+1];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    int no_atoms;
    Atom atom[MAX_NO_ATOMS];
    int interface;
} Residue;

# define ALMT_NAME_LENGTH 30

typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int * seq_gaps;
    int * column_gaps;
}  Alignment;

Residue * sequence1, *sequence2;
int no_res_1, no_res_2;
Alignment alignment;


int main ( int argc, char * argv[]) {

    double **x, **y;
    double x_mp[N], y_mp[N];
    char pdbname1 [150] = "\0";
    char pdbname2 [150] = "\0";
    char pdbid1   [50]  = "\0"; /* this might be too restrictive */
    char pdbid2   [50]  = "\0";
    char filename [150] = "\0";
    
    int  component, ctr;
    int  no_vectors;
    int  h, i, j, k;
    double sum = 0;
    Residue * sequence_new;

    double ATA     [4][4] = {{0.0}};
    double prev_ATA[4][4] = {{0.0}};
    double ATA_sum [4][4] = {{0.0}};
    double a[3] = {0.0}, b[3] = {0.0};
    double q[4] = {0.0};
    double R[3][3], T[3];
    
    int add_matrices  (double matrix1[4][4],double matrix2[4][4],
		       double result[4][4]);
    int construct_ATA (double ATA[4][4], double a[3], double  b[3]);
    int quat_to_R (double quat[4], double R[3][3]);
    int read_pdb  (char * pdbname, Residue ** sequence, int *no_res);
    
    /* note how we pass the matrix: pointer to the first element in the block */
    void dsyev_ (char * jobz, char *uplo,  int *n,
		  double *A, int * lda, double * w, double * work, int * lwork, int *info);

    int calphas_to_XY ( double *** xptr,  double ***yptr, char * name_x, char *name_y, int  no_res);
    int transform (double tfm_matrix[][N], double * transl_vector,
		   Residue * seqeunce, int no_res, Residue * seqeunce_new);
    int pdb_output ( char *filename, Residue * sequence_new, int no_res);

    if ( argc < 2 ) {
	fprintf ( stderr, "Usage: %s <pdbname1> <pdbname2>.\n", argv[0] );
	fprintf ( stderr, "(To transform <pdbname1> onto <pdbname2>;");
	fprintf ( stderr, " the alignment btw the two pdb files assumed).\n");
	exit (1);
    } 
    sprintf ( pdbname1, "%s", argv[1]);
    sprintf ( pdbname2, "%s", argv[2]);
   
    memcpy (pdbid1, pdbname1, strlen (pdbname1) - 4); /* get rid of the pdb extension */
    memcpy (pdbid2, pdbname2,  strlen (pdbname2) - 4);
    
    /* input two pdbs */
    if ( read_pdb ( pdbname1, &sequence1, &no_res_1)) exit (1);
    if ( read_pdb ( pdbname2, &sequence2, &no_res_2))  exit (1);

    if ( no_res_1 != no_res_2) {
	fprintf ( stderr, "Error: number of residues in %s and in %s is not the same.\n",
		  pdbid1, pdbid2);
	exit (1);
    }
    /* will need it below */
    sequence_new = emalloc ( no_res_1 * sizeof(Residue));
      
    /* turn the matching atoms into vectors x and y - use only c-alphas*/
    calphas_to_XY ( &x, &y, pdbid1, pdbid2, no_res_1);
    no_vectors = no_res_1;
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
    double rmsd;
    
    if ( !( A=dmatrix(4,4) ) ) exit (1);
    memcpy (A[0], ATA_sum[0], 4*4*sizeof(double));

# if 0
    /* test: */
    /* the e-vals should be -2.0531 -0.5146 -0.2943 12.8621 */
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    if ( i >= j) {
		A[i][j] = i+1;
	    } else {
		A[i][j] = j+1;
		
	    }
	}
    }
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    printf ("%8.3lf ", A[i][j]);
	}
	printf ("\n");
    }
    
# endif
   /* note how we pass the matrix: */
    dsyev_ ( &jobz, &uplo,  &n, A[0], &lda, w, work, &lwork, &info);
    if (  ! info) {
	rmsd = sqrt (w[0]);
	for (i=0; i<4; i++ ) q[i] = A[0][i];
	if (0) {
	    /* w contains the eigenvalues */
	    printf ("\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", w[i]);
	    printf ("\nrmsd: %8.3lf \n", rmsd);
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
   
    for ( i =0; i < 3; i++ ) {
	for ( j =0; j < 3; j++ ) {
	    printf ("%8.3lf ",  R[i][j]);
	}
	printf ("%8.3lf \n", T[i]);
    }
    

# if 0
    /* rotate and translate the first chain */
    memcpy ( sequence_new, sequence1,  no_res_1 * sizeof(Residue));
    transform ( R, T, sequence1, no_res_1, sequence_new );
    /* output the transformed chain */
    sprintf (filename, "%s", "rotated.pdb");
    pdb_output ( filename, sequence_new, no_res_1);

    double aux;
    double xx[3], yy[3];
    int resctr,atomctr;

    rmsd = 0;
    for ( resctr =0; resctr < no_vectors; resctr++ ) {
	
	for (atomctr=0; atomctr < sequence_new[resctr].no_atoms; atomctr++) {
	    xx[0] = sequence_new[resctr].atom[atomctr].x;
	    xx[1] = sequence_new[resctr].atom[atomctr].y;
	    xx[2] = sequence_new[resctr].atom[atomctr].z;
	    
	    yy[0] = sequence2[resctr].atom[atomctr].x;
	    yy[1] = sequence2[resctr].atom[atomctr].y;
	    yy[2] = sequence2[resctr].atom[atomctr].z;

	    for ( i =0; i < 3; i++ ) {
		aux = xx[i] - yy[i];
		rmsd += aux*aux;
	    }
	}
    }
    
    rmsd = sqrt(rmsd/no_vectors);
    printf ("rmsd %8.3lf \n\n",  rmsd);
    
# endif
    

    return 0;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int transform (double  tfm_matrix[][N], double * transl_vector,
	       Residue * sequence, int no_res, Residue * sequence_new) {
    
    int resctr, atomctr;
    int i,j;
    double  x[N], new_x[N];
    
    for ( resctr= 0; resctr < no_res; resctr ++ ) {
	for (atomctr=0; atomctr < sequence[resctr].no_atoms; atomctr++) {
	    x[0] = sequence[resctr].atom[atomctr].x;
	    x[1] = sequence[resctr].atom[atomctr].y;
	    x[2] = sequence[resctr].atom[atomctr].z;
	    for (i=0; i <N; i++ ) {
		new_x[i] = 0.0;
		for (j=0; j <N; j++ ) {
		    new_x[i] += tfm_matrix[i][j]*x[j];
		}
		new_x[i] += transl_vector[i];
	    }
	    sequence_new[resctr].atom[atomctr].x = new_x[0];
	    sequence_new[resctr].atom[atomctr].y = new_x[1];
	    sequence_new[resctr].atom[atomctr].z = new_x[2];
	}
    }
    return 0;
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int pdb_output ( char *filename, Residue * sequence, int no_res){
    
    int resctr, atomctr, serial;
    FILE * fptr;
    Atom * atomptr;
    
    /* open file */
    fptr = fopen ( filename, "w");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", filename);
	return 1;
    }

    
    serial = 0;
    for ( resctr= 0; resctr < no_res; resctr ++ ) {
	for (atomctr=0; atomctr < sequence[resctr].no_atoms; atomctr++) {
	    serial ++;
	    atomptr = sequence[resctr].atom + atomctr;
	    fprintf (fptr,  "%-6s%5d  %-3s%1s%-3s %1s%4s%1s   %8.3lf%8.3lf%8.3lf\n", 
		"ATOM",  serial ,  atomptr->type,   " ",   sequence[resctr].res_type,
		"Z", sequence[resctr].pdb_id,   " " ,   atomptr->x,  atomptr->y,   atomptr->z);
	}
    }

    fclose (fptr);
    
    return 0;
}




#define BUFFLEN 250
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/* calphas_to_XY expects pdb's and alignment on the stack*/
int  calphas_to_XY ( double *** xptr,  double ***yptr, char * name_x, char *name_y, int  no_res){

    int  resctr;
    double **x, **y;
    int find_calpha ( Residue * sequence, int no_res, int resctr, int x_ctr, double ** x);
    

    
     
    /* allocate */
    if ( ! (x = dmatrix (N, no_res)))  exit(1);
    if ( ! (y = dmatrix (N,   no_res)))  exit(1);

	    
    for (resctr=0; resctr<no_res; resctr++) {
	
	if ( find_calpha ( sequence1, no_res, resctr, resctr, x) ) {
	    fprintf ( stderr, "Error constructing x for %s.\n", name_x);
	    exit (1);
	}
	if (find_calpha ( sequence2, no_res, resctr, resctr,  y) ) {
	    fprintf ( stderr, "Error constructing y for %s.\n", name_y);
	    exit (1);
	}
    }

    *xptr = x;
    *yptr = y;
    
    return 0;
}
/******/
int find_calpha ( Residue * sequence, int no_res, int resctr, int x_ctr, double ** x) {
    int atomctr;
    int found_calpha;
    
    found_calpha = 0;
    for (atomctr=0; atomctr < sequence[resctr].no_atoms; atomctr++) {
	if (sequence[resctr].atom[atomctr].backbone &&
	    ! strncmp ( sequence[resctr].atom[atomctr].type, "CA", 2)  ) {
	    found_calpha = 1;
	    x[0][x_ctr] = sequence[resctr].atom[atomctr].x;
	    x[1][x_ctr] = sequence[resctr].atom[atomctr].y;
	    x[2][x_ctr] = sequence[resctr].atom[atomctr].z;
	}
    }
    if ( ! found_calpha) {
	fprintf ( stderr, "C-alpha not found for residue %d.\n", resctr);
	return 1;
    }
   
    return 0;
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int read_pdb ( char * pdbname, Residue ** sequence_ptr, int * no_res_ptr) {

    Residue * sequence;
    FILE * fptr = NULL;
    char line[BUFFLEN];
    char oldresno[PDB_ATOM_RES_NO_LEN+1];
    char oldrestype [PDB_ATOM_RES_NAME_LEN+1];
    char tmp[PDB_ATOM_X_LEN+1], *auxptr;
    int atomctr, resctr, no_res, ctr, nonblank;

    
    /* open file */
    fptr = fopen ( pdbname, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", pdbname);
	return 1;
    }

    /* count residues */
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+1);
    resctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if( ! strncmp(line,"ATOM", 4) ||  ! strncmp(line,"HETATM", 6)){
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN) ) {
		
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN);
		oldresno[PDB_ATOM_RES_NO_LEN] = '\0';
		/* printf ( "New residue number:  %s \n", oldresno); */
		resctr ++;
	    }
	}
    }
    no_res = resctr;
    *no_res_ptr = no_res;
    /* printf ("no residues: %d\n", no_res); */

    /* allocate space */
    sequence = NULL;
    sequence = calloc ( no_res, sizeof (Residue));
    if ( ! sequence ) {
	fprintf ( stderr, "Error allcating sequence space.\n");
	exit (0);
    }
    *sequence_ptr = sequence;

    /* read in the atom */
    rewind ( fptr);
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+1);
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+1);
    resctr= -1;
    atomctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if( ! strncmp(line,"ATOM", 4)  ||  ! strncmp(line,"HETATM", 6)){
	    /* adjust the counters */ 
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN) ) {
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN);
		strncpy (oldrestype, line+PDB_ATOM_RES_NAME, PDB_ATOM_RES_NAME_LEN);
		oldresno[PDB_ATOM_RES_NO_LEN] = '\0';
		oldrestype[PDB_ATOM_RES_NAME_LEN] = '\0';
		resctr ++;
		atomctr = 0;
		
		sequence[resctr].no_atoms = 1;
		strncpy ( sequence[resctr].pdb_id,oldresno, PDB_ATOM_RES_NO_LEN);
		sequence[resctr].pdb_id[PDB_ATOM_RES_NO_LEN] = '\0';
		strncpy ( sequence[resctr].res_type,oldrestype, PDB_ATOM_RES_NAME_LEN);
		sequence[resctr].res_type[PDB_ATOM_RES_NAME_LEN] = '\0';
	   
	    } else {
		atomctr ++;
		sequence[resctr].no_atoms = atomctr + 1;
		if ( atomctr >= MAX_NO_ATOMS ) {
		    fprintf ( stderr, "Error: I thought every aa has < %d atoms.\n",
			      MAX_NO_ATOMS );
		    exit (1);
		}
	    }
	    /* read in atom info */
	    
	    auxptr = line+ PDB_ATOM_ATOM_NAME;
	    memset ( tmp, 0, PDB_ATOM_ATOM_NAME_LEN+1);
	    /* skip initial blanks*/
	    ctr  = 0;
	    while ( !(isalpha (*(auxptr + ctr))) &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) ctr++;
	    /* copy alphanum info */
	    nonblank = 0;
	    while (  isalpha (*(auxptr +ctr))  &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) {
		tmp[nonblank] =  *(auxptr +ctr);
		nonblank ++;
		ctr++;
	    }

	    strncpy ( sequence[resctr].atom[atomctr].type, tmp, PDB_ATOM_ATOM_NAME_LEN );

	    /* is this a backbone atom?*/
	    sequence[resctr].atom[atomctr].backbone = 0;
	    if ( nonblank == 1) {
		  sequence[resctr].atom[atomctr].backbone =
		      !(  strcmp ( tmp, "N") && strcmp ( tmp, "C") && strcmp ( tmp, "O")  );
	    } else if (  nonblank == 2) {
		  sequence[resctr].atom[atomctr].backbone = ! strcmp ( tmp, "CA" );
	    }
	    /* printf ( " %4d %4d %4s is backbone: %1d \n", resctr, atomctr, */
	    /* sequence[resctr].atom[atomctr].type, sequence[resctr].atom[atomctr].backbone); */
	    strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
	    tmp[PDB_ATOM_X_LEN] = '\0';
	    sequence[resctr].atom[atomctr].x=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
	    tmp[PDB_ATOM_Y_LEN] = '\0';
	    sequence[resctr].atom[atomctr].y=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
	    tmp[PDB_ATOM_Z_LEN] = '\0';
	    sequence[resctr].atom[atomctr].z=atof(tmp);
	   
	}
	
    }

    /* close file */
    fclose (fptr);
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

    
/**********************************************************/

int quat_to_R ( double quat[4], double R[3][3] ){

    
    double q0 = quat[0];
    double q1 = quat[1];
    double q2 = quat[2];
    double q3 = quat[3];


    R[0][0] = 1 -2*q2*q2 -2*q3*q3;
    R[1][1] = 1 -2*q3*q3 -2*q1*q1 ;
    R[2][2] = 1 -2*q1*q1 -2*q2*q2 ;

    R[0][1] = 2*q1*q2 - 2*q0*q3;
    R[1][0] = 2*q1*q2 + 2*q0*q3;
    
    R[0][2] = 2*q1*q3 + 2*q0*q2;
    R[2][0] = 2*q1*q3 - 2*q0*q2;
    
    R[1][2] = 2*q2*q3 - 2*q0*q1;
    R[2][1] = 2*q2*q3 + 2*q0*q1;

    return 0;

}
