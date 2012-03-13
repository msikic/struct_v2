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
    double R[N][N], T[N];
    double Q[N][N], tau[N];
    double rQ[N][N], rQnew[N][N], H[N][N], v[N];
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

    int  read_pdb ( char * pdbname, Residue ** sequence, int *no_res);
    void dgels_ (char * trans, int * no_rows, int * no_columns,  int * ,
		 double ** scratch, int *,
		  double **A, int *, double * work, int * lwork, int *info);
    void dgeqrf_ (int *M, int *, double **A, int *LDA, double * TAU,
		  double * WORK, int * LWORK, int *INFO );
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
    
    /* make the fourth component of x equal to 1 - a trick to incorporate translation into A */
    for (ctr =0; ctr < no_vectors; ctr++)  x[3][ctr] = 1.0;
  
    
    /* solve the least squares problem  - use the nomenclature from dgels */
    char trans= 'N';
    int  info;
    int  lwork = 2*no_vectors;
    double work [2*no_vectors];
    int  n = N;
   
    double A[N+1][no_vectors];
    int no_rows = no_vectors, no_columns = N+1;
    int lead_dim_A = no_vectors;
    double B[N][no_vectors];
    int lead_dim_B = no_vectors;
    
    memcpy (A[0], x[0], (N+1)*no_vectors*sizeof(double));
    memcpy (B[0], y[0], N*no_vectors*sizeof(double));

    dgels_ ( &trans, &no_rows, &no_columns,  &n,  &A, &lead_dim_A,
	     &B, &lead_dim_B, work, &lwork, &info);


# if 0 
    printf ("******************************************************\n");
    printf (" solution: \n" );
    for ( ctr =0; ctr < N; ctr++ ) {
	for ( component=0; component<N+1; component++) {
	    printf ("%10.3lf", B[ctr][component]);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
# endif

    /* rotation and translation parts*/
    for ( i =0; i < N; i++ ) {
	for ( j =0; j < N; j++ ) {
	    R[i][j] = B[i][j];
	}
	T[i] = B[i][N];
    }
    /* rotate and translate the first chain */
    sequence_new = emalloc ( no_res_1 * sizeof(Residue));
    memcpy ( sequence_new, sequence1,  no_res_1 * sizeof(Residue));
    transform ( R, T, sequence1, no_res_1, sequence_new );
    /* output the transformed chain */
    sprintf (filename, "%s", "transformed.pdb");
    pdb_output ( filename, sequence_new, no_res_1);
    
# if 0   
    /* is the solution orthogonal? */
    printf (" orthogonal?\n" );
    for ( i =0; i < N; i++ ) {
	for ( j =0; j < N; j++ ) {
	    sum = 0.0;
	    for ( component=0; component<N; component++) {
		sum +=  R[i][component]*R[j][component];
	    }
	    printf ("%10.3lf", sum);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
# endif
    
    /* no reason to believe it will be orthogonal, so orthogonalize using QR decomp: */
    /* find decomposition: */
    for ( i =0; i < N; i++ ) {
	for ( j =0; j < N; j++ ) {
	    Q[i][j] = B[j][i];
	}
    }
    n = N;
    dgeqrf_ ( &n, &n, &Q, &n, tau, work, &lwork, &info);
    if ( info ) {
	fprintf ( stderr, "Error running dgeqrf. Info: %d.\n", info);
	exit (1);
    }
    
    /* reconstruct Q: */
    /*extract R*/
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < i; j++ ) {
	    R[i][j] = 0.0;
	}
	for ( j =i; j < n; j++ ) {
	    R[i][j] = Q[j][i];
	}
    }
    /* reconstruct Q (I could not get the orginal LAPACK function to work: */
    
    memset( rQ[0], 0, n*n*sizeof(double));
    rQ[0][0] = rQ[1][1] = rQ[2][2] = 1.0;
    for ( h =0; h < n; h++ ) {
	/* find vh*/
	for ( i=0; i<h; i++ ) v[i] = 0.0;
	v[h] = 1.0;
	for ( i=h+1; i<n; i++ ) v[i] = Q[h][i];
	
	/* find Hh */
	for ( i =0; i < n; i++ ) {
	    H[i][i] = 1.0 -tau[h]*v[i]*v[i];
	    for ( j =i+1; j < n; j++ ) {
		H[i][j] = H[j][i] = -tau[h]*v[i]*v[j];
	    }
	}
	
	/* multiply rQ by Hi */
	for ( i =0; i < n; i++ ) {
	    for ( j =0; j < n; j++ ) {
		rQnew[i][j] = 0.0;
		for ( k =0; k < n; k++ ) {
		    rQnew[i][j] += rQ[i][k]*H[k][j];
		}
		
	    }
	}
	memcpy ( rQ[0], rQnew[0], n*n*sizeof(double));
    }

    /* to get as close as possible to the original matrix, require that diagonals
       in R be positive (in the limiting case when the input matrix is already
       orthogonal, R should be I */
    
    for ( i =0; i < n; i++ ) {
	if ( R[i][i] < 0 ) {
	    for ( j =0; j < n; j++ ) {
		rQ[j][i] *= -1;
		R [i][j] *= -1;
	    }
	} 
    }
  

# if 0   
    printf ("Q reconstructed    \n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    printf ("%10.3lf",  rQ[i][j]);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    printf ("R:   \n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    printf ("%10.3lf",  R[i][j]);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    printf ("final orthogonality\n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    sum = 0.0;
	    for ( component=0; component<n; component++) {
		sum +=  rQ[component][i]*rQ[component][j];
	    }
	    printf ("%10.3lf", sum);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");
    printf ("QRproduct\n");
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    sum = 0.0;
	    for ( component=0; component<n; component++) {
		sum +=  rQ[i][component]*R[component][j];
	    }
	    printf ("%10.3lf", sum);
	}
	printf ("\n");
    }
    printf ("\n");
    printf ("******************************************************\n");


# endif
    /* rotate and translate the first chain */
    // memcpy ( sequence_new, sequence1,  no_res_1 * sizeof(Residue));
    transform ( rQ, T, sequence1, no_res_1, sequence_new );
    /* output the transformed chain */
    sprintf (filename, "%s", "rotated.pdb");
    pdb_output ( filename, sequence_new, no_res_1);

    double aux, rmsd = 0;
    double xx[3], yy[3];
    int resctr,atomctr;
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
    
    for ( i =0; i < n; i++ ) {
	for ( j =0; j < n; j++ ) {
	    printf ("%8.3lf ",  rQ[i][j]);
	}
	printf ("%8.3lf \n", T[i]);
    }


    

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
    if ( ! (x = dmatrix (N+1, no_res)))  exit(1);
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
