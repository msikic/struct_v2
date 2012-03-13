# include "struct.h"


int  geometric_descriptors (Protein * protein) {
    
    int helix_ctr, strand_ctr, res_ctr, atom_ctr, ctr, no_of_points;
    double **point;
    double  p[3]; /* the direction vector */
    double cm[3]; /* the CM vector */
    int sheet[100][100] = {{0}}; /* TODO I can have up to a 100 sheets, 100 starnds each */
    char sheet_id [100][PDB_SHEET_SHEET_LEN+1] = {{'\0'}}; /* TODO id for each sheet */
    int sheet_ctr,  no_of_sheets = 0;
    int n_params = 7;
    double *params;
    Residue * residue;
    int lm_fit (double **points, int n_points, double *params, int n_params);
    int line_fit (double **point, int no_of_points,
	      double p[], double center_of_mass[]);
    int perp_line_fit (double **point, int no_of_points,
		       double p[], double cm[], double avg[]);
    double normal_to_catenoid (double **point, int no_of_points,
			       double *params, double normal[3], double foot[3]);
    
    if ( ! (params = emalloc (n_params*sizeof(double) ) ) ) return 1;
    /* allocate 3 times the protein length (= numbr of residues),
       bcs for the fiting purposes
       we might extract C,N and perhaps O backbone atoms */
    if ( ! (point = dmatrix (3*protein->length, 3) ) ) return 1;
     
    /************************************************************************/
    /************************************************************************/
    /************************************************************************/
    /*****               HELIX   PROCESSING  ********************************/
    /************************************************************************/
    /************************************************************************/
    /*for each helix, extact C-alpha in an array, and fit a line through it */
     
    for (helix_ctr=0; helix_ctr<protein->no_helices; helix_ctr++) {
	ctr = 0;

	/* extract C-alphas */
	for (res_ctr = 	protein->helix[helix_ctr].begin;
	     res_ctr <= protein->helix[helix_ctr].end; res_ctr++) {
	    residue = protein->sequence+res_ctr;
	    for (atom_ctr = 0; atom_ctr < residue->no_atoms; atom_ctr++) {
		if ( ! strcmp (residue->atom[atom_ctr].type, "CA") ) {
		    point[ctr][0] = residue->atom[atom_ctr].x;
		    point[ctr][1] = residue->atom[atom_ctr].y;
		    point[ctr][2] = residue->atom[atom_ctr].z;
		    ctr++;
		    break;
		}
	    }
	}
	no_of_points = ctr;
	
	/* fit a line through them */
	line_fit (point, no_of_points, p, cm);
	if ( ! (protein->helix[helix_ctr].p  = dmatrix (1,3)) ) return 1;
	if ( ! (protein->helix[helix_ctr].cm = dmatrix (1,3)) ) return 1;
	memcpy (protein->helix[helix_ctr].p[0],   p, 3*sizeof(double) );
	memcpy (protein->helix[helix_ctr].cm[0], cm, 3*sizeof(double) );
	protein->helix[helix_ctr].no_of_vectors = 1;

    }
    
    /************************************************************************/
    /************************************************************************/
    /************************************************************************/
    /*****              SHEET    PROCESSING  ********************************/
    /************************************************************************/
    /************************************************************************/
    /*for the moment, do the same for each strand */
    /*for each strand, extact C-alpha in an array, and fit a line through it */
    for (strand_ctr=0; strand_ctr<protein->no_strands; strand_ctr++) {

	/* which sheet does this strand belong to? */
	for ( sheet_ctr=0; sheet_ctr < no_of_sheets; sheet_ctr++ ) {
	    if ( ! strcmp (sheet_id[sheet_ctr],
			   protein->strand[strand_ctr].sheet_id) ) {
		break;
	    }
	}
	
	/* some sheet we have not seen so far: */
	if ( sheet_ctr == no_of_sheets ) {
	    memcpy (sheet_id[sheet_ctr],  protein->strand[strand_ctr].sheet_id,
		    PDB_SHEET_SHEET_LEN); 
	    no_of_sheets ++; 
	}
	sheet[sheet_ctr][0] ++; 
	sheet[sheet_ctr][ sheet[sheet_ctr][0] ] = strand_ctr; 
	protein->strand[strand_ctr].sheet_number = sheet_ctr+1;

	/**** associate a direction with each strand ****/
	/* extract C's and N's */
	ctr = 0;
	for (res_ctr = 	protein->strand[strand_ctr].begin;
	     res_ctr <= protein->strand[strand_ctr].end; res_ctr++) {
	    residue = protein->sequence+res_ctr;
	    for (atom_ctr = 0; atom_ctr < residue->no_atoms; atom_ctr++) {
		if ( ! strcmp (residue->atom[atom_ctr].type, "C")      
		     ||  ! strcmp (residue->atom[atom_ctr].type, "N") ) {
		    point[ctr][0] = residue->atom[atom_ctr].x;
		    point[ctr][1] = residue->atom[atom_ctr].y;
		    point[ctr][2] = residue->atom[atom_ctr].z;
		    ctr++;
		}
	    }
	}
	no_of_points = ctr;
	
	/* allocate */
	if ( ! (protein->strand[strand_ctr].p  = dmatrix (1,3)) ) return 1;
	if ( ! (protein->strand[strand_ctr].cm = dmatrix (1,3)) ) return 1;
	/***********/
	
	/* fit a line through C's and N's ,--- cannot use CA trace*/
	line_fit (point, no_of_points, p, cm);
	memcpy (protein->strand[strand_ctr].p[0],   p, 3*sizeof(double) );
	memcpy (protein->strand[strand_ctr].cm[0], cm, 3*sizeof(double) );
	protein->strand[strand_ctr].no_of_vectors = 1;
	/***********/


    }
    protein->no_sheets  = no_of_sheets;
    if ( ! no_of_sheets) return 0;

    
    free (params);
    free_dmatrix (point);
    
    return 0;
}

/************************************************/
/************************************************/
/************************************************/
/************************************************/
int line_fit (double **point, int no_of_points,
	      double p[], double center_of_mass[]) {
    
    double I[3][3] = {{0.0}}; /* the "moments of inertia" */
    double my_point[no_of_points][3];
    int i, x, y;
    double normalize (double *p);

    /***************************/
    /* find the center of mass */
    /***************************/
    for ( x=0; x<3; x++)  center_of_mass[x] = 0.0;
    for (i=0; i< no_of_points; i++ ) {
	for ( x=0; x<3; x++) {
	    center_of_mass[x] += point[i][x];
	}
    }
    for ( x=0; x<3; x++) {
	center_of_mass[x] /= no_of_points;
    }
    
    /***********************************/
    /* move the points to the cm frame */
    /***********************************/
    for (i=0; i< no_of_points; i++ ) {
	for ( x=0; x<3; x++) {
	    my_point[i][x]  = point[i][x] - center_of_mass[x];
	}
    }

    /**********************************/
    /* find the "moments of inertia"  */
    /**********************************/
    for (i=0; i< no_of_points; i++ ) {
	for ( x=0; x<3; x++) {  /* modulo = circular permutation */
	    I[x][x] += my_point[i][(x+1)%3]*my_point[i][(x+1)%3] +
		my_point[i][(x+2)%3]*my_point[i][(x+2)%3];
	    for ( y=x+1; y<3; y++) { /* off diag elements */
		I[x][y] -= my_point[i][x]*my_point[i][y];
	    }
	}
    }
    for ( x=0; x<3; x++) { 
	for ( y=x+1; y<3; y++) {
	    I[y][x] =  I[x][y];
	}
    }

    /*****************************************/
    /* diagonalize I[][], pick the direction
       with the smallest moment of inertia,
       and rotate back to the initial frame */
    /*****************************************/
    void dsyev_ ( char * jobz, char * uplo, int* N, double * A, int * leading_dim,
		  double * eigenvalues, double *workspace, int *workspace_size, int * retval);
    char jobz = 'V'; /* find evalues and evectors */
    char uplo = 'L'; /* amtrix is stored as lower (fortran convention) */
    int  N = 3; /* the order of matrix */
    int leading_dim = N;
    int retval;
    double A[N*N];
    double eigenvalues[N];
    double workspace[3*N];
    int workspace_size = 3*N;

    for ( x=0; x<3; x++) {
	for ( y=0; y<3; y++) {
	    A[x*3+y] = I[x][y];
	}
    }
   
    dsyev_ ( &jobz, &uplo, &N, A,  &leading_dim, eigenvalues,
	     workspace, &workspace_size, &retval);

    if ( retval ) {
	fprintf ( stderr, "Dsyev  error: %d.\n", retval);
	exit (1);
    }

    /* the eigenvalues are returned in ascending order, so the first guy is mine: */
    x = 0;
    for ( y=0; y<3; y++) {
	p[y] = A[x*3+y];/*this is  p, the direction vector   */
    }

    /* is it pointing toward C-terminal of my helix? */
    /* scalar product between the vector from the
       first to the last point in the helix and p -
       if negative, change the sign of p */
    {
	double *pt_last  = my_point[no_of_points-1];
	double *pt_first = my_point[0];
	double scp = 0.0;
	for ( y=0; y<3; y++) {
	    scp += ( pt_last[y]-pt_first[y] )*p[y];
	}
	if ( scp < 0 )  for ( y=0; y<3; y++) p[y] = - p[y];
	
    }
    return 0;
    
}

/************************************************/
/************************************************/
/************************************************/
/************************************************/
int perp_line_fit (double **point, int no_of_points,
		   double p[], double cm[], double avg[]) {

    /* p and the center of the mass define a line */
    /* we need a new set of lines perpencidular through this one,
       going through the points given in the array */
    int i, j;
    double dp, d[3], p_perp[3], p_perp_0[3], cosine, norm;
    double normalize (double *p);

    memset (avg, 0, 3*sizeof(double));
    
    /* if d is the vector from CM to the atom O,
       p_perp = d -(dp)p (all vectors) */
    for (i=0; i< no_of_points; i++) {
	dp = 0;
	for (j=0; j<3; j++) {
	    d[j] = point[i][j] - cm[j];
	    dp += d[j]*p[j];
	}
	norm = 0;
	for (j=0; j<3; j++) {
	    p_perp[j] = d[j] - dp*p[j];
	    norm +=  p_perp[j]*p_perp[j];
	}
	norm = sqrt (norm);
	for (j=0; j<3; j++) p_perp[j] /= norm;

	if ( i== 0 ) {
	    for (j=0; j<3; j++) p_perp_0[j] = p_perp[j];
	}
	unnorm_dot (p_perp_0, p_perp, &cosine);
	if ( cosine < 0 ) {
	    for (j=0; j<3; j++)  p_perp[j] = -p_perp[j] ;
	    cosine = -cosine;
	}
	for (j=0; j<3; j++) avg[j] +=  p_perp[j];
# if 0
	/* TODO - this might be the place to start
	   cleaning up the line fit - all vectors should
	   be parallel */ 
	for (j=0; j<3; j++)  printf (" %8.2lf", p_perp[j]);
	printf ("   %8.2lf\n", cosine);
# endif
    }

    norm = 0;
    for (j=0; j<3; j++) {
	avg[j] /= no_of_points;
	norm +=  avg[j]*avg[j];
    }
    norm = sqrt (norm);
    for (j=0; j<3; j++) avg[j] /= norm;
    
     
    return 0;      
}





