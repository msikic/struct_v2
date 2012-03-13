# ifndef _BC_H
# define _BC_H
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include <time.h>
# include <assert.h>
# include "struct_utils.h"
# include "struct_pdb.h"
# include "struct_geometry.h"
# include "hungarian.h"
# include <omp.h>
# ifdef DMALLOC
#   include "dmalloc.h"
# endif


# define  BUFFLEN 150

/* used in preprocessor */
# define MIN_HELIX_LENGTH 8
# define MIN_STRAND_LENGTH 3
/* pdb parser errors */
# define ERR_SSE_NONE    2
# define ERR_CA_TRACE    4
# define ERR_NONSENSE    8
# define ERR_CONTAINER  16
# define ERR_OVERLAP    32
# define ERR_MAX_ATOMS  64
# define ERR_NO_FILE_OR_CHAIN  128

/* used in alignment functions */
# define FAR_FAR_AWAY -1000

/* for matching elements of different type */
#define REWARD   1
#define PENALTY  0

typedef struct {
    double sheet_cosine;  /* min cos angle between two strands
			     in a sheet to be represented by the
			     same vector */ 
    double alpha;         /* gaussian width for the scoring fn F */
                          /* note that it is set from the input table */
                          /* rather than from the cmd file        */
    double tol;           /* ? */
    double z_max_store;   /* z-score cutoff to store a map */
    double z_max_out;     /* z-score cutoff for map output  */
    double F_guess_max;   /* max  value of F for an intial guess    */
    double F_eff_max;     /* max  value of F after gradient descent */
    double z_max_corr;    /* max value of z for which two maps are */
			  /* considered correlated */
    double z_min_compl;   /* min value of z for which two maps are */
			  /* considered complementary */
    double grad_step_size;/* step size in gradient descent */
    double grad_stop_tol; /* stopping precision in gradient descent */
    double far_far_away;  /* nonsense distance for the pairwise alignment*/
    double gap_open;      /* gap penalties for the pairwise alignment */
    double gap_extend;
    double endgap;
    
    double far_away_cosine; /* minimum cosine for F_effective estimate*/

    double H_length_mismatch_tol; /* tolerance for the difference in length
				     in matched helices */
    double S_length_mismatch_tol; /* tolerance for the difference in length
				     in matched strands */
    
    int min_no_SSEs;      /* min number of SSEs we are willing to consider as a hit */
    int smith_waterman;   /* use Smith-Waterman rather than Needleman-Wunsch */
    int use_endgap;
    int grid_size;        /* minimal number of points for the sphere grid */
    int number_maps_cpl;  /* number of top scoring maps among which to   */
                          /* look for a complement */ 
    int number_maps_out;  /* number of top scoring maps to output     */
    int grad_max_step;    /* max number of steps in  gradient descent */
    int exp_table_size;   /* size of the lookup table for the exp funcion */
    int verbose;          /* toggles verbose output */
    int use_perp;         /* use direction perp to beta sheet? */
    int use_length;       /* use length of SSEs in the structure matching */
    int print_header;     /* print header in the short output file  */
    int report_no_sse_overlap; /*produce short output line even when
				 there is no overlap in the SSE type */
    int  postprocess;     /* produce an output for postprocessing */ 
    char pdbf_tgt[BUFFLEN];  /* for postprocessing, we'll need the full set of
			     coordinates - the path to PDB file */
    char pdbf_qry[BUFFLEN];
    char chain_tgt;
    char chain_qry;
  
    char outname[BUFFLEN];/* name for the output file   */
    char path[BUFFLEN];   /* path to the integral table */
    
    int sw_score; // use smith-waterman score in database search ?
    int score_out; // 0 - sw_score, 1 - hungarian, 2 - both
    
} Options;

typedef struct {
    double theta;
    double phi;
    void * subgrid;
} Point;



typedef Strand Element; /* pseudonym */

typedef struct {
    char name[SHORTSTRING];
    int no_of_residues, no_of_elements; 
    int no_of_helices, no_of_p_sheets, no_of_a_sheets; 
    int no_of_strands;
    Element  * element;  /* strand has all the fields that I need,
			   but whether it is actually a strand or helix
			   is stored in type array:
			   1 for helix, -1 for strand */
    int * type;
    int * sheet_id; /* sheet id for each strand */
    int * length;
} Descr;


typedef struct{
    double ** full;
    double ** cm; 
    double ** translation; /* this one will change
			      with changing CM of the subset */
    double *  transl_norm; /* this one will change
			      with changing CM of the subset */
    double origin[3];      /* origin point that traslation vectors refer to */
    int * full_type;
    int * length;
    int N_full;
    int full_no_of_strands;
    int *full_sheet_id;
    
    double ** compact;
    int *compact_type;
    int N_compact;
    int **represents;
    int * is_rep_by;
    
} Representation;


# define MAX_MULT_MATCH 100

# define MAP_MAX  30 /*number of maps to keep */
// # define MAP_MAX   1000 /*number of maps to keep */
typedef struct {
    int size;
    int matches;
    int *x2y, *y2x;
    int x2y_size, y2x_size;
    int *x2y_residue_level, *y2x_residue_level;
    int x2y_residue_l_size, y2x_residue_l_size;
    int *submatch_best;
    int res_almt_length;
    double assigned_score;
    double avg_length_mismatch;
    double score_with_children;
    double z_score;
    double compl_z_score;
    double F;
    double aln_score;
    double avg, avg_sq;
    double **cosine;
    double **image;
    double q[4]; /* rotation (quaternion) */
    double T[3]; /* translation (once we believe we know it) */
    double rmsd; /* rmsd for the centers of matched SSEs */
    double res_rmsd;
    double res_almt_score;
} Map;

typedef struct {
    double total_assigned_score;
    double z_score;
    double w_submap_z_score;
    double gap_score;
    double rmsd;
    double conformer;
    double avg_length_mismatch;
    double submap_avg_length_mismatch;
    /* scores that appear in potprocessing */
    double res_almt_score;
    double res_rmsd; /* rmsd for aligned CAs */
    int res_almt_length;
    int number_of_maps;
    double q_rmsd[3];
} Score;

/*************************************/
/* penalty scheme for the alignment  */
/*************************************/
typedef struct {
    int endgap_special_treatment;
    double gap_opening;
    double gap_extension;
    double endgap;
    double *custom_gap_penalty_x; /* special, per-element gap penalties*/
    double *custom_gap_penalty_y;

} Penalty_parametrization;


/******************************/
/* global data:               */
/******************************/
extern Options options;  /* defined in struct.c */

/******************************/
/* integral lookup table :    */
/******************************/



# define MAX_EXP_VALUE 10
# define NR_POINTS 500
# define TABLE_SIZE NR_POINTS+1

# define INTEGRAL_TABLE "/home/ivanam/baubles/struct/08_data/tmp.1000.table"


extern double int_table [TABLE_SIZE][TABLE_SIZE];
extern double exp_table [TABLE_SIZE];

/******************************/
/* function declarations :    */
/******************************/

int check_gap_lengths  (Map * map, double *gap_score );

int complement_match_old (Representation* X_rep, Representation* Y_rep,
		      Map * map, int map_max,
		      int * map_ctr, int * map_best, int best_max, int parent_map);
int construct_translation_vecs ( Representation *X_rep,  Representation *Y_rep,
				 Map *map );
int descr_out (FILE * fptr, Descr *descr);
int descr_init ( Descr * description);
int descr_shutdown ( Descr * description );

double F (double **X, int * x_type,  int NX,
	  double **Y, int * y_type, int NY, double alpha);
int F_moments (double **x, int *x_type, int NX,
	       double **y, int *y_type, int NY, double alpha,
	       double *avg_ptr, double *avg_sq_ptr,
	       double *std_ptr);
int find_next_pair (double **X, double **Y, 
		    int *x_type, int *y_type,
		    int NX, int NY, int nbr_range,
		    int x_ctr, int y_ctr,
		    int *x_next_ptr, int *y_next_ptr);
int find_quat (double ** Xin, double **Yin, int no_vectors,
		 double * quat, double *qof);
int find_quat_exp (double ** X, int NX, double **Y, int NY,
		     double alpha, double * quat, double *qof);
int free_map (Map *map);
int geometric_descriptors (Protein * protein);
int get_next_descr (FILE * fptr,  Descr * description);
int initialize_map (Map *map, int NX, int NY );
int input  (FILE * fptr, Descr * description);
double lookup ( double alpha, double beta);
int map_assigned_score ( Representation *X_rep,  Map* map);
int map_complementarity (Map *map1,Map *map2,  double *z);
//int map_consistence ( int NX, int NY, Map *map1, Map *map2,
//		     double *total_ptr, double * gap_score, FILE *fptr);
int map_consistence (  int NX, int NY,  Map * combined_map, Map *map1, Map *map2,
		       double *total_ptr,  double *gap_score, FILE *fptr) ;
int mat_out (double A[4][4], char *name);
int mat_mult (double new [4][4], double  A[4][4], double  B[4][4]);
int mat_sum (double sum[4][4], double  new_term[4][4]);
int mat_diag ( double B[4][4], double eval[4], double evect[4][4] );
int mat_exp (double expB[4][4], double B[4][4]);
int match_length (int N, int *x2y);
int multiply (double *quat1, double *quat2_in,
		int conjugate,  double *product);
int needleman_wunsch (int max_i, int max_j, double **distance,
		      int *map_i2j, int * map_j2i, double *aln_score);

int smith_waterman (Penalty_parametrization *params, int max_i, int max_j, double **similarity,
		    int *map_i2j, int * map_j2i, double * aln_score);
int scoring (int NX, int NY, double **similarity, int * x2y, int * y2x, double * alignment );
void similarity_to_scoring(int NX, int NY, int multiplier, double** similarity, int ** scoring );

int normalized_cross (double *x, double *y, double * v, double *norm_ptr);
int output (FILE * fptr,char *name, char chain,  Protein * protein);
int postprocess (Descr *descr1, Protein * structure1, Representation *rep1, 
		 Descr *descr2, Protein * structure2, Representation *rep2, 
		 Map *map, Score * score);

int print_map (FILE *fptr, Map * map, Descr * descr1, Descr * descr2,
	       Protein *protein1, Protein *protein2,  int tab);
    
double qnorm (double *quat);
double q_dotprod (double *quat1, double *quat2);
int quat_to_R ( double quat[4], double **R);
double quat_rmsd (double parent_map_q[4], double  current_map_q[4]);
int random_q ( double exp_s[4],  double theta_step );
int read_integral_table (char * file_name );
int read_pdb (char * pdbname,  char chain, Protein * protein, int postp_flag);
int rep_initialize (Representation * rep,  Descr * descr  );
int rep_shutdown  (Representation * rep);
int rotate(double **Ynew, int NY, double **R, double ** Y);
int set_match_algebra ();
int set_up_exp_table ();
char single_letter ( char code[]);
int store_image (Representation *X_rep,  Representation *Y_rep, 
		 double **R, double alpha,  Map *map);
int unnorm_dot (double *x, double *y, double * dot);
int vec_out (double *vec, int dim,  char * name );

int find_Calpha (Protein *protein, int  resctr, double ca[3] );
double two_point_distance (double point1[3], double point2[3]);
int point_rot_tr (double point_in[3], double **R, double T[3],double point_out[3]); 
int print_score(FILE * digest, Descr * qry_descr, Descr * tgt_descr, Score * score, double fraction_assigned, int overlap) ;

int find_best_triples(Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
        double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
        double **best_quat);
int find_best_triples_old(Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
        double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
        double **best_quat);
int assign_score(Map *map, int *map_best, Representation *X_rep, Representation *Y_rep, Score *score);

int find_submaps(int NX, int NY, Map * map, int * map_best);
Map * init_map(int NX, int NY, int map_max);
int clean_map(Map *map, int map_max) ;
# endif
