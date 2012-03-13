#include "struct.h"

/* ss elements: */
# define HELIX 1
# define PARALLEL 2
# define PERP  4
/* maximal allowed overlap between structural elements */
# define MAX_ALLOWED_OVERLAP 3

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 100 /* so I can handle things like heme */

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+2];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    int no_atoms;
    Atom atom[MAX_NO_ATOMS];
    Atom *Ca;
    int interface;
    int solvent_accessible;
    int belongs_to_helix;
    int belongs_to_strand;
    int alt_belongs_to_helix; /* we allow a couple of residues overlap in SSEs */
    int alt_belongs_to_strand;
} Residue;



typedef struct {
    char begin_id[PDB_HELIX_BEGIN_LEN+2]; /* this will make reading in easier*/
    char end_id[PDB_HELIX_END_LEN+2]; /* this will make reading in easier*/
    int length;
    int begin, end; /* this may be added as post-processing step */
    int no_of_vectors;
    double **p, **cm;
} Helix;

typedef struct {
    char sheet_id[PDB_SHEET_SHEET_LEN+1];
    char begin_id[PDB_SHEET_BEGIN_LEN+2]; /* this will make reading in easier*/
    char end_id[PDB_SHEET_END_LEN+2]; /* this will make reading in easier*/
    int sheet_number;
    int length;
    int begin, end;
    int no_of_vectors;
    double **p, **cm, **perp, **foot; /* if the strand is curved, it will be
			 associated with an array of vectors (well, 2 for now*/ 
} Strand;

typedef struct {
    char sheet_id[PDB_SHEET_SHEET_LEN+1];
    int parallel;
    double p[3];
    int number_of_strands;
    Strand * first_strand;
} Sheet;


typedef struct {
    int length;
    Residue * sequence;
    int no_helices;
    Helix   * helix;
    int no_strands;
    Strand   * strand;
    int no_sheets;
    Sheet * sheet;
    int *sse_sequence;
} Protein;


int protein_shutdown (Protein * protein);
