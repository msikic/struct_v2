
/*#define MIN_HELIX_LENGTH  8*/
/*#define MIN_HELIX_LENGTH  3*/
#define MIN_HELIX_LENGTH  5
#define MIN_STRAND_LENGTH 3

# define NO_BB_ATOMS 4         /* for use in line_data */
# define MAX_ALLOWED_OVERLAP 3 /* maximal allowed overlap between SSEs */
# define MAX_NO_ATOMS 100      /* so I can handle things like heme */

/* ss elements: */
typedef enum  {HELIX, STRAND, NO_SSE_TYPES} SSE_type;

#define N_HELIX_TYPES 4
typedef enum {
  HELIX_TYPE_UNDEF=0,
  HELIX_TYPE_ALPHA=1,
  HELIX_TYPE_310=2,
  HELIX_TYPE_PI=3
} HelixType;

typedef struct {
  char type [PDB_ATOM_ATOM_NAME_LEN+1];
  double x,y,z;
  int backbone;
} Atom;


typedef struct {
  char pdb_id[PDB_ATOM_RES_NO_LEN+2];
  char res_type[PDB_ATOM_RES_NAME_LEN+1];
  char res_type_short;
  int  no_atoms;
  Atom  atom[MAX_NO_ATOMS];
  int interface;
  int solvent_accessible;
  int belongs_to_helix;
  int belongs_to_strand;
  int alt_belongs_to_helix; /* allow a couple of residues overlap in SSEs */
  int alt_belongs_to_strand;
  HelixType helix_type;
} Residue;

typedef struct {
    double **p, **cm;
    int *from, *to;
    int *no_residues;
    int no_atoms;
    Atom ** atom;    /*for each line: the set of bb atoms*/
} Line_data;
    
typedef struct {
    int length;
    int begin, end; /* this may be added as post-processing step */
    int n_lines;
    void* lines;
} Helix;

typedef struct {
    int length;
    int begin, end;
    int n_lines;
    void* lines;
} Strand;

typedef struct {
    int length;
    Residue * sequence;
    int no_helices;
    Helix   * helix;
    int no_strands;
    Strand   * strand;
    int *sse_sequence;
} Protein;
