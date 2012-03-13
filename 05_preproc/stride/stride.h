
#define STRIDE_RESNUM  16
#define STRIDE_RESTYPE  3

typedef struct StrideSseStruct {
  char subtype;
  char chain_id;
  char bgn_resnum[STRIDE_RESNUM+1];
  char bgn_restype[STRIDE_RESTYPE+1];
  char end_resnum[STRIDE_RESNUM+1];
  char end_restype[STRIDE_RESTYPE+1];
} StrideSse;

typedef struct StrideDataStruct {
  const char * errmsg;
  int n_helices;
  StrideSse* helices;
  int n_strands;
  StrideSse* strands;
} StrideData;

void strideInit(StrideData * sd);
void strideRead(StrideData * sd, const char * chains,
		FILE * fptr, const char * pdbfilename);
void strideDone(StrideData * sd);
