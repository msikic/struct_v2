
/* pdb parser errors */
#define ERR_SSE_NONE    2
#define ERR_CA_TRACE    4
#define ERR_NONSENSE    8
#define ERR_CONTAINER  16
#define ERR_OVERLAP    32
#define ERR_MAX_ATOMS  64
#define ERR_NO_FILE_OR_CHAIN  128

#ifdef __cplusplus
extern "C" 
#endif
int readPDB(const char *pdbname,  char chain, Protein *protein);
