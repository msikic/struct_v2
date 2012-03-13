

/* escape.c */
int escape(int RetVal, char *format, ... );

/* memory.c */
void *ckalloc(size_t bytes);

/* splitstr.c */
int SplitString(char *Buffer, char **Fields, int MaxField);

/* strutil.c */
char SpaceToDash(char Id);
BOOLEAN ChInStr(char *String, char Char);
void ExtractAsn(CHAIN **Chain, int Cn, char *Asn);
BOOLEAN ExistsSecStr(CHAIN **Chain, int NChain);
int PdbN2SeqN(CHAIN *Chain, char *PdbN, int *SeqN);
int FindAtom(CHAIN *Chain, int ResNumb, char *Atom, int *AtNumb);
int Boundaries(char *Asn, int L, char SecondStr, int (*Bound)[2]);
int FindChain(CHAIN **Chain, int NChain, char ChainId);
/*
float SecStrContent(CHAIN *Chain, int *HelAlp, int *HelPI, int *Hel310,
		    int *Sheet, int *Turn);
*/

