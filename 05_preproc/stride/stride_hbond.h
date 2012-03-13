
/* place_h.c */
int PlaceHydrogens(CHAIN *Chain);

/* hydrbond.c */
int FindPolInt(HBOND **HBond, RESIDUE *Res1, RESIDUE *Res2);
int FindBnd(HBOND **HBond, RESIDUE *Res1, RESIDUE *Res2);
int NoDoubleHBond(HBOND **HBond, int NHBond);
const char * FindHydrogenBonds(CHAIN **Chain, int NChain,
			       HBOND **HBond, int * NHBond,
			       StrideCmd *Cmd);
