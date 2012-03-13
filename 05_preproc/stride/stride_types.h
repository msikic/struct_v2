
#define Pi                        3.1415927
#define Eps                       0.000001
#define RAD(x)                    (x)*Pi/180.0
#define DEG(x)                    (x)*180.0/Pi
#define RADDEG                    57.2958
#define BREAKDIST                 2.5

#define SUCCESS                   1
#define FAILURE                   0
#define YES                       1
#define NO                        0
#define ERR                      -1

#define BUFSZ                     1024
#define MAX_FIELD                 50

#define MAX_CHAIN                 100
#define MAX_RES                   20000
#define MAX_HELIX                 500
#define MAX_SHEET                 500
#define MAX_STRAND_IN_SHEET       20
#define MAX_INFO                  1000
#define MAX_AT_IN_RES             75
#define MAXRESDNR                 6
#define MAXRESACC                 6 

#define RES_FIELD                 6
#define AT_FIELD                  5

#define MAXCONDITIONS             20

#define MAXHYDRBOND               50000

#define MINPHIPSI                -180.0
#define MAXPHIPSI                 180.0


enum METHOD {XRay, NMR, Model};
enum HYBRID {Nsp2, Nsp3, Osp2, Osp3, Ssp3};
enum GROUP  {Peptide, Trp, Asn, Gln, Arg, His, Lys, Ser, Thr,
	     Tyr, Asp, Glu, Met, Cys};

typedef char BUFFER[BUFSZ+1];

typedef int BOOLEAN;


typedef struct {
  float Phi, Psi, Omega;
  int PhiZn, PsiZn;
  /*float Solv, DsspSolv;*/
  char Asn; /*, PdbAsn, DsspAsn;*/
} PROPERTY;

typedef struct {
  int HBondDnr[MAXRESDNR];
  int HBondAcc[MAXRESACC];
  int NBondDnr, NBondAcc;
  BOOLEAN InterchainHBonds;
} INVOLVED;


typedef struct {
  int NAtom;
  char PDB_ResNumb[RES_FIELD];
  char ResType[RES_FIELD];
  char AtomType[MAX_AT_IN_RES][AT_FIELD];
  float Coord[MAX_AT_IN_RES][3];
  float Occupancy[MAX_AT_IN_RES];
  float TempFactor[MAX_AT_IN_RES];
  PROPERTY *Prop;
  INVOLVED *Inv;
} RESIDUE;


typedef struct {
  char Res1[RES_FIELD];
  char Res2[RES_FIELD];
  char PDB_ResNumb1[RES_FIELD], PDB_ResNumb2[RES_FIELD];
  char InsCode1, InsCode2;
  int Class;
} HELIX;

typedef struct {
  int  NStrand;
  char SheetId[RES_FIELD];
  char ResType1[MAX_STRAND_IN_SHEET][RES_FIELD];
  char ResType2[MAX_STRAND_IN_SHEET][RES_FIELD];
  char PDB_ResNumb1[MAX_STRAND_IN_SHEET][RES_FIELD];
  char PDB_ResNumb2[MAX_STRAND_IN_SHEET][RES_FIELD];
  char InsCode1[MAX_STRAND_IN_SHEET];
  char InsCode2[MAX_STRAND_IN_SHEET];
  int  Sence[MAX_STRAND_IN_SHEET];
  int  RegYN[MAX_STRAND_IN_SHEET];
  char AtomNameReg1[MAX_STRAND_IN_SHEET][AT_FIELD];
  char AtomNameReg2[MAX_STRAND_IN_SHEET][AT_FIELD];
  char ResTypeReg1[MAX_STRAND_IN_SHEET][RES_FIELD];
  char ResTypeReg2[MAX_STRAND_IN_SHEET][RES_FIELD];
  char PDB_ResNumbReg1[MAX_STRAND_IN_SHEET][RES_FIELD];
  char PDB_ResNumbReg2[MAX_STRAND_IN_SHEET][RES_FIELD];
  char InsCodeReg1[MAX_STRAND_IN_SHEET];
  char InsCodeReg2[MAX_STRAND_IN_SHEET];
} SHEET;

typedef struct {
  char Res1[RES_FIELD];
  char Res2[RES_FIELD];
  char PDB_ResNumb1[RES_FIELD], PDB_ResNumb2[RES_FIELD];
  char InsCode1, InsCode2;
  char TurnType;
} TURN;

typedef struct {
  int NRes, NHetRes, NonStandRes, Ter;
  int NHet, NAtom, NonStandAtom, NHelix, NSheet;
  int NTurn, NAssignedTurn, NBond, NHydrBond, NHydrBondInterchain, NHydrBondTotal, NInfo;
  char Id, *File;
  float Resolution;
  enum METHOD Method;
  BOOLEAN Valid, Published, DsspAssigned;

  RESIDUE **Rsd;
  /*
  HETERORESIDUE **HetRsd;
  HET     **Het;
  */
  HELIX   **Helix;
  SHEET   **Sheet;
  TURN    **Turn;
  /*
  TURN    **AssignedTurn;
  SSBOND  **SSbond;
  */
  char    **Info;
  char    PdbIdent[5];
} CHAIN;

typedef struct {
  CHAIN *Chain;
  int D_Res, DD_Res, DDI_Res;
  int D_At, DD_At, DDI_At, H;
  enum HYBRID Hybrid;
  enum GROUP Group;
  float HB_Radius;
} DONOR;

                 
typedef struct {
  CHAIN *Chain;
  int A_Res, AA_Res, AA2_Res;
  int A_At, AA_At, AA2_At;
  enum HYBRID Hybrid;
  enum GROUP Group;
  float HB_Radius;
} ACCEPTOR;


typedef struct {
  DONOR *Dnr;
  ACCEPTOR *Acc;
  BOOLEAN ExistPolarInter, ExistHydrBondRose, ExistHydrBondBaker;
  float Energy, Er, Et, Ep, ti, to, p;
  float AccDonDist, OHDist, AngNHO, AngCOH;
  float AccAng, DonAng, AccDonAng, DonAccAng;
} HBOND;


typedef struct PAT {
  HBOND *Hb1, *Hb2;
  struct PAT *Nei1, *Nei2;
  BOOLEAN ExistPattern;
  BUFFER Type;
} PATTERN;
