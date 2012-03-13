
typedef struct StrideCmdStruct {
  BUFFER InputFile;

  int NActive;
  char Active[MAX_CHAIN+1]; 

  int NProcessed;
  char Processed[MAX_CHAIN+1];

  char Cond[MAXCONDITIONS];

  BOOLEAN SideChainHBond, MainChainHBond, MainChainPolarInt;
  BOOLEAN Info, Truncate;
  BOOLEAN Stringent;

  char EnergyType;

  float DistCutOff, PhiPsiStep;

  float C1_H, C2_H, C1_E, C2_E;
  float Treshold_H1, Treshold_H2, Treshold_H3, Treshold_H4;
  float Treshold_E1, Treshold_E2, Treshold_E3, Treshold_E4;

  float MinResolution, MaxResolution;
  int MinLength, MaxLength;

  int NPixel;

  /*
  BUFFER InputFile, OutFile, SeqFile;
  BUFFER MapFileHelix, MapFileSheet;
  BUFFER MolScriptFile, DsspFile;
  char EnergyType, Active[MAX_CHAIN+1]; 
  char Processed[MAX_CHAIN+1], Cond[MAXCONDITIONS];
  char FirstResidue[RES_FIELD], LastResidue[RES_FIELD];
  
  int NPixel, NActive, NProcessed;
  int MinLength, MaxLength;

  float PhiPsiStep, DistCutOff;
  float Treshold_H1, Treshold_H2, Treshold_H3, Treshold_H4;
  float Treshold_E1, Treshold_E2, Treshold_E3, Treshold_E4;
  float MinResolution, MaxResolution;
  float C1_H, C2_H, C1_E, C2_E;

  BOOLEAN SideChainHBond, MainChainHBond, MainChainPolarInt;
  BOOLEAN Published, DsspAssigned, UseResolution, Info, Truncate;
  BOOLEAN ExposedArea, ReportSummaryOnly, ReportBonds, BrookhavenAsn, DsspAsn;
  BOOLEAN MolScript, OutSeq, Stringent, Measure, ContactOrder, ContactMap;
  */

  /* Not used by STRIDE */
  /*
  BOOLEAN Shrink;
  BUFFER CarteFile, MapFile, MathFile;
  char AsnSource, Mode, SecStrType;
  int FilterOrder;
  int NStepA, NStepB, NStepC, NStepD;
  float Treshold_From_A, Treshold_To_A, StepA;
  float Treshold_From_B, Treshold_To_B, StepB;
  float Treshold_From_C, Treshold_To_C, StepC;
  float Treshold_From_D, Treshold_To_D, StepD;
  */
} StrideCmd;

const char * strideCmdInit(StrideCmd * cmd,
			   const char * chains,
			   const char * pdbfilename);
