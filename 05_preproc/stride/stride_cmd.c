
#include <string.h>

#include "stride_types.h"
#include "stride_cmd.h"

const char * strideCmdInit(StrideCmd * cmd,
			   const char * chains,
			   const char * pdbfilename) {

  if (strlen(pdbfilename) > BUFSZ)
    return "strideCmdInit(): pdbfilename too long";
  strncpy(cmd->InputFile, pdbfilename, BUFSZ);

  cmd->NActive = strlen(chains); 
  if (cmd->NActive > MAX_CHAIN)
    return "strideCmdInit(): too many chains";
  strcpy(cmd->Active, chains);

  cmd->NProcessed = 0;
  strcpy(cmd->Processed,"");

  strcpy(cmd->Cond,"");

  cmd->SideChainHBond    = NO;
  cmd->MainChainHBond    = YES;
  cmd->MainChainPolarInt = YES;
  cmd->Info              = NO;
  cmd->Truncate          = YES;
  cmd->Stringent         = NO;

  cmd->EnergyType        = 'G';

  cmd->DistCutOff        =  6.0;
  cmd->PhiPsiStep        =  0.0;

  cmd->C1_H              = -1.0;
  cmd->C2_H              =  1.0;
  cmd->C1_E              = -0.2;
  cmd->C2_E              =  0.2;

  cmd->Treshold_H1       = -230.0;
  cmd->Treshold_H3       =  0.12;
  cmd->Treshold_H4       =  0.06;
  cmd->Treshold_E1       = -240.0;
  cmd->Treshold_E2       = -310.0;

  cmd->MinResolution     =  0.1;
  cmd->MaxResolution     =  100.0;

  cmd->MinLength         = 0;
  cmd->MaxLength         = MAX_RES;

  cmd->NPixel            = 0;

  return 0;
  
  /*
    Cmd->SideChainHBond    = NO;
  Cmd->MainChainHBond    = YES;
  Cmd->MainChainPolarInt = YES;
  Cmd->Published         = NO;
  Cmd->DsspAssigned      = NO;
  Cmd->UseResolution     = NO;
  Cmd->Info              = NO;
  Cmd->Truncate          = YES;
  Cmd->ExposedArea       = YES;
  Cmd->ReportSummaryOnly = NO;
  Cmd->ReportBonds       = NO;
  Cmd->BrookhavenAsn     = NO;
  Cmd->DsspAsn           = NO;
  Cmd->MolScript         = NO;
  Cmd->OutSeq            = NO;
  Cmd->Stringent         = NO;
  Cmd->Measure           = NO;

  Cmd->EnergyType        = 'G';

  Cmd->DistCutOff        =  6.0;
  Cmd->PhiPsiStep        =  0.0;

  Cmd->C1_H              = -1.0;
  Cmd->C2_H              =  1.0;
  Cmd->C1_E              = -0.2;
  Cmd->C2_E              =  0.2;

  Cmd->Treshold_H1       = -230.0;
  Cmd->Treshold_H3       =  0.12;
  Cmd->Treshold_H4       =  0.06;
  Cmd->Treshold_E1       = -240.0;
  Cmd->Treshold_E2       = -310.0;

  Cmd->MinResolution     =  0.1;
  Cmd->MaxResolution     =  100.0;

  Cmd->MinLength         = 0;
  Cmd->MaxLength         = MAX_RES;
  
  Cmd->NPixel            = 0;
  Cmd->NActive           = 0;
  Cmd->NProcessed        = 0;

  strcpy(Cmd->FirstResidue,"");
  strcpy(Cmd->LastResidue,"");

  strcpy(Cmd->MapFileHelix,""); 
  strcpy(Cmd->MapFileSheet,""); 
  strcpy(Cmd->OutFile,"");
  strcpy(Cmd->Active,"");
  strcpy(Cmd->Processed,"");
  strcpy(Cmd->Cond,"");
  */
}
