
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stride.h"

#include "stride_types.h"
#include "stride_util.h"
#include "stride_cmd.h"
#include "stride_pdb.h"
#include "stride_chain.h"
#include "stride_phipsi.h"
#include "stride_map.h"
#include "stride_hbond.h"
#include "stride_helix.h"
#include "stride_sheet.h"

static char * strideSseAdd(StrideSse** sses, int *n_sses,
			   char subtype, char chain_id,
			   char * bgn_resnum, char * bgn_restype,
			   char * end_resnum, char * end_restype) {
  StrideSse* sse;

  (*n_sses)++;
  *sses = realloc(*sses, (*n_sses) * sizeof(StrideSse));
  if (!*sses)
    return "strideSseAdd(): memory allocation failed";
  sse = &((*sses)[(*n_sses) - 1]);
  
  sse->subtype = subtype;
  sse->chain_id = chain_id;

  if (strlen(bgn_resnum) > STRIDE_RESNUM)
    return "strideSseAdd(): bgn_resnum too big";
  strcpy(sse->bgn_resnum, bgn_resnum);

  if (strlen(bgn_restype) > STRIDE_RESTYPE)
    return "strideSseAdd(): bgn_restype too big";
  strcpy(sse->bgn_restype, bgn_restype);

  if (strlen(end_resnum) > STRIDE_RESNUM)
    return "strideSseAdd(): end_resnum too big";
  strcpy(sse->end_resnum, end_resnum);

  if (strlen(end_restype) > STRIDE_RESTYPE)
    return "strideSseAdd(): end_restype too big";
  strcpy(sse->end_restype, end_restype);

  return 0;
}

void strideInit(StrideData * sd) {
  sd->errmsg = 0;
  sd->n_helices = 0;
  sd->helices = 0;
  sd->n_strands = 0;
  sd->strands = 0;
}

#define MAX_ASSIGN 500

void strideRead(StrideData * sd, const char * chains,
		FILE * fptr, const char * pdbfilename) {
  StrideCmd cmd;
  const char * errmsg;
  int NChain;
  CHAIN ** Chain;
  int NHBond;
  HBOND **HBond;
  int Cn, ValidChain, i;
  float **PhiPsiMapHelix, **PhiPsiMapSheet;
  /* Working vars for filling StrideData */
  char * Asn;
  int NStr;
  int Bound[MAX_ASSIGN][2];
  CHAIN* chnptr;
  char chain_id;
  RESIDUE* rsd_bgn;
  RESIDUE* rsd_end;
  char * helix_types = "HGI";
  char * ht;

  if (strchr(chains, ' ') != 0) {
    sd->errmsg = "strideRead(): pass chains=\"\" to process chain=' '";
    return;
  }

  NChain = 0;
  Chain = (CHAIN  **)ckalloc(MAX_CHAIN*sizeof(CHAIN *));
  NHBond = 0;
  HBond = (HBOND  **)ckalloc(MAXHYDRBOND*sizeof(HBOND *));

  errmsg = strideCmdInit(&cmd, chains, pdbfilename);
  if (errmsg) {
    sd->errmsg = errmsg;
    return;
  }

  NChain = 0;
  if( !ReadPDBFile(fptr, Chain, &NChain, &cmd) || !NChain ) {
    sd->errmsg = "ReadPDBFile(): failed";
    return;
  }

  ValidChain = 0;
  for(Cn = 0; Cn < NChain; Cn++)
    ValidChain += CheckChain(Chain[Cn], &cmd);
  if( !ValidChain ) {
    sd->errmsg = "strideRead(): no valid chains";
    return;
  }

  BackboneAngles(Chain, NChain);

  PhiPsiMapHelix = DefaultHelixMap(&cmd);
  PhiPsiMapSheet = DefaultSheetMap(&cmd);

  for(Cn = 0; Cn < NChain; Cn++)
    PlaceHydrogens(Chain[Cn]);

  errmsg =  FindHydrogenBonds(Chain, NChain, HBond, &NHBond, &cmd);
  if (errmsg) {
    sd->errmsg = errmsg;
    return;
  }
  if (NHBond == 0) {
    sd->errmsg = "strideRead(): no hydrogen bonds found";
    return;
  }

  NoDoubleHBond(HBond, NHBond);

  DiscrPhiPsi(Chain, NChain, &cmd);

  /* Area(Chain,NChain,Cmd); ?Add eventually? */


  for(Cn = 0; Cn < NChain; Cn++) {

    if(Chain[Cn]->Valid) {
    
      Helix(Chain, Cn, HBond, &cmd, PhiPsiMapHelix);
      
      for(i = 0; i < NChain; i++ )
        if(Chain[i]->Valid )
          Sheet(Chain, Cn, i, HBond, &cmd, PhiPsiMapSheet);    
      
      /*BetaTurn(Chain,Cn);
      GammaTurn(Chain,Cn,HBond);
      */
    }
  }

  /* unnecessary? */
  if(!ExistsSecStr(Chain, NChain)) {
    sd->errmsg = "strideRead(): no secondary structure found";
    return;
  }


  for(Cn = 0; Cn < NChain; Cn++) {
    if(!Chain[Cn]->Valid)
      continue;

    chnptr = Chain[Cn];
    chain_id = chnptr->Id;

    Asn = (char *) ckalloc(chnptr->NRes * sizeof(char));
    ExtractAsn(Chain, Cn, Asn);

    /* alpha-helix = 'H' */
    /* 310-helix   = 'G' */
    /* pi-helix    = 'I' */

    ht = helix_types;
    while (*ht != '\0') {

      NStr = Boundaries(Asn, chnptr->NRes, *ht, Bound);
      for(i = 0; i < NStr; i++) {
	rsd_bgn = chnptr->Rsd[Bound[i][0]];
	rsd_end = chnptr->Rsd[Bound[i][1]];
	errmsg = strideSseAdd(&(sd->helices), &(sd->n_helices),
			      *ht, chain_id,
			      rsd_bgn->PDB_ResNumb, rsd_bgn->ResType,
			      rsd_end->PDB_ResNumb, rsd_end->ResType);
	if (errmsg) {
	  sd->errmsg = errmsg;
	  return;
	}
      }
      ht++;
    }

    NStr = Boundaries(Asn, chnptr->NRes, 'E', Bound);
    for(i = 0; i < NStr; i++) {
      rsd_bgn = chnptr->Rsd[Bound[i][0]];
      rsd_end = chnptr->Rsd[Bound[i][1]];
      errmsg = strideSseAdd(&(sd->strands), &(sd->n_strands),
			    'E', chain_id,
			    rsd_bgn->PDB_ResNumb, rsd_bgn->ResType,
			    rsd_end->PDB_ResNumb, rsd_end->ResType);
      if (errmsg) {
	sd->errmsg = errmsg;
	return;
      }
    }

    free(Asn);

  }


  if (Chain != 0)
    free(Chain);
  if (HBond != 0)
    free(HBond);
}

void strideDone(StrideData * sd) {
  sd->errmsg = 0;

  sd->n_helices = 0;
  if (sd->helices != 0)
    free(sd->helices);

  sd->n_strands = 0;
  if (sd->strands != 0)
    free(sd->strands);
}
