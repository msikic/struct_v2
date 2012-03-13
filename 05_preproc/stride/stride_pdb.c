
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stride_types.h"
#include "stride_util.h"
#include "stride_cmd.h"
#include "stride_chain.h"

#include "stride_pdb.h"

static int Process_MODEL(enum METHOD *Method);
static int Process_ENDMDL(BUFFER Buffer, CHAIN **Chain, int *ChainNumber);
static int Process_ATOM(BUFFER Buffer, CHAIN **Chain, int *ChainNumber, 
			BOOLEAN *First_ATOM, StrideCmd *Cmd);

int ReadPDBFile(FILE * pdb, CHAIN **Chain, int *Cn, StrideCmd *Cmd)
{

  int ChainCnt, InfoCnt, i;
  enum METHOD Method = XRay;
  BOOLEAN First_ATOM, Published=YES, DsspAssigned=NO;
  float Resolution = 0.0;
  /*FILE *pdb;*/
  BUFFER Buffer;
  char *Info[MAX_INFO], PdbIdent[5];
  RESIDUE *r;
  CHAIN *c;

  *Cn= 0;
  InfoCnt = 0;
  strcpy(PdbIdent,"~~~~");

  /*
  if( !(pdb = fopen(Cmd->InputFile,"r")) )
    return(FAILURE); 
  */
     
  First_ATOM = YES;
  
  while( fgets(Buffer,BUFSZ,pdb) ) {
    
    if(!strncmp(Buffer,"HEADER",6)) {
      Info[InfoCnt] = (char *)ckalloc(BUFSZ*sizeof(char));
      strcpy(Info[InfoCnt],"HDR  ");
      strcat(Info[InfoCnt++],Buffer+10);
      strncpy(PdbIdent,Buffer+62,4);
      PdbIdent[4] = '\0';
    }
    /*
    else
    if(!strncmp(Buffer,"AUTHOR",6)) {
      Info[InfoCnt] = (char *)ckalloc(BUFSZ*sizeof(char));
      strcpy(Info[InfoCnt],"AUT  ");
      strcat(Info[InfoCnt++],Buffer+10);
    }
    else
    if(!strncmp(Buffer,"SOURCE",6)) {
      Info[InfoCnt] = (char *)ckalloc(BUFSZ*sizeof(char));
      strcpy(Info[InfoCnt],"SRC  ");
      strcat(Info[InfoCnt++],Buffer+10);
    }
    else
    if(!strncmp(Buffer,"COMPND",6)) {
      if( !Process_COMPND(Buffer,&Method) ) 
	return(FAILURE);
      else {
	Info[InfoCnt] = (char *)ckalloc(BUFSZ*sizeof(char));
	strcpy(Info[InfoCnt],"CMP  ");
	strcat(Info[InfoCnt++],Buffer+10);
      }
    }
    else if(!strncmp(Buffer,"JRNL",4) && !Process_JRNL(Buffer,&Published)) 
	return(FAILURE);
    else if(!strncmp(Buffer,"REMARK",6) && !Process_REMARK(Buffer,&Method,&Resolution,
							   &DsspAssigned)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"EXPDTA",6) && !Process_EXPDTA(Buffer,&Method)) 
      return(FAILURE);
    */
    else if(!strncmp(Buffer,"MODEL",5) && !Process_MODEL(&Method)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"ENDMDL",6)) {
      Process_ENDMDL(Buffer,Chain,Cn);
      break;
    }
    /*
    else if(!strncmp(Buffer,"HELIX",5) && !Process_HELIX(Buffer,Chain,Cn,Cmd)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"SHEET",5) && !Process_SHEET(Buffer,Chain,Cn,Cmd)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"TURN",4) && !Process_TURN(Buffer,Chain,Cn,Cmd)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"SSBOND",6) && !Process_SSBOND(Buffer,Chain,Cn,Cmd)) 
      return(FAILURE);
    */
    else if(!strncmp(Buffer,"ATOM",4) && !Process_ATOM(Buffer,Chain,Cn,&First_ATOM,Cmd)) 
      return(FAILURE);
  }
  /*
  fclose(pdb);
  */

  for( ChainCnt=0; ChainCnt< *Cn; ChainCnt++ ) {
    c = Chain[ChainCnt];
    if( c->NRes != 0  && !FindAtom(c,c->NRes,"CA",&i) )
      c->NRes--;
    strcpy(c->File,Cmd->InputFile);

    strcpy(c->PdbIdent,PdbIdent);
    if( c->NRes != 0 )  c->NRes++;
    if( c->NSheet != -1 ) c->NSheet++;
    c->Resolution = Resolution;
    c->Method = Method;
    c->Published = Published;
    c->DsspAssigned = DsspAssigned;
    c->NInfo = InfoCnt;
    for(i=0; i<InfoCnt; i++) {
      c->Info[i] = (char *)ckalloc(BUFSZ*sizeof(char));
      strcpy(c->Info[i],Info[i]);
      c->Info[i][71] = '\0';
    }
    for( i=0; i<c->NRes; i++ ) {
      r = c->Rsd[i];
      r->Inv =  (INVOLVED *)ckalloc(sizeof(INVOLVED));
      r->Prop = (PROPERTY *)ckalloc(sizeof(PROPERTY));
      r->Inv->NBondDnr = 0;
      r->Inv->NBondAcc = 0;
      r->Inv->InterchainHBonds = NO;
      r->Prop->Asn     = 'C';
      /*
      r->Prop->PdbAsn  = 'C';
      r->Prop->DsspAsn = 'C';
      r->Prop->Solv    = 0.0;
      */
      r->Prop->Phi     = 360.0;
      r->Prop->Psi     = 360.0;
    }
  }
  
  for(i=0; i<InfoCnt; i++)
    free(Info[i]);
  
  return(SUCCESS);
}


static int Process_MODEL(enum METHOD *Method) {
  if( *Method == XRay ) 
    *Method = Model;
  return(SUCCESS);
}

static int Process_ENDMDL(BUFFER Buffer, CHAIN **Chain, int *ChainNumber) {
  int CC;
  for( CC=0; CC < *ChainNumber; CC++ )
    Chain[CC]->Ter = 1;
  return(SUCCESS);
}

static int Process_ATOM(BUFFER Buffer, CHAIN **Chain, int *ChainNumber, 
			BOOLEAN *First_ATOM, StrideCmd *Cmd) {
  char *Field[MAX_FIELD];
  BUFFER Tmp;
  int CC, NR, NA;
  static char LastRes[MAX_CHAIN][RES_FIELD];
  RESIDUE *r;

  if( Cmd->NActive && !ChInStr(Cmd->Active,SpaceToDash(Buffer[21])) )
     return(SUCCESS);

  if( Buffer[16] != 'A' && Buffer[16] != ' ' && Buffer[16] != '1' ) 
    return(SUCCESS);

  if( *First_ATOM ) {
    for( CC=0; CC<MAX_CHAIN; CC++ ) 
      strcpy(LastRes[CC],"XXXX");
    *First_ATOM = NO;
  }
  
  for( CC=0; CC < *ChainNumber && Chain[CC]->Id != Buffer[21] ; CC++ );
  
  if( CC == *ChainNumber ) {
    InitChain(&Chain[CC]); 
    Chain[CC]->Id = Buffer[21];
    (*ChainNumber)++;
  }
  else
  if( Chain[CC]->Ter == 1 ) 
    return(SUCCESS);

  if( Buffer[34] != '.' || Buffer[42] != '.' || Buffer[50] != '.' )
    return(escape(FAILURE,"File %s has no coordinates\n",Cmd->InputFile));

  
  if( Cmd->Stringent && Buffer[63] != '.')
    return(escape(FAILURE,"File %s has no temperature factor\n",Cmd->InputFile));


  SplitString(Buffer+22,Field,1);
  if( strcmp(Field[0],LastRes[CC]) ) {
    if( strcmp(LastRes[CC],"XXXX") && !FindAtom(Chain[CC],Chain[CC]->NRes,"CA",&NA) ) {
      free(Chain[CC]->Rsd[Chain[CC]->NRes]);
      Chain[CC]->NRes--;
    }
    if( strcmp(LastRes[CC],"XXXX") ) Chain[CC]->NRes++;
    NR = Chain[CC]->NRes;
    strcpy(LastRes[CC],Field[0]);
    Chain[CC]->Rsd[NR] = (RESIDUE *)ckalloc(sizeof(RESIDUE));
    strcpy(Chain[CC]->Rsd[NR]->PDB_ResNumb,LastRes[CC]);
    Chain[CC]->Rsd[NR]->NAtom = 0;
    SplitString(Buffer+17,Field,1);
    strcpy(Chain[CC]->Rsd[NR]->ResType,Field[0]);
  }
  else 
    NR = Chain[CC]->NRes;
  
  NA = Chain[CC]->Rsd[NR]->NAtom;

  if( Buffer[16] != ' ' ) {
    strcpy(Tmp,Buffer);
    Tmp[16] = ' ';
    SplitString(Tmp+12,Field,1);
  }
  else
    SplitString(Buffer+12,Field,1);
  
  r = Chain[CC]->Rsd[NR];
  strcpy(r->AtomType[NA],Field[0]);


  strcpy(Tmp,Buffer);
  Buffer[38] = ' ';
  SplitString(Tmp+30,Field,1);
  r->Coord[NA][0] = atof(Field[0]);

  strcpy(Tmp,Buffer);
  Buffer[46] = ' ';
  SplitString(Tmp+38,Field,1);
  r->Coord[NA][1] = atof(Field[0]);

  strcpy(Tmp,Buffer);
  Buffer[54] = ' ';
  SplitString(Tmp+46,Field,1);
  r->Coord[NA][2] = atof(Field[0]);

  if( Buffer[57] == '.' ) {
    strcpy(Tmp,Buffer);
    Tmp[60] = ' ';
    SplitString(Tmp+54,Field,1);
    r->Occupancy[NA] = atof(Field[0]);
  }
  else 
    r->Occupancy[NA] = -1.00;
  
  SplitString(Buffer+63,Field,1);
  r->TempFactor[NA] = atof(Field[0]);

  r->NAtom++;

  if( r->NAtom > MAX_AT_IN_RES-1 )
    return(escape(FAILURE,"File %s has too many atoms in residue %s %s %c\n",
		  Cmd->InputFile,r->ResType,r->PDB_ResNumb,Chain[CC]->Id));

  return(SUCCESS);
}
