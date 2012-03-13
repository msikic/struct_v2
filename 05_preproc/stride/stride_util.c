
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "stride_types.h"
#include "stride_util.h"


int escape(int RetVal, char *format, ... ) {
  va_list ptr;
  va_start(ptr,format);
  vfprintf(stderr,format,ptr);
  va_end(ptr);
  return(RetVal);
}

void *ckalloc(size_t bytes) {
  register void *ret;
  if( !(ret = malloc(bytes)) ) {
    fprintf(stderr, "stride: ckalloc(): unable to malloc() %lu bytes: exiting", (unsigned long) bytes);
    exit(1);
  }
  return ret;   
}


int SplitString(char *Buffer, char **Fields, int MaxField) {
  int FieldCnt, SymbCnt, FieldFlag, BuffLen;
  static char LocalBuffer[BUFSZ];

  FieldCnt =0; FieldFlag = 0;
  BuffLen = (int)strlen(Buffer) - 1;

  strcpy(LocalBuffer,Buffer);

  for(SymbCnt=0; SymbCnt<BuffLen; SymbCnt++) {
    if( (isspace(LocalBuffer[SymbCnt])) &&
	FieldFlag == 0 &&
	SymbCnt != BuffLen-1 ) continue;
    if( (!isspace(LocalBuffer[SymbCnt])) &&
	FieldFlag == 1 &&
	SymbCnt == BuffLen-1 ) {
      LocalBuffer[SymbCnt+1] = '\0';
      return(FieldCnt);
    }
    else
    if( (isspace(LocalBuffer[SymbCnt])) && FieldFlag == 1 ) {
      LocalBuffer[SymbCnt] = '\0';
      FieldFlag = 0;
      if( FieldCnt == MaxField ) return(FieldCnt);
    }
    else
    if( (!isspace(LocalBuffer[SymbCnt])) && FieldFlag == 0 ) {
      FieldFlag = 1;
      Fields[FieldCnt] = LocalBuffer+SymbCnt;
      FieldCnt++;
    }
  }
  
  return(FieldCnt);
}


char SpaceToDash(char Id)
{
  static char NewId;

  if( Id == ' ' )
    NewId = '-';
  else
    NewId = Id;

  return(NewId);
}

BOOLEAN ChInStr(char *String, char Char)
{

  if( strchr(String,toupper(Char)) || 
      strchr(String,Char) ||
      strchr(String,tolower(Char)) )
    return(YES);
  
  return(NO);
}

void ExtractAsn(CHAIN **Chain, int Cn, char *Asn)
{
  register int Res;

  for( Res=0; Res<Chain[Cn]->NRes; Res++ )
    Asn[Res] = Chain[Cn]->Rsd[Res]->Prop->Asn;
}

BOOLEAN ExistsSecStr(CHAIN **Chain, int NChain)
{
  register int i, Cn;
  
  for( Cn=0; Cn<NChain; Cn++ )
    for( i=0; i<Chain[Cn]->NRes; i++ )
      if( Chain[Cn]->Rsd[i]->Prop->Asn != 'C' )
        return(YES);

  return(NO);
}


/*************************************************************************
**                                                                      **
** Find sequential residue number for a PDB residue number              **
**                                                                      **
** INPUT:  *Chain    Pointer to a protein chain                         **
**         *PdbN     String containing PDB residue number               **
**                                                                      **
** OUTPUT: *SeqN     Pointer to the sequential residue number           **
**                                                                      **
*************************************************************************/
int PdbN2SeqN(CHAIN *Chain, char *PdbN, int *SeqN)
{

  for( (*SeqN)=0; (*SeqN)<Chain->NRes; (*SeqN)++ )
    if( !strcmp(Chain->Rsd[(*SeqN)]->PDB_ResNumb,PdbN) ) 
      return(SUCCESS);

  return(FAILURE);
}

/*************************************************************************
**                                                                      **
** Find atom of specified type in a residue                             **
**                                                                      **
** INPUT:  *Chain    Pointer to a protein chain                         **
**         ResNumb   Number of residue in the protein chain             **
**         *Atom     String containing atom name                        **
**                                                                      **
** OUTPUT: *AtNumb   Pointer to the atom number in the residue          **
**                                                                      **
*************************************************************************/
int FindAtom(CHAIN *Chain, int ResNumb, char *Atom, int *AtNumb)
{

  for( (*AtNumb)=0; (*AtNumb)<Chain->Rsd[ResNumb]->NAtom; (*AtNumb)++ )
    if( !strcmp(Atom,Chain->Rsd[ResNumb]->AtomType[(*AtNumb)]) )
       return(SUCCESS);

  *AtNumb = ERR;
  return(FAILURE);
}

/*************************************************************************
**                                                                      **
** Find beginning and end residues of each secondary structure element  **
** in a secondary structure assignment                                  **
**                                                                      **
** INPUT:   *Asn       One letter secondary structure assignment        **
**          L          Length of the protein chain                      **
**          SecondStr  Secondary structure type                         **
**                                                                      **
** OUTPUT:  *(Bound)[2] Two dimensional array containing numbers of     **
**                      first and last residue of each secondary        **
**                      structure element                               **
**                                                                      **
** RETURNS: Number of the Secondary structure elements                  **
**                                                                      **
*************************************************************************/
int Boundaries(char *Asn, int L, char SecondStr, int (*Bound)[2])
{
  register int Res;
  int NStr = 0, Flag = 0;
  
  for( Res=0; Res<L; Res++ ) {
    if( Asn[Res] == SecondStr && Flag == 0 ) {
      Flag = 1; 
      Bound[NStr][0] = Res;
    }
    else
    if( Asn[Res] != SecondStr && Flag == 1 ) {
      Flag = 0; 
      Bound[NStr++][1] = Res-1;
    }
  }
  return(NStr);
}

/*************************************************************************
**                                                                      **
** Find protein chain with a given chain identifier number              **
**                                                                      **
** INPUT:  **Chain   Array of protein chains                            **
**         NChain    Number of chains                                   **
**         ChainId   Chain identifier                                   **
**                                                                      **
** RETURNS: Number of the protein chain with identifier ChainId         **
**                                                                      **
*************************************************************************/
int FindChain(CHAIN **Chain, int NChain, char ChainId)
{
  register int i;

  for( i=0; i<NChain; i++ )
    if( Chain[i]->Id == ChainId )
      return(i);

  return(ERR);
}

/*************************************************************************
**                                                                      **
** Calculate the number of residues in helical or beta sheet state and  **
** what percent they constitute from the total number of residues in a  **
** protein chain                                                        **
**                                                                      **
** INPUT:   *Chain     Pointer to a protein chain                       **
**                                                                      **
** OUTPUT:  *HelAlp    Number of alpha-helical residues                 **
**          *HelPI     Number of residues in pi-helices                 **
**          *Hel310    Number of residues in 3-10 helices               **
**          *Sheet     Number of residues in beta sheet                 **
**                                                                      **
** RETURNS: Secondary structure content                                 **
**                                                                      **
*************************************************************************/
/*
float SecStrContent(CHAIN *Chain, int *HelAlp, int *HelPI, int *Hel310,
		    int *Sheet, int *Turn) {

  int Res;
  float Content = 0.0;

  *HelAlp = 0; *HelPI = 0; *Hel310 = 0; *Sheet = 0; *Turn = 0;
  
  for( Res=0; Res<Chain->NRes; Res++ ) {
    switch( Chain->Rsd[Res]->Prop->PdbAsn ) {
    case 'H' : (*HelAlp)++;
               break;
    case 'G' : (*Hel310)++;
               break;
    case 'I' : (*HelPI)++;
               break;
    case 'E' : (*Sheet)++;
               break;
    case 'T' : (*Turn)++;
               break;
    }
  }

  Content = ( (float)( (*HelAlp)+(*HelPI)+(*Hel310)+(*Sheet)+(*Turn) ) )/(float)Chain->NRes;

  return(Content);
}
*/
