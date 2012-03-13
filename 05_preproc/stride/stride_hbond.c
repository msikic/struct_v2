
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "stride_types.h"
#include "stride_util.h"
#include "stride_cmd.h"
#include "stride_geom.h"

#include "stride_hbond.h"

#define MAXDONOR     MAX_RES
#define MAXACCEPTOR  MAX_RES

#define DIST_N_H 1.0
#define RmGRID                    3.0
#define EmGRID                   -2.8
#define CGRID                    -3.0*EmGRID*pow(RmGRID,8.0)
#define DGRID                    -4.0*EmGRID*pow(RmGRID,6.0)
#define K1GRID                    0.9/pow(cos(RAD(110.0)),6.0)
#define K2GRID                    pow(cos(RAD(110.0)),2.0)

#define MINACCANG_SP2   90.0
#define MAXACCANG_SP2   180.0
#define MINACCANG_SP3   60.0
#define MAXACCANG_SP3   180.0
#define MINDONANG_SP2   90.0
#define MAXDONANG_SP2   180.0
#define MINDONANG_SP3   90.0
#define MAXDONANG_SP3   180.0
#define ACCDONANG       60.0
#define DONACCANG       90.0


static void GRID_Energy(float *CA2, float *C, float *O, float *H, float *N,
			StrideCmd *Cmd, HBOND *HBond);

static int FindDnr(CHAIN *Chain, DONOR **Dnr, int *NDnr, StrideCmd *Cmd);
static int DefineDnr(CHAIN *Chain, DONOR **Dnr, int *dc, int Res,
		     enum HYBRID Hybrid, enum GROUP Group,
		     float HB_Radius, int N);
static int FindAcc(CHAIN *Chain, ACCEPTOR **Acc, int *NAcc, StrideCmd *Cmd);
static int DefineAcceptor(CHAIN *Chain, ACCEPTOR **Acc, int *ac, int Res,
			  enum HYBRID Hybrid, enum GROUP Group,
			  float HB_Radius, int N);
static void PrintHydrBond(char *Text, HBOND *HBond);


int FindPolInt(HBOND **HBond, RESIDUE *Res1, RESIDUE *Res2) {
  register int i, j, hb;
  INVOLVED *p1, *p2;

  p1 = Res1->Inv;
  p2 = Res2->Inv;

  if( p1->NBondDnr && p2->NBondAcc ) {
    for( i=0; i<p1->NBondDnr; i++ ) {
      hb = p1->HBondDnr[i];
      for( j=0; j<p2->NBondAcc; j++ )
	if( hb == p2->HBondAcc[j] && HBond[hb]->ExistPolarInter ) 
	  return(hb);
    }
  }

  return(ERR);
}

int FindBnd(HBOND **HBond, RESIDUE *Res1, RESIDUE *Res2) {
  register int i, j, hb;
  INVOLVED *p1, *p2;

  p1 = Res1->Inv;
  p2 = Res2->Inv;

  if( p1->NBondDnr && p2->NBondAcc ) {
    for( i=0; i<p1->NBondDnr; i++ ) {
      hb = p1->HBondDnr[i];
      for( j=0; j<p2->NBondAcc; j++ )
	if( hb == p2->HBondAcc[j] && HBond[hb]->ExistHydrBondRose ) 
	  return(hb);
    }
  }

  return(ERR);
}

int PlaceHydrogens(CHAIN *Chain) {

  int Res, i, N, C, CA, H, PlacedCnt=0;
  float Length_N_C, Length_N_CA, Length_N_H;
  RESIDUE *r, *rr;
  
  for( Res=1; Res<Chain->NRes; Res++ ) {
    
    r  = Chain->Rsd[Res];
    rr = Chain->Rsd[Res-1];

      if( !strcmp(r->ResType,"PRO") ) continue;

      /* Replace deiterium atoms by hydrogens */
      if( FindAtom(Chain,Res,"D",&H) ) 
	strcpy(r->AtomType[H],"H");

      if( !FindAtom(Chain,Res,"H",&H)  && FindAtom(Chain,Res,"N",&N)   &&
	   FindAtom(Chain,Res-1,"C",&C) && FindAtom(Chain,Res,"CA",&CA) ) {

	H = r->NAtom;
	
	Length_N_C   = Dist(r->Coord[N],rr->Coord[C]);
	Length_N_CA  = Dist(r->Coord[N],r->Coord[CA]);
	
	for( i=0; i<3; i++ )
	  r->Coord[H][i] = r->Coord[N][i] - 
	      ( (rr->Coord[C][i] -  r->Coord[N][i])/Length_N_C +
	        (r->Coord[CA][i]  -  r->Coord[N][i])/Length_N_CA );

	Length_N_H = Dist(r->Coord[N],r->Coord[H]);

	for( i=0; i<3; i++ )
	  r->Coord[H][i] = r->Coord[N][i] + 
	    DIST_N_H*(r->Coord[H][i]-r->Coord[N][i])/Length_N_H;

	strcpy(r->AtomType[H],"H");
	r->NAtom++;
	PlacedCnt++;
      }
    }
  return(PlacedCnt);
}


int NoDoubleHBond(HBOND **HBond, int NHBond)
{

  int i, j, NExcl=0;

  for( i=0; i<NHBond-1; i++ )
    for( j=i+1; j<NHBond; j++ ) 
      if( HBond[i]->Dnr->D_Res == HBond[j]->Dnr->D_Res &&
	  HBond[i]->Dnr->Chain->Id == HBond[j]->Dnr->Chain->Id &&
	  HBond[i]->ExistPolarInter && HBond[j]->ExistPolarInter ) {
	if( HBond[i]->Energy < 5.0*HBond[j]->Energy ) {
	  HBond[j]->ExistPolarInter = NO;
	  NExcl++;
	}
	else 
	if( HBond[j]->Energy < 5.0*HBond[i]->Energy ) {
	  HBond[i]->ExistPolarInter = NO;
	  NExcl++;
	}
      }
  
  return(NExcl);
}



/*************************************************
 Calculate the hydrogen bond energy as defined by
 Boobbyer et al., 1989
**************************************************/

static void GRID_Energy(float *CA2, float *C, float *O, float *H, float *N,
			StrideCmd *Cmd, HBOND *HBond)
{

  float ProjH[3];

 /***** Distance dependence ( 8-6 potential ) ****/

  if( Cmd->Truncate && HBond->AccDonDist < RmGRID ) 
    HBond->AccDonDist = RmGRID;
  HBond->Er = CGRID/pow(HBond->AccDonDist,8.0) - DGRID/pow(HBond->AccDonDist,6.0);

 /************** Angular dependance ****************/

 /* Find projection of the hydrogen on the O-C-CA plane */
  Project4_123(O,C,CA2,H,ProjH); 


 /* Three angles determining the direction of the hydrogen bond */
  HBond->ti = fabs(180.0 - Ang(ProjH,O,C)); 
  HBond->to = Ang(H,O,ProjH);             
  HBond->p  = Ang(N,H,O);

 /* Calculate both angle-dependent HB energy components Et and Ep */ 
  if( HBond->ti >= 0.0 && HBond->ti < 90.0 )   
    HBond->Et = cos(RAD(HBond->to))*(0.9+0.1*sin(RAD(2*HBond->ti)));
  else
  if( HBond->ti >= 90.0 && HBond->ti < 110.0 ) 
    HBond->Et = K1GRID*cos(RAD(HBond->to))*
      (pow((K2GRID-pow(cos(RAD(HBond->ti)),2.0)),3.0));
  else
    HBond->Et = 0.0;

  if( HBond->p > 90.0 && HBond->p < 270.0 )
    HBond->Ep = pow(cos(RAD(HBond->p)),2.0);
  else
    HBond->Ep = 0.0;

    /******** Full hydrogen bond energy *********************/
  HBond->Energy = 1000.0*HBond->Er*HBond->Et*HBond->Ep;
}


static int FindDnr(CHAIN *Chain, DONOR **Dnr, int *NDnr, StrideCmd *Cmd)
{

  int Res, dc;
  char Rsd[RES_FIELD];
  
  dc = *NDnr;

  for( Res=0; Res<Chain->NRes; Res++ ) {

    strcpy(Rsd,Chain->Rsd[Res]->ResType);

    DefineDnr(Chain,Dnr,&dc,Res,Nsp2,Peptide,1.90,0);

    if( !Cmd->SideChainHBond ) continue;

    if( !strcmp(Rsd,"TRP") ) 
      DefineDnr(Chain,Dnr,&dc,Res,Nsp2,Trp,1.90,0); 
    else if( !strcmp(Rsd,"ASN") ) DefineDnr(Chain,Dnr,&dc,Res,Nsp2,Asn,1.90,0); 
    else if( !strcmp(Rsd,"GLN") ) DefineDnr(Chain,Dnr,&dc,Res,Nsp2,Gln,1.90,0); 
    else if( !strcmp(Rsd,"ARG") ) {
      DefineDnr(Chain,Dnr,&dc,Res,Nsp2,Arg,1.90,1);
      DefineDnr(Chain,Dnr,&dc,Res,Nsp2,Arg,1.90,2);
      DefineDnr(Chain,Dnr,&dc,Res,Nsp2,Arg,1.90,3);
    }
    else if( !strcmp(Rsd,"HIS") ) {
      DefineDnr(Chain,Dnr,&dc,Res,Nsp2,His,1.90,1);
      DefineDnr(Chain,Dnr,&dc,Res,Nsp2,His,1.90,2);
    }
    else if( !strcmp(Rsd,"LYS") ) DefineDnr(Chain,Dnr,&dc,Res,Nsp3,Lys,2.10,0); 
    else if( !strcmp(Rsd,"SER") ) DefineDnr(Chain,Dnr,&dc,Res,Osp3,Ser,1.70,0);
    else if( !strcmp(Rsd,"THR") ) DefineDnr(Chain,Dnr,&dc,Res,Osp3,Thr,1.70,0);
    else if( !strcmp(Rsd,"TYR") ) DefineDnr(Chain,Dnr,&dc,Res,Osp2,Tyr,1.70,0); 
  }

  *NDnr = dc;
  return(dc);
}

static int DefineDnr(CHAIN *Chain, DONOR **Dnr, int *dc, int Res,
		     enum HYBRID Hybrid, enum GROUP Group,
		     float HB_Radius, int N)
{

  Dnr[*dc] = (DONOR *)ckalloc(sizeof(DONOR));
  
  Dnr[*dc]->Chain = Chain;
  Dnr[*dc]->D_Res = Res;
  if( Group != Peptide )
    Dnr[*dc]->DD_Res = Res;
  else
    Dnr[*dc]->DD_Res = Res-1;
  Dnr[*dc]->DDI_Res = Res;
  Dnr[*dc]->Hybrid = Hybrid;
  Dnr[*dc]->Group = Group;
  Dnr[*dc]->HB_Radius = HB_Radius;
  
  if( Group == Peptide ) {
    if( Res != 0 ) {
      FindAtom(Chain,Res,"N",&Dnr[*dc]->D_At);
      FindAtom(Chain,Res-1,"C",&Dnr[*dc]->DD_At);
    }
    else {
      Dnr[*dc]->D_At  = ERR;
      Dnr[*dc]->DD_At = ERR;
    }
    FindAtom(Chain,Res,"CA",&Dnr[*dc]->DDI_At);
    FindAtom(Chain,Res,"H",&Dnr[*dc]->H);
  }
  else if( Group == Trp ) {
    FindAtom(Chain,Res,"NE1",&Dnr[*dc]->D_At);
    FindAtom(Chain,Res,"CE2",&Dnr[*dc]->DD_At);
    FindAtom(Chain,Res,"CD1",&Dnr[*dc]->DDI_At);
  }
  else if( Group == Asn ) {
    FindAtom(Chain,Res,"ND1",&Dnr[*dc]->D_At);
    FindAtom(Chain,Res,"CG",&Dnr[*dc]->DD_At);
    FindAtom(Chain,Res,"CB",&Dnr[*dc]->DDI_At);
  }
  else if( Group == Gln ) {
    FindAtom(Chain,Res,"NE2",&Dnr[*dc]->D_At);
    FindAtom(Chain,Res,"CD",&Dnr[*dc]->DD_At);
    FindAtom(Chain,Res,"CG",&Dnr[*dc]->DDI_At);
  }
  else if( Group == Arg ) {
    if( N == 1 ) {
      FindAtom(Chain,Res,"NE",&Dnr[*dc]->D_At);
      FindAtom(Chain,Res,"CZ",&Dnr[*dc]->DD_At);
      FindAtom(Chain,Res,"CD",&Dnr[*dc]->DDI_At);
    }
    else
    if( N == 2 ) {
      FindAtom(Chain,Res,"NH1",&Dnr[*dc]->D_At);
      FindAtom(Chain,Res,"CZ",&Dnr[*dc]->DD_At);
      FindAtom(Chain,Res,"NE",&Dnr[*dc]->DDI_At);
    }
    else
    if( N == 3 ) {
      FindAtom(Chain,Res,"NH2",&Dnr[*dc]->D_At);
      FindAtom(Chain,Res,"CZ",&Dnr[*dc]->DD_At);
      FindAtom(Chain,Res,"NE",&Dnr[*dc]->DDI_At);
    }
  }
  else if( Group == His ) {
      if( N == 1 ) {
	FindAtom(Chain,Res,"ND1",&Dnr[*dc]->D_At);
	FindAtom(Chain,Res,"CG",&Dnr[*dc]->DD_At);
	FindAtom(Chain,Res,"CE1",&Dnr[*dc]->DDI_At);
      }
      else if( N == 2 ) {
	FindAtom(Chain,Res,"NE2",&Dnr[*dc]->D_At);
	FindAtom(Chain,Res,"CE1",&Dnr[*dc]->DD_At);
	FindAtom(Chain,Res,"CD2",&Dnr[*dc]->DDI_At);
      }
    }
  else if( Group == Tyr ) {
    FindAtom(Chain,Res,"OH",&Dnr[*dc]->D_At);
    FindAtom(Chain,Res,"CZ",&Dnr[*dc]->DD_At);
    FindAtom(Chain,Res,"CE1",&Dnr[*dc]->DDI_At);
  }
  else if( Group == Lys ) {
    FindAtom(Chain,Res,"NZ",&Dnr[*dc]->D_At);
    FindAtom(Chain,Res,"CE",&Dnr[*dc]->DD_At);
  }
  else if( Group == Ser ) {
    FindAtom(Chain,Res,"OG",&Dnr[*dc]->D_At);
    FindAtom(Chain,Res,"CB",&Dnr[*dc]->DD_At);
  }
  else if( Group == Thr ) {
    FindAtom(Chain,Res,"OG1",&Dnr[*dc]->D_At);
    FindAtom(Chain,Res,"CB",&Dnr[*dc]->DD_At);
  }      
  
  if( Dnr[*dc]->H == ERR || Dnr[*dc]->D_At   == ERR || Dnr[*dc]->DD_At  == ERR ||
      (Dnr[*dc]->DDI_At == ERR && (Hybrid == Nsp2 || Hybrid == Osp2 )) ) {
    free(Dnr[*dc]); return(FAILURE);
  }
  else (*dc)++;
  return(SUCCESS);
}


static int FindAcc(CHAIN *Chain, ACCEPTOR **Acc, int *NAcc, StrideCmd *Cmd)
{

  int Res, ac;
  char Rsd[RES_FIELD];
  
  ac = *NAcc;

  for( Res=0; Res<Chain->NRes; Res++ ) {
    strcpy(Rsd,Chain->Rsd[Res]->ResType);

    DefineAcceptor(Chain,Acc,&ac,Res,Osp2,Peptide,1.60,0);

    if( !Cmd->SideChainHBond ) continue;

    if( !strcmp(Rsd,"HIS") ) {
      DefineAcceptor(Chain,Acc,&ac,Res,Nsp2,His,1.60,0);
      DefineAcceptor(Chain,Acc,&ac,Res,Nsp2,His,1.60,0);
    }
    else if( !strcmp(Rsd,"SER") ) DefineAcceptor(Chain,Acc,&ac,Res,Osp3,Ser,1.70,0);
    else if( !strcmp(Rsd,"THR") ) DefineAcceptor(Chain,Acc,&ac,Res,Osp3,Thr,1.70,0);
    else if( !strcmp(Rsd,"ASN") ) DefineAcceptor(Chain,Acc,&ac,Res,Osp2,Asn,1.60,0);
    else if( !strcmp(Rsd,"GLN") ) DefineAcceptor(Chain,Acc,&ac,Res,Osp2,Gln,1.60,0);
    else if( !strcmp(Rsd,"ASP") ) {
      DefineAcceptor(Chain,Acc,&ac,Res,Osp2,Asp,1.60,1);
      DefineAcceptor(Chain,Acc,&ac,Res,Osp2,Asp,1.60,2);
    }
    else if( !strcmp(Rsd,"GLU") ) {
      DefineAcceptor(Chain,Acc,&ac,Res,Osp2,Glu,1.60,1);
      DefineAcceptor(Chain,Acc,&ac,Res,Osp2,Glu,1.60,2);
    }
    else if( !strcmp(Rsd,"TYR") ) DefineAcceptor(Chain,Acc,&ac,Res,Osp2,Tyr,1.70,0);
    else if( !strcmp(Rsd,"MET") ) DefineAcceptor(Chain,Acc,&ac,Res,Ssp3,Met,1.95,0);
    else if( !strcmp(Rsd,"CYS") ) DefineAcceptor(Chain,Acc,&ac,Res,Ssp3,Cys,1.70,0);
    }

  *NAcc = ac;
  return(ac);
}


static int DefineAcceptor(CHAIN *Chain, ACCEPTOR **Acc, int *ac, int Res,
			  enum HYBRID Hybrid, enum GROUP Group,
			  float HB_Radius, int N)
{

  Acc[*ac] = (ACCEPTOR *)ckalloc(sizeof(ACCEPTOR));

  Acc[*ac]->Chain = Chain;
  Acc[*ac]->A_Res    = Res;
  Acc[*ac]->AA_Res   = Res;
  Acc[*ac]->AA2_Res   = Res;
  Acc[*ac]->Hybrid    = Hybrid;
  Acc[*ac]->Group     = Group;
  Acc[*ac]->HB_Radius = HB_Radius;
  
  if( Group == Peptide ) {
    if( Res != Chain->NRes-1 ) {
      FindAtom(Chain,Res,"O",&Acc[*ac]->A_At);
      FindAtom(Chain,Res,"C",&Acc[*ac]->AA_At);
    }
    else {
      Acc[*ac]->A_At = ERR;
      Acc[*ac]->AA_At = ERR;
    }
    FindAtom(Chain,Res,"CA",&Acc[*ac]->AA2_At);
  }
  else if( Group == His ) {
    if( N == 1 ) {
      FindAtom(Chain,Res,"ND1",&Acc[*ac]->A_At);
      FindAtom(Chain,Res,"CG",&Acc[*ac]->AA_At);
      FindAtom(Chain,Res,"CE1",&Acc[*ac]->AA2_At);
    }
    else if( N == 2 ) {
      FindAtom(Chain,Res,"NE2",&Acc[*ac]->A_At);
      FindAtom(Chain,Res,"CE1",&Acc[*ac]->AA_At);
      FindAtom(Chain,Res,"CD2",&Acc[*ac]->AA2_At);
    }
  }
  else if( Group == Asn ) {
    FindAtom(Chain,Res,"OD1",&Acc[*ac]->A_At);
    FindAtom(Chain,Res,"CG",&Acc[*ac]->AA_At);
    FindAtom(Chain,Res,"CB",&Acc[*ac]->AA2_At);
  }
  else if( Group == Gln ) {
    FindAtom(Chain,Res,"OE1",&Acc[*ac]->A_At);
    FindAtom(Chain,Res,"CD",&Acc[*ac]->AA_At);
    FindAtom(Chain,Res,"CG",&Acc[*ac]->AA2_At);
  }
  else if( Group == Asp ) {
    if( N == 1 ) {
      FindAtom(Chain,Res,"OD1",&Acc[*ac]->A_At);
      FindAtom(Chain,Res,"CG",&Acc[*ac]->AA_At);
      FindAtom(Chain,Res,"CB",&Acc[*ac]->AA2_At);
    }
    else if( N == 2 ) {
      FindAtom(Chain,Res,"ND2",&Acc[*ac]->A_At);
      FindAtom(Chain,Res,"CG",&Acc[*ac]->AA_At);
      FindAtom(Chain,Res,"CB",&Acc[*ac]->AA2_At);
    }
  }
  else if( Group == Glu ) {
    if( N == 1 ) {
      FindAtom(Chain,Res,"OE1",&Acc[*ac]->A_At);
      FindAtom(Chain,Res,"CD",&Acc[*ac]->AA_At);
      FindAtom(Chain,Res,"CG",&Acc[*ac]->AA2_At);
    }
    else if( N == 2 ) {
      FindAtom(Chain,Res,"NE2",&Acc[*ac]->A_At);
      FindAtom(Chain,Res,"CD",&Acc[*ac]->AA_At);
      FindAtom(Chain,Res,"CG",&Acc[*ac]->AA2_At);
    }
  }
  else if( Group == Tyr ) {
    FindAtom(Chain,Res,"OH",&Acc[*ac]->A_At);
    FindAtom(Chain,Res,"CZ",&Acc[*ac]->AA_At);
    FindAtom(Chain,Res,"CE1",&Acc[*ac]->AA2_At);
  }
  else if( Group == Ser ) {
    FindAtom(Chain,Res,"OG",&Acc[*ac]->A_At);
    FindAtom(Chain,Res,"CB",&Acc[*ac]->AA_At);
  }
  else if( Group == Thr ) {
    FindAtom(Chain,Res,"OG1",&Acc[*ac]->A_At);
    FindAtom(Chain,Res,"CB",&Acc[*ac]->AA_At);
  }
  else if( Group == Met ) {
    FindAtom(Chain,Res,"SD",&Acc[*ac]->A_At);
    FindAtom(Chain,Res,"CG",&Acc[*ac]->AA_At);
  }
  else if( Group == Cys ) {
    FindAtom(Chain,Res,"SG",&Acc[*ac]->A_At);
    FindAtom(Chain,Res,"CB",&Acc[*ac]->AA_At);
  }
  
  if( Acc[*ac]->A_At   == ERR || Acc[*ac]->AA_At  == ERR ||
      (Acc[*ac]->AA2_At == ERR && (Hybrid == Nsp2 || Hybrid == Osp2 )) ) {
    free(Acc[*ac]); return(FAILURE); 
  }
  else (*ac)++;
  return(SUCCESS);
}


const char * FindHydrogenBonds(CHAIN **Chain, int NChain,
			       HBOND **HBond, int * NHBond,
			       StrideCmd *Cmd) {
  DONOR **Dnr;
  ACCEPTOR **Acc;
  BOOLEAN *BondedDonor, *BondedAcceptor;
  int NDnr=0, NAcc=0;
  int dc, ac, ccd, cca, cc, hc=0, i;
  void (*HBOND_Energy)();
  BUFFER Text;

  Dnr = (DONOR **)ckalloc(MAXDONOR*sizeof(DONOR *));
  Acc = (ACCEPTOR **)ckalloc(MAXACCEPTOR*sizeof(ACCEPTOR *));
  
  for( cc=0; cc<NChain; cc++ ) { 
    FindDnr(Chain[cc],Dnr,&NDnr,Cmd); 
    FindAcc(Chain[cc],Acc,&NAcc,Cmd);
  }

  BondedDonor    = (BOOLEAN *)ckalloc(NDnr*sizeof(BOOLEAN));
  BondedAcceptor = (BOOLEAN *)ckalloc(NAcc*sizeof(BOOLEAN));
  for( i=0; i<NDnr; i++ )
    BondedDonor[i] = NO;
  for( i=0; i<NAcc; i++ )
    BondedAcceptor[i] = NO;

  /*
  if( Cmd->EnergyType == 'D' ) 
    HBOND_Energy = DSSP_Energy; 
  else 
    HBOND_Energy = GRID_Energy;
  */
  if( Cmd->EnergyType != 'G' )
    return "FindHydrogenBonds(): EnergyType must be 'G'";
  HBOND_Energy = GRID_Energy;


  for( dc=0; dc<NDnr; dc++ ) {

    if( Dnr[dc]->Group != Peptide && !Cmd->SideChainHBond ) continue;
    
    for( ac=0; ac<NAcc; ac++ ) {
      
      if( abs(Acc[ac]->A_Res - Dnr[dc]->D_Res) < 2 && Acc[ac]->Chain->Id == Dnr[dc]->Chain->Id ) 
	continue;
      
      if( Acc[ac]->Group != Peptide && !Cmd->SideChainHBond ) continue;

      if( hc == MAXHYDRBOND )
	return "FindHydrogenBonds(): Number of hydrogen bonds exceeds current limit (MAXHYDRBOND)";
      HBond[hc] = (HBOND *)ckalloc(sizeof(HBOND));

      HBond[hc]->ExistHydrBondRose = NO;
      HBond[hc]->ExistHydrBondBaker = NO;
      HBond[hc]->ExistPolarInter = NO;
      
      if( (HBond[hc]->AccDonDist = 
	   Dist(Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->D_At],
		Acc[ac]->Chain->Rsd[Acc[ac]->A_Res]->Coord[Acc[ac]->A_At]) ) <=
	 Cmd->DistCutOff ) {
	
	
	if( Cmd->MainChainPolarInt && Dnr[dc]->Group == Peptide && 
	   Acc[ac]->Group == Peptide && Dnr[dc]->H != ERR) {
	  HBOND_Energy(Acc[ac]->Chain->Rsd[Acc[ac]->AA2_Res]->Coord[Acc[ac]->AA2_At],
		       Acc[ac]->Chain->Rsd[Acc[ac]->AA_Res]->Coord[Acc[ac]->AA_At],
		       Acc[ac]->Chain->Rsd[Acc[ac]->A_Res]->Coord[Acc[ac]->A_At],
		       Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->H],
		       Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->D_At],
		       Cmd,HBond[hc]);
	  
	  if( HBond[hc]->Energy < -10.0 &&
	     ( (Cmd->EnergyType == 'G' && fabs(HBond[hc]->Et) > Eps && 
		fabs(HBond[hc]->Ep) > Eps ) || Cmd->EnergyType != 'G' ) )
	    HBond[hc]->ExistPolarInter = YES;
	}

	if( Cmd->MainChainHBond && 
	    (HBond[hc]->OHDist = 
	     Dist(Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->H],
		  Acc[ac]->Chain->Rsd[Acc[ac]->A_Res]->Coord[Acc[ac]->A_At])) <= 2.5 &&
	    (HBond[hc]->AngNHO = 
	     Ang(Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->D_At],
		 Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->H],
		 Acc[ac]->Chain->Rsd[Acc[ac]->A_Res]->Coord[Acc[ac]->A_At])) >= 90.0 &&
	     HBond[hc]->AngNHO <= 180.0 &&
	    (HBond[hc]->AngCOH = 
	     Ang(Acc[ac]->Chain->Rsd[Acc[ac]->AA_Res]->Coord[Acc[ac]->AA_At],
		 Acc[ac]->Chain->Rsd[Acc[ac]->A_Res]->Coord[Acc[ac]->A_At],
		 Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->H])) >= 90.0 &&
	     
	     HBond[hc]->AngCOH <= 180.0 )
	  HBond[hc]->ExistHydrBondBaker = YES;

	if( Cmd->MainChainHBond && 
	   HBond[hc]->AccDonDist <= Dnr[dc]->HB_Radius+Acc[ac]->HB_Radius ) {
	  
	  HBond[hc]->AccAng = 
	    Ang(Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->D_At],
		Acc[ac]->Chain->Rsd[Acc[ac]->A_Res]->Coord[Acc[ac]->A_At],
		Acc[ac]->Chain->Rsd[Acc[ac]->AA_Res]->Coord[Acc[ac]->AA_At]);
	  
	  if( ( ( Acc[ac]->Hybrid == Nsp2 || Acc[ac]->Hybrid == Osp2 ) && 
	       ( HBond[hc]->AccAng >= MINACCANG_SP2 && 
		HBond[hc]->AccAng <= MAXACCANG_SP2 ) ) || 
	     ( ( Acc[ac]->Hybrid == Ssp3 ||  Acc[ac]->Hybrid == Osp3 ) && 
	      ( HBond[hc]->AccAng >= MINACCANG_SP3 && 
	       HBond[hc]->AccAng <= MAXACCANG_SP3 ) ) ) {
	    
	    HBond[hc]->DonAng = 
	      Ang(Acc[ac]->Chain->Rsd[Acc[ac]->A_Res]->Coord[Acc[ac]->A_At],
		  Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->D_At],
		  Dnr[dc]->Chain->Rsd[Dnr[dc]->DD_Res]->Coord[Dnr[dc]->DD_At]);
	    
	    if( ( ( Dnr[dc]->Hybrid == Nsp2 || Dnr[dc]->Hybrid == Osp2 ) && 
		 ( HBond[hc]->DonAng >= MINDONANG_SP2 && 
		  HBond[hc]->DonAng <= MAXDONANG_SP2 ) ) || 
	       ( ( Dnr[dc]->Hybrid == Nsp3 || Dnr[dc]->Hybrid == Osp3 ) && 
		( HBond[hc]->DonAng >= MINDONANG_SP3 &&
		 HBond[hc]->DonAng <= MAXDONANG_SP3 ) ) ) {
	      
	      if( Dnr[dc]->Hybrid == Nsp2 || Dnr[dc]->Hybrid == Osp2 ) {
		HBond[hc]->AccDonAng = 
		  fabs(Torsion(Dnr[dc]->Chain->Rsd[Dnr[dc]->DDI_Res]->Coord[Dnr[dc]->DDI_At],
			       Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->D_At],
			       Dnr[dc]->Chain->Rsd[Dnr[dc]->DD_Res]->Coord[Dnr[dc]->DD_At],
			       Acc[ac]->Chain->Rsd[Acc[ac]->A_Res]->Coord[Acc[ac]->A_At]));
		
		if( HBond[hc]->AccDonAng > 90.0 && HBond[hc]->AccDonAng < 270.0 )
		  HBond[hc]->AccDonAng = fabs(180.0 - HBond[hc]->AccDonAng);
		
	      }
	      
	      if( Acc[ac]->Hybrid == Nsp2 || Acc[ac]->Hybrid == Osp2 ) {
		HBond[hc]->DonAccAng = 
		  fabs(Torsion(Dnr[dc]->Chain->Rsd[Dnr[dc]->D_Res]->Coord[Dnr[dc]->D_At],
			       Acc[ac]->Chain->Rsd[Acc[ac]->A_Res]->Coord[Acc[ac]->A_At],
			       Acc[ac]->Chain->Rsd[Acc[ac]->AA_Res]->Coord[Acc[ac]->AA_At],
			       Acc[ac]->Chain->Rsd[Acc[ac]->AA2_Res]->Coord[Acc[ac]->AA2_At]));
		
		if(HBond[hc]->DonAccAng > 90.0 && HBond[hc]->DonAccAng < 270.0)
		  HBond[hc]->DonAccAng = fabs(180.0 - HBond[hc]->DonAccAng);
		
	      }
	      
	      if( ( Dnr[dc]->Hybrid != Nsp2 && Dnr[dc]->Hybrid != Osp2 && 
		   Acc[ac]->Hybrid != Nsp2 && Acc[ac]->Hybrid != Osp2 ) ||
		 ( Acc[ac]->Hybrid != Nsp2 && Acc[ac]->Hybrid != Osp2 &&
		  ( Dnr[dc]->Hybrid == Nsp2 || Dnr[dc]->Hybrid == Osp2 ) &&
		  HBond[hc]->AccDonAng <= ACCDONANG ) ||
		 ( Dnr[dc]->Hybrid != Nsp2 && Dnr[dc]->Hybrid != Osp2 &&
		  ( Acc[ac]->Hybrid == Nsp2 || Acc[ac]->Hybrid == Osp2 ) &&
		  HBond[hc]->DonAccAng <= DONACCANG ) ||
		 ( ( Dnr[dc]->Hybrid == Nsp2 || Dnr[dc]->Hybrid == Osp2 ) && 
		  ( Acc[ac]->Hybrid == Nsp2 || Acc[ac]->Hybrid == Osp2 ) &&
		  HBond[hc]->AccDonAng <= ACCDONANG &&
		  HBond[hc]->DonAccAng <= DONACCANG ) )
		HBond[hc]->ExistHydrBondRose = YES;
	    }
	  }
	}
	
      }

      if( (HBond[hc]->ExistPolarInter && HBond[hc]->Energy < 0.0) 
	  || HBond[hc]->ExistHydrBondRose || HBond[hc]->ExistHydrBondBaker ) {
	HBond[hc]->Dnr = Dnr[dc];
	HBond[hc]->Acc = Acc[ac];
	BondedDonor[dc] = YES;
	BondedAcceptor[ac] = YES;
	if( (ccd = FindChain(Chain,NChain,Dnr[dc]->Chain->Id)) != ERR ) {
	  if( Chain[ccd]->Rsd[Dnr[dc]->D_Res]->Inv->NBondDnr < MAXRESDNR )
	    Chain[ccd]->Rsd[Dnr[dc]->D_Res]->Inv->
	      HBondDnr[Chain[ccd]->Rsd[Dnr[dc]->D_Res]->Inv->NBondDnr++] = hc;
	  else
	    fprintf(stderr,"Residue %s %s of chain %s%c is involved in %d hydrogen bonds (%d are allowed)\n",
		Chain[ccd]->Rsd[Dnr[dc]->D_Res]->ResType,
		Chain[ccd]->Rsd[Dnr[dc]->D_Res]->PDB_ResNumb,
		Chain[ccd]->File,SpaceToDash(Chain[ccd]->Id),
		Chain[ccd]->Rsd[Dnr[dc]->D_Res]->Inv->NBondDnr,MAXRESDNR-1);
	}
	if( (cca  = FindChain(Chain,NChain,Acc[ac]->Chain->Id)) != ERR ) {
	  if( Chain[cca]->Rsd[Acc[ac]->A_Res]->Inv->NBondAcc < MAXRESACC )
	    Chain[cca]->Rsd[Acc[ac]->A_Res]->Inv->
	      HBondAcc[Chain[cca]->Rsd[Acc[ac]->A_Res]->Inv->NBondAcc++] = hc;
	  else
	    fprintf(stderr,"Residue %s %s of chain %s%c is involved in %d hydrogen bonds (%d are allowed)\n",
		Chain[cca]->Rsd[Acc[ac]->A_Res]->ResType,
		Chain[cca]->Rsd[Acc[ac]->A_Res]->PDB_ResNumb,
		Chain[cca]->File,SpaceToDash(Chain[cca]->Id),
		Chain[cca]->Rsd[Acc[ac]->A_Res]->Inv->NBondAcc,MAXRESACC-1);
	}
	if( ccd != cca && ccd != ERR ) {
	  Chain[ccd]->Rsd[Dnr[dc]->D_Res]->Inv->InterchainHBonds = YES;
	  Chain[cca]->Rsd[Acc[ac]->A_Res]->Inv->InterchainHBonds = YES;
	  if( HBond[hc]->ExistHydrBondRose ) {
	    Chain[0]->NHydrBondInterchain++;
	    Chain[0]->NHydrBondTotal++;
	  }
	}
	else
	if( ccd == cca && ccd != ERR && HBond[hc]->ExistHydrBondRose ) {
	  Chain[ccd]->NHydrBond++;
	  Chain[0]->NHydrBondTotal++;
	}
	hc++;
      }
      else
	free(HBond[hc]);
    }
  }
  if( Cmd->Info )
    for( i=0; i<hc; i++ ) { 
      if( HBond[i]->Energy < 0.0 ) {
	sprintf(Text,"%3d ",i); 
	PrintHydrBond(Text,HBond[i]); 
      }
    } 
    
  for( i=0; i<NDnr; i++ )
    if( !BondedDonor[i] )
      free(Dnr[i]);
  for( i=0; i<NAcc; i++ )
    if( !BondedAcceptor[i] )
      free(Acc[i]);

  if( NDnr )
    free(BondedDonor);
  if( NAcc )
    free(BondedAcceptor);

  *NHBond = hc;
  return 0;
}

static void PrintHydrBond(char *Text, HBOND *HBond)
{

  fprintf(stdout,"HB %s %20s %3s %4s %4d %c <> %3s %4s %4d %c ",Text,
	  HBond->Dnr->Chain->File,
	  HBond->Dnr->Chain->Rsd[HBond->Dnr->D_Res]->ResType,
	  HBond->Dnr->Chain->Rsd[HBond->Dnr->D_Res]->PDB_ResNumb,
	  HBond->Dnr->D_Res,HBond->Dnr->Chain->Id,
	  HBond->Acc->Chain->Rsd[HBond->Acc->A_Res]->ResType,
	  HBond->Acc->Chain->Rsd[HBond->Acc->A_Res]->PDB_ResNumb,
	  HBond->Acc->A_Res,HBond->Acc->Chain->Id);

  fprintf(stdout," %7.1f ",HBond->AccDonDist);
  if( HBond->ExistPolarInter ) 
    fprintf(stdout,"%7.1f ",HBond->Energy);
  else
    fprintf(stdout,"XXXXXXX ");

  if( HBond->ExistHydrBondRose ) 
    fprintf(stdout,"YES ");
  else
    fprintf(stdout,"NO ");

  if( HBond->ExistHydrBondBaker ) 
    fprintf(stdout,"YES\n");
  else
    fprintf(stdout,"NO\n");

}
