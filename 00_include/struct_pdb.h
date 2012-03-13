#ifndef PDB_H
# define	PDB_H



/* field start and length in atom record: */
    
/*  1 -  6        Record name     "ATOM  " */
# define        PDB_ATOM_REC_NAME          0
# define        PDB_ATOM_REC_NAME_LEN      6
/*     7 - 11        Integer         serial        Atom serial number. */
#define         PDB_ATOM_ATOM_NO           6       
#define         PDB_ATOM_ATOM_NO_LEN       5
  
/* 13 - 16        Atom            name          Atom name. */
#define         PDB_ATOM_ATOM_NAME         12
#define         PDB_ATOM_ATOM_NAME_LEN     4     
    

/* 17             Character       altLoc        Alternate location indicator. */
#define         PDB_ATOM_ALTLOC            16  
#define         PDB_ATOM_ALTLOC_LEN        1  

/* 18 - 20        Residue name    resName       Residue name. */
# define        PDB_ATOM_RES_NAME          17
# define        PDB_ATOM_RES_NAME_LEN      3

/* 22             Character       chainID       Chain identifier. */
#define        PDB_ATOM_CHAINID            21
#define        PDB_ATOM_CHAINID_LEN         1

/* 23 - 26        Integer         resSeq        Residue sequence number */
#define        PDB_ATOM_RES_NO             22
#define        PDB_ATOM_RES_NO_LEN         4

/* 27             AChar           iCode         Code for insertion of residues. */
#define        PDB_ATOM_INS_CODE           26
#define        PDB_ATOM_INS_CODE_LEN       1
    
/* 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in */
/*                                              Angstroms. */
#define        PDB_ATOM_X                  30
#define        PDB_ATOM_X_LEN               8
    
/* 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in */
/*                                              Angstroms. */
#define        PDB_ATOM_Y                  38
#define        PDB_ATOM_Y_LEN               8
    

/* 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in */
/*                                              Angstroms. */
#define        PDB_ATOM_Z                  46
#define        PDB_ATOM_Z_LEN               8
    


/* 55 - 60        Real(6.2)       occupancy     Occupancy. */
#define        PDB_ATOM_OCC                54
#define        PDB_ATOM_OCC_LEN             6
    

/* 61 - 66        Real(6.2)       tempFactor    Temperature factor. */
#define        PDB_ATOM_TEMP              60
#define        PDB_ATOM_TEMP_LEN           6  
    

/* 73 - 76        LString(4)      segID         Segment identifier, left-justified. */
#define        PDB_ATOM_SEGID             72
#define        PDB_ATOM_SEGID_LEN          4
    

/* 77 - 78        LString(2)      element       Element symbol, right-justified. */
#define        PDB_ATOM_ELMT              76
#define        PDB_ATOM_ELMT_LEN           2
     

/* 79 - 80        LString(2)      charge        Charge on the atom. */
#define        PDB_ATOM_CHARGE            78
#define        PDB_ATOM__CHARGE_LEN        2  
    

/*22 - 25       Integer          initSeqNum   Sequence number of the initial
                                            residue.
26            AChar            initICode    Insertion code of the initial residue.*/
# define       PDB_HELIX_BEGIN_ID  21
# define       PDB_HELIX_BEGIN_LEN  5

/*20            Character        initChainID  Chain identifier for the chain */
# define       PDB_HELIX_BEGIN_CHAIN 19
# define       PDB_HELIX_BEGIN_CHAIN_LEN 1


/*34 - 37       Integer          endSeqNum    Sequence number of the terminal
                                            residue.
38            AChar            endICode     Insertion code of the terminal
residue.*/

# define       PDB_HELIX_END_ID  33
# define       PDB_HELIX_END_LEN  4

/*32            Character        endChainID   Chain identifier for the chain */
# define       PDB_HELIX_END_CHAIN 31
# define       PDB_HELIX_END_CHAIN_LEN 1

/*72 - 76       Integer          length       Length of this helix.*/
# define       PDB_HELIX_LENGTH 71
# define       PDB_HELIX_LENGTH_LEN 4

/*********************************/
/*12 - 14     LString(3)       sheetID      Sheet identifier.  */

# define      PDB_SHEET_SHEET_ID   11
# define      PDB_SHEET_SHEET_LEN   3

/*22          Character        initChainID  Chain identifier of initial 
                                          residue in strand..*/
# define       PDB_SHEET_BEGIN_CHAIN 21
# define       PDB_SHEET_BEGIN_CHAIN_LEN 1

/*23 - 26     Integer          initSeqNum   Sequence number of initial 
                                          residue in strand..*/
//SHEET    2   A 6 TYR A 106  GLU A 111 -1  N  VAL A 110   O  ALA A 301
//0123456789012345678901234567890123456789
# define       PDB_SHEET_BEGIN_ID  22
# define       PDB_SHEET_BEGIN_LEN  4

/*34 - 37     Integer          endSeqNum    Sequence number of terminal
                                          residue.*/
# define       PDB_SHEET_END_ID  33
# define       PDB_SHEET_END_LEN  4


/*  39 - 40     Integer          sense        Sense of strand with respect to
                                          previous strand in the sheet. 0
                                          if first strand, 1 if parallel,
                                          -1 if anti-parallel. */
# define      PDB_SHEET_SENSE  39
# define      PDB_SHEET_SENSE_LEN  2

#endif /* PDB_H */
