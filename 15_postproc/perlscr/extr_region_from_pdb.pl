#! /usr/bin/perl -w
use IO::Handle;         #autoflush
# FH -> autoflush(1);

 
(@ARGV > 2  && !($#ARGV%2)  )  ||
    die "Usage: extr_region_from_pdb.pl <pdb_file> [-c <chain>]  <count from>  <count to> ".
" [<count from> <count to> ...].\n"; 

$pdbfile  = shift  @ARGV;
$chain    = "";

if ( $ARGV[0] eq "-c" ) {
    shift  @ARGV;
    $chain = shift  @ARGV;
}
@from_to  =  @ARGV;

open ( PDB, "<$pdbfile" ) ||
    die "Cno $pdbfile: $! .\n";

while ( <PDB> ) {
    last if ( /^ENDMDL/);
    next if ( ! (/^ATOM/ || /^HETATM/) ) ;

    $res_seq  = substr $_, 22, 4;  $res_seq=~ s/\s//g;
    $chain_id = substr $_, 21, 1;
    next if ( $chain && ($chain ne $chain_id) );

    for ( $ctr=0; $ctr < $#from_to; $ctr+=2 ) {
	if ( $from_to[$ctr] <= $res_seq  && 
	     $res_seq <= $from_to[$ctr+1] ) {
	    print;
	    last;
	}
    }
     
}

close PDB;
