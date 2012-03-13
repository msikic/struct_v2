#!/usr/bin/perl -w

use strict;
use Switch;
sub needleman_wunsch (@);


(@ARGV >= 3) ||
    die "Usage: $0  <db search struct_out>  <params>  <top number>\n";

my $pymol     = "/home/ivanam/downloads/pymol/pymol.exe";
my $struct    = "/home/ivanam/baubles/struct/struct";
my $dbdir     = "/home/ivanam/databases/struct_dbfiles"; 
my $pdbdir    = "/home/ivanam/databases/pdbfiles";
my $pdbaffine = "/home/ivanam/perlscr/pdb_manip/pdb_affine_tfm.pl";

my ( $db_search, $params_file,$top_number) = @ARGV;
my ($struct_out, $struct_tfm);

my $filename;
my ($ignore, $target, $query);
my ($reading_map, $reading_rotation, $i, $j);
my @rotation;
my ($region, %region);
my $reading;
my %map = ();
my @qry_sse_ids = ();
my ($sse1, $sse2);
my ($tgt_pdb_id, $qry_pdb_id);
my $file;
my $cmd;

my ($kwd, $outname);
my @similarity;
my (@map_i2j, @map_j2i);
my ( $tgt_db, $qry_db);
my ($hit, @top_list);
my ($qry_name, $tgt_name);
my ($tgt_pdb_name, $qry_pdb_name, $tgt_pdb, $qry_pdb);
my ($tgt_chain, $qry_chain);
my $new_params;

############################################
# do we have everything we need?
foreach ( $struct,  $params_file, $dbdir, $pymol,$pdbaffine, "tfms") {
    (-e $_) || die "$_ not found\n";
}


############################################
@top_list = split "\n", `sort -grk 5 $db_search | head -n $top_number`;



(-e "postp_out") && `rm postp_out`;
`touch postp_out`;

############################################
foreach $hit (@top_list ) {

    ($qry_name, $tgt_name) = split " ", $hit;

    # locate dbfiles
    $tgt_db = `ls $dbdir/*/$tgt_name.db`; chomp $tgt_db;
    $qry_db = `ls $dbdir/*/$qry_name.db`; chomp $qry_db;

    # locate pdbfiles    <<<<<<
    $tgt_pdb_name = substr $tgt_name, 0, 4;
    $qry_pdb_name = substr $qry_name, 0, 4;

    $tgt_pdb = "$pdbdir/$tgt_pdb_name.pdb";
    $qry_pdb = "$pdbdir/$qry_pdb_name.pdb";


    foreach ($tgt_db, $qry_db) {
	(-e $_) || die "$_ not found\n";
    }

    # make sure that params file correspods
    # to tgt and $qry
    $tgt_chain = substr $tgt_name, 4, 1;
    $qry_chain = substr $qry_name, 4, 1;

    $new_params = $params_file.".new";
    $cmd = "grep -v pdbf_ $params_file  | grep -v chain_ | ".
	"grep -v postp | grep -v verbose > $new_params";
    `$cmd`;
    `echo postp >> $new_params`;

    `echo pdbf_tgt  $tgt_pdb  >> $new_params`;
    `echo chain_tgt $tgt_chain  >> $new_params`;

    `echo pdbf_qry  $qry_pdb  >> $new_params`;
    `echo chain_qry $qry_chain  >> $new_params`;

    `echo verbose >> $new_params`;
    
 
    ############################################
    # run struct

    $outname = $qry_name;

    ($struct_out, $struct_tfm) = ("$outname.long_out", "$outname.struct_out.for_postp");
    $cmd = "$struct  $tgt_db  $qry_db $new_params > $struct_out ";
    (system $cmd) && die "Error running $cmd.\n";

    `grep -v done $outname.struct_out >>  postp_out`;
    `mv $outname.struct_out.for_postp tfms/$tgt_name\_to_$qry_name.tfm`;
}

############################
# sort the output by the aligned bb score
`sort -grk 9 postp_out > postp_out.sorted`;
