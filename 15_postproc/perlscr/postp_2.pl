#!/usr/bin/perl -w

use strict;
use Switch;
sub needleman_wunsch (@);


(@ARGV >= 3) ||
    die "Usage: $0  <db search struct_out>  <params>  <top number>\n";

my $pymol              = "/home/ivanam/downloads/pymol/pymol.exe";
my $struct             = "/home/ivanam/baubles/struct/struct";


my ( $db_search, $params_file,$top_number) = @ARGV;
my ($struct_out, $struct_tfm);
my ($tgt_pdb, $tgt_chain, $qry_pdb, $qry_chain);
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
my $qry_name;
my ($kwd, $outname);
my ($tgt_1chain_pdb, $qry_1chain_pdb);
my @similarity;
my (@map_i2j, @map_j2i);
my ( $tgt_db, $qry_db);


############################################
# do we have everything we need?
foreach ( $struct,  $params_file,
	  $pdbaffine, $pymol) {
    (-e $_) || die "$_ not found\n";
}


@top_list = split "\n", `sort -grk 5 $db_search | head -n $top_number`;


############################################
# run struct
# note: here we are assuming the disabled 
# outn option in params file (see the grepping above)
($ignore, $qry_name) = split " ", `grep name $qry_db`;
chomp $qry_name;

$outname = $qry_name;

($struct_out, $struct_tfm) = ("$outname.long_out", "$outname.struct_out.for_postp");
$cmd = "$struct  $tgt_db  $qry_db $params_file > $struct_out ";
(system $cmd) && die "Error running $cmd.\n";

############################################
# read in the names and the mapped structures
# from the long output of struct
$filename = $struct_out;
open (IF, "<$filename") ||  die "Cno $filename: $!.\n";
($reading_map, $reading_rotation, $i) = (0,0,0);
%region = ();
while ( <IF>) {

    print;
 
    if (/^[^\w\d\s\.]{2}\s/) {
	($ignore, $tgt_pdb_id, $qry_pdb_id) = split " ";
	print "$tgt_pdb_id, $qry_pdb_id\n"; 

   } elsif  (/^map/ ) {
	$reading_map = 1;

    } elsif ( /^total/ ) {
	$reading_map = 0;

    } elsif ( /^rotation/ ) {
	$reading_rotation = 1;

    } elsif ( $reading_map  ) {

	/(\d+)\s*\-\-\>\s*(\d+)/;
	($sse2, $sse1) = ($1, $2); 
	push @qry_sse_ids, $sse1;
	$map{$sse1} = $sse2;

	/\s+(\d+)\s*\-\s*(\d+)\s+(\d+)\s*\-\s*(\d+)/;
	$region{"2_$sse2"} = "$1 $2";
	$region{"1_$sse1"} = "$3 $4";
	
    } elsif ( $reading_rotation  ) {
	chomp;
	@{$rotation[$i]} = split " ", $_;
	$i++;
	last if ($i == 3)
    }

}




######################################################################
# transform the original pdb file
my $tfmd_pdb = $tgt_1chain_pdb;
$tfmd_pdb =~ s/\.pdb/\.tfmd.pdb/;
$cmd = "$pdbaffine $tgt_1chain_pdb $struct_tfm > $tfmd_pdb";
;
(system $cmd) && die "Error running $cmd.\n";

######################################################################
# create and run pymol script --> create a pymol session
my $pml_name    = "$tgt_pdb_id\_$qry_pdb_id.pml";
my $pml_pse     = "$tgt_pdb_id\_$qry_pdb_id.pse";
my $pml_string  = "";
my $matched_set = "";

$pml_string .= "load $tfmd_pdb, $tgt_pdb_id\n";
$pml_string .= "load $qry_1chain_pdb, $qry_pdb_id\n";
$pml_string .= "hide all\n";
$pml_string .= "bg_color white\n";


$pml_string .= "create $qry_pdb_id\_outline, $qry_pdb_id\n";
$pml_string .= "create $tgt_pdb_id\_outline, $tgt_pdb_id\n";

$matched_set = "";
foreach $sse1 ( @qry_sse_ids ){
    my $region_name   = "qry_$sse1";
    my $from_to = $region{"1_$sse1"};
    $from_to =~ s/\s/-/;
    $pml_string .= "select $region_name, resi $from_to and $qry_pdb_id\n";
    if ( ! $matched_set ) {
	$matched_set = $region_name;
    } else {
	$matched_set = "matched_set or $region_name";
    }
    $pml_string .= "select matched_set, $matched_set\n";

    $sse2 = $map{$sse1};
    $region_name   = "tgt_$sse2";
    $from_to = $region{"2_$sse2"};
    $from_to =~ s/\s/-/;
    $pml_string .= "select $region_name, resi $from_to and $tgt_pdb_id\n";
    $matched_set = "matched_set or $region_name";
    $pml_string .= "select matched_set, $matched_set\n";

}
$pml_string .= "select matched_compl, ($qry_pdb_id or $tgt_pdb_id) and not matched_set\n";
$pml_string .= "remove matched_compl\n";
$pml_string .= "delete matched_compl\n";

$pml_string .= "color tv_blue, $qry_pdb_id\n";
$pml_string .= "color raspberry, $tgt_pdb_id\n";

$pml_string .= "color bluewhite, $qry_pdb_id\_outline\n";
$pml_string .= "color lightpink, $tgt_pdb_id\_outline\n";

$pml_string .= "set cartoon_transparency, 0.5, $qry_pdb_id\_outline\n";
$pml_string .= "set cartoon_transparency, 0.5, $tgt_pdb_id\_outline\n";

$pml_string .= "show cartoon\n";
$pml_string .= "deselect\n";

$pml_string .= "save $pml_pse\n";

$file = $pml_name;
open ( POF, ">$file") || die "Cno $file: $!.\n";
printf POF $pml_string;
close POF;

$cmd = "$pymol -qc -u $pml_name > /dev/null";
(system $cmd) && die "Error running $cmd.\n";



