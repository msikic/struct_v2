#!/usr/bin/perl -w

use strict;
use Switch;
sub needleman_wunsch (@);


(@ARGV >= 5) ||
    die "Usage: $0 <tgt single chain PDB>  <tgt db> ".
    " <qry  single chain PDB>  <qry db> <params> \n";

my $extract_region     = "/home/ivanam/perlscr/extractions/extr_region_from_pdb.pl";
my $pdb_extract_chain  = "/home/ivanam/perlscr/pdb_manip/pdb_extract_chain.pl";
my $pdbtfm             = "/home/ivanam/baubles/struct/15_postproc/pdbtfm";
my $pdbaffine          = "/home/ivanam/perlscr/pdb_manip/pdb_affine_tfm.pl";
my $pymol              = "/home/ivanam/downloads/pymol/pymol.exe";
my $struct             = "/home/ivanam/baubles/struct/struct";

my $FAR_FAR_AWAY =  -10000;

my ($tgt_1chain_pdb, $tgt_db, $qry_1chain_pdb, $qry_db, $params_file) = @ARGV;
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

############################################
# do we have everything we need?
foreach ( $tgt_1chain_pdb, $tgt_db, $qry_1chain_pdb, $qry_db, "params",
	  $extract_region, $pdb_extract_chain,  $struct,  $params_file,
	  $pdbtfm, $pdbaffine, $pymol) {
    (-e $_) || die "$_ not found\n";
}

############################################
# check that params file has potprocessing and verbose options
# in case they are commented - just make a new one
#`grep -v postp $params_file | grep -v verbose | grep -v outn > $params_file.new`;
`grep -v postp $params_file | grep -v verbose | grep -v header > $params_file.new`;
foreach my $option ( "postp", "verbose") {
    `echo $option >> $params_file.new `;
}

#`mv $params_file.tmp  $params_file`;
$params_file .= ".new";

############################################
# run struct
# note: here we are assuming the disabled 
# outn option in params file (see the grepping above)
($ignore, $qry_name) = split " ", `grep name $qry_db`;
chomp $qry_name;

($kwd, $outname) = split " ", `grep outname $params_file`;
chomp $outname;


($struct_out, $struct_tfm) = ("$outname.long_out", "$outname.for_postp");
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

############################################
# do we have the pdbs? 
# since the user may be uploading the structure(s)
# this is not to be decidede here;
# here, for the reference, is how one
# would extract the needed chain
#
#$tgt_pdb   = (substr $tgt_pdb_id, 0, 4).".pdb";
#$tgt_chain =  substr $tgt_pdb_id, 4, 1;
#
#$qry_pdb   = (substr $qry_pdb_id, 0, 4).".pdb";
#$qry_chain =  substr $qry_pdb_id, 4, 1;
#
#foreach ( $qry_pdb, $tgt_pdb) {
#    (-e $_) || die "$_ not found\n";
#}
#
#$tgt_1chain_pdb = $tgt_pdb_id.".pdb";
#$qry_1chain_pdb = $qry_pdb_id.".pdb";
#
#if ( ! -e $tgt_1chain_pdb) {
#    $cmd = "$pdb_extract_chain $tgt_pdb $tgt_chain > $tgt_1chain_pdb ";
#    (system $cmd) && die "Error running $cmd.\n";
#}
#if ( !-e $qry_1chain_pdb) {
#    $cmd = "$pdb_extract_chain $qry_pdb $qry_chain > $qry_1chain_pdb ";
#    (system $cmd) && die "Error running $cmd.\n";
#}


############################################
# extract CAs
my (%CAs);
my ($ret, $key);
my ($num, $sse_id, $from, $to, $pdb);
my $line;

%CAs = ();


while ( ($key, $region) = each (%region) ) {

    ($num, $sse_id) = split "_", $key;

    if ( $num == 1 ) {
	$pdb  = $qry_1chain_pdb
    } else {
	$pdb  = $tgt_1chain_pdb;
    }

    ($from, $to) = split " ", $region;

    $cmd = "$extract_region  $pdb  $from  $to ".
	" | awk \'\$3==\"CA\"\'";
    $ret = `$cmd`;
    $CAs{$key} =  $ret;
}


############################################
# read in the 
# preliminary translation from struct
my %transl =  ();
my @aux;
my %cm = ();
my %origin = ();
my $cm_ctr = 0;

$filename = $struct_tfm;
open (IF, "<$filename") ||  die "Cno $filename: $!.\n";
$reading = "";

while ( <IF> ) {
    if ( /^(\w+)/) {
	$reading = $1;
	if ( $reading eq "origin1") {
	    chomp;
	    @aux = split " ";
	    shift @aux;
	    @{$origin{$tgt_pdb_id}} = @aux;
	} elsif ( $reading eq "origin2") {
	    chomp;
	    @aux = split " ";
	    shift @aux;
	    @{$origin{$qry_pdb_id}} = @aux;
	} 
    } elsif ( $reading eq "translations") {
	chomp;
	@aux = split " ";
	$sse1 = $aux[4];
	$sse2 = $aux[0];
	@{$transl{"1_$sse1"}} = @aux[5..7];
	@{$transl{"2_$sse2"}} = @aux[1..3];

	( $map{$sse1} eq $sse2 )  ||
	    die "mismatch btw $struct_out and $struct_tfm\n";

    } elsif ( $reading eq "cms") {
	chomp;
	@aux = split " ";
	$sse1 = $aux[4];
	$sse2 = $aux[0];
	@{$cm{"1_$sse1"}} = @aux[5..7];
	@{$cm{"2_$sse2"}} = @aux[1..3];

	( $map{$sse1} eq $sse2 )  ||
	    die "mismatch btw $struct_out and $struct_tfm\n";

    } elsif ( $reading eq "rotation") {
	# use as a check
    }
}
close IF;

for $i (0 .. 2) {
    printf "%8.3f %8.3f %8.3f\n", @{$rotation[$i]};
}

############################################
# construct two sets of vectors
# qry has to be translated
# and tgt has to be translated and rotated
my @qry_coords;
my @tgt_coords;
my (@r, @t, @T, @r_aux);
my (@qry_block_size, @tgt_block_size);
my ($qry_block_ctr, $tgt_block_ctr) = (0, 0);
my ($qry_ctr, $tgt_ctr, $new_line); 
my ($number_of_qry_vecs, $number_of_tgt_vecs);
my (@qry_orig_coords,  @tgt_orig_coords);
$qry_ctr = 0;
$tgt_ctr = 0;
printf "translating $qry_pdb_id by  %8.3f %8.3f %8.3f\n",  @{$origin{$qry_pdb_id}};
printf "translating $tgt_pdb_id by  %8.3f %8.3f %8.3f\n",  @{$origin{$tgt_pdb_id}};

foreach $sse1 ( @qry_sse_ids ) {

    #############################
    # query vectors
    # extract
    $key = "1_$sse1";
    # translation which places tthe cm at the origin
    @T = @{$cm{$key}};
    $qry_block_size[$qry_block_ctr] = 0;
    foreach ( split "\n", $CAs{$key} ) {
	# store the original coordinates
	$qry_orig_coords[$qry_ctr] = $_;

	# extract the  vectors
	$r[0] = substr $_, 30, 8;  $r[0] =~ s/\s//g;
	$r[1] = substr $_, 38, 8;  $r[1] =~ s/\s//g;
	$r[2] = substr $_, 46, 8;  $r[2] =~ s/\s//g;

	# translate to origin
	for $i (0 .. 2) {
	    $r[$i] -= $T[$i]; 
	}
	# translate out
	@t = @{$transl{$key}};
	for $i (0 .. 2) {
	    $r[$i] += $t[$i];
	}
	@{$qry_coords[$qry_ctr]} = @r;
	$qry_ctr++;
	$qry_block_size[$qry_block_ctr] ++;
    }
    $qry_block_ctr++;

    #############################
    # target vectors
    $sse2 = $map{$sse1};
    # extract
    $key = "2_$sse2";
    # for each SSE: translation which places the cm  at the origin
    # cm - t = 0
    @T = @{$cm{$key}};
    
    $tgt_block_size[$tgt_block_ctr] = 0.0;
    foreach ( split "\n", $CAs{$key} ) {
	# store the original coordinates
	$tgt_orig_coords[$tgt_ctr] = $_;

	# extract the  vectors
	$r[0] = substr $_, 30, 8;  $r[0] =~ s/\s//g;
	$r[1] = substr $_, 38, 8;  $r[1] =~ s/\s//g;
	$r[2] = substr $_, 46, 8;  $r[2] =~ s/\s//g;

	# translate to cm
	for $i (0 .. 2) {
	    $r[$i] -= $T[$i]; 
	}
	# rotate 
	for $i (0 .. 2) {
	    $r_aux[$i] = 0;
	    for $j (0 .. 2) {
		$r_aux[$i] += $rotation[$i][$j]*$r[$j];
	    }
	}
	@r = @r_aux;

	# rotate the translation vector
	for $i (0 .. 2) {
	    $t[$i] = 0;
	    for $j (0 .. 2) {
		$t[$i] += $rotation[$i][$j]*$transl{$key}[$j];
	    }
	}

	# translate out
	for $i (0 .. 2) {
	    $r[$i] += $t[$i];
	}

	@{$tgt_coords[$tgt_ctr]} = @r;
	$tgt_ctr++;
	$tgt_block_size[$tgt_block_ctr] ++;
    }
    $tgt_block_ctr++;
}
($number_of_qry_vecs, $number_of_tgt_vecs) = ($qry_ctr, $tgt_ctr);


###########################################
# find the mapping between individual CAs
# implied by these tfmtions:
# construct similarity (proximity) matrix
my @similarity = ();
my @distance = ();
my ($dim1, $dim2) = ($number_of_qry_vecs, $number_of_tgt_vecs);
my $FAR = 1000;
my (@rq, @rt);

foreach $qry_ctr ( 0 .. $dim1-1) {
    foreach $tgt_ctr ( 0 .. $dim2-1) {
	$similarity[$qry_ctr][$tgt_ctr] = -1.0;
    }
}

my $block_ctr = 0;
my $qry_block_end = $qry_block_size[0]-1;
my ($tgt_block_begin, $tgt_block_end) = (0, $tgt_block_size[0]-1);
my (@map_i2j, @map_j2i);
my $avg_sim = 0;
my $aux_ctr = 0;
my $max = -100;
my $min = 1000;
foreach $qry_ctr ( 0 .. $dim1-1) {

    if ( $qry_ctr > $qry_block_end) {
	$block_ctr  ++;
	$qry_block_end   += $qry_block_size[$block_ctr];
	$tgt_block_begin  = $tgt_block_end + 1;
	$tgt_block_end   += $tgt_block_size[$block_ctr];
    }
    @rq = @{$qry_coords[$qry_ctr]};


    ($tgt_block_end < $dim2 ) || die "Error calculating block bounds.\n";
    foreach $tgt_ctr ( $tgt_block_begin .. $tgt_block_end) {
	@rt= @{$tgt_coords[$tgt_ctr]};
	my $d = 0.0;
	for $i (0 ..2) {
	    my $aux = $rt[$i] - $rq[$i];
	    $d += $aux*$aux;
	}
	$d = sqrt ($d);
	$distance[$qry_ctr][$tgt_ctr] = $d;

	$similarity[$qry_ctr][$tgt_ctr] = exp(-$d/5);

	$avg_sim += $similarity[$qry_ctr][$tgt_ctr];
	($max < $similarity[$qry_ctr][$tgt_ctr]) && ($max = $similarity[$qry_ctr][$tgt_ctr]);
	($min > $similarity[$qry_ctr][$tgt_ctr]) && ($min = $similarity[$qry_ctr][$tgt_ctr]);
	$aux_ctr++;
   }
}

$avg_sim /= $aux_ctr;

#printf "$qry_ctr x $tgt_ctr ... map size: %5d    ", $aux_ctr;
#printf "sim:  min %15.8f   max %15.8f   avg  %15.8f\n", $min, $max,  $avg_sim;


my $retval = needleman_wunsch ($dim1, $dim2, 1);
( $retval  == $FAR_FAR_AWAY ) && exit 10;

######################################################################
# output - to be used by the c routine doing the linear algebra
# I want the original coordinates here, to get finally the
# transformation which gets me from one original pdb to the other
$file = "$qry_pdb_id.out.pdb";
open ( QRY_OF, ">$file") || die "Cno $file: $!.\n";
$file = "$tgt_pdb_id.out.pdb";
open ( TGT_OF, ">$file") || die "Cno $file: $!.\n";

=pod
foreach $qry_ctr ( 0 .. $dim1-1) {
    $tgt_ctr = $map_i2j[$qry_ctr];
    printf "  $qry_ctr    $tgt_ctr\n";
}
exit;
=cut

foreach $qry_ctr ( 0 .. $dim1-1) {
    $tgt_ctr = $map_i2j[$qry_ctr];
    next if ($tgt_ctr<0);
    next if  $distance[$qry_ctr][$tgt_ctr] > 6;
    print QRY_OF $qry_orig_coords[$qry_ctr]."\n";
    print TGT_OF $tgt_orig_coords[$tgt_ctr]."\n";
}

close QRY_OF;
close TGT_OF;



######################################################################
# find the transformation mapping this particular set
# of atoms one onto another
my $tfmfile = "$tgt_pdb_id\_onto_$qry_pdb_id.tfm";
$cmd = "$pdbtfm $tgt_pdb_id.out.pdb $qry_pdb_id.out.pdb ".
    " > $tfmfile  ";
(system $cmd) && die "Error running $cmd.\n";

######################################################################
# transform the original pdb file
my $tfmd_pdb = $tgt_1chain_pdb;
$tfmd_pdb =~ s/\.pdb/\.tfmd.pdb/;
$cmd = "$pdbaffine $tgt_1chain_pdb $tfmfile > $tfmd_pdb";
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




######################################################################
######################################################################
######################################################################
sub needleman_wunsch (@) {

    my  ($max_i, $max_j, $use_endgap) = @_;

    my @F;  #alignment_scoring table
    my @direction;
    my $gap_opening   = -0.2;
    my $gap_extension = -0.1;
    #my $gap_opening   =  -1.0;
    #my $gap_extension =  -0.2;
    my $endgap        =  -0.0;
    my ($i_sim, $j_sim, $diag_sim, $max_sim) = (0.0, 0.0, 0.0, 0.0) ;
    my $score = 0.0;
    my ($i,$j, $matched);
    my $penalty = 0.0;
    # inititalize F and direction
    foreach $i ( 0 .. $max_i) {
	foreach $j ( 0 .. $max_j) {
	    $F[$i][$j] = 0.0;
	    $direction[$i][$j] = 0;
	}
    }
    #  fill the table 
    foreach $i ( 0 .. $max_i) {
	foreach $j ( 0 .. $max_j) {

	    if ( !$i && !$j ) { # upper left corner 
		$F[0][0] = 0;
		$direction[$i][$j] = 'd';
		next;
	    }
	    
	    if ( $i && $j ){ 
		if ( $direction[$i-1][$j] eq 'i' ) {
		    #  gap extension  */
		    $penalty = ($use_endgap&&$j==$max_j) ? $endgap : $gap_extension;		    
		} else {
		    #  gap opening  */
		    $penalty = ($use_endgap&&$j==$max_j) ? $endgap : $gap_opening;
		}
		$i_sim =  $F[$i-1][$j] + $penalty;

		
		if ( $direction[$i][$j-1] eq 'j' ) {
		    $penalty = ($use_endgap && $i==$max_i) ? $endgap : $gap_extension;		    
		} else {
		    $penalty = ($use_endgap && $i==$max_i) ? $endgap : $gap_opening;		    
		}
		$j_sim = $F[$i][$j-1] +  $penalty;
		
		
		$diag_sim =  $F[$i-1][$j-1] + $similarity [$i-1][$j-1] ;
		
	    } elsif ( $j ) {
		
		if ( $use_endgap) {
		    $penalty = $endgap;
		} else {
		    if ( $direction[$i][$j-1] eq 'j' ) {
			$penalty =  $gap_extension;
		    } else {
			$penalty =  $gap_opening;
		    }
		}
		$j_sim = $F[$i][$j-1] + $penalty;
		
		$i_sim = $diag_sim =  $FAR_FAR_AWAY;

	    } elsif ( $i ) {
		if ( $use_endgap) {
		    $penalty = $endgap;
		} else {
		    if ( $direction[$i-1][$j] eq 'i' ) {
			$penalty =  $gap_extension;
		    } else {
		        $penalty =  $gap_opening;
		    }
		}
		$i_sim = $F[$i-1][$j] + $penalty;
		
		$j_sim = $diag_sim =  $FAR_FAR_AWAY;
		
	    } 

	    $max_sim = $diag_sim;
	    $direction[$i][$j] = 'd';
	    if ( $i_sim > $max_sim ){
		$max_sim = $i_sim;
		$direction[$i][$j] = 'i';
	    }
	    if ( $j_sim > $max_sim ) {
		$max_sim = $j_sim;
		$direction[$i][$j] = 'j';
	    }

	    $F[$i][$j] = $max_sim;
	    
	}
    }

    $score = $F[$max_i][$max_j];

    #retrace*/
    $i = $max_i;
    $j = $max_j;
    $matched = 0;
    while ( $i>0 ||  $j >0 ) {
	#printf (" %4d  %4d  %8.3f  \n", i, j, $F[i][j]);
	switch ( $direction[$i][$j] ) {
	    case 'd' {
		#printf ( " %4d  %4d \n",  i, j);	    
		$map_i2j [$i-1] = $j-1;
		$map_j2i [$j-1] = $i-1;
		$i--;
		$j--; 
		$matched++;
	    }
	    case 'i' {
		#printf ( " %4d  %4d \n",  i, -1);
		$map_i2j [$i-1] = $FAR_FAR_AWAY;
		$i--; 
	    } 
	    case 'j' {
		#printf ( " %4d  %4d \n",  -1, j);
		$map_j2i [$j-1] = $FAR_FAR_AWAY;
		$j--; 
	    } 
	    else {
		die "Retracing error.\n";
	    } 
	}
    }

    $matched ||  return $FAR_FAR_AWAY;
    return $score; 
   
    
}

# find the translation which best matches
# the extracted CA


# apply to the whole chain


# produce pymol session

=pod
###########################################
foreach $sse1 ( @qry_sse_ids ) {
    $sse2 = $map{$sse1};
    print "\n********************************\n";
    printf " %2d (%s) maps to %2d  (%s)\n", 
    $sse1, $region{"1_$sse1"}, $sse2, $region{"2_$sse2"};

    print  $CAs{"1_$sse1"},"\n";   
    print  "***\n";
    print  $CAs{"2_$sse2"},"\n";

}
=cut
=pod
######################################################
# output the CAs - for checking purposes

$file = "$qry_pdb_id.in.pdb";
open (QRY_OF, ">$file") || die "Cno $file: $!.\n";

$file = "$tgt_pdb_id.in.pdb";
open (TGT_OF, ">$file") || die "Cno $file: $!.\n";

foreach $sse1 ( @qry_sse_ids ) {

    $key = "1_$sse1";
    print QRY_OF $CAs{$key};
    $sse2 = $map{$sse1};
    $key = "2_$sse2";
    print TGT_OF $CAs{$key};
}
 
close QRY_OF;
close TGT_OF;
=cut


=pod
######################################################
# visualization - replace coords in CA blocks

$file = "$qry_pdb_id.out.pdb";
open ( QRY_OF, ">$file") || die "Cno $file: $!.\n";
$file = "$tgt_pdb_id.out.pdb";
open ( TGT_OF, ">$file") || die "Cno $file: $!.\n";

$qry_ctr = 0;
$tgt_ctr = 0;

foreach $sse1 ( @qry_sse_ids ) {

    $key = "1_$sse1";
    foreach  $new_line ( split "\n", $CAs{$key} ) {
	@r = @{$qry_coords[$qry_ctr]};
	for $i ( 0 .. 2 ) {
	    my $tmp = sprintf "%8.3f",  $r[$i];
	    (substr $new_line, 30+$i*8, 8) = $tmp;
	}
	$qry_ctr ++;
	print QRY_OF "$new_line\n";
    }

    $sse2 = $map{$sse1};
    $key = "2_$sse2";
    foreach  $new_line ( split "\n", $CAs{$key} ) {
	@r = @{$tgt_coords[$tgt_ctr]};
	for $i ( 0 .. 2 ) {
	    my $tmp = sprintf "%8.3f",  $r[$i];
	    (substr $new_line, 30+$i*8, 8) = $tmp;
	}
	$tgt_ctr ++;
	print TGT_OF "$new_line\n";
    }
}

close QRY_OF;
close TGT_OF;
=cut
