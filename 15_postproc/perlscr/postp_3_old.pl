#!/usr/bin/perl -w

use strict;
use Switch;
sub needleman_wunsch (@);


(@ARGV >= 3) ||
    die "Usage: $0  <db search struct_out>  <params>  <top number>\n";

my $pymol     = "/home/ivanam/downloads/pymol/pymol.exe";
my $struct    = "/home/ivanam/kode/03_struct/struct";
my $dbdir     = "cathdbdir"; 
my $pdbdir    = "cathpdb";
my $pdbaffine = "/home/ivanam/perlscr/pdb_manip/pdb_affine_tfm.pl";
my $pdb_extr_chain = "/home/ivanam/perlscr/pdb_manip/pdb_extract_chain.pl";

foreach ($pymol, $struct, $dbdir, $pdbdir, $pdbaffine, $pdb_extr_chain) {
    (-e $_) || die "$_ not found\n";
} 



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

my  $tgt_chain_pdb;
my  $qry_chain_pdb;
my $tfmd_pdb;

############################################
# do we have everything we need?
foreach ( $struct,  $params_file, $dbdir, $pymol,$pdbaffine, 
	  "tfms", "maps", "pdbs") {
    (! -e $_) && `mkdir $_` ;
}


############################################
@top_list = split "\n", `sort -grk 5 $db_search| grep -v done | head -n $top_number`;

`mv $db_search $db_search.full`;

(-e "postp_out") && `rm postp_out`;
`touch postp_out`;

############################################
foreach $hit (@top_list ) {

    ($qry_name, $tgt_name) = split " ", $hit;

    # locate dbfiles
    $tgt_db = `ls $dbdir/$tgt_name.db`; chomp $tgt_db;
    $qry_db = `ls $dbdir/$qry_name.db`; chomp $qry_db;

    # locate pdbfiles    <<<<<<
    #$tgt_pdb_name = substr $tgt_name, 0, 4;
    #$qry_pdb_name = substr $qry_name, 0, 4;
    $tgt_pdb_name = substr $tgt_name, 0, 7;
    $qry_pdb_name = substr $qry_name, 0, 7;

    $tgt_pdb = "$pdbdir/$tgt_pdb_name.pdb";
    $qry_pdb = "$pdbdir/$qry_pdb_name.pdb";

    if (  -z  $tgt_pdb) {
	`rm $tgt_pdb`;
	$cmd = "/home/ivanam/perlscr/downloading/pdbdownload.pl  $tgt_pdb_name";
	(system $cmd) && die "error running $cmd\n";
    }
    if (  -z  $qry_pdb) {
	`rm $qry_pdb`;
	$cmd = "/home/ivanam/perlscr/downloading/pdbdownload.pl  $qry_pdb_name";
	(system $cmd) && die "error running $cmd\n";
    }

    foreach ($tgt_db, $qry_db, $tgt_pdb, $qry_pdb) {
	(-e $_) || die "$_ not found\n";
 	(! -z $_) || die "$_ is empty\n";
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

    `echo pdbf_tgt  $tgt_pdb   >> $new_params`;
    `echo chain_tgt $tgt_chain >> $new_params`;

    `echo pdbf_qry  $qry_pdb   >> $new_params`;
    `echo chain_qry $qry_chain >> $new_params`;

    `echo verbose >> $new_params`;
    
 
    ############################################
    # run struct

    $outname = $qry_name;

    ($struct_out, $struct_tfm) = ("$outname.long_out",
				  "$outname.struct_out.for_postp");
    $cmd = "$struct  $tgt_db  $qry_db $new_params > $struct_out ";
    if (system $cmd) {
	warn "Error running $cmd.\n";
    } else {
	`grep -v done $outname.struct_out >>  postp_out`;
	`mv $outname.struct_out.for_postp tfms/$tgt_name\_to_$qry_name.tfm`;
	`mv $struct_out maps/$tgt_name\_to_$qry_name.map`;
    }

}

############################
# sort the output by the aligned bb score
`sort -grk 9 postp_out > postp_out.sorted`;

@top_list = split "\n", `head -n 30 postp_out.sorted`;

my $loop_ctr = 0;


foreach $hit (@top_list ) {

    print $hit."\n";

    ($qry_name, $tgt_name) = split " ", $hit;

    ($qry_name eq $tgt_name) && next;

    # locate dbfiles
    $tgt_db = `ls $dbdir/$tgt_name.db`; chomp $tgt_db;
    $qry_db = `ls $dbdir/$qry_name.db`; chomp $qry_db;

    # locate pdbfiles    <<<<<<
    $tgt_pdb_name = substr $tgt_name, 0, 7;
    $qry_pdb_name = substr $qry_name, 0, 7;
    $tgt_chain = substr $tgt_name, 4, 1;
    $qry_chain = substr $qry_name, 4, 1;

    $tgt_pdb = "$pdbdir/$tgt_pdb_name.pdb";
    $qry_pdb = "$pdbdir/$qry_pdb_name.pdb";
    $struct_tfm = "tfms/$tgt_name\_to_$qry_name.tfm";
    $struct_out = "maps/$tgt_name\_to_$qry_name.map";
    foreach ($tgt_db, $qry_db, $struct_tfm, $struct_out) {
	(-e $_) || die "$_ not found\n";
	(! -z $_) || die "$_ is empty\n";
    }
    

    ############################################
    # read in the names and the mapped structures
    # from the long output of struct
    $filename = $struct_out;
    open (IF, "<$filename") ||  die "Cno $filename: $!.\n";
    ($reading_map, $reading_rotation, $i) = (0,0,0);
    %region = ();
    @qry_sse_ids = ();
    %map = ();
    @rotation = ();
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
    # extract single chains


    $tgt_chain_pdb = "pdbs/$tgt_name.pdb";
    $cmd = "$pdb_extr_chain $tgt_pdb  $tgt_chain > $tgt_chain_pdb";
    (system $cmd) && die "Error running $cmd.\n";

    $qry_chain_pdb = "pdbs/$qry_name.pdb";
    $cmd = "$pdb_extr_chain $qry_pdb  $qry_chain > $qry_chain_pdb";
    (system $cmd) && die "Error running $cmd.\n";




    ######################################################################
    # transform the original pdb file
    $tfmd_pdb = "$tgt_name\_to_$qry_name.tfmd.pdb";
    $cmd = "$pdbaffine $tgt_chain_pdb $struct_tfm > $tfmd_pdb";
    
    (system $cmd) && die "Error running $cmd.\n";

=pod
    ######################################################################
    # create and run pymol script --> create a pymol session


    (! -e "pymol") && `mkdir pymol` ;

    my $pml_name    = "$tgt_pdb_id\_$qry_pdb_id.pml";
    my $pml_pse     = "$tgt_pdb_id\_$qry_pdb_id.pse";
    my $pml_string  = "";
    my $matched_set = "";

    $pml_string .= "load $tfmd_pdb, $tgt_pdb_id\n";
    $pml_string .= "load $qry_chain_pdb, $qry_pdb_id\n";
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

    `mv $pml_name $pml_pse pymol`;
=cut
    `mv $tfmd_pdb pdbs`;

    $loop_ctr++;

=pod
    ######################################################################
    # create chimera script --> create a chimera session
    my $outdir = "chimera";
    my $chi_script = "$tgt_pdb_id\_$qry_pdb_id.com";
    my $chi_session = "$tgt_pdb_id\_$qry_pdb_id.py";
    my $chi_string  = "";

    (! -e $outdir) && `mkdir $outdir` ;

    # Consolidate data
    my @chi_data_qry = ();
    my @chi_data_tgt = ();
    for my $sse ( @qry_sse_ids ) {
	my %qry = ();
	my %tgt = ();
	$qry{id} = $sse;
	$tgt{id} = $map{$sse};
	$qry{fromto} = join('-', split(/\s/, $region{"1_$qry{id}"}));
	$tgt{fromto} = join('-', split(/\s/, $region{"2_$tgt{id}"}));
	push(@chi_data_qry, \%qry);
	push(@chi_data_tgt, \%tgt);
    }
    my %chi_data = ('qry' => {'file_path' => $qry_chain_pdb,
                          'pdb_id' => $qry_pdb_id,
                          'model_id' => 0,
                          'sse_list' => \@chi_data_qry},
                'tgt' => {'file_path' => $tfmd_pdb,
                          'pdb_id' => $tgt_pdb_id,
                          'model_id' => 1,
                          'sse_list' => \@chi_data_tgt});

    # load files
    for my $role ('qry', 'tgt') {
	$chi_string .=
	    "open noprefs $chi_data{$role}{model_id}" .
        " $chi_data{$role}{file_path}\n";
    }

    # general set up
    $chi_string .=
	"set bg_color white\n" .
	"~show\n" .
	"\n" .
	"colordef qry_matched   0.0 0.0 1.0 1.0\n" .
	"colordef qry_unmatched 0.0 0.0 1.0 0.2\n" .
	"\n" .
	"colordef tgt_matched   1.0 0.0 0.0 1.0\n" .
	"colordef tgt_unmatched 1.0 0.0 0.0 0.2\n" .
	"\n";
    
    # name and highlight matched sse's
    for my $role ('qry', 'tgt') {
	my $model_id = $chi_data{$role}{model_id};
	$chi_string .=
	    "select #$model_id\n" .
	    "namesel $chi_data{$role}{pdb_id}_full\n" .
	    "color ${role}_unmatched selected\n" .
	    "\n";
	for my $sse (@{$chi_data{$role}{sse_list}}) {
	    $chi_string .=
		"select #$model_id:$sse->{fromto}\n" .
		"namesel qry_$sse->{id}\n" .
		"color ${role}_matched selected\n" .
		"\n";
	}
    }

    # name entire matched set
    sub chi_sse_string {
	my ($role) = @_;
	my @fromto_list = ();
	for my $sse (@{$chi_data{$role}{sse_list}}) {
	    push(@fromto_list, $sse->{fromto});
	}
	"#$chi_data{$role}{model_id}:" . join(',', @fromto_list);
    }
    $chi_string .=
	"select " . chi_sse_string('qry') . " | " . chi_sse_string('tgt') . "\n" .
	"namesel matched_set\n" .
	"\n";

    # display as ribbons,
    # save session
    $chi_string .=
	"~select\n" .
	"ribbon\n" .
	"\n" .
	"save $chi_session\n";

    my $chi_fd;
    open($chi_fd, ">$chi_script") || die "Cno $chi_script: $!.\n";
    print($chi_fd $chi_string);
    close $chi_fd;

    `mv $chi_script $outdir`;
    `mv $chi_session $outdir`;
=cut
}



