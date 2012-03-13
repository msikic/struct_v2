#! /usr/bin/perl -w

defined ( $ARGV[0]  ) ||
    die "Usage: $0   <pdb_file>   [<chain_name>].\n";

$pdbfile = $ARGV[0];
if ( defined $ARGV[1] ) {
    $query_chain_name =$ARGV[1] ;
} else {
    $query_chain_name ="" ;
}

open ( IF, "<$pdbfile") ||
    die "Cno $pdbfile: $!.\n";

$exceptions = "ACE_PCA_DIP_NH2_LTA_MSE";
    
while ( <IF> ) {

    last if ( /^ENDMDL/);
    #last if ( /^TER/);
    next if ( ! /^ATOM/ && ! /^HETATM/ );
    $line = $_;
    $chain_name = substr ( $line,  21, 1) ;
    next if ( $query_chain_name &&   $chain_name ne " "  && 
	      $chain_name ne $query_chain_name );
    $res_name = substr $_,  17, 3; $res_name =~ s/\s//g;
    next if ( /^HETATM/ && $exceptions !~ $res_name );
    print $line;
    
}

close IF;
