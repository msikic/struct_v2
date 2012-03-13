#!/usr/bin/perl -wT

use strict;

if (@ARGV == 0) {
    print("Usage: $0 pdfiles...\n");
    exit 1;
}
my @pdbfiles = @ARGV;

for my $pdbfile (@pdbfiles) {
    open(my $fh, $pdbfile) or die("Couldn't open: '$pdbfile'\n");
    while (<$fh>) {
	next if (!/^HELIX/);
	my @c = split(//);
	my $type = "$c[38]$c[39]";
	next if ($type eq ' 1');
	if ($type eq '  ') {
	    print(STDERR "blank helix id: $pdbfile\n");
	    next;
	}
	if ($type eq ' 5') {
	    print("'$type': $pdbfile\n");
	    next;
	}
	print(STDERR "Unrecognized type: '$type'\n");
    }
    close($fh);
}
