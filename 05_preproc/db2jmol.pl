#!/usr/bin/perl -w

use strict;

if (!((@ARGV == 1) || (@ARGV == 2))) {
    print(STDERR "Usage: $0 dbfile [pdbfile]\n");
    exit(1);
}
my $dbfile  = $ARGV[0];
my $pdbfile = (@ARGV == 2) ? $ARGV[1] : undef;


if (defined $pdbfile) {
    print("load $pdbfile\n" .
	  "spacefill OFF\n" .
	  "wireframe OFF\n" .
	  #"color GROUP\n" .
	  "color blue\n" .
	  #"trace 0.5\n");
	  "backbone\n");
}

my $id = 0;
open(my $fh, $dbfile) or die("Unable to open '$dbfile': $!\n");
while (<$fh>) {
    if (/^(HELIX|STRAND)/) {
	my @f = split(' ');
	(@f == 12) or die("Wrong number of fields: " . scalar(@f) . "\n");
	my $sse = $f[0];
	my $len = $f[3];
	my @dir = @f[6,7,8];
	my @pos = @f[9,10,11];

	my $name = lc($sse) . 'X' . $id;

	if ($sse eq 'HELIX') {
	    $len *= 1.5 / 2.0;
	}
	if ($sse eq 'STRAND') {
	    $len *= 3.0 / 2.0;
	}
	@dir = map({$len * $_} @dir);
	#@dir = ($len * $dir[0], $len * $dir[1], $len * $dir[2]);

	my @bgn = ($pos[0] - $dir[0],
		   $pos[1] - $dir[1],
		   $pos[2] - $dir[2]);
	my @end = ($pos[0] + $dir[0],
		   $pos[1] + $dir[1],
		   $pos[2] + $dir[2]);

	print("draw $name arrow" .
	      " {$bgn[0] $bgn[1] $bgn[2]}" .
	      " {$end[0] $end[1] $end[2]}" .
	      " width 0.5; color \$$name gold\n");

	$id++;
    }
}
close($fh);#unnecesary
