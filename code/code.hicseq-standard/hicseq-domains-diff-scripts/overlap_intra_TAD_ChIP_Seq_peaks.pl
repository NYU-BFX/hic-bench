#perl

use strict;
use warnings;

# BED1 has to be the TSV from TAD disruptions, BED2 has to be diffBind BED file
my $tads = $ARGV[0];
my $bed_1 = $ARGV[1];
my $extension = $ARGV[2];
my $out = $ARGV[3];
#my $out2 = $ARGV[4]
#print "out=$out\n";
#exit;

my $ARGC = @ARGV;

if ($ARGC == 0) {
	die("Usage: perl SCRIPT bed_1 bed_2 extension-bp[INT] out");
}

my $bin_size = 40000;
my %bed_1;
open(BED, "<$bed_1");
while (my $line = <BED>) {
	chomp($line);
	my @splitted_line = split(/\t/, $line);
	unless (defined($bed_1{$splitted_line[0]})) {
		my @array = ();
		$bed_1{$splitted_line[0]} = \@array;
	}
	push(@{$bed_1{$splitted_line[0]}}, \@splitted_line);
#	push(@bed_1, \@splitted_line);
}
close(BED);


my %seen;

my $boundary_ext=160000;
my $boundary_ext_inner=40000;

open(TADS, "<$tads");
open(OUT, ">$out");
#open(OUT2, ">$out2")
my $header = <TADS>;
while (my $line = <TADS>) {
	chomp($line);
	my @splitted_line = split(/\t/, $line);
#	my $found = 0;

	my $bed_1_hits = 0;
	my $bed_2_hits = 0;

	my $sum_bed_1 = 0;
	my $sum_bed_2 = 1;
	foreach my $bed_entry (@{$bed_1{$splitted_line[0]}}) {
#		if ((($splitted_line[2] + $boundary_ext) >= ($$bed_entry[1] - $extension) && ($splitted_line[2] - $boundary_ext) <= ($$bed_entry[2] + $extension)) || 
#				(($splitted_line[4] + $boundary_ext) >= ($$bed_entry[1] - $extension) && ($splitted_line[4] - $boundary_ext) <= ($$bed_entry[2] + $extension))) {
		if ((($splitted_line[2] - $extension) <= $$bed_entry[1]) && (($splitted_line[4] + $extension) >= $$bed_entry[2])) {
#		if ((($splitted_line[2]) >= ($$bed_entry[1] - $extension) && ($splitted_line[2] - $boundary_ext) <= ($$bed_entry[2] + $extension)) || 
#				(($splitted_line[4] + $boundary_ext) >= ($$bed_entry[1] - $extension) && ($splitted_line[4]) <= ($$bed_entry[2] + $extension))) {
			$sum_bed_1 += $$bed_entry[5];
			$bed_1_hits++;
			#print OUT2 "$splitted_line[0]\t.$splitted_line[2]\t.$splitted_line[4]\t.$$bed_entry[1]\t.$$bed_entry[2]"."\n"
		}
	}
#	foreach my $bed_entry (@{$bed_2{$splitted_line[0]}}) {
#		if (($splitted_line[3]*$bin_size) >= ($$bed_entry[1] - $extension) && ($splitted_line[2]*$bin_size) <= ($$bed_entry[2] + $extension)) {
#			$sum_bed_2 += $$bed_entry[3];
#			$bed_2_hits++;
#		}
#	}

#	print "bed_1_hits=$bed_1_hits vs bed_2_hits=$bed_2_hits\n";

#	my $sumFC = ($sum_bed_2 - $sum_bed_1);

	# OR TAKE 11 FOR SUM
	print OUT "$sum_bed_1\t".$splitted_line[8]."\n";
#	print "sumFC=$sum_bed_1\n";
	
}
close(TADS);
close(OUT);
#close(OUT2);

sub min() {
	my $val1 = shift();
	my $val2 = shift();
	if ($val1 <= $val2) {
		return $val1;
	} else {
		return $val2;
	}
}

sub max() {
	my $val1 = shift();
	my $val2 = shift();
	if ($val1 >= $val2) {
		return $val1;
	} else {
		return $val2;
	}
}




