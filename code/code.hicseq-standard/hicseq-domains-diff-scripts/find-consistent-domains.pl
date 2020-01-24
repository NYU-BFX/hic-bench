#perl

use strict;
use warnings;

my $sample_1_tads = $ARGV[0];
my $sample_2_tads = $ARGV[1];
my $max_diff = $ARGV[2];
my $bin_size = $ARGV[3];
my $out_sample_1 = $ARGV[4];
my $out_sample_2 = $ARGV[5];
# Added code to print a file with the common TADs
my $out_common = $ARGV[6];

my %sample_1;
open(sample_1, "<$sample_1_tads");
while (my $line = <sample_1>) {
	chomp($line);
	my @splitted_line = split(/\t/, $line);
	$sample_1{$splitted_line[0] . ":" . $splitted_line[1] . "-" . $splitted_line[2]} = \@splitted_line;
}
close(sample_1);

open(sample_2, "<$sample_2_tads");
open(OUT_1, ">$out_sample_1");
open(OUT_2, ">$out_sample_2");
open(OUT_common, ">$out_common");

while (my $line = <sample_2>) {
	chomp($line);
	my @splitted_line = split(/\t/, $line);

	for (my $i = (-1)*$max_diff; $i < $max_diff; $i++) {
		for (my $j = (-1)*$max_diff; $j < $max_diff; $j++) {
			my $offset_i = $i * $bin_size;
			my $offset_j = $j * $bin_size;
			if (defined($sample_1{$splitted_line[0] . ":" . ($splitted_line[1] + $offset_i) . "-" . ($splitted_line[2] + $offset_j)})) {
				print "FOUND MATCH with offsets $offset_i and $offset_j\n";
				print OUT_1 join("\t", @{$sample_1{$splitted_line[0] . ":" . ($splitted_line[1] + $offset_i) . "-" . ($splitted_line[2] + $offset_j)}}) . "\n";
				print OUT_2 join("\t", @splitted_line) . "\n";
				print OUT_common $splitted_line[0] . "\t" . &min($splitted_line[1], ($splitted_line[1]+$offset_i)) . "\t" . &max($splitted_line[2], ($splitted_line[2] + $offset_j)) . "\n"; 
			} 
		}
	}

}
close(OUT_2);
close(OUT_1);
close(OUT_common);
close(sample_2);

sub max {
	my $val1 = shift();
	my $val2 = shift();
	if ($val1 > $val2) {
		return($val1);
	}
	return($val2);
}

sub min {
	my $val1 = shift();
	my $val2 = shift();
	if ($val1 < $val2) {
		return($val1);
	}
	return($val2);
}


