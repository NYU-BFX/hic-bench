#peeeeeeerl

use strict;
use warnings;

my $bed_file = $ARGV[0];
my $gene_anno = $ARGV[1];
my $out_promoter = $ARGV[2];
my $out_others = $ARGV[3];

# maximum distance between outer BED entry in bp
my $max_range = 1000;

my %gene_anno;
open(GENEANNO, "<$gene_anno");
while (my $line = <GENEANNO>) {
	chomp($line);

	my @splitted_line = split(/\t/, $line);
	unless (defined($gene_anno{$splitted_line[0]})) {
		my @array = ();
		push(@array, \@splitted_line);
		$gene_anno{$splitted_line[0]} = \@array;
	} else {
		push(@{$gene_anno{$splitted_line[0]}}, \@splitted_line);
	}
}
close(GENEANNO);

my $counter = 0;
open(BED, "<$bed_file");
open(OUT, ">$out_others");
open(OUTP, ">$out_promoter");
#my $header = <BED>;
while (my $line = <BED>) {
	chomp($line);
	my @splitted_line = split(/[\t,]/, $line);
	my $array_ref = $gene_anno{$splitted_line[0]};
	my $hit = 0;
	my $peak_id = "NA";
	my $start = "NA";
	my $end = "NA";
#	my @splitted_id = split(/\./, $splitted_line[3]);
	my $gene_TSS;
#	if ($splitted_id[0] eq "ENSG00000165806") {
#		print "gene_TSS=$gene_TSS\n";
#	}
	foreach my $gene_anno_ref (@{$array_ref}) {

		if ($$gene_anno_ref[5] eq "+") {
			$gene_TSS = $$gene_anno_ref[1];
		} elsif ($$gene_anno_ref[5] eq "-") {
			$gene_TSS = $$gene_anno_ref[2];
		} else {
#			print "Couldn't load strand for gene $splitted_id[0]. Next\n";
			next;
		}


		if ($splitted_line[1] <= ($gene_TSS + $max_range) && ($splitted_line[2] + $max_range) > $gene_TSS) {
			$hit = 1;
			print OUTP $line . "\n";
			last;
		}
	}
	unless ($hit) {
		print OUT $line . "\n";
	}
	$counter++;
	if ($counter % 10000 == 0) {
		print "total of $counter lines passed...\n";
	}
}
close(OUTP);
close(OUT);
close(BED);





