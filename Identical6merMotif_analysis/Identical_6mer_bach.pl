use strict;
my $list=shift;
open(IN, $list)or die $!;
while(my $line=<IN>){
	chomp $line;
	my @temp=split /\s/,$line;
	my $codes="perl identical6mer_enrichment.pl $temp[0] $temp [1] $temp[2]";
	print `$codes`;
}

close(IN);