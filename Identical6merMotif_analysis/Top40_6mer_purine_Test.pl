use strict;
use warnings;
use List::Util qw(sum shuffle);




my $file=shift; # identical 6mer of same strand fragments of a sample group
my $list=shift; # 6mer motifs in generated background fragments of  a sample group

my @expriment_group=();
my @backgroud_group=();

open(IN,$file)or die $!;
my %K6mer=();
while(my $line=<IN>){
	chomp $line;
	my @temp=split /\t/,$line;
	if(!exists($K6mer{$temp[0]})){
		$K6mer{$temp[0]}=$temp[1];
	}
}
close(IN);

open(LI,$list)or die $!;
my %BG_K6mer=();
while(my $BG_file=<LI>){
	chomp $BG_file;
	open(BG,$BG_file)or die $!;
		while(my $line=<BG>){
			chomp $line;
			my @temp=split /\t/,$line;
			if(!exists($BG_K6mer{$temp[0]})){
				$BG_K6mer{$temp[0]}=calculate_purine_content($temp[0]);
			}
			if(!exists($BG_K6mer{$temp[1]})){
				$BG_K6mer{$temp[1]}=calculate_purine_content($temp[1]);
			}
		}
	
	close(BG);
}
close(LI);
my @BMer=keys(%BG_K6mer);

my @Mer=sort{$K6mer{$b}<=>$K6mer{$a}}keys(%K6mer);
for (my $i=0;$i<=39;$i++){
	my $GC=calculate_gc_content($Mer[$i]);
	my $Purine=calculate_purine_content($Mer[$i]);
	push @expriment_group,$Purine;
	generate_gc_matched_background_motif($GC,10);
}


my ($obs, $p, $null_ref)=permutation_test(\@expriment_group,\@backgroud_group,10000,'greater');
my $ci_ref = bootstrap_mean_diff(\@expriment_group,\@backgroud_group);

open(OUT,">>Purine_Permutation_Test.txt")or die $!;
printf OUT "\n=== $file ===Permutation Test (Perl) ===\n";
printf OUT "Observed mean difference = %.4f\n", $obs;
printf OUT "Permutation test one tailed p-value < 1e%d\n", $p;
printf OUT "95%% CI of mean difference = [%.3f, %.3f]\n", @$ci_ref;
close(OUT);

sub calculate_gc_content {
    my $seq = shift;
    my $gc_count = ($seq =~ tr/GCgc//);
    return $gc_count / 6;
}

sub calculate_purine_content {
    my $seq = shift;
    my $purine_count = ($seq =~ tr/AGag//);
    return $purine_count / 6;
}

sub generate_gc_matched_background_motif{
	my ($GC,$num_backgrounds)=@_;
	my %dif_GC=();
	#my %gc_matched_background_motif=();
	for my $bg_motif(@BMer){
		$dif_GC{$bg_motif}=abs($BG_K6mer{$bg_motif}-$GC);
	}
	my @dif_BMer=sort{$dif_GC{$a}<=>$dif_GC{$b}}keys(%dif_GC);
	

	my %value_count=();
	my $count=0;
	for my $bg_motif(@dif_BMer){
		if($count<$num_backgrounds){
			#$gc_matched_background_motif{$bg_motif}=$BG_K6mer{$bg_motif};
			push @backgroud_group, $BG_K6mer{$bg_motif};
			$count++;
		}else{
			last;
		}
	}
}

sub permutation_test {
    my ($g1_ref, $g2_ref, $n_perm, $alternative) = @_;
    $alternative //= 'greater';

    my @group1 = @$g1_ref;
    my @group2 = @$g2_ref;

    # Observed effect size: mean difference
    my $mean1 = sum(@group1) / @group1;
    my $mean2 = sum(@group2) / @group2;
    my $obs_diff = $mean1 - $mean2;

    # Permutation after merging.
    my @combined = (@group1, @group2);
    my $n1 = @group1;
    my @null_diffs;

    for (1 .. $n_perm) {
        my @shuffled = shuffle @combined;
        my @perm1 = @shuffled[0 .. $n1-1];
        my @perm2 = @shuffled[$n1 .. $#shuffled];

        my $m1 = sum(@perm1) / @perm1;
        my $m2 = sum(@perm2) / @perm2;
        push @null_diffs, $m1 - $m2;
    }
    # compute p-value
    my $power;
    if ($alternative eq 'greater') {
        my $count = grep { $_ >= $obs_diff } @null_diffs;
        my $p_adj = ($count + 1) / (@null_diffs + 1);
		$power = int(log($p_adj) / log(10)) - 1;  
        
    }
    elsif ($alternative eq 'two-sided') {
        my $count = grep { abs($_) >= abs($obs_diff) } @null_diffs;
        my $p_adj = ($count + 1) / (@null_diffs + 1);
		$power = int(log($p_adj) / log(10)) - 1;  
    }
    else {
        die "Unknown alternative: $alternative\n";
    }

    return ($obs_diff, $power, \@null_diffs);
}

#----------------------------------------------------------
# Bootstrap mean difference confidence interval
# Return the [2.5%, 97.5%] percentiles
#----------------------------------------------------------
sub bootstrap_mean_diff {
    my ($g1_ref, $g2_ref, $n_boot) = @_;
    $n_boot //= 10000;
    my @group1 = @$g1_ref;
    my @group2 = @$g2_ref;
    my $n1 = @group1;
    my $n2 = @group2;

    my @boot_diffs;
    for (1 .. $n_boot) {
        my @s1 = map { $group1[int rand $n1] } 1 .. $n1;
        my @s2 = map { $group2[int rand $n2] } 1 .. $n2;
        push @boot_diffs, (sum(@s1)/@s1) - (sum(@s2)/@s2);
    }

    @boot_diffs = sort { $a <=> $b } @boot_diffs;
    my $lo = $boot_diffs[ int(0.025 * $#boot_diffs) ];
    my $hi = $boot_diffs[ int(0.975 * $#boot_diffs) ];
    return [$lo, $hi];
}
