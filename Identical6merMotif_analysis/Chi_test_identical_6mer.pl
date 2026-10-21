#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum min max);
use POSIX qw(log10);
use Statistics::Distributions;

# ============================================
# statistical functions
# ============================================

# Computing the chi-square test statistic (with Yates' continuity correction)
sub chi2_test {
    my ($a, $b, $c, $d) = @_;
    
    my $n = $a + $b + $c + $d;
    return (0, 1, 0) if $n == 0;  # empty table
    
    # Computing the expected frequencies
    my $e11 = ($a + $b) * ($a + $c) / $n;
    my $e12 = ($a + $b) * ($b + $d) / $n;
    my $e21 = ($c + $d) * ($a + $c) / $n;
    my $e22 = ($c + $d) * ($b + $d) / $n;
    
    # Yates' continuity correction
    my $chi2 = 0;
    $chi2 += (abs($a - $e11) - 0.5)**2 / $e11 if $e11 > 0;
    $chi2 += (abs($b - $e12) - 0.5)**2 / $e12 if $e12 > 0;
    $chi2 += (abs($c - $e21) - 0.5)**2 / $e21 if $e21 > 0;
    $chi2 += (abs($d - $e22) - 0.5)**2 / $e22 if $e22 > 0;
    
    # free degree = 1 (2x2 table)
    my $df = 1;
    
    # compute p-value
    my $p_value = chi2_pvalue($chi2, $df);
    
    return ($chi2, $p_value, $df);
}

sub chi2_pvalue {
    my ($chi2, $df) = @_;
    
    if ($df == 1 && $chi2 > 100) {
        my $z = sqrt($chi2);
        my $z2 = $z * $z;
        my $log_p = log(2) - 0.5*$z2 - log($z) - 0.5*log(2*3.141592653589793);
        $log_p += log(1 - 1/$z2 + 3/($z2*$z2));
        
        return $log_p;  # returen log p-value
    }
    
    my $p = Statistics::Distributions::chisqrprob($df, $chi2);
    return log($p);  # returen log p-value
}



# Computing the OR (Odds Ratio) and confidence interval
sub calculate_or_ci {
    my ($a, $b, $c, $d) = @_;
    
    # Handling division by zero
    if ($b == 0 || $c == 0) {
        return ('Inf', 'Inf', 'Inf');
    }
    
    # Computing the OR 
    my $or = ($a * $d) / ($b * $c);
    
    # Computing the STD 
    my $se = sqrt(1/$a + 1/$b + 1/$c + 1/$d);
    
    # 95% confidence interval 
    my $z = 1.96;
    my $log_or = log($or);
    my $log_lower = $log_or - $z * $se;
    my $log_upper = $log_or + $z * $se;
    
    my $ci_lower = exp($log_lower);
    my $ci_upper = exp($log_upper);
    
    return ($or, $ci_lower, $ci_upper);
}

# Formatting numerical output


sub format_number {
    my ($num) = @_;
    
    return 'Inf' if $num eq 'Inf' || $num > 1e10;
    
    if ($num < 0.001 && $num > 0) {
        return sprintf("%.2e", $num);
    } else {
        return sprintf("%.4f", $num);
    }
}

sub fdr_correction_log {
    my @log_p_values = @_;
    my $num_tests = scalar(@log_p_values);
    
   
    return @log_p_values if $num_tests <= 1;
    
    
    my @tests;
    for my $i (0 .. $num_tests-1) {
        push @tests, {
            original_index => $i,
            log_p => $log_p_values[$i]
        };
    }
    
  
    @tests = sort { $a->{log_p} <=> $b->{log_p} } @tests;
    
    
   
    my $log_m = log($num_tests);
    my $last_log_q = 0;  # log(1) = 0
   
    for my $idx (reverse 0 .. $num_tests-1) {
        my $rank = $idx + 1;  
        my $test=$tests[$idx];
        # BH formula: q = p * m / rank
        my $log_q = $test->{log_p} + $log_m - log($rank);
        
        # constrains: q ≤ 1
        $log_q = 0 if $log_q > 0;
        
       
        $log_q = $last_log_q if $log_q > $last_log_q;
        $test->{log_p}= $log_q;
        $last_log_q = $log_q;
    }
     
    my @corrected = (0) x $num_tests;
    foreach my $test (@tests) {
        $corrected[$test->{original_index}] = $test->{log_p};
    }
    
    return @corrected;
}

# Formatting p-values
sub format_log_pvalue {
    my ($log_p, $type) = @_;
    $type ||= 'pvalue';
    
    #
    return "1.0000" if $log_p >= 0;  # p ≥ 1
    

    my $exponent = int($log_p / log(10));
    
    if ($exponent < -100) {
       
        return "< 10^$exponent";
    } elsif ($exponent < -10) {
        my $mantissa = exp($log_p - $exponent * log(10));
        return sprintf("%.2e", $mantissa * (10 ** $exponent));
    } else {
        my $p = exp($log_p);
        return $type eq 'qvalue' ? sprintf("%.6f", $p) : sprintf("%.4f", $p);
    }
}

sub is_significant_log {
    my ($log_p, $alpha) = @_;
    $alpha ||= 0.05; 
    return $log_p < log($alpha);
}


# ============================================
# Inputfile process
# ============================================

#parse_input_file
sub parse_input_file {
    my ($filename) = @_;
    my %sample_data;
    
    open(my $fh, '<', $filename) or die "Error: Cannot open file '$filename'\n";
    
    my $line_num = 0;
    while (my $line = <$fh>) {
        $line_num++;
       chomp $line;
        my @parts = split(/\t/, $line);
        
        
        my $sample_name = $parts[0];
        my $a=$parts[2];
        my $b=$parts[1]-$parts[2];
        my $c=$parts[3];
        my $d=$parts[1]-$parts[3];
        $sample_data{$sample_name} = [$a, $b, $c, $d];
    }
    
    close($fh);
    
    die "Error: No valid data found in file\n" unless %sample_data;
    
    return %sample_data;
}

# ============================================
# main 
# ============================================

# Command-line argument checking
if (@ARGV != 1) {
    print "Usage: perl $0 <input_file>\n";
    print "\nInput file format: sample_name N True_Candidate True_background\n";
    print "Example: Sample1	213601	197000	4553\n";
    exit 1;
}

my $input_file = $ARGV[0];

print "Processing file: $input_file\n";

# 1. read data
my %sample_data = parse_input_file($input_file);
my @sample_names = keys %sample_data;
my $num_samples = scalar(@sample_names);

print "Found $num_samples samples\n\n";

# 2. analyze each sample 
print "INDIVIDUAL SAMPLE ANALYSIS\n";
print "=" x 60 . "\n";

my @results;
my @p_values;

foreach my $sample (@sample_names) {
    my ($a, $b, $c, $d) = @{$sample_data{$sample}};
    
    print "\nSample: $sample\n";
    print "  Table: $a $b $c $d\n";
    
    # chi-square test 
    my ($chi2, $p_value, $df) = chi2_test($a, $b, $c, $d);
    print "  Chi2($df) = " . format_number($chi2) . ", p = " . format_log_pvalue($p_value) . "\n";
    
    # OR (Odds Ratio) and confidence interval
    my ($or, $ci_lower, $ci_upper) = calculate_or_ci($a, $b, $c, $d);
    print "  OR = " . format_number($or) . ", 95% CI: [" . 
          format_number($ci_lower) . ", " . format_number($ci_upper) . "]\n";
    
    push @results, {
        sample => $sample,
        a => $a, b => $b, c => $c, d => $d,
        chi2 => $chi2,
        p_value => $p_value,
        df => $df,
        or => $or,
        ci_lower => $ci_lower,
        ci_upper => $ci_upper,
        significant_raw => is_significant_log($p_value) ? 1 : 0
    };
    
    push @p_values, $p_value;
}

# 3. fdr_correction
print "\n\nMULTIPLE TESTING CORRECTION\n";
print "=" x 60 . "\n";

if ($num_samples <= 1) {
    print "Only one sample, skipping FDR correction\n";
} else {
	print "pvalue:   \n@p_values\n\n";
	
    my @corrected_p = fdr_correction_log(@p_values);
    
    print "corrected_p:   \n@corrected_p\n\n";
    
    my $raw_sig = 0;
    my $fdr_sig = 0;
    
    for my $i (0 .. $#results) {
        $results[$i]{p_fdr} = $corrected_p[$i];
        $results[$i]{significant_fdr} = is_significant_log($corrected_p[$i]) ? 1 : 0;
        
        $raw_sig++ if $results[$i]{significant_raw};
        $fdr_sig++ if $results[$i]{significant_fdr};
    }
    
    print "Significant samples (raw p<0.05): $raw_sig/$num_samples\n";
    print "Significant samples (FDR q<0.05): $fdr_sig/$num_samples\n";
    
    print "\nFDR-corrected results:\n";
    foreach my $res (@results) {
        printf "  %-10s: p=%-8s q=%-8s %s\n",
            $res->{sample},
            format_log_pvalue($res->{p_value}),
            format_log_pvalue($res->{p_fdr}),
            ($res->{significant_fdr} ? "SIGNIFICANT" : "not significant");
    }
}

# 4. Combined analysis of a group of samples 
print "\n\nCOMBINED ANALYSIS\n";
print "=" x 60 . "\n";
my %total_a=();
my %total_b=();
my %total_c=();
my %total_d=();

foreach my $sample (@sample_names) {
	my @temp= split /\//, $sample;
	if(! exists($total_a{$temp[2]})){
		$total_a{$temp[2]}=0;
		$total_b{$temp[2]}=0;
		$total_c{$temp[2]}=0;
		$total_d{$temp[2]}=0;
	}
    $total_a{$temp[2]} += $sample_data{$sample}[0];
    $total_b{$temp[2]} += $sample_data{$sample}[1];
    $total_c{$temp[2]} += $sample_data{$sample}[2];
    $total_d{$temp[2]} += $sample_data{$sample}[3];
    print "Combined table: $total_a{$temp[2]} $total_b{$temp[2]} $total_c{$temp[2]} $total_d{$temp[2]}\n";
}



my @groups=keys(%total_a);
my %chi2_total=();
my %p_total=();
my %df_total=();
my %or_total=();
my %ci_lower_total=();
my %ci_upper_total=();
my %n_total=();
my %phi=();
foreach my $group (@groups){
	# chi-square test 
	($chi2_total{$group}, $p_total{$group}, $df_total{$group}) = chi2_test($total_a{$group}, $total_b{$group}, $total_c{$group}, $total_d{$group});
	print "Chi2($df_total{$group}) = " . format_number($chi2_total{$group}) . ", p = " . format_log_pvalue($p_total{$group}) . "\n";

	# OR (Odds Ratio) and confidence interval
	($or_total{$group}, $ci_lower_total{$group}, $ci_upper_total{$group}) = calculate_or_ci($total_a{$group}, $total_b{$group}, $total_c{$group}, $total_d{$group});
	print "OR = " . format_number($or_total{$group}) . ", 95% CI: [" . 
      format_number($ci_lower_total{$group}) . ", " . format_number($ci_upper_total{$group}) . "]\n";

	# effect size Phi 
	$n_total{$group} = $total_a{$group} + $total_b{$group} + $total_c{$group} + $total_d{$group};
	$phi{$group} = sqrt($chi2_total{$group} / $n_total{$group}) if $n_total{$group} > 0;
	print "Effect size (Phi) = " . format_number($phi{$group}) . "\n";

	# effect size explantation
	if ($phi{$group} < 0.1) {
   	 print "Effect interpretation: Small\n";
	} elsif ($phi{$group} < 0.3) {
   	 print "Effect interpretation: Medium\n";
	} else {
  	  print "Effect interpretation: Large\n";
	}

}



# 5. save results
print "\nSaving results...\n";

# save the result of each sample 
my $individual_file = "Enrich_identical_6mer.tsv";
open(my $out, '>', $individual_file) or die $!;
print $out "INDIVIDUAL SAMPLE ANALYSIS RESULTS\n";
print $out "=" x 40 . "\n";
print $out "Group\tSample\tProportion of Same-Strand 5' Ends fragments with Identical 6mer ends\tProportion of background fragments with Identical 6mer ends\tChi2\tP_value\tOR\t95%CI_lower\t95%CI_upper\tFDR P_value\n";


my %group_sample_out=();
foreach my $res (@results) {
	my @temp=split /\//, $res->{sample};
	$temp[3]=~s/_spc.*$//;
	my @sample_id=split /_/,$temp[3];
	if(!exists($group_sample_out{$temp[2]})){
		$group_sample_out{$temp[2]}="";
	}
	$group_sample_out{$temp[2]}.=join("\t",
        $temp[2],$sample_id[1],
        $res->{a}/($res->{a}+$res->{b}), $res->{c}/($res->{a}+$res->{b}),
        format_number($res->{chi2}),
        format_log_pvalue($res->{p_value}),
        format_number($res->{or}),
        format_number($res->{ci_lower}),
        format_number($res->{ci_upper}),
        format_log_pvalue(($res->{p_fdr}))
   	)."\n";
		
    

}

foreach my $group (@groups){
	print $out $group_sample_out{$group};
}
print $out "=" x 40 . "\n";
print $out "\n";


# save the result of a group of samples

print $out "COMBINED ANALYSIS RESULTS\n";
print $out "=" x 40 . "\n";
print $out "Group\tProportion of Same-Strand 5' Ends fragments with Identical 6mer ends\tProportion of background fragments with Identical 6mer ends\tChi2\tP-value\tOR\t95% CI_lower\t95% CI_upper\tEffect size (Phi)\n";
foreach my $group (@groups){
	print $out join("\t",
	$group,$total_a{$group}/($total_a{$group}+$total_b{$group}),$total_c{$group}/($total_a{$group}+$total_b{$group}),
	format_number($chi2_total{$group}),
	format_log_pvalue($p_total{$group}),
	format_number($or_total{$group}),
	format_number($ci_lower_total{$group}),
	format_number($ci_upper_total{$group}),
	format_number($phi{$group})
	)."\n";
	
}
print $out "=" x 40 . "\n";
print $out "\n";
close($out);


print "\nAnalysis completed successfully!\n";