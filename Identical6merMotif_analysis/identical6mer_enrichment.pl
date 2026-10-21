#!/usr/bin/perl
use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use List::Util qw(shuffle sum);

# ============================================
# 参数parameters
# ============================================
if (@ARGV != 3) {
    die "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <Target_Fragment_File>\n" .
        "Target_Fragment_Format：pos5end_file\n" .
        "BXX:5:1207:10318:15328	R1	chr3	176646228	+	CCCAACTTGATCTATAGATTTAATGCAATCCCAA	R2	chr1	233804790	+	CCCAACTTGATCTATAGATTTAATGCACTCCCAATCAAAATCTCAGAAAATTCGTTTGTGGATATGAATAAAATG	3	179	R11>..>R20;	chr3:176646228-176646261;chr1:233804790-233804938;
K00306:55:HGW35B" .
        "Output: Background_Fragment_File（Each Line: left_6mer\\tright_6mer\\tleft_gc\\tright_gc）\n";
}

my ($r1_file, $r2_file, $target_file) = @ARGV;

# ============================================
# functions
# ============================================



# ============================================
# main 
# ============================================

# Step 1: 
print STDERR "Step 1: Read the target fragment file...\n";
my %target_ids;
my @target_fragments;

my $N1=0;
my $N2=0;
my $N3=0;
my $N_t=10;
open(IN,$target_file)or die $!;
	while(my $line=<IN>){
		chomp $line;
		my @temp=split /\t/,$line;
		#extract same-strand fragments: same chromosome, same strand, 5'ends have valid coordinates
		if(($temp[2] eq  $temp[7]) && ($temp[4] eq $temp[9]) && ($temp[3] ne "-") && ($temp[8] ne "-")){
			my $left_6mer ="";
			my $right_6mer = "";
			my $left_gc =0;
			my $right_gc = 0;
			if($temp[4] eq "-"){
				$right_6mer =reverse_complement(substr($temp[10],-6));
				$left_6mer =reverse_complement(substr($temp[5],-6));
			}else{
				$left_6mer = substr($temp[5], 0, 6);
				$right_6mer = substr($temp[10], 0, 6);
				
			}
			$left_gc = calculate_gc_content($left_6mer);
			$right_gc = calculate_gc_content($right_6mer);
			
			if($temp[0]=~/#/){
				$temp[0]=~ s/#.*$//;
				my @temp2=split /_/,$temp[0];
				$temp[0]=$temp2[1].":".$temp2[2].":".$temp2[3];
				
			}
			
			if(length($left_6mer) == 6 && length($right_6mer) == 6){
				push @target_fragments, {
				     id => $temp[0],
				     strand => $temp[4],
					 end_dis=> abs($temp[3]-$temp[8]),
       			     left_6mer => $left_6mer,
       			     left_gc => calculate_gc_content($left_6mer),
       			     left_purine => calculate_purine_content($left_6mer),
       			     right_6mer => $right_6mer,
        		     right_gc => calculate_gc_content($right_6mer),
        		     right_purine => calculate_purine_content($right_6mer)
    		   };
    		   
			   $target_ids{$temp[0]}="";
			   if($N_t>0){
			   	print STDERR $temp[0]."\n";
			   	$N_t--;
			   }
			}
			
		}
	}
close(IN);
# Step 2: 	
print STDERR "\nStep 2: construct background 6-mer pool(left and right)...\n";

my %left_pool;  # $left_pool{gc_bucket} = [@entries]
my %right_pool;

my $left_count = extract_6mers_from_fastq($r1_file, \%left_pool);
my $right_count = extract_6mers_from_fastq($r2_file, \%right_pool);
	

# Step 3: generate background for each target fragment 
print STDERR "\nStep 3: generate background fragment...\n";
my @base_file=split /\//, $target_file;
$base_file[-1]=~/pos5end_(.+)_spc/;
#$target_file=~/\/(\w+)$/;
my $output_file = "BG_".$1;
my $output_file2 = "TS_".$1;
my $statistic_file = "Identical_6mer_Enrichment.txt";
open(OUT1, ">$output_file") or die $!;
print OUT1 "left_6mer\tright_6mer\tleft_gc\tright_gc\n";
open(OUT2, ">$output_file2") or die $!;



my $target_identical_6mer=0;
my $background_identical_6mer=0;
my $string1="";
my $string2="";
my $string3="";

my $processed = 0;
foreach my $target (@target_fragments) {
	
    # GC matched 6-mer of left 5'end
    my $background_left = find_gc_matched_6mer(
        $target->{left_gc},
        \%left_pool,
    );
    my $background_left_gc = calculate_gc_content($background_left);  
    # GC matched 6-mer of right 5'end  
    my $background_right = find_gc_matched_6mer(
        $target->{right_gc},
        \%right_pool,
    );
    my $background_right_gc = calculate_gc_content($background_right);
    my $identical_flag=&identical_6mer( $target->{left_6mer},$target->{right_6mer});
    $target_identical_6mer+=$identical_flag;
    $background_identical_6mer+=&identical_6mer($background_left,$background_right );
    
    # 
    $string1.=join("\t",
        $background_left,
        $background_right,
        sprintf("%.3f", $background_left_gc),
        sprintf("%.3f", $background_right_gc)) . "\n";
    $string2.=$target->{id}."\t".$target_ids{$target->{id}}."\t".$target->{strand}."\t".$target->{left_6mer}."\t".$target->{right_6mer}."\tGC:".sprintf("%.3f", $target->{left_gc}).",".sprintf("%.3f", $target->{right_gc})."\tPurine:".sprintf("%.3f", $target->{left_purine}).",".sprintf("%.3f", $target->{right_purine})."\t".$identical_flag."\n";
}

$string3=$target_file."\t".@target_fragments."\t".$target_identical_6mer."\t".$background_identical_6mer."\n";
print OUT1 $string1;
print OUT2 $string2;

print $N1."\t".$N2."\t".$N3."\n";
open(OUT3, ">>$statistic_file") or die $!;
print OUT3 $string3;
close(OUT3);
close();



# compute the GC content of 6-mer motif
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


# extract 6-mer from fastq file
sub extract_6mers_from_fastq {
    my ($fastq_file, $pool_ref) = @_;
    
    print STDERR "Processing $fastq_file ...\n";
    
    my $gz = IO::Uncompress::Gunzip->new($fastq_file) 
        or die "Can not open $fastq_file: $GunzipError\n";
    
    my $line_count = 0;
    my $read_count = 0;
    my $skipped_target = 0;
    my $total_added = 0;
    my $N_b=10;
    
    while (my $read_id = <$gz>) {
        # Four line of a read in fastqfile
        my $seq = <$gz>;
        <$gz>;
        <$gz>;
       
        
        
        
        # check whether the read_id in target fragments, if yes, skip it）
        chomp($read_id);
        my @temp=split /\s+/,$read_id;
        $temp[0]=~/^@(.+)$/;
        $read_id=$1;
        chomp($seq);
        
        $read_count++;
        if($N_b>0){
			   	print STDERR $read_id."\n";
			   	$N_b--;
			   }
      #  print $read_id."\n";
        if (exists $target_ids{$read_id}) {
        	if($target_ids{$read_id} eq ""){
        		$target_ids{$read_id}=$seq;
        	}else{
        		$target_ids{$read_id}.="\t".$seq;
        	}
        	
            $skipped_target++;
            next;
        }
        
        # extract 6-mer of 5' end
        my $six_mer = substr($seq, 0, 6);
        
        # valid 6-mer motif containing 6 bases
        if ($six_mer =~ /^[ATCG]{6}$/i) {
          
            # GC content
            my $gc_content = calculate_gc_content($six_mer);
            
            # GC bins  
            my $gc_bucket = int($gc_content * 100);
            
            # store 6-mer and GC content  
            if (!exists $pool_ref->{$gc_bucket}) {
                $pool_ref->{$gc_bucket} = [];
            }
            
            push @{$pool_ref->{$gc_bucket}}, {
                seq => $six_mer,
                gc => $gc_content
            };
            $total_added++;
              
        }
        if ($read_count % 100000 == 0) {
       		 print "$read_count reads have been processed\n";
        }
    
    }
    
    $gz->close();
    print STDERR "Finished $fastq_file: $read_count reads are processed,  $total_added 6-mer motifs are added, $skipped_target targed reads are skipped\n";
    return $total_added;
}

sub reverse_complement {
    my $seq = shift;
    $seq = reverse($seq);
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}

sub find_gc_matched_6mer {
    my ($target_gc, $pool_ref) = @_;
    
    # GC bin（0-100）
    my $target_gc_bucket = int($target_gc * 100);
    
    # try to find a random read in the same GC bin  
    if (exists $pool_ref->{$target_gc_bucket} && @{$pool_ref->{$target_gc_bucket}} > 0) {
        my $candidates = $pool_ref->{$target_gc_bucket};
        my $random_index = int(rand(@$candidates));
        $N1++;
        return $candidates->[$random_index]{seq};
    }
    
    # if not find in the same GC bin, extend the scop to  ±5% GC content
    for my $delta (1..5) {
        my $lower_bucket = $target_gc_bucket - $delta;
        my $upper_bucket = $target_gc_bucket + $delta;
        
        # collect all candidates
        my @all_candidates;
        foreach my $bucket ($lower_bucket..$upper_bucket) {
            next if $bucket < 0 || $bucket > 100;
            next unless exists $pool_ref->{$bucket};
            
            push @all_candidates, @{$pool_ref->{$bucket}};
        }
        
        if (@all_candidates > 0) {
            my $random_index = int(rand(@all_candidates));
            $N2++;
            return $all_candidates[$random_index]{seq};
        }
        
    }
    
    # if not find in the extended bins, select from the total pool
    my @all_sequences;
    foreach my $bucket (keys %$pool_ref) {
        push @all_sequences, @{$pool_ref->{$bucket}};
    }
    
    if (@all_sequences > 0) {
        my $random_index = int(rand(@all_sequences));
        $N3++;
        return $all_sequences[$random_index]{seq};
    }
    
}

sub identical_6mer{
	my ($S6mer1,$S6mer2) = @_ ;
	
    my $k5_0_0=substr($S6mer1,0,5);
    my $k5_0_1=substr($S6mer1,-5);
    my $k5_1_0=substr($S6mer2,0,5);
    my $k5_1_1=substr($S6mer2,-5);
     	   
    if($S6mer1 eq $S6mer2 or $k5_0_0 eq $k5_1_1 or $k5_0_1 eq $k5_1_0){
		return 1
     }
    return 0;
}
