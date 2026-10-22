=cut
	 
	Function: 
		1. Extract two ends of each fragment, following the sequencing direction 5'-3'
		
		2. Perform statistic analysis on the split-aligned fragments from each sample:
    	(1) Proportions of fragments with 1-junction, 2-junction count, or 3-junction count; 
    	(2) Proportions of fragments with 1-junction derived from one chromosome, or >1 chromosomes;
		(3) Proportions of each type of fragments from one chromosome (Intra-cis Determinate, Intra-cis Circular; Intra-cis Complex); 
		(4) Proportions of each type of fragments from two chromosome (Inter-trans Determinate; Inter-trans Complex); 
		(5) Proportions of two type of fragments (Same-Strand 5' Ends fragments; Opposite-Strand 5' Ends fragments);
		(6) Proportion of fragments with Identical 6-mer end motifs in each type (Same-Strand 5' Ends fragments; Opposite-Strand 5' Ends fragments);
		
		3. Perform length analysis on the split-aligned fragments by group of samples
		(1) Fragment length distribution based on Intra-cis Determinate and Inter-trans Determinate fragments;
		(2) eccDNA molecules size distribution based on Intra-cis Circular;
		
		4. Perform lengths, endpoint distances, and identical 6-mer end motifs analysis on the split-aligned fragments from each group:
		(1) Endpoint distance distribution of Same-Strand 5' Ends fragments
        (2) Endpoint distance distribution of Opposite-Strand 5' Ends fragments
        (3) Identical 6-mer end motifs profile of Same-Strand 5' Ends fragments
        (4) Identical 6-mer end motifs profile of split-aligned fragments


	Input: a file containing a list of file names which corresponding to split-aligned fragments from a group of samples  
		  For example: list_HCC
		      which contains three lines and there is a filename on each line, such as following:
				         HCC1_spc_paired.sam
				         HCC2_spc_paired.sam
				         HCC3_spc_paired.sam
				
	Output:
		1. statics_proportion_length_list_HCC; # corresponding to  function 2
		(optional)
		2. File of Fragment length distribution
		3. File of eccDNA molecules size distribution
		4. File of endpoint distance distribution of Same-Strand 5' Ends fragments
		5. File of endpoint distance distribution of Opposite-Strand 5' Ends fragments
		6. File of Identical 6-mer end motifs profile of Same-Strand 5' Ends fragments
		7. File of Identical 6-mer end motifs profile of split-aligned fragments fragments
		(optional)
		8. pos5end_HCC1_spc_paired.sam
		9. pos5end_HCC2_spc_paired.sam
		10.pos5end_HCC3_spc_paired.sam
		...
	         
=cut

use strict;
my $list=shift;
my $pos5end_output_flag=1; # 1: output the 5'end information file for each sample; 0: not output
my $distribution_output_flag=1; # 1: output the length distributions, endpoint distance distributions, and Identical 6-mer end motifs profiles for a group; 0: not output
my $ReadLenth=75;
my %split_count=();
my %read=();


my %lenth_dis=(); #Fragment length distribution based on Intra-cis Determinate and Inter-trans Determinate fragments
my %lenth_ecc=(); #eccDNA molecules size distribution based on Intra-cis Circular

my %dis_op=();  #Endpoint distance distribution of Opposite-Strand 5' Ends fragments
my %dis_same=(); #Endpoint distance distribution of Same-Strand 5' Ends fragments


my $J1=0;   #Count of fragments with One Junction
my $J2=0;   #Count of fragments with Two Junction
my $J3=0;   #Count of fragments with Three Junction
my $J=0;    #Count of fragments with Junctions

#Fragments from one chromosome with one junction
my $N_S_T=0; #Count of Total
my $N_S_1=0; #Count of Intra-cis Determinate Type 
my $N_S_2=0; #Count of Intra-cis Circular Type
my $N_S_3=0; #Count of Intra-cis Complex Type

#Fragments from two or more chromosome with one junction
my $N_T_T=0; #Count of Total
my $N_T_1=0; #Count of Inter-trans Determinate Type
my $N_T_2=0; #Count of Inter-trans Complex Type


my $N_same5end=0;   #Same-Strand 5' Ends fragments
my $N_opposite5end=0;   #Opposite-Strand 5' Ends fragments



my %K6_motif=(); #K6_motif in Same-Strand 5' Ends fragments
my %K6_motif_split=(); #K6_motif in Split-aligned fragments

my $Ends_equal_k6=0;     #Number of Identical K6_motif in Split-aligned fragments
my $Ends_equal_k6_same=0;   #Number of Identical K6_motif in Same-Strand 5' Ends fragments
my $Ends_equal_k6_op=0;      #Number of Identical K6_motif in Opposite-Strand 5' Ends fragments


my $header="";
my $N7_reads="";

my $statics_f="statics_proportion_length_".$list;
open(OUT3,">$statics_f")or die $!;
	print OUT3 "Filename\t#junction=1\t#junction=2\t#junction>2\tFragments from one chromosome\tFragments from more than one chromosomes\tIntra-cis Determinate Type\tIntra-cis Circular Type\tIntra-cis Complex Type\tInter-trans Determinate Type\tInter-trans Complex Type\tSame-Strand 5' Ends fragments\tOpposite-Strand 5' Ends fragments\tIdentical 6-mer in Same-Strand 5' Ends\tIdentical 6-mer in Opposite-Strand 5' Ends\n";

open(LI,$list)or die $!;
while(my $file=<LI>){
	chomp $file;
    %split_count=();
	%read=();
	#	%lenth_dis=();
	
	#Fragments from one chromosome with one junction
    $N_S_T=0; #Count of Total
    $N_S_1=0; #Count of Intra-cis Determinate Type 
    $N_S_2=0; #Count of Intra-cis Circular Type
    $N_S_3=0; #Count of Intra-cis Complex Type

   #Fragments from two or more chromosome with one junction
    $N_T_T=0; #Count of Total
    $N_T_1=0; #Count of Inter-trans Determinate Type
    $N_T_2=0; #Count of Inter-trans Complex Type
	
	
    $N_same5end=0; 
    $N_opposite5end=0; 
	
    $Ends_equal_k6=0;
    $Ends_equal_k6_same=0;
    $Ends_equal_k6_op=0;
   
    $header="";
    
    $N7_reads="";
	# %dis_total=();
	# %dis_op=();
	#	%dis_same=();
	print "Processing....$file\n";
	open(IN,$file)or die $!;
	while(my $line=<IN>){
		if($line=~/^\@/){
			$header.=$line;	
		}else{
			my @temp=split /\t/,$line;
			if(!exists($read{$temp[0]})){
				$split_count{$temp[0]}=1;
				$read{$temp[0]}=$line;
			}else{
				$split_count{$temp[0]}++;
				$read{$temp[0]}.=$line;
			}
		}
	}
	close(IN);
	
	
	
#	my $N7reads_f="N7reads_".$file;
#	open(OUT7,">$N7reads_f")or die $!;
    my $cluster_string="";
    my $detail_string="";
    my $end_string="";
    my @read_id=keys(%read);
	my %statics=();
	for(my $i=0;$i<@read_id;$i++){
		if(exists($statics{$split_count{$read_id[$i]}})){
	    	  $statics{$split_count{$read_id[$i]}}++;
		}else{
        	 $statics{$split_count{$read_id[$i]}}=1;
		}
		$cluster_string.=$read{$read_id[$i]};
    	my ($details,$end_inf)=&call_5endto3end($read{$read_id[$i]});
    	$detail_string.=$details;
		$end_string.=$end_inf;
	
	}
	if($pos5end_output_flag==1){
		my $pos5end_f="pos5end_".$file;
		open(OUT1,">$pos5end_f")or die $!;
		print OUT1 $end_string;
		close(OUT1);
	}																																																																										
#	print OUT7  $header.$N7_reads;
#	close(OUT7);
	
	$J1=0;
	$J2=0;
	$J3=0;
	$J=0;
	
	my @count=sort{$a<=>$b}keys(%statics);
	foreach(@count){
		$J+=$statics{$_};
		if($_==3){
			$J1+=$statics{$_};
		}elsif($_==4){
			$J2+=$statics{$_};
		}else{
		    $J3+=$statics{$_};
		}
		
	}
	
	$N_S_3=$N_S_T-$N_S_1-$N_S_2;
	my $output_string=$file;
	$output_string.="\t".($J1/$J)."\t".($J2/$J)."\t".($J3/$J)."\t".($N_S_T/$J)."\t".($N_T_T/$J)."\t".($N_S_1/$J)."\t".($N_S_2/$J)."\t".($N_S_3/$J)."\t".($N_T_1/$J)."\t".($N_T_2/$J)."\t".($N_same5end/$J)."\t".($N_opposite5end/$J)."\t".($Ends_equal_k6_same/$N_same5end)."\t".($Ends_equal_k6_op/$N_opposite5end)."\n";
	print OUT3 $output_string;
	

}
close();
if($distribution_output_flag==1){
	&print_distributions();
}




sub print_distributions{
	
	my $len_ecc_string="";
	my @len=sort{$a<=>$b}keys(%lenth_ecc);
	foreach(@len){
		$len_ecc_string.=$_."\t".$lenth_ecc{$_}."\n";
	}
	my $out="len_ecc_dis_".$list;
	open(OUTC,">$out")or die $!;
	print OUTC $len_ecc_string;
	close(OUTC);
	

	my $len_string="";
	my @len=sort{$a<=>$b}keys(%lenth_dis);
	foreach(@len){
		$len_string.=$_."\t".$lenth_dis{$_}."\n";
	}
	my $out="len_dis_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $len_string;
	close(OUTT);
	
	
	
	my $end_dis_op="";
	my @end_dis=sort{$a<=>$b}keys(%dis_op);
	foreach(@end_dis){   
   		$end_dis_op.=$_."\t".$dis_op{$_}."\n";
   		print  $_."\t".$dis_op{$_}."\n";
	}
	my $out="end_dis_op_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $end_dis_op;
	close(OUTT);
	

	my $end_dis_sd="";
	my @end_dis=sort{$a<=>$b}keys(%dis_same);
	foreach(@end_dis){   
    	$end_dis_sd.=$_."\t".$dis_same{$_}."\n"; 
	}
	my $out="end_dis_sd_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $end_dis_sd;
	close(OUTT);

	

	my $kmer="";
	my @k6mer=sort{$a<=>$b}keys(%K6_motif);
	foreach(@k6mer){   
   		$kmer.=$_."\t".$K6_motif{$_}."\n"; 
	}
	my $out="k6mer_ SameStrand5Ends_equal_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $kmer;
	close(OUTT);
	
	
	
	my $kmer="";
	my @k6mer=sort{$a<=>$b}keys(%K6_motif_split);
	foreach(@k6mer){   
   		$kmer.=$_."\t".$K6_motif_split{$_}."\n"; 
	}
	my $out="k6mer_split_aligned_equal_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $kmer;
	close(OUTT);

}



sub call_5endto3end{
	my $pair=shift;
	my $end5="";
	my %start=();
    my @temp=split /\n/, $pair;
	for(@temp){
		if($_ ne ""){
			my @temp1=split /\t/,$_;
			my $read_flag=$temp1[1] & 192;
			if(!exists($start{$read_flag})){
				$start{$read_flag}=$_;
			}else{
				$start{$read_flag}.="\n".$_;    #the alignments from the same read
			}
		}
	}
	my $read_order=0;
	my $start_inf="";
	my $read_id="";
	my @read=keys(%start);
	for(@read){
	    $read_order++;
		my @temp_fragment=split /\n/,$start{$_};
		my $length=0; #length mapped to reference
		my $length2=0; #length of segment
		if(@temp_fragment==1){    # A read has only one alignment
		    my @temp2=split /\t/,$temp_fragment[0];
		    
		    if(($temp2[1] & 4) ==4){
		    	$start_inf.=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t*\t*\t*\t$ReadLenth\t*\t$temp2[9];\n";
		    	$end5.="\tR$read_order\t$temp2[2]\t*\t*\t$temp2[9]";
		    	next;
		    }
			
			my @CIGAR=split /[A-Z|=]/,$temp2[5];
			my @CIGAR2=split /[0-9]+/,$temp2[5];
			if($CIGAR2[0] eq ""){
				shift @CIGAR2;
			}	
			$read_id=$temp2[0];
			if(@CIGAR == @CIGAR2){
				my $start_base=0;
				if((($temp2[1] & 16) ==0) && $CIGAR2[0] ne "M"){
					$start_base=$CIGAR[0];
				}
				if((($temp2[1] & 16) ==16) && $CIGAR2[-1] ne "M"){
					$start_base=$CIGAR[-1];
				}

				for (my $i=0;$i<@CIGAR2;$i++){
					if($CIGAR2[$i] eq "M"){
						$length+=$CIGAR[$i];
						$length2+=$CIGAR[$i];
					}
					if($CIGAR2[$i] eq "I"){
						$length2+=$CIGAR[$i];
					}
					if($CIGAR2[$i] eq "D"){
						$length+=$CIGAR[$i];
					}
				}
				 
				my $test=$temp2[1] & 16;
				my $binary_number = sprintf("%b", $temp2[1]);
				 
				if($test == 0){
			    	$start_inf.=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t".$temp2[3]."\t".($temp2[3]+$length-1)."\t".$start_base."\t".($length2)."\t+\t$temp2[9];\n";
					if($start_base==0){
						$end5.="\tR$read_order\t$temp2[2]\t$temp2[3]\t+\t$temp2[9]";       # The mapping of 5'end without hard clip
					}else{
						$end5.="\tR$read_order\t$temp2[2]\t-\t+\t$temp2[9]";               # The mapping of 5'end has hard clip
					}
					 
		   		}elsif($test == 16){
			   		$start_inf.=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t".($temp2[3]+$length-1)."\t".$temp2[3]."\t".$start_base."\t".($length2)."\t-\t$temp2[9];\n";
					if($start_base==0){
						$end5.="\tR$read_order\t$temp2[2]\t".($temp2[3]+$length)."\t-\t$temp2[9]";    # The mapping of 5'end without hard clip
					}else{
						$end5.="\tR$read_order\t$temp2[2]\t-\t-\t$temp2[9]";                          # The mapping of 5'end has hard clip
					}
				}else{
					print "error1...\n";
				}
				 
			}else{
					print "error...@CIGAR...@CIGAR2...$temp2[5]\n";
			}
			
		}elsif(@temp_fragment >= 2){          # A read has two alignments (split-aligned)
			my %start_frag=();
			my %start_base=();
			for my $each_frag(@temp_fragment){
			    my @temp2=split /\t/,$each_frag;
				if(($temp2[1] & 4) ==4){
					$length=length($temp2[9]);
		    		$start_inf.=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t*\t*\t*\t$length\t*\t$temp2[9];";
		    		next;
		        }
			    $length=0;
			    $length2=0;
				
                my @CIGAR2=split /[0-9]+/,$temp2[5];
			    my @CIGAR=split /[A-Z|=]/,$temp2[5];
		     	if($CIGAR2[0] eq ""){
			    	shift @CIGAR2;
			    }	
				my $strand="";
				$read_id=$temp2[0];
			   if(@CIGAR == @CIGAR2){
					$start_base{$each_frag}=0;
					my $i=0;
					my $j=0;
					
		   			if(($temp2[1] & 16) ==0){
						$i=0;
						$j=1;
						$strand="+";
               		 }elsif(($temp2[1] & 16) == 16){
						$i=-1;
						$j=-1;
						$strand="-";
					}
					my $test=$temp2[1] & 16;
                    #  my $binary_number = sprintf("%b", $temp2[1]);
					#   print "test...$binary_number\t".$test."\n";
					while (($i<0 && abs($i)<=@CIGAR2)||($i>=0 && $i<@CIGAR2)){
                        if($CIGAR2[$i] eq "M"){
                           $length+=$CIGAR[$i];
                           $length2+=$CIGAR[$i];
						   $i+=$j;
						   last;
                        }else{
	                       $start_base{$each_frag}+=$CIGAR[$i];  # record the slipped bases at the 5'end in this alignment
						}
						$i+=$j;
						 
					}
					for ($i;($i<0 && abs($i)<=@CIGAR2)||($i>=0 && $i<@CIGAR2);$i+=$j){
						if($CIGAR2[$i] eq "M"){
							$length+=$CIGAR[$i];
							$length2+=$CIGAR[$i];
						}elsif($CIGAR2[$i] eq "I"){
							$length2+=$CIGAR[$i];
                    	}elsif($CIGAR2[$i] eq "D"){
                        	$length+=$CIGAR[$i];
                    	}
						 
					}
               
					if(($temp2[1] & 16) ==0){
						$start_frag{$each_frag}=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t".$temp2[3]."\t".($temp2[3]+$length-1)."\t".$start_base{$each_frag}."\t".($length2)."\t+\t$temp2[9];";
						
					}elsif(($temp2[1] & 16) == 16){
						$start_frag{$each_frag}=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t".($temp2[3]+$length-1)."\t".$temp2[3]."\t".$start_base{$each_frag}."\t".($length2)."\t-\t$temp2[9];";
					}
				
               }else{
                        print "error...@CIGAR...@CIGAR2...$temp2[5]\n";
               }
		   }
		   
		   # two alignments of the split-aligned read are arranged following the sequencing direction 5' to 3'. 
		   my @sort_frag=sort{$start_base{$a}<=>$start_base{$b}}keys(%start_base);  
		   my $k=0;
		   for my $each_frag(@sort_frag){
			   if($k==0){
			   	   my @temp3=split /\t/,$start_frag{$each_frag};
				   my $strand="";
				   if($temp3[-2]=~/(.+)/){
				   		$strand=$1;
			   		}else{
						$strand=$temp3[-2];
					}
					#	print $strand."\n";
					$temp3[-1]=~/(.+);/;
				   if($temp3[-4]==0){
				   		$end5.="\tR".$read_order."\t".$temp3[2]."\t".$temp3[3]."\t".$strand."\t".$1;
			   		}else{
						$end5.="\tR".$read_order."\t".$temp3[2]."\t-\t".$strand."\t".$1;
					}
			   }
			   $k++;
			   $start_inf.=$start_frag{$each_frag};
	  		 }
			$start_inf.="\n";

		}
	}
	
	$end5=$read_id.$end5;
	
	my($flag,$len,$con,$cnv)=&estimate_length($start_inf);
	$end5.="\t".$flag."\t".$len."\t".$con."\t".$cnv."\n";
	if($flag==2 && $len<0){
		print $end5;
	}
	$start_inf=$start_inf.$end5;
	&ends5_analysis($end5,$pair);
	return ($start_inf,$end5);
}



sub ends5_analysis{
	my $fragment_ends=shift;
	my $pair_origin=shift;
	my @temp_end=split /\t/,$fragment_ends;
	
	# extract k-mer motif and count the identical k-mer motif
	my $flag_equal_k6=0;
	my @temp_k6=();
   	my $k6="";
	if(($temp_end[3] ne "*" && $temp_end[8] ne "*") && ($temp_end[3] ne "-") && ($temp_end[8] ne "-")){
		#### k-mer at 5'end of the former read
     	my $k6="";
           if($temp_end[4] eq "+"){
            	$k6=substr($temp_end[5],0,6);
            }else{
            	$k6=reverse(substr($temp_end[5],-6));
            	my @tempkmer=split //,$k6;
            	my $tempk6="";
            	for(@tempkmer){
            		if($_ eq "A"){
            			$tempk6.="T"
            		}elsif($_ eq "T"){
            			$tempk6.="A"
            		}elsif($_ eq "G"){
            			$tempk6.="C"
            		}elsif($_ eq "C"){
            			$tempk6.="G"
            		}else{  # for N
						$tempk6+=$_;
            			#print "error...$k6....$tempk6....$temp_end[5]...\n";
    
            	 
            		}
            	}
            	$k6=$tempk6;	
            }
   	 	   
     	    push @temp_k6,$k6;
          #### k-mer at 5'end of the last read ()
            if($temp_end[9] eq "+"){
            	$k6=substr($temp_end[10],0,6);
            }else{
            	$k6=reverse(substr($temp_end[10],-6));
            	my @tempkmer=split //,$k6;
            	my $tempk6="";
            	for(@tempkmer){
            		if($_ eq "A"){
            			$tempk6.="T"
            		}elsif($_ eq "T"){
            			$tempk6.="A"
            		}elsif($_ eq "G"){
            			$tempk6.="C"
            		}elsif($_ eq "C"){
            			$tempk6.="G"
            		}else{  # for N
						$tempk6+=$_;
            		}
            	}
            	$k6=$tempk6;	
            }
     	    push @temp_k6,$k6;
     	   
		  
           if(@temp_k6==2){
     	      my $k5_0_0=substr($temp_k6[0],0,5);
     	      my $k5_0_1=substr($temp_k6[0],-5);
     	      my $k5_1_0=substr($temp_k6[1],0,5);
     	      my $k5_1_1=substr($temp_k6[1],-5);
     	   
     	      if($temp_k6[0] eq $temp_k6[1] or $k5_0_0 eq $k5_1_1 or $k5_0_1 eq $k5_1_0){
			     $flag_equal_k6=1; 
			      for my $k6(@temp_k6){
            	    if(exists($K6_motif_split{$k6})){
                       $K6_motif_split{$k6}++;
                    }else{
                      $K6_motif_split{$k6}=1; 
                   }
           		 }  
     	      }
     	
           }
	
	}
	
	# two ends of a fragment aligned to the same chromsome, without hard clip at both 5'end
	if(($temp_end[2] eq $temp_end[7]) && ($temp_end[3] ne "*" && $temp_end[8] ne "*") && ($temp_end[3] ne "-") && ($temp_end[8] ne "-")){
		my $dis=abs($temp_end[3]-$temp_end[8]);
		if( $flag_equal_k6==1) {
			$Ends_equal_k6++;
			   
     		if($temp_end[4] eq $temp_end[9]){
     		     $Ends_equal_k6_same++;
     		 }else{
     		     $Ends_equal_k6_op++;
     		 }
		}
		########same orientation from different molecular
		if($temp_end[4] eq $temp_end[9]){
			$N_same5end++;
			if(exists($dis_same{$dis})){
                $dis_same{$dis}++;
            }else{
                $dis_same{$dis}=1; 
            }
            
			if($flag_equal_k6==1){
            	for my $k6(@temp_k6){
            	 	if(exists($K6_motif{$k6})){
                     	$K6_motif{$k6}++;
                	 }else{
                    	 $K6_motif{$k6}=1; 
                	 }
            	}
		    }
            
		}else{
			$N_opposite5end++;
			$N7_reads.=$pair_origin;
			if($temp_end[4] eq "+"){
				$dis=$temp_end[8]-$temp_end[3];
			}else{
				$dis=$temp_end[3]-$temp_end[8];
			}
			
			if(exists($dis_op{$dis})){
                $dis_op{$dis}++;
            }else{
                $dis_op{$dis}=1;
		    }
		    
		}
	}
	
	
}


sub estimate_length{
	#fragments has different aligned chromosome, return 150+
	#two or more fragments have the same chromosome, try to estimate length
	my $read_pair=shift;
	my $con_inf="";
	my $cnv_inf="";
	
	#print $read_pair;
	my @temp_read=split /\n/,$read_pair;
	my @read_1=split /;/,$temp_read[0];
	my @read_2=split /;/,$temp_read[1];
	if ($read_1[-1] eq ""){
		pop @read_1;
	}
	if ($read_2[-1] eq ""){
        pop @read_2;
    } 
    my %chr=();
    my @r1h=split /\t/,$read_1[0];
	my @r1t=split /\t/,$read_1[-1];
	my @r2t=split /\t/,$read_2[-1];
	my @r2h=split /\t/,$read_2[0];
    my $flag=0;
    my $length=0;
    my $length_con=0;
    
  	# fragments with one junction
#$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t".($temp2[3]+$length-1)."\t".$temp2[3]."\t".$start_base{$each_frag}."\t".($length2)."\t-\t$temp2[9];"
	if(@read_1+@read_2==3 && $r1h[3] ne "*" && $r1t[3] ne "*" && $r2h[3] ne "*" && $r2t[3] ne "*"){
      if(($r1h[2] eq $r2h[2] && $r1h[2] eq $r1t[2] && $r1h[2] eq $r2t[2]) && $r1h[2] ne "*" ){    # from one chromsome
    	$N_S_T++;
    	if(@read_1==2 && $r2h[-3]>=($ReadLenth-5) && ($r1h[-2] ne "*") && ($r1t[-2] ne "*") && ($r2h[-2] ne "*")){  #Read1 has two alignments R10 and R11, Read 2 has one alignment R20
    			if(($r1h[-2] ne $r1t[-2]) && $r1t[-2] ne $r2h[-2]){     # The aligned strands of R11 and R20 are opposite
    				if($r1t[-2] eq "+" && $r2h[4]>=$r1t[4]){          #The coordinates of R11 and R20 are linearly increasing or decreasing.
    					#my $min=($r2h[4]<$r1t[3])?$r2h[4]:$r1t[3];
						$length+=$ReadLenth+$r2h[3]-$r1t[4];
						$length_con=$r2h[3]-$r1t[3]+1;	
						$cnv_inf.="$r1h[2]:$r1h[4]-$r1h[3];$r1t[2]:$r1t[3]-$r2h[3];";
						$flag=1;
						$con_inf.="R11<..<R20;";
						$N_S_1++;                   #Intra-cis Determinate type   
    				}elsif($r1t[-2] eq "-" && $r2h[4]<=$r1t[4]){	   #The coordinates of R11 and R20 are linearly decreasing.
				      #  my $max=($r2h[4]<$r1t[3])?$r1t[3]:$r2h[4];
				        $length+=$ReadLenth+$r1t[4]-$r2h[3];
				        $length_con=$r1t[3]-$r2h[3]+1;
				        $cnv_inf.="$r1h[2]:$r1h[3]-$r1h[4];$r1t[2]:$r2h[3]-$r1t[3];";
				        $flag=1;
				        $con_inf.="R11>..>R20;";
				        $N_S_1++;              #Intra-cis Determinate type  
			       }	
    			}elsif(($r1h[-2] eq $r1t[-2]) && ($r1t[-2] ne $r2h[-2])){
    				if($r1t[-2] eq "+" && $r1h[3]>$r2h[3] && $r1t[3]<$r2h[3]){
    					my $min=($r2h[4]<$r1t[3])?$r2h[4]:$r1t[3];
						$length+=$r1h[4]-$min+1;	
						$cnv_inf.="$r1t[2]:$min-$r1h[4];";
						$flag=2;
						$con_inf.="R11<..R20<..<R10;";
						$N_S_2++;           #Intra-cis Circular  
    				}elsif($r1t[-2] eq "-" && $r1h[3]<$r2h[3] && $r2h[3]<$r1t[3]){	 
				        my $max=($r2h[4]<$r1t[3])?$r1t[3]:$r2h[4];
				        $length+=$max-$r1h[4]+1;
				        $cnv_inf.="$r1t[2]:$r1h[4]-$max;";
				        $flag=2;
				        $con_inf.="R10<..R20<..<R11;";
				        $N_S_2++;           #Intra-cis Circular
			       }else{
			            if($r1t[-2] eq "+" && $r2h[4]>=$r1t[4]){
			                 $length+=$ReadLenth+$r2h[3]-$r1t[4];
						     $length_con=$r2h[3]-$r1t[3]+1;	
						     $cnv_inf.="$r1h[2]:$r1h[3]-$r1h[4];$r1t[2]:$r1t[4]-$r2h[3];";
						     $flag=4;
						     $con_inf.="R11<..<R20;";
						     $N_S_1++;      #Intra-cis Determinate type
			            }elsif($r1t[-2] eq "-" && $r2h[4]<=$r1t[4]){
			                 $length+=$ReadLenth+$r1t[4]-$r2h[3];
				             $length_con=$r1t[3]-$r2h[3]+1;
				             $cnv_inf.="$r1h[2]:$r1h[4]-$r1h[3];$r1t[2]:$r2h[3]-$r1t[4];";
				             $flag=4;
				             $con_inf.="R11>..>R20;";
				             $N_S_1++;     #Intra-cis Determinate type
			            }
			       }
    			
    		}}elsif(@read_2==2 &&  $r1h[-3]>=($ReadLenth-5) && ($r2t[-2] ne "*") && ($r1h[-2] ne "*") && ($r2h[-2] ne "*") ){ #Read1 has one alignment R10, Read 2 has two alignments R20 and R21
    				if(($r2t[-2] ne $r1h[-2]) && ($r2t[-2] ne $r2h[-2])){
    				   if($r2t[-2] eq "+" && $r2t[4] <= $r1h[4]){
    				     #	my $max=($r1h[4]>$r2t[3])?$r1h[4]:$r2t[4];
						  $length+=$ReadLenth+$r1h[3]-$r2t[4];
						  $length_con=$r1h[3]-$r2t[3]+1;
						  $cnv_inf.="$r2t[2]:$r2t[3]-$r1h[3];$r2h[2]:$r2h[4]-$r2h[3];";
						  $flag=1;
						  $con_inf.="R10>..>R21;";
					      $N_S_1++;        #Intra-cis Determinate type; 
    				   }elsif($r2t[-2] eq "-" && $r2t[4]>=$r1h[4]){	 
				          #  my $min=($r1h[4]>$r2t[3])?$r2t[3]:$r1h[4];
				          $length+=$ReadLenth+$r2t[4]-$r1h[3];
				          $length_con=$r2t[3]-$r1h[3]+1;
				          $cnv_inf.="$r1t[2]:$r1h[3]-$r2t[3];$r2h[2]:$r2h[3]-$r2h[4];";
				          $flag=1;
				          $con_inf.="R10<..<R21;";
				          $N_S_1++;       #Intra-cis Determinate type; 
			          }
			          	
   		         }elsif(($r2t[-2] ne $r1h[-2]) && ($r2t[-2] eq $r2h[-2])){
   		
    				if($r2t[-2] eq "+" && $r1h[3]<$r2h[3] && $r1h[3]>$r2t[3]){
    					my $min=($r1h[4]>$r2t[3])?$r2t[3]:$r1h[4];
				        $length+=$r2h[4]-$min+1;
				        $cnv_inf.="$r1t[2]:$min-$r2h[4];";
				        $flag=2;
				        $con_inf.="R21<..R10<..<R20;";
					    $N_S_2++;       #Intra-cis Circular
    				}elsif($r2t[-2] eq "-" && $r1h[3]>$r2h[3] && $r2t[3]>$r1h[3]){	 
				        my $max=($r1h[4]>$r2t[3])?$r1h[4]:$r2t[3];
						$length+=$max-$r2h[4]+1;	
						$cnv_inf.="$r2t[2]:$r2h[4]-$max;";
						$flag=2;
						$con_inf.="R20<..R10<..<R21;";
				        $N_S_2++;       #Intra-cis Circular
			       }else{
			         if($r2t[-2] eq "+" && $r2t[4] <= $r1h[4]){
						$length+=$ReadLenth+$r1h[3]-$r2t[4];
				        $length_con=$r1h[3]-$r2t[3]+1;
				        $cnv_inf.="$r1t[2]:$r2t[3]-$r1h[3];$r2h[2]:$r2h[3]-$r2h[4];";
						$flag=4;
						 $con_inf.="R10>..>R21;";
					    $N_S_1++;      #Intra-cis Determinate type;
    				  }elsif($r2t[-2] eq "-" && $r2t[4]>=$r1h[4]){	 
				        $length+=$ReadLenth+$r2t[4]-$r1h[3];
						$length_con=$r2t[3]-$r1h[3]+1;
						$cnv_inf.="$r2t[2]:$r1h[3]-$r2t[3];$r2h[2]:$r2h[4]-$r2h[3];";
				        $flag=4;
				        $con_inf.="R10<..<R21;";
				        $N_S_1++;      #Intra-cis Determinate type;
			          }
			       
			       
			       }
    			}
            }
 
      }else{      # from two or more chromsomes
      		$N_T_T++;
    	    if(@read_1==2 && $r2h[-3]>=($ReadLenth-5) && ($r1h[2] ne $r1t[2]) &&  ($r1t[2] eq $r2h[2]) && ($r1t[-2] ne $r2h[-2])){  #Read1 has two alignments R10 and R11, Read 2 has one alignment R20; R11 and R20 from a same chromosome
    	  	       if($r1t[-2] eq "+"  && $r1t[4]<=$r2h[4]){
    					#my $min=($r2h[4]<$r1t[3])?$r2h[4]:$r1t[3];
						$length+=$ReadLenth+$r2h[3]-$r1t[4];
						$length_con=$r2h[3]-$r1t[3];
						if($r2h[-2] eq "+"){
						     $cnv_inf.="$r1h[2]:$r1h[3]-$r1h[4];";
						}elsif($r2h[-2] eq "-"){
						     $cnv_inf.="$r1h[2]:$r1h[4]-$r1h[3];";
						}
						$cnv_inf.="$r1t[2]:$r1t[3]-$r2h[3];";
						$flag=3;
						$con_inf.="R11<..<R20;";
						$N_T_1++;           #Inter-cis Determinate type;
    				}elsif($r1t[-2] eq "-" && $r1t[4]>=$r2h[4]){	 
				        #my $max=($r2h[4]<$r1t[3])?$r1t[3]:$r2h[4];
				        $length+=$ReadLenth+$r1t[4]-$r2h[3];
				        $length_con=$r1t[3]-$r2h[3];
				        if($r2h[-2] eq "+"){
						     $cnv_inf.="$r1h[2]:$r1h[3]-$r1h[4];";
						}elsif($r2h[-2] eq "-"){
						     $cnv_inf.="$r1h[2]:$r1h[4]-$r1h[3];";
						}
				        $cnv_inf.="$r1t[2]:$r2h[3]-$r1t[3];";
				        $flag=3;
				        $con_inf.="R11>..>R20;";
				        $N_T_1++;     #Inter-cis Determinate type;
			       }else{
			       	    $N_T_2++;   #Inter-trans Complex type;
			       }
    	  }elsif(@read_2==2 && $r1h[-3]>=($ReadLenth-5) && ($r2h[2] ne $r2t[2]) &&  ($r1h[2] eq $r2t[2]) && ($r1h[-2] ne $r2t[-2])){   #Read2 has two alignments R20 and R21, Read 2 has one alignment R20; R10 and R21 from a same chromosome
    	   		   
    	           if($r2t[-2] eq "+" && $r2t[4]<=$r1h[4]){
    					#my $max=($r1h[4]>$r2t[3])?$r1h[4]:$r2t[3];
						$length+=$ReadLenth+$r1h[3]-$r2t[4];
						$length_con=$r1h[3]-$r2t[3];	
						$cnv_inf.="$r2t[2]:$r2t[4]-$r1h[3];";
						if($r2h[-2] eq "+"){
						     $cnv_inf.="$r2h[2]:$r2h[4]-$r2h[3];";
						}elsif($r2h[-2] eq "-"){
						     $cnv_inf.="$r2h[2]:$r2h[3]-$r2h[4];";
						}
						$flag=3;
						$con_inf.="R10>..>R21;";
					    $N_T_1++;      #Inter-cis Determinate type;
    				}elsif($r2t[-2] eq "-" && $r2t[4]>=$r1h[4]){	 
				       # my $min=($r1h[4]>$r2t[3])?$r2t[3]:$r1h[4];
				        $length+=$ReadLenth+$r2t[4]-$r1h[3];
				        $length_con=$r2t[3]-$r1h[3];
				        $cnv_inf.="$r2t[2]:$r1h[3]-$r2t[4];";
				        if($r2h[-2] eq "+"){
						     $cnv_inf.="$r2h[2]:$r2h[4]-$r2h[3];";
						}elsif($r2h[-2] eq "-"){
						     $cnv_inf.="$r2h[2]:$r2h[3]-$r2h[4];";
						}
				        $flag=3;
				        $con_inf.="R10<..<R21;";
				        $N_T_1++;       #Inter-cis Determinate type;
			         }else{
			       	    $N_T_2++;     #Inter-trans Complex type;
			         }
			    }else{    # from two or more chromsomes
			    	$N_T_2++; 
			    	
			    }
    	   }
     
        }


    
  
   if($flag==0){
        $length=$ReadLenth*2;
        my @read_all=(@read_1,@read_2);
        for my $frag_temp(@read_all){
        	my @r1=split /\t/,$frag_temp;
        	if($r1[3]<$r1[4]){
				 $cnv_inf.="$r1[2]:$r1[3]-$r1[4];";
		 	}else{
				$cnv_inf.="$r1[2]:$r1[4]-$r1[3];";
			}
        }
    	
			$con_inf.="R10-R1".(@read_1-1)."--R20-R2".(@read_2-1).";";
   
   }else{
	  if($flag==1 or $flag ==3 or $flag==4){    #Inter-cis Determinate type and  #Intra-cis Determinate type;
       
         if(!exists($lenth_dis{$length})){     #fragment length distribution based on Inter-cis Determinate type and  #Intra-cis Determinate type
			  $lenth_dis{$length}=1;
		 }else{
			 $lenth_dis{$length}++;
		 }
	
		 if($length<$ReadLenth-5){
		 	print $read_pair;
		 }

     }elsif($flag ==2){                         #Intra-cis Circular
         if(!exists($lenth_ecc{$length})){      #eccDNA molecular size distribution based on Intra-cis Circular
		     $lenth_ecc{$length}=1;
		 }else{
			 $lenth_ecc{$length}++;
	 	 }
     }
    }
	return ($flag,$length,$con_inf,$cnv_inf);

}

