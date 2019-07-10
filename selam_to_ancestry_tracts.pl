#perl! -w

use Math::Round;
use Math::Random;

if(@ARGV<8){
    print "perl selam_to_ancestry_tracts.pl selam_tracts use_macs_0_1 genome1 genome2 simulation_id num_macs_haps_par1 num_macs_haps_par2 use_map map_or_global_rate program_path\n";
}#print usage

my $selam=shift(@ARGV); chomp $selam;
open TRACTS, $selam or die "cannot open SELAM tracts\n";

my $use_macs=shift(@ARGV); chomp $use_macs;

my $genome1=shift(@ARGV); chomp $genome1;

my $genome2=shift(@ARGV); chomp $genome2;

my $id=shift(@ARGV); chomp $id;

my $par1total=shift(@ARGV); chomp $par1total;
$par1total=$par1total-1;

my $par2total=shift(@ARGV); chomp $par2total;
$par2total=$par2total;

my $use_map=shift(@ARGV); chomp $use_map;

my $map_or_rate=shift(@ARGV); chomp $map_or_rate;

my $rec_map_file=shift(@ARGV); chomp $rec_map_file;

my $program_path=shift(@ARGV); chomp $program_path;

if($use_map eq 0){$map_or_rate=$map_or_rate;}

my $bothpar=$par1total+$par2total;

my $counter=0;  my $chr_length=""; my $chrom_name=""; my $morganl="";

$chrom_name=qx(head -n 1 $genome1 | perl -p -e 's/>//g'); chomp $chrom_name;
$chr_length=qx(cat $genome1 | tail -n +2 | perl -p -e 's/\n//g' | wc -c); chomp $chr_length;

$morganl=($chr_length/1000)*$map_or_rate;

print "length in morgans is $morganl\n";

if($morganl<0.5){
    print "WARNING: length of the chromosome in Morgans is less than 0.5 ("."$morganl".")\n";
}#print warning

open OUT, ">indiv"."$id".".fa";
open LOG, ">indiv"."$id"."_log";

#now read in SELAM tracts for haplotype 1 and 2 for all individuals
#then generate tracts depending on whether macs or genomes are being used
my $seq=""; my $par1=""; my $par2="";
my $start_bed=1; my $stop_bed=1; my $hap_counter=1; my $indiv_counter=1;
open BED1, ">indiv"."$id"."_hap1.bed";
open BED2, ">indiv"."$id"."_hap2.bed";

my $last_start=1; my $flag=0; my $doubleflag=0; my $start=""; my $stop="";
while(my $line2=<TRACTS>){
    chomp $line2;
    #print "$line2\n";
    my @elements=split(/\t/,$line2);

    #data line
    if(scalar(@elements) eq 9){

    my $par_curr=$elements[6]; chomp $par_curr;
    my $start_raw=$elements[7]; chomp $start_raw;
    my $stop_raw=$elements[8]; chomp $stop_raw;

    if($start_raw == 0){$flag=0; $start_raw_prev=$start_raw; $stop_raw_prev=$stop_raw; $stop=1; $last_start=1}
    if($flag == 1){$stop_raw=$stop_raw_prev; $stop=$chr_length+1;}

    #convert from coordinates in morgans to tracts in basepairs 
    if($use_map eq 0){
	$start=round($start_raw*(1/$map_or_rate)*1000)+1; if($start eq 0){$start=1;}
	$stop=round($stop_raw*(1/$map_or_rate)*1000);
    }#use global rate to convert
    if(($use_map eq 1) && ($hap_counter le 2)){
		
	if(($start_raw_prev eq 0) && ($stop_raw_prev eq 0)){$stop_raw=0;}
	#print "$start_raw_prev\t$stop_raw_prev\t$stop_raw\n";

	if(($flag eq 0) && ($stop_raw != 0) && ($stop < $chr_length)){
	my $results=qx(Rscript $program_path/morgans_to_bp_recmap.R $rec_map_file $start_raw $stop_raw $map_or_rate $chr_length $last_start $start_raw_prev); chomp $results;
        my @mapinfo=split(/\t/,$results);
	#print "Rscript morgans_to_bp_recmap.R $rec_map_file $start_raw $stop_raw $map_or_rate $chr_length $last_start $start_raw_prev\n";
        $start=$mapinfo[0]; $stop=$mapinfo[1]; $flag=$mapinfo[2];
	$start_raw_prev=$stop_raw; $stop_raw_prev=$stop_raw;
	$last_start=$stop+1;
	} elsif($stop le $chr_length){
	    $stop=$chr_length+1; $start_raw_prev=0; $flag=1; $stop_raw_prev=0; $last_start=$stop;
	} #last window must finish the chromosome

	#print "morgans $start_raw\t$stop_raw to bp $start\t$stop\t$hap_counter\t$flag\n";

    }#use user-provided recombination map

    if(($use_macs eq 0) && ($hap_counter le 2)){

    if($par_curr eq 0){
	my $raw_seq=qx(fastahack $genome1 -r $chrom_name:$start..$stop); chomp $raw_seq;
	$seq="$seq"."$raw_seq";
	if($hap_counter eq 1){
	print BED1 "$chrom_name\t$start\t$stop\tpar1\n";
	} else{
	    print BED2 "$chrom_name\t$start\t$stop\tpar1\n";
	}#print to correct file
    }#collect par1 sequence
    if($par_curr eq 1){
	my $raw_seq=qx(fastahack $genome2 -r $chrom_name:$start..$stop); chomp $raw_seq;
	$seq="$seq"."$raw_seq";
	if($hap_counter eq 1){
	    print BED1 "$chrom_name\t$start\t$stop\tpar2\n";
	} else{
            print BED2 "$chrom_name\t$start\t$stop\tpar2\n";
        }#print to correct file  
    }#collect par2 seq

    } elsif(($use_macs eq 1) && ($hap_counter le 2) && ($stop_raw_prev != 0) && ($start < $chr_length) && ($stop <= $chr_length)){

	if($par_curr eq 0){
	    my $par1hap=random_uniform_integer(1,0,$par1total);
	    my $par1_chr=qx(grep -w $par1hap macs_simulation_results_trees.phy); chomp $par1_hap;
	    print LOG "using parent1 hap $par1hap for $hap_counter from $start to $stop\n";
	    $par1_chr=~ s/ +/\t/g;

	    my @par1elements=split(/\t/,$par1_chr);

	    my $tmp_par1="parent1_indiv"."$id".".fa";

	    open CHR, ">"."$tmp_par1";
	    print CHR ">chr1\n$par1elements[1]\n";

	    my $raw_seq=qx(fastahack $tmp_par1 -r $chrom_name:$start..$stop); chomp $raw_seq;
	    $seq="$seq"."$raw_seq";

	    if($hap_counter eq 1){
		print BED1 "chr1\t$start\t$stop\tpar1\n";
	    } else{
		print BED2 "chr1\t$start\t$stop\tpar1\n";
	    }#print to correct file  
 
	}#collect par1 sequence                                                                                                                                                                                                                                                                                       
	if($par_curr eq 1){

	    my $par2hap=random_uniform_integer(1,$par2total,$bothpar);
	    my $par2_chr=qx(grep -w $par2hap macs_simulation_results_trees.phy); chomp $par2_hap;
	    print LOG "using parent2 hap $par2hap for $hap_counter from $start to $stop\n";
	    $par2_chr=~ s/ +/\t/g;

	    my @par2elements=split(/\t/,$par2_chr);

	    my $tmp_par2="parent2_indiv"."$id".".fa";

	    open CHR, ">"."$tmp_par2";
	    print CHR ">chr1\n$par2elements[1]\n";

	    my $raw_seq=qx(fastahack $tmp_par2 -r $chrom_name:$start..$stop); chomp $raw_seq;
	    $seq="$seq"."$raw_seq";

	    if($hap_counter eq 1){
                print BED1 "chr1\t$start\t$stop\tpar2\n";
            } else{
                print BED2 "chr1\t$start\t$stop\tpar2\n";
            }#print to correct file  
	}#collect par2 seq                                         

    }#use macs or don't


    }#don't reset counter, same chrom and hap
    else{
    #transition line

	if(($hap_counter eq 1) or ($hap_counter <= 2)){
	print OUT ">"."indiv"."$id"."_hap"."$hap_counter\n";
	print OUT "$seq\n";
	}#only print non-sex chromosomes

	$seq=""; 

	$hap_counter++;

    }#print individual and reset counter

  }#read in all the ancestry tract info                                                                  

my $rm1="parent1_indiv"."$id".".fa"."*";
my $rm2="parent2_indiv"."$id".".fa"."*";

system("rm $rm1 $rm2");
