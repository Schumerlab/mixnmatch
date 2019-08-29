#perl! -w

if(@ARGV<5){
    print "perl post_hmm_accuracy_shell.pl ancestry-par1.tsv ancestry-par2.tsv simulation_folder posterior_thresh path_to_Simulator_install\n"; exit;
}

my $file1=shift(@ARGV); chomp $file1;
my $file2=shift(@ARGV); chomp $file2;

my $bed_list=shift(@ARGV); chomp $bed_list;

my $pp=shift(@ARGV); chomp $pp;

my $path=shift(@ARGV); chomp $path;

my $genos="genotypes"."_"."simulation_accuracy";
system("perl $path/parsetsv_to_genotypes_v2.pl $file1 $file2 $pp $genos");

my $num_indiv=qx(ls ./$bed_list/indiv*bed | wc -l | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $num_indiv;

my $outfile="results_summary_"."$bed_list"; 
open OUT, ">results_summary_"."$bed_list";
for my $i (1..$num_indiv){

    my $indiv_name="indiv"."$i";
    my $bed="$bed_list"."/"."$indiv_name".".bed";

    #print "Rscript Determine_accuracy.R $bed $indiv_name $genos\n";
    my $current=qx(Rscript $path/Determine_accuracy.R $bed $indiv_name $genos $path); chomp $current;
    print OUT "$current\n";

}#for all individuals

system("rm accuracy_indiv*");

#run summary plots:

system("Rscript $path/plot_mixnmatch_accuracy.R $outfile");
