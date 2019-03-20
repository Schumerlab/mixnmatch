#perl! -w

if(@ARGV<3){

    print "perl generate_reads_hybrid_genomes_v3.pl folder num_indiv num_reads\n"; exit;

}

my $folder=shift(@ARGV); chomp $folder;

my $num_indiv=shift(@ARGV); chomp $num_indiv;
$num_indiv=$num_indiv-1;

my $num_reads=shift(@ARGV); chomp $num_reads;

for my $i (0..$num_indiv){

    my $hap1="./$folder"."/"."$i"."_genome_hap1.txt";
    my $hap2="./$folder"."/"."$i"."_genome_hap2.txt";


    my $id1read1="indiv"."$i"."hap1_read1.fq";
    my $id1read2="indiv"."$i"."hap1_read2.fq";

    my $id2read1="indiv"."$i"."hap2_read1.fq";
    my $id2read2="indiv"."$i"."hap2_read2.fq";

    my $seed=rand(1000);
    system("/home/groups/schumer/shared_bin/wgsim -N$num_reads -X0.85 -1150 -2150 -d0 -S$seed -e0.01 -r0.0012 -R1 $hap1 $id1read1 $id1read2");

    my $seed=rand(1000);
    system("/home/groups/schumer/shared_bin/wgsim -N$num_reads -X0.85 -1150 -2150 -d0 -S$seed -e0.01 -r0.0012 -R1 $hap2 $id2read1 $id2read2");

    system("perl -pi -e 's/5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/g' $id1read1 $id1read2 $id2read1 $id2read2");

    my $combined1="$folder"."/"."indiv"."$i"."combined_read1.fq.gz";
    my $combined2="$folder"."/"."indiv"."$i"."combined_read2.fq.gz";
    system("cat $id1read1 $id2read1 | gzip > $combined1");
    system("cat $id1read2 $id2read2 | gzip > $combined2");

    system("rm $id1read1 $id1read2 $id2read1 $id2read2");

}
