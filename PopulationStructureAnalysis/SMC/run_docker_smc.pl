#!/usr/bin/perl -w
use strict;

my $in0         = $ARGV[0]; ##- 21SBR.JZS.snp.miss01MAF005.plink.vcf.gz
my $sampleArray = $ARGV[1]; ##- CC:Chiifu,CCB,21SB0218,21SB0132
my $prefix      = $ARGV[2]; ##- CC 


##- create vcf2smc folder
system("rm  $prefix"), if(-d  $prefix);
system("mkdir $prefix");
my $vcf2smcCMD = $prefix.".vcf2smc.cmds";
my @chrs = ("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10");

open OUT, ">$vcf2smcCMD";
for my $eachChr(@chrs){
    print OUT "docker   run --rm -v \$PWD:/mnt terhorst/smcpp:latest   vcf2smc  $in0   $prefix/$eachChr.smc.gz  $eachChr  $sampleArray", "\n";
}
close OUT;
system("/10t/caix/Bra_2023/Bra_ECD04/snpsMAF005/SMC/ParaFly -c $vcf2smcCMD  -CPU 10");

##- analysis
my $analysisFolder = $prefix."_analysis";
system("rm -rf  $analysisFolder"), if(-d  $analysisFolder);
system("mkdir $analysisFolder");
system("docker   run --rm -v \$PWD:/mnt terhorst/smcpp:latest   estimate   -o  $analysisFolder  9e-9  $prefix/*.smc.gz");

##- plot Ne
system("docker   run --rm -v \$PWD:/mnt terhorst/smcpp:latest   plot  plot.$prefix.pdf  $analysisFolder/model.final.json");
