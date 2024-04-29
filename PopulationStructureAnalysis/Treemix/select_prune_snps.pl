#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- BraV3.snps.miss0.5.maf0.05.addName.ftLD.prune.in
my $in1 = $ARGV[1]; ##- BraV3.snps.miss0.5.maf0.05.addName.vcf.gz

my $out = "ECD04.snps.maf0.05.fm.prune.vcf";

my %pruneSNP = ();
open IN0, $in0;
while(<IN0>){
  chomp;
  my $pos = $_;
     $pruneSNP{$pos} = "Y"; 
}
close IN0;

open IN1, "pigz -dc $in1 | ";
open OUT, ">$out";
while(<IN1>){
  if(/^##/){
     print OUT $_;
  }
  elsif(/^#CHROM/){
     chomp;
     my @temp = split(/\t/, $_);
     for(my $m=0; $m<=$#temp; $m++){
         if($temp[$m] eq "A03_D1607281"){
            $temp[$m] = "A03";
         }
     }
     print OUT join("\t", @temp), "\n"; 
  }
  else{
     chomp;
     my @temp = split(/\t/, $_);
     if(exists $pruneSNP{$temp[2]}){
        print OUT $_, "\n";
     }
  }
}
close IN1;
close OUT;

system("bgzip  $out");
system("tabix -p vcf $out.gz");
