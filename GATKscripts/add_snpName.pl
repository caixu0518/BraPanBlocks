#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- BraV3.snps.miss0.5.maf0.05.vcf.gz

my $out = "BolJZS.v2.snps.maf0.05.vcf.recode.addname.vcf";

my $count = 0;
open IN0, "pigz -p 10 -dc $in0 | ";
open OUT, ">$out";
while(<IN0>){
  if(/^#/){
     print OUT $_;
  }
  else{
     chomp;
     $count += 1;
     my @temp = split(/\t/, $_);
     $temp[2] = $temp[0].":".$temp[1];
     print OUT join("\t", @temp), "\n";
  }
}
close IN0;
close OUT;

`bgzip  $out`;
`tabix -p vcf  BolJZS.v2.snps.maf0.05.vcf.recode.addname.vcf.gz`;
