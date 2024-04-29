#!/usr/bin/perl -w
use strict;
use Cwd; 

my $in0 = $ARGV[0]; ##- Z1v2.bed
my $in1 = $ARGV[1]; ##- gvcfs.fm.lst


my $dir = getcwd;
my $gatk = "/data/mg1/caix/src/NGSpipe/gatk";

system("rm -rf tmpDir"), if(-d "tmpDir");
system("mkdir tmpDir");

system("rm -rf GATKdatabase"), if(-d "GATKdatabase");
system("mkdir GATKdatabase");

system("rm -rf LogDir"), if(-d "LogDir");
system("mkdir LogDir");

my $cmd = "batch.DBImport.cmds";

open OUT, ">$cmd";
open IN0, $in0;
while(<IN0>){
  chomp;
  my @temp = split(/\t/, $_);
  my ($chr, $start, $end) = ($temp[0], $temp[1], $temp[2]);

  my $interval = $chr.":".$start."-".$end;
  
  my $eachDtabase = $interval.".database";
  my $eachtmpDir = $interval.".tmp";   

  system("rm -rf $dir/tmpDir/$eachtmpDir"), if(-d "$dir/tmpDir/$eachtmpDir");
  system("mkdir $dir/tmpDir/$eachtmpDir");

  my $cmdStrings = "$gatk --java-options -Xmx120G  GenomicsDBImport -R  $dir/Linked.path.fasta  -L  $interval  --sample-name-map   $in1  --genomicsdb-workspace-path  $dir/GATKdatabase/$eachDtabase  --tmp-dir  $dir/tmpDir/$eachtmpDir  --reader-threads  20  --max-num-intervals-to-import-in-parallel 30";
 
  print OUT $cmdStrings, "   &>  LogDir/$interval.log", "\n"; 

}
close IN0;
close OUT;
