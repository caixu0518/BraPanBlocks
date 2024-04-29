#!/usr/bin/perl -w
use strict;

my $in0    = $ARGV[0]; ##- ref.fasta
my $in1    = $ARGV[1]; ##- query.fasta
my $prefix = $ARGV[2]; ##- out prefix

my $refName   =  $in0;
   $refName   =~ s/\.gz$//;
my $queryName =  $in1;
   $queryName =~ s/\.gz$//;
  
   if(not -e $refName){
      &fastaExtract($in0, $refName);
   }
   if(not -e $queryName){
      &fastaExtract($in1, $queryName);
   }

   ##- nucmer
   system("/120t/caix/src/mummer-4.0.0beta2/nucmer -t 20  --mum  $refName  $queryName  --prefix  $prefix");

   system("/120t/caix/src/mummer-4.0.0beta2/delta-filter -l 10000 -r -q   $prefix.delta   >  $prefix.filtered_L10kb_rq.delta");

   system("/120t/caix/src/mummer-4.0.0beta2/show-coords  -TrHcl  $prefix.filtered_L10kb_rq.delta   >  $prefix.filtered.coords");

   system("rm  $prefix.delta  $prefix.filtered_L10kb_rq.delta");

   
sub fastaExtract {

    my ($fasta, $fastaFm) = @_;
 
    if($fasta =~ /gz$/){
       `pigz -dc $fasta > $fastaFm`;
    }

}
