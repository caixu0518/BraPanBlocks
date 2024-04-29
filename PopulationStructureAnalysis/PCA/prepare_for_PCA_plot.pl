#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- pca.results.eigenvec
my $in1 = $ARGV[1]; ##- accessions.groups.lst
my $out = "PCA.results.format.lst";


my %sam2pc1 = ();
my %sam2pc2 = ();
   &readPCAresults($in0, \%sam2pc1, \%sam2pc2);

   &output($in1, \%sam2pc1, \%sam2pc2, $out);

sub output {

    my ($group, $sam2pc1, $sam2pc2, $output) = @_;

    open IN1, $group;
    <IN1>;
    open OUT, ">$output";
    print OUT join("\t", "Sample", "Group", "pca1", "pca2"), "\n";
    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_) ;
      print OUT join("\t", @temp, $sam2pc1 ->{$temp[0]}, $sam2pc2 ->{$temp[0]}), "\n";
    }
    close OUT;
    close IN1;

}

sub readPCAresults {

    my ($PCAfile, $sam2pc1, $sam2pc2) = @_;

    open IN0, $PCAfile;
    while(<IN0>){
      chomp;
      my @temp = split(/\s+/, $_);
         $sam2pc1 ->{$temp[1]} = $temp[2];
         $sam2pc2 ->{$temp[1]} = $temp[3];
    }
    close IN0;

}
