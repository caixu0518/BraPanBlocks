#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- A06.selected.snp.vcf.gz
my $in1 = $ARGV[1]; ##- ordered sam
my $in2 = $ARGV[2]; ##- control sample id

my @ordersam = ();
   &readSam($in1, \@ordersam);

my %genotypeindex = ();
   &readVcf($in0, \%genotypeindex);

my $out = $in2.".snp.similarity.txt";
   &output($in2, \@ordersam, \%genotypeindex, $out);

##------- all subs ------------------------------------------------------------------
sub output {

    my ($controlid, $ordersam, $genotypeindex, $output) = @_;

    open OUT, ">$output";
    my $total = 0;
    my %identicalLoci = ();
    for my $chr(keys %{$genotypeindex}){
        for my $pos(keys %{$genotypeindex ->{$chr}}){
            my $controlGt = $controlid;
            die "cannot find : $chr\t$pos\t$controlid\n", if(not exists $genotypeindex ->{$chr} ->{$pos} ->{$controlid});
            $controlGt = $genotypeindex ->{$chr} ->{$pos} ->{$controlid};
            next, if($controlGt eq './.');   ##- filter missing pos

            $total += 1;

            for my $each(@{$ordersam}){
                my $eachgt = $each;
                   die "cannot find : $chr\t$pos\t$each\n", if(not exists $genotypeindex ->{$chr} ->{$pos} ->{$each});
                   $eachgt = $genotypeindex ->{$chr} ->{$pos} ->{$each};
                if($eachgt eq $controlGt){
                   $identicalLoci{$each} += 1;
                }
            }
        }
    }  
 
    for my $eachtmp(@{$ordersam}){
        print  OUT  join("\t", $eachtmp, $total, $identicalLoci{$eachtmp}, ($identicalLoci{$eachtmp}/$total)), "\n";
    }
    close OUT;

}

sub readVcf {

    my ($fileIn, $genotypeindex) = @_;

    my %index2sam = ();
    open IN1, "pigz -dc $fileIn | ";
    while(<IN1>){
      if(/^##/){
         next;
      }
      elsif(/^#CHROM/){
        chomp;
        my @title = split(/\t/, $_); 
        for(my $m=9; $m<=$#title; $m++){
            $index2sam{$m} = $title[$m];
        }
      }
      else{
        chomp;
        my @temp = split(/\t/, $_);
        for(my $m=9; $m<=$#temp; $m++){
            my $samid = $index2sam{$m};
            $genotypeindex ->{$temp[0]} ->{$temp[1]} ->{$samid} = $temp[$m];
        }
      }
    }
    close IN1;

}

sub readSam {

    my ($fileIn, $samArray) = @_;
   
    open IN0, $fileIn;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
      push(@{$samArray}, $temp[0]);
    }
    close IN0;

}
