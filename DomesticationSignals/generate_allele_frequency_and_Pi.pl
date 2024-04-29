#!/usr/bin/perl -w
use strict;
use threads;

my $in0 = $ARGV[0]; ##- a file contain diff group informations
my $in1 = $ARGV[1]; ##- a vcf file

my @indexNames = ();
my @groupsFile = ();
   &readGroups($in0, \@indexNames, \@groupsFile);

my %AlleleFrequency = ();
   &generatevcfs(\@indexNames, \@groupsFile, $in1, \%AlleleFrequency);

   &output(\@indexNames, $in1, \%AlleleFrequency);

##-------------------------------------------------------------------------------
sub output {
  
    my ($indexNames, $vcfFile, $AlleleFrequency) = @_;

    my $prefix = $vcfFile;
       $prefix    =~ s/\.vcf\.gz//;
    
    my $outPi        = $prefix.".groups.pi";
    my $outFrequency = $prefix.".groups.alleles.freq";
        
    ##- merge Pi
    my %Pirecord = ();
    for my $eachgroup(@{$indexNames}){
        my $PiFile = $prefix.".".$eachgroup.".sites.pi";
   
        if(-e $PiFile){ 
           #die "cannot find Pi file for the group $eachgroup: $PiFile.\n", if(not -e $PiFile);
           open IN5, $PiFile;
           <IN5>;
           while(<IN5>){
             chomp;
             my @temp = split(/\t/, $_);
                $Pirecord{$temp[0]}{$temp[1]}{$eachgroup} = $temp[2];
           }
           close IN5;  
           system("rm $PiFile -rf"), if(-e "$PiFile");
        }
    } 
   
    ##- output merged Pi   
    open OUTPI, ">$outPi";
    print OUTPI join("\t", "Chr", "Pos", @{$indexNames}), "\n";
    for my $chr (sort keys %Pirecord){
        for my $pos(sort {$a<=>$b} keys %{$Pirecord{$chr}}){
            my @pivalues = ();
            my $tmppi    = "NaN";
            for my $each(@{$indexNames}){
                if(exists $Pirecord{$chr}{$pos}{$each}){
                   $tmppi = $Pirecord{$chr}{$pos}{$each};
                }
                push(@pivalues, $tmppi);
            }
            print OUTPI join("\t", $chr, $pos, @pivalues), "\n"; 
        }
    }
    close OUTPI;

    ##- output merged frequency
    open OUTFreq, ">$outFrequency";       
    print OUTFreq join("\t", "Chr", "Pos", @{$indexNames}), "\n";
    for my $chr(sort keys %{$AlleleFrequency}){
        for my $pos(sort {$a<=>$b} keys %{$AlleleFrequency ->{$chr}}){
            my @info = ();
            my $tmpInfo = "NaN";
            for my $eachtmp(@{$indexNames}){
                if(exists $AlleleFrequency ->{$chr} ->{$pos} ->{$eachtmp}){
                   $tmpInfo = $AlleleFrequency ->{$chr} ->{$pos} ->{$eachtmp};
                }
                push(@info, $tmpInfo); 
            }
            print OUTFreq join("\t", $chr, $pos, @info), "\n"; 
        }       
    }
    close OUTFreq;

}

sub generatevcfs {

    my ($indexNames, $groupsFile, $vcfFile, $AlleleFrequency) = @_;
    
    my $count = scalar(@{$indexNames});

    my @thr = ();
    my $i;
    for($i=0;$i<=($count-1);$i+=1){
        $thr[$i] = threads->create(\&generatesubVcf, $vcfFile, $indexNames ->[$i], $groupsFile ->[$i]);

    }
    for($i=0;$i<=($count-1);$i+=1){
        $thr[$i]->join;
    } 
   
    ##- record frequencies
    my $prefix    =  $vcfFile;
       $prefix    =~ s/\.vcf\.gz//;   

    for my $each(@{$indexNames}){
        
        my $eachvcfFile       = $prefix.".$each.gz"; 
        my $pahsedVcfFile     =  $prefix.".$each.phased.gz";;

        my %Posindex = (); 
           if(-e $eachvcfFile){  
              &recordPos($eachvcfFile, \%Posindex);
              &calAlleleFrequency($pahsedVcfFile, \%Posindex, $each, $AlleleFrequency);      
           }
           %Posindex = ();  
        
        ##- clean
        system("rm -rf $eachvcfFile  $eachvcfFile.tbi"); 
        system("rm -rf $pahsedVcfFile  $pahsedVcfFile.tbi"); 
    } 

}

sub generatesubVcf {

    my ($vcfFile, $indexNames, $groupsFile) = @_;
   
    my $prefix    =  $vcfFile;
       $prefix    =~ s/\.vcf\.gz//;
    my $outFile   =  $prefix.".$indexNames.gz";
       system("vcftools --gzvcf  $vcfFile  --recode --recode-INFO-all --stdout  --keep  $groupsFile   | bgzip  > $outFile");
       system("tabix -p vcf $outFile");
       ##- calculate Pi
       if(-e $outFile){
          my $Piout  = $prefix.".".$indexNames;
          system("vcftools --gzvcf  $outFile  --site-pi   --out  $Piout");
          system("rm  $Piout.log"), if(-e "$Piout.log");

          my $phasedOut = $prefix.".$indexNames.phased.gz";
          system("python3.9  /120t/caix/src/genomics_general/VCF_processing/parseVCF.py  -i  $outFile --addRefTrack   | bgzip > $phasedOut");
          system("tabix -p vcf $phasedOut");
       }
}

sub calAlleleFrequency {

    my ($phasedVcf, $posindex, $indexName, $alleleIndex) = @_;

    open IN3, "pigz -dc $phasedVcf | ";
    my $firstLine = <IN3>;
       chomp($firstLine);
    my @firstLine = split(/\t/, $firstLine);

    while(<IN3>){
      chomp;
      my @temp = split(/\t/, $_);
      my ($refAllele, $altAllele) = split(/\t/, $posindex ->{$temp[0]} ->{$temp[1]});
       
      my $refCount = 0;
      my $altCout  = 0;
      my $heter    = 0;
      my $missLoci = 0;
             
      for(my $m=3; $m<=$#temp; $m++){
          my $info = $temp[$m];
             $info =~ s/\///;
             $info =~ s/\|//;
             if($info eq $refAllele.$refAllele){
                $refCount += 1;
             }
             elsif($info eq $altAllele.$altAllele){
                $altCout  += 1;
             }
             elsif(($info eq $refAllele.$altAllele) || ($info eq $altAllele.$refAllele)){
                $heter    += 1;
             }
             elsif($info eq "NN"){
                $missLoci += 1;
             }
             else{
                print "warining(s): cannot calculate such pahsed genotype: $info.\n";
             }
      } 
      $alleleIndex ->{$temp[0]} ->{$temp[1]}  ->{$indexName} = "$refAllele/$altAllele;$refCount:$altCout:$heter:$missLoci";     
    }
    close IN3;

}

sub recordPos {

    my ($vcfIn, $posindex) = @_;

    open IN2, "pigz -dc $vcfIn | ";
    while(<IN2>){
      next, if(/^#/);
      chomp;
      my @temp = split(/\t/, $_);
         $posindex ->{$temp[0]} ->{$temp[1]} = $temp[3]."\t".$temp[4];
    }
    close IN2;

}

sub readGroups {

    my ($filein, $arrayIndexName, $groupFile) = @_;

    open IN0, $filein;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
         push(@{$arrayIndexName}, $temp[0]);
         push(@{$groupFile}, $temp[1]);
    }
    close IN0;

}
