#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);
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
    
    my $outFrequency = $prefix.".groups.alleles.freq";
        
    ##- output merged frequency
    open OUTFreq, ">$outFrequency";      
    my $indexAgroup = $indexNames ->[0];
    my $indexBgroup = $indexNames ->[1];
       $indexBgroup =~ s/r1$//;
    
    print OUTFreq join("\t", "Chr", "Pos", $indexAgroup, $indexBgroup), "\n";

    for my $chr(sort keys %{$AlleleFrequency}){
        for my $pos(sort {$a<=>$b} keys %{$AlleleFrequency ->{$chr}}){


               if(exists $AlleleFrequency ->{$chr} ->{$pos} ->{$indexAgroup}){
                  my $groupAalleleinfo = $AlleleFrequency ->{$chr} ->{$pos} ->{$indexAgroup};

                  my @Ainfo = split(/;/, $groupAalleleinfo);

                  my ($refcount, $altcount, $hybridcount, $missingcount) = (0, 0, 0, 0, 0);
                  my $round = 0;
                  for my $eachtmp(@{$indexNames}){
                      next, if($eachtmp eq $indexAgroup);
                      $round += 1;
                      if(exists $AlleleFrequency ->{$chr} ->{$pos} ->{$eachtmp}){
                         my @tmp = split(/;/, $AlleleFrequency ->{$chr} ->{$pos} ->{$eachtmp});
                         my @tmpcount = split(/:/, $tmp[1]);
                            $refcount += $tmpcount[0];
                            $altcount += $tmpcount[1];
                            $hybridcount  += $tmpcount[2];
                            $missingcount += $tmpcount[3];
                      }
                  }

                  my ($refcounta, $altcounta, $hybridcounta, $missingcounta) = (int($refcount/$round), int($altcount/$round), int($hybridcount/$round), int($missingcount/$round));
                  my $groupBa = join(":", $refcounta, $altcounta, $hybridcounta, $missingcounta);
                  my $groupBinfo = $Ainfo[0].";".$groupBa;
                  print OUTFreq join("\t", $chr, $pos, $groupAalleleinfo, $groupBinfo), "\n";
               }
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
           if(-e $eachvcfFile &&  -e $pahsedVcfFile){  
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

       if(-e $outFile){
          my $phasedOut = $prefix.".$indexNames.phased.gz";
          system("python3.9  /120t/caix/src/genomics_general/VCF_processing/parseVCF.py  -i  $outFile --addRefTrack   | bgzip > $phasedOut");
          system("tabix -p vcf $phasedOut");
          
          ##- clean
          system("rm -rf $groupsFile");
  
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


    my ($groupAindex, $groupASam);
    my ($groupBindex, $groupBSam);

    my $count = 0;
    open IN0, $filein;
    while(<IN0>){
      chomp;
      $count += 1;
      die "Error:  Only surport two groups, see example: 
           Line 1: SCHN	SCHN.sam.lst; 
           Line 2: NCHN	NCHN.sam.lst.\n", if($count > 2);

      my @temp = split(/\t/, $_);

         if($count == 1){
            ($groupAindex, $groupASam) = ($temp[0], $temp[1]);
         }
         if($count == 2){
            ($groupBindex, $groupBSam) = ($temp[0], $temp[1]);
         }

    }
    close IN0;
    
    my @Asam = ();
       &sam2Array($groupASam, \@Asam);
    my @Bsam = ();
       &sam2Array($groupBSam, \@Bsam);
    
    my $mincount = scalar(@Asam); 
    my @smallGroup = @Asam;
    my $smallindex = $groupAindex;
    my @bigGroup   = @Bsam;
    my $bigindex   = $groupBindex; 

       if(scalar(@Bsam) < scalar(@Asam)){
          @smallGroup = @Bsam;
          @bigGroup   = @Asam;    
          $mincount   = scalar(@Bsam); 
          $smallindex = $groupBindex;
          $bigindex   = $groupAindex;
       }       
       
       open  OUTSamllGroup, ">$smallindex.$in1.samLst";
       print OUTSamllGroup join("\n", @smallGroup), "\n";
       close OUTSamllGroup;
       push(@{$arrayIndexName}, $smallindex);
       push(@{$groupFile},  "$smallindex.$in1.samLst");

       for(my $r=1; $r<=10; $r++){                     ##- it depends__ default: 10 
           my @rand_data = shuffle(@bigGroup); 
           my @tmp = @rand_data[1..$mincount];
           
           my $fileOut = $bigindex.".r$r.$in1.samLst";  
           open OUTX, ">$fileOut";
           print OUTX join("\n", @tmp), "\n";
           close OUTX;

           push(@{$arrayIndexName}, $bigindex."r$r");
           push(@{$groupFile}, $fileOut);

       }         
      
}


sub sam2Array {

    my ($fileIn, $Samarray) = @_;
 
    open IN0X, $fileIn;
    while(<IN0X>){
      chomp;
      my @temp = split(/\t/, $_);
      push(@{$Samarray}, $temp[0]);
    }
    close IN0X;

}


