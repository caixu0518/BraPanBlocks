#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- sample chr chrLen
my $in1 = $ARGV[1]; ##- A01

my $coordDir = "/10t/caix/Pan-genome/coords/";

my %sam2chrLen = ();
my @orderedSam = ();
   &chr2Len($in0, \@orderedSam, \%sam2chrLen);

   &calculateHaps(\@orderedSam, \%sam2chrLen, $in1);

sub calculateHaps {

    my ($orderedSam, $sam2chrLen, $chrname) = @_;

    my %sam2Diffhaps = ();

    for(my $m=0; ($m+1)<=$#orderedSam; $m++){
        my $refId = $orderedSam[$m];
        
        for(my $n=$m+1; $n<=$#orderedSam; $n++){
            my $queryId = $orderedSam[$n];
            
            my $AlignCoords = $coordDir.$refId."_".$queryId.".filtered.coords";
            if(-e $AlignCoords){

               open IN1, $AlignCoords;
               while(<IN1>){
                 chomp;
                 my @temp        = split(/\t/, $_);
                 next, if($temp[11] ne $chrname);   ##- get one-2-one align between same chromosome
                 next, if($temp[12] ne $chrname);   ##- get one-2-one align between same chromosome

                 my @refcoords   = sort {$a<=>$b}($temp[0], $temp[1]);
                 my $refMatch    = $temp[11].":".$refcoords[0]."-".$refcoords[1]."#".$refId;
                 my @queryCoords = sort {$a<=>$b}($temp[2], $temp[3]);
                 my $queryMatch  = $temp[12].":".$queryCoords[0]."-".$queryCoords[1]."#".$refId;

                    ##- record ref                 
                    if(not exists $sam2Diffhaps{$refId}){
                       $sam2Diffhaps{$refId} = $refMatch;
                    }
                    else{
                       $sam2Diffhaps{$refId} .= "\t".$refMatch;
                    }

                    ##- record query
                    if(not exists $sam2Diffhaps{$queryId}){
                       $sam2Diffhaps{$queryId} = $queryMatch;
                    }
                    else{
                       $sam2Diffhaps{$queryId} .= "\t".$queryMatch;
                    }
               }            
               close IN1;

            }
            else{
               print "ERROR: cannot find: $AlignCoords.\n";
            } 
         
        }

    }

=pop
    ##- generate aligned out
    for my $eachSam(@{$orderedSam}){

        my $fileOut = $eachSam.".out";         
        open OUTX, ">$fileOut";
        my @info = split(/\t/, $sam2Diffhaps{$eachSam});
        print OUTX join("\n", @info), "\n";
        close OUTX;
 
    }
=cut
   ##- find new aligned blocks
   for(my $m=0; $m<=$#orderedSam; $m++){
       my $currentRef = $orderedSam[$m];
       my @info = split(/\t/, $sam2Diffhaps{$currentRef});
      
       my %samregions = ();
       for my $each(@info){
           my $element = $each;
           if($element =~ /(\S+?):(\S+?)-(\S+?)#(\S+)/){
              if(not exists $samregions{$4}){
                 $samregions{$4} = $1."\t".$2."\t".$3;
              }
              else{
                 $samregions{$4} .= "\n".$1."\t".$2."\t".$3;
              }           
           }  
       }
       
       my @alignedSams = ();
       for my $eachsamTmp(@orderedSam){

           if(exists $samregions{$eachsamTmp}){
               
              push(@alignedSams, $eachsamTmp);
              my $fileOut = $currentRef.".part.".$eachsamTmp;
              open  OUTfileOut, ">$fileOut.tmp"; 
              print OUTfileOut  $samregions{$eachsamTmp}, "\n";          
              close OUTfileOut;         

              system("sort -k1,1 -k2,2n -k3,3n  $fileOut.tmp |  bedtools merge  > $fileOut");
              system("rm  $fileOut.tmp");

           }
       }
        
       ##- run bedtools subtract ...
       if(scalar(@alignedSams) >= 2){
          my %record = ();
          for(my $i=0; $i<=$#alignedSams; $i++){
              if($i == 0){
                 my $tmpFile = "$currentRef.part.$alignedSams[$i]";
               
                 open IN3, $tmpFile || die "cannot find $tmpFile"; 
                 while(<IN3>){
                   chomp;
                   my @temp = split(/\t/, $_);
                      $record{$temp[0]}{$temp[1]}{$temp[2]}{$alignedSams[$i]} = "Y";
                 }        
                 close IN3;
                  
              } 
              else{
                ##- run bedtools subtract ...
                my $currentBed = "$currentRef.part.$alignedSams[$i]";
                my $mergedBed  = "$currentBed.covered";
                system("rm $mergedBed.tmp"), if(-e "$mergedBed.tmp");
                for(my $n=0; $n<=($i-1); $n++){
                    system("cat  $currentRef.part.$alignedSams[$n]  >> $mergedBed.tmp");
                }

                system("sort -k1,1 -k2,2n -k3,3n $mergedBed.tmp | bedtools merge  > $mergedBed");
                system("bedtools subtract -a $currentBed  -b $mergedBed  > $currentBed.subtract"); 

                open IN4, "$currentBed.subtract";
                while(<IN4>){
                      chomp;
                      my @temp = split(/\t/, $_);
                         $record{$temp[0]}{$temp[1]}{$temp[2]}{$alignedSams[$i]} = "Y";
                }
                close IN4;
                system("rm  $mergedBed.tmp  $currentBed.subtract  $mergedBed");

              }
          }
          ##- records to file
          my $subtractOut = "$currentRef.filtered";
          open OUTsubtractOut, ">$subtractOut";
          for my $key1(sort keys %record){
              for my $key2(sort {$a<=>$b} keys %{$record{$key1}}){
                  for my $key3(sort {$a<=>$b} keys %{$record{$key1}{$key2}}){
                      for my $key4(sort keys %{$record{$key1}{$key2}{$key3}}){
                          print OUTsubtractOut join("\t", $key1, $key2, $key3, $key4), "\n";
                      }
                  }
              }
          }
          close OUTsubtractOut;

       }
       else{
         ##- only record reference aligned
         open IN5, "$currentRef.part.$alignedSams[0]";
         open OUT5, ">$currentRef.filtered";
         while(<IN5>){
           chomp;
           my @temp = split(/\t/, $_);
           print OUT5 join("\t", $temp[0], $temp[1], $temp[2], $alignedSams[0]), "\n";
         }
         close IN5;
         close OUT5;
       }
       ##- clean
       for my $tmp(@alignedSams){
           system("rm  $currentRef.part.$tmp");
       }  

   }

}

sub chr2Len {

    my ($in0, $orderedSam, $sam2chrLen) = @_;

    open IN0, $in0;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
         push(@{$orderedSam}, $temp[0]);
         $sam2chrLen ->{$temp[0]} = $temp[2];
    }
    close IN0;

}
