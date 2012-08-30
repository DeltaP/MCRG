#!/usr/bin/perl
use warnings;
use strict;

#This program performs MCRG matching for SU(3) with Nf numbers of fermions

#adds packages
use lib '..';
use lib::stat_mod qw(:DEFAULT);                                  #tells perl to use the jackknife package
use PDL;
use PDL::Fit::Polynomial;
use Math::Polynomial::Solve qw(poly_roots);
use Chart::Gnuplot;
use Proc::ParallelLoop;
my $N_core=1;                                               #finds the number of cores on the machine (check for linux)
if(`cat /proc/cpuinfo | grep processor`=~m/\d/){$N_core=`cat /proc/cpuinfo | grep processor | wc -l`;}    
if(`system_profiler SPHardwareDataType | grep Total`=~m/(\d)/){$N_core=$1;}
print"Number of Cores found:  $N_core\n";
#static variables
my %Nb;
$Nb{1224}=2;
$Nb{2448}=3;
my @Vol= (1224,2448);
my $Max_block=3;
#user input
my $Alphain=$ARGV[0];
my $N_chunks=$ARGV[1];
#dynamic variables
my (%Betavalues, %Mass, %T_start, %T_max, %T_corr, %Start)=();
my (%Data, %Avg, %Err, %B_avg, %B_err, %Full_delta_beta, %Delta_beta, %Delta_beta_avg, %Delta_beta_err, %Alpha_optimal, %Delta_beta_optimal)=();
my (@Alphavalues)=();
#-----

#gathers input from files
open (ARG, '<', 'alpha.txt')||die"$!\n";                    #reads in possible alpha values from file
chomp(my @arg=<ARG>);
close ARG;
@Alphavalues=grep(/^$Alphain/,@arg);                        #looks for a match
@Alphavalues=split/\s+/,shift(@Alphavalues);                #splits the match by spaces
shift(@Alphavalues);                                        #pulls off the first item which is used for matching
print"\n*\t*\t*\nRead in alpha values:  @Alphavalues\n";    
open (ARG, '<', 'param.txt')||die"$!\n";                    #reads in parameters from text file
chomp(@arg=<ARG>);
close ARG;
my ($gate, $m, $v)=0;
foreach my $l (@arg) {                                      #loops through param.txt
  my @line=();
  if ($l=~m/^#/){next;}
  else{
    @line=split/\s+/,$l;
    if (($line[0]=~'begin')&&(grep {$_ eq $line[1]} @Vol)){
      $m=$line[2];
      $v=$line[1];
      $Mass{$v}=$m;                                         #when there is one mass at a given volume
    }
    elsif(($line[0]=~'beta')){
      push(@{$Betavalues{$v}{$m}}, $line[1]);
      $T_start{$v}{$m}{$line[1]}=$line[2];
      $T_max{$v}{$m}{$line[1]}=$line[3];
      $T_corr{$v}{$m}{$line[1]}=$line[4];
      $Start{$v}{$m}{$line[1]}=$line[5];
    }
    else{next;}
  }
}
#-----

#main program calls
sub fetch();
sub full();
sub chunk($);
sub find_deltabeta($);
sub combine();
sub alpha_optimal();

fetch();
full();
foreach my $chunk (0..($N_chunks-1)){                       #repeates the matching looping over
  chunk($chunk);                                            #jack knife chunks
  find_deltabeta($chunk);
}
combine();
alpha_optimal();
#-----

#fetches the data and stores it in memory
sub fetch() {
  foreach my $v (keys %Mass){                             
    my $m=$Mass{$v};
    foreach my $b (@{$Betavalues{$v}{$m}}){                 #next line globs for all matching files
    my $s=$Start{$v}{$m}{$b};
    my $base="${Alphain}/mcrg_${s}_${v}_${b}_${m}";
      chomp(my @files=<$base*>);
      foreach my $f (@files){                               #loops over slelected files
        my($base, $n) = $f =~ m/(.*\.)(.*)$/;               #splits into the file name and the number
        if($n<$T_start{$v}{$m}{$b}){next;}                  #skips runs prior to thermalization
        open (IN, "<$f") || die "Can't open $f:  $!\n";     #opens the file
        chomp(my @in=<IN>);
        close IN;
        foreach my $l (@in){                                #parses the file
          if($l=~/^LOOPS/){
            my ($t0, $label, $t1, $nb, $a, $value)=split/\s+/,$l;
            push(@{$Data{$v}{$m}{$b}{$label}{$nb}{$a}},$value);
          }
          if($l=~/^POLYA ORIG/){
            my ($t0, $t1, $real, $imag, $t2, $tag)=split/\s+/,$l;
            push(@{$Data{$v}{$m}{$b}{'POLY'}{0}{0}},sqrt($real**2+$imag**2));
          }
          if($l=~/^POLYA NHYP/){
            my ($t0, $t1, $nb, $a, $real, $imag, $t2, $t3)=split/\s+/,$l;
            push(@{$Data{$v}{$m}{$b}{'POLY'}{$nb}{$a}},sqrt($real**2+$imag**2));
          }
        }
      }
    }
  }
  print"*\t*\t*\nDone\n";
}
#-----

#computes the average and standard deviation of the whole data set
sub full() {
  print"*\t*\t*\nFinding statistics on whole set\n";
  foreach my $v (keys %Mass){ 
  my $m=$Mass{$v};         
    foreach my $b (@{$Betavalues{$v}{$m}}){
      foreach my $nb (1..$Nb{$v}){
        foreach my $a (@Alphavalues){
          foreach my $label (qw/0 1 2 3 4 5 POLY/){
            if(($nb==0)&&($a!=0)){next;}                    #weeds out flags without
            if(($nb!=0)&&($a==0)){next;}                    #data
            $Avg{$v}{$m}{$b}{$label}{$nb}{$a}               #computes the average
             =stat_mod::avg(@{$Data{$v}{$m}{$b}{$label}{$nb}{$a}});
            $Err{$v}{$m}{$b}{$label}{$nb}{$a}               #computes the error
             =stat_mod::stdev(@{$Data{$v}{$m}{$b}{$label}{$nb}{$a}});
          }
        }
      }
    }
  }
  my $bm=$Mass{2448};                                       #finds delta beta for the whole set
  my $sm=$Mass{1224};
  foreach my $match (2..$Max_block){                        #loops over matchings will need to be changed for finite volume
    my $sblock=$match-1;
    my $lblock=$match;
    foreach my $alpha (@Alphavalues){                       #find delta beta for each alpha and loop
      foreach my $label (qw/0 1 2 3 4 5 POLY/){
        foreach my $bb (@{$Betavalues{2448}{$bm}}){         #loops over the large volume beta values
          my $large=$Avg{2448}{$bm}{$bb}{$label}{$lblock}{$alpha};
          my @x1=@{$Betavalues{1224}{$sm}};
          my @y1=();
          my @e1=();
          my @x2=@{$Betavalues{2448}{$bm}};
          my @y2=$Avg{2448}{$bm}{$bb}{$label}{$lblock}{$alpha};
          my @e2=$Err{2448}{$bm}{$bb}{$label}{$lblock}{$alpha};
          foreach my $sb (@{$Betavalues{1224}{$sm}}){       #loops over the small volume beta values
            push(@y1,$Avg{1224}{$sm}{$sb}{$label}{$sblock}{$alpha});
            push(@e1,$Err{1224}{$sm}{$sb}{$label}{$sblock}{$alpha});
          }
          my $x=pdl(@x1);                                   #puts data for small volumes into a piddle for fitting
          my $y=pdl(@y1);
          my $e=pdl(@e1);
          (my $fit,my $coeffs)=fitpoly1d $x, $y, 4;         #fits the small volumes
          my $a=$coeffs->at(3);                             #extracts out the coefficients
          my $b=$coeffs->at(2);
          my $c=$coeffs->at(1);
          my $d=$coeffs->at(0);
          my @roots=poly_roots($a,$b,$c,($d-$large));       #solves for the difference between the fit and the large mass value
          my $i=0;
          my @r=();
          foreach my $r (@roots){
            if ($r =~ /i$/){next;}                          #skips imaginary roots
            push(@r,$r);
          }
          if(@r==0){print"no intersection for deltabeta found\n";}
          if(@r==1){$Full_delta_beta{$match}{$bb}{$label}{$alpha}=($bb-pop(@r));}
          if(@r>2) {$Full_delta_beta{$match}{$bb}{$label}{$alpha}=($bb-stat_mod::avg(@r));}
          my $chart = Chart::Gnuplot->new(                  #Create chart object 
            output => "${Alphain}/deltabeta/${match}_${bb}_${alpha}_${label}_full.png",
            title  => "Deltabeta for matching ${match} given alpha ${alpha} for the full data set",
            xlabel => "Beta",
            ylabel => "Expectation Value",
          );
          $chart->command("set label 1 \"Delta Beta:  $Full_delta_beta{$match}{$bb}{$label}{$alpha}\"");
          $chart->command("set label 1 at graph 0.02, 0.85 tc lt 3");
          my $dataSet0 = Chart::Gnuplot::DataSet->new(      #Create dataset object for small volumes
            xdata => \@x1,
            ydata => [\@y1, \@e1],
            title => "Small Volume: ${label}",
            style => "yerrorbars",
          );
          my $dataSet1 = Chart::Gnuplot::DataSet->new(      #Create dataset object for large volume
            xdata => \@x2,
            ydata => [\@y2, \@e2],
            title => "Large Volume: ${label}",
            style => "yerrorbars",
          );
          my $dataSet2 = Chart::Gnuplot::DataSet->new(      #Create dataset object for the fit
            func => "$a*x**3+$b*x**2+$c*x+$d",
            title => "Fit: ${label}",
          );
          $chart->plot2d($dataSet0, $dataSet1, $dataSet2);  #plots the chart
        }
      }
    }
  }

  print"*\t*\t*\nDone\n";
}
#-----

#breaks the data into jackknife chunks
sub chunk($) {
  my $chunk=shift(@_);
  print"*\t*\t*\nChunking block number $chunk\n";
  foreach my $v (keys %Mass){                              
    my $m=$Mass{$v};                               
    foreach my $b (@{$Betavalues{$v}{$m}}){
      my $chunksize=int(($T_max{$v}{$m}{$b}-$T_start{$v}{$m}{$b})/$N_chunks);
      foreach my $nb (1..$Nb{$v}){
        foreach my $a (@Alphavalues){
          foreach my $label (qw/0 1 2 3 4 5 POLY/){
            if(($nb==0)&&($a!=0)){next;}                    #weeds out flags without
            if(($nb!=0)&&($a==0)){next;}                    #data
            my @temp=();                                    #temp array throws out some of the data to form Jackknife chunks
            if ($chunk==0) 
             {@temp=@{$Data{$v}{$m}{$b}{$label}{$nb}{$a}}[(($chunk+1)*$chunksize)..(@{$Data{$v}{$m}{$b}{$label}{$nb}{$a}}-1)];}
            if ($chunk>0)  
             {@temp=@{$Data{$v}{$m}{$b}{$label}{$nb}{$a}}[0..($chunk*$chunksize-1),(($chunk+1)*$chunksize)..(@{$Data{$v}{$m}{$b}{$label}{$nb}{$a}}-1)];}
            push(@{$B_avg{$v}{$m}{$b}{$label}{$nb}{$a}},stat_mod::avg(@temp));
            push(@{$B_err{$v}{$m}{$b}{$label}{$nb}{$a}},stat_mod::stdev(@temp));
          }
        }
      }
    }
  }
  print"*\t*\t*\nDone\n";
}

#fits the curves for observables and finds the difference delta beta for each chunk
sub find_deltabeta($) {
  my $chunk=shift(@_);
  print"*\t*\t*\nFinding delta beta for block number $chunk\n";
  my $bm=$Mass{2448};
  my $sm=$Mass{1224};
  foreach my $match (2..$Max_block){                        #loops over matchings will need to be changed for finite volume
    my $sblock=$match-1;
    my $lblock=$match;
    foreach my $alpha (@Alphavalues){                       #find delta beta for each alpha and loop
      foreach my $label (qw/0 1 2 3 4 5 POLY/){
        foreach my $bb (@{$Betavalues{2448}{$bm}}){         #loops over the large volume beta values
          my $large=$B_avg{2448}{$bm}{$bb}{$label}{$lblock}{$alpha}[$chunk];
          my @x1=@{$Betavalues{1224}{$sm}};
          my @y1=();
          my @e1=();
          my @x2=@{$Betavalues{2448}{$bm}};
          my @y2=$B_avg{2448}{$bm}{$bb}{$label}{$lblock}{$alpha}[$chunk];
          my @e2=$B_err{2448}{$bm}{$bb}{$label}{$lblock}{$alpha}[$chunk];
          foreach my $sb (@{$Betavalues{1224}{$sm}}){       #loops over the small volume beta values
            push(@y1,$B_avg{1224}{$sm}{$sb}{$label}{$sblock}{$alpha}[$chunk]);
            push(@e1,$B_err{1224}{$sm}{$sb}{$label}{$sblock}{$alpha}[$chunk]);
          }
          my $x=pdl(@x1);                                   #puts data for small volumes into a piddle for fitting
          my $y=pdl(@y1);
          my $e=pdl(@e1);
          (my $fit,my $coeffs)=fitpoly1d $x, $y, 4;         #fits the small volumes
          my $a=$coeffs->at(3);                             #extracts out the coefficients
          my $b=$coeffs->at(2);
          my $c=$coeffs->at(1);
          my $d=$coeffs->at(0);
          my @roots=poly_roots($a,$b,$c,($d-$large));         #solves for the difference between the fit and the large mass value
          my $i=0;
          my @r=();
          foreach my $r (@roots){
            if ($r =~ /i$/){next;}                          #skips imaginary roots
            push(@r,$r);
          }
          if(@r==0){print"no intersection for deltabeta found\n";}
          if(@r==1){push(@{$Delta_beta{$match}{$bb}{$label}{$alpha}},($bb-pop(@r)));}
          if(@r>2) {push(@{$Delta_beta{$match}{$bb}{$label}{$alpha}},($bb-stat_mod::avg(@r)));}
        }
      }
    }
  }
  print"*\t*\t*\nDone\n";
}
#-----

#jackknife average the delta beta values
sub combine(){
  print"*\t*\t*\nCalculating jackknife error for delta beta fits.\n";
  foreach my $match (2..$Max_block){                        #loops over matchings will need to be changed for finite volume
    foreach my $b (@{$Betavalues{2448}{$Mass{2448}}}){
      foreach my $a (@Alphavalues){ 
        foreach my $label (qw/0 1 2 3 4 5 POLY/){
          $Delta_beta_avg{$match}{$b}{$label}{$a}           #finds the average value of delta beta
           =stat_mod::avg(@{$Delta_beta{$match}{$b}{$label}{$a}});
          $Delta_beta_err{$match}{$b}{$label}{$a}           #finds the error on delta beta
           =stat_mod::jack_error2($Full_delta_beta{$match}{$b}{$label}{$a},@{$Delta_beta{$match}{$b}{$label}{$a}});
        }
      }
    }
  }
  print"*\t*\t*\nDone\n";
}
#-----

#finds alpha optimal
sub alpha_optimal(){
  print"*\t*\t*\nFinding alpha optimal\n";
  foreach my $b (@{$Betavalues{2448}{$Mass{2448}}}){
    foreach my $label (qw/0 1 2 3 4 5 POLY/){
      my @x=@Alphavalues;
      my @y1=();
      my @e1=();
      my @y2=();
      my @e2=();
      foreach my $a (@Alphavalues){                         #collects data into arrays
        push(@y1,$Delta_beta_avg{2}{$b}{$label}{$a});
        push(@e1,$Delta_beta_err{2}{$b}{$label}{$a});
        push(@y2,$Delta_beta_avg{3}{$b}{$label}{$a});
        push(@e2,$Delta_beta_err{3}{$b}{$label}{$a});
      }
      my $x1=pdl(@x);                                       #puts data into a piddle for fitting
      my $y1=pdl(@y1);
      my $e1=pdl(@e1);
      my $x2=pdl(@x);                                       #puts data into a piddle for fitting
      my $y2=pdl(@y2);
      my $e2=pdl(@e2);
      (my $fit1,my $coeffs1)=fitpoly1d $x1, $y1, 2;         #fits the data
      (my $fit2,my $coeffs2)=fitpoly1d $x2, $y2, 2;
      my $b1=$coeffs1->at(0);                               #extracts out the coefficients
      my $a1=$coeffs1->at(1);
      my $b2=$coeffs2->at(0);      
      my $a2=$coeffs2->at(1);

      my @roots=poly_roots($a1-$a2,$b1-$b2);                #finds intersection i.e. alpha optimal
      foreach my $r (@roots){
        if ($r =~ /i$/){next;}                              #skips imaginary roots
        elsif(($r>=$x[0])&&($r<=$x[$#x])){                  #looks for roots in the fit range
          $Alpha_optimal{$label}{$b}=$r;
          $Delta_beta_optimal{$label}{$b}=$a1*$r+$b1;
          print"$label\t$r\t$Delta_beta_optimal{$label}{$b}\n";
        }
        my $chart = Chart::Gnuplot->new(                    #Create chart object 
          output => "${Alphain}/alphaoptimal/${b}_${label}.png",
          title  => "Alpha optimal given beta ${b}",
          xlabel => "Alpha",
          ylabel => "Delta Beta",
        );
        $chart->command("set label 1 \"Alpha Optimal:  $Alpha_optimal{$label}{$b}\"");
        $chart->command("set label 1 at graph 0.02, 0.85 tc lt 3");
        my $dataSet0 = Chart::Gnuplot::DataSet->new(        #Create dataset object for small volumes
          xdata => \@x,
          ydata => [\@y1, \@e1],
          title => "Blocked Twice: ${label}",
          style => "yerrorbars",
        );
        my $dataSet1 = Chart::Gnuplot::DataSet->new(        #Create dataset object for large volume
          xdata => \@x,
          ydata => [\@y2, \@e2],
          title => "Blocked Thrice: ${label}",
          style => "yerrorbars",
        );
        my $dataSet2 = Chart::Gnuplot::DataSet->new(        #Create dataset object for the fit
          func => "$a1*x+$b1",
          title => "Fit: Blocked Twice",
        );
        my $dataSet3 = Chart::Gnuplot::DataSet->new(        #Create dataset object for the fit
          func => "$a2*x+$b2",
          title => "Fit: Blocked Thrice",
        );
        $chart->plot2d($dataSet0, $dataSet1, $dataSet2, $dataSet3);
      }
    }
  }
  print"*\t*\t*\nDone\n";
}
#-----
