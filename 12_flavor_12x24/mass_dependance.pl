#!/usr/bin/perl
use strict;
use warnings;

use stat_mod qw(:DEFAULT);
use Chart::Gnuplot;

my $start="high";                                                                             #run paramaters
my $nt=24;
my $ns=12;
my $nf=12;
my @mass=(0.005,0.01,0.015);
my @beta=(qw/2.9 3.0 3.1 3.2/);
my @alpha=(qw/0 0.400 0.500 0.600/);
my $start_t=20;
my $end_t=100;
my $stub="Out/mcrg_${start}_${ns}${nt}";

my %data=();                                                                                  #the main event
my %avg=();
my %err=();

foreach my $b (@beta){                                                                        #loops over beta values
  print"* * * * *\nCurrently Parsing Beta = $b\n* * * * *\n";
  foreach my $m (@mass){                                                                      #loops over masses
    foreach my $t ($start_t..$end_t){                                                         #selects out thermalized
      my $file="${stub}_${b}_${m}.${t}";                                                      #configurations
      open (IN, "<$file")  || die "Can't open $file: $!\n";
      chomp(my @in=<IN>);
      close IN;
      foreach my $l (@in){                                                                    #parses the output file
        if($l=~/^LOOPS/){
          my ($tag, $label, $nb, $a, $value);
          ($tag, $label, $tag, $nb, $a, $value)=split/\s+/,$l;
          push(@{$data{$b}{$m}{$label}{$nb}{$a}},$value);
        }
        if($l=~/^POLYA ORIG/){
          my ($tag, $real, $imag);
          ($tag, $tag, $real, $imag, $tag, $tag)=split/\s+/,$l;
          push(@{$data{$b}{$m}{'POLY'}{0}{0}},sqrt($real**2+$imag**2));
        }
        if($l=~/^POLYA NHYP/){
          my ($tag, $nb, $a, $real, $imag);
          ($tag, $tag, $nb, $a, $real, $imag, $tag, $tag)=split/\s+/,$l;
          push(@{$data{$b}{$m}{'POLY'}{$nb}{$a}},sqrt($real**2+$imag**2));
        }

      }
    }
    foreach my $nb (0..2){
      foreach my $a (@alpha){
        foreach my $label (qw/0 1 2 3 4 5 POLY/){
          if(($nb==0)&&($a!=0)){next;}                                                        #weeds out flags without
          if(($nb!=0)&&($a==0)){next;}                                                        #data
          $avg{$b}{$m}{$label}{$nb}{$a}=stat_mod::avg(@{$data{$b}{$m}{$label}{$nb}{$a}});     #computes the average
          $err{$b}{$m}{$label}{$nb}{$a}=stat_mod::stdev(@{$data{$b}{$m}{$label}{$nb}{$a}});   #computes the error
        }
      }
    }
  }
  foreach my $nb (0..2){
    foreach my $a (@alpha){
      if(($nb==0)&&($a!=0)){next;}
      if(($nb!=0)&&($a==0)){next;}
      # Create chart object and specify the properties of the chart
      my $chart = Chart::Gnuplot->new(
        output => "fig/${b}_${nb}_${a}.png",
        title  => "Mass dependance at beta=$b blocked $nb with paramater $a",
        xlabel => "Mass",
        ylabel => "Expectation Value",
      );

      my @y0=();                                                                              #makes arrays for y
      my @err0=();                                                                            #values and errors
      foreach my $m (@mass){
        push (@y0,$avg{$b}{$m}{0}{$nb}{$a});
        push (@err0,$err{$b}{$m}{0}{$nb}{$a});
      }
      # Create dataset object and specify the properties of the dataset
      my $dataSet0 = Chart::Gnuplot::DataSet->new(
        xdata => \@mass,
        ydata => [\@y0, \@err0],
        title => "Observable: 0",
        style => "yerrorlines",
      );

      my @y1=();
      my @err1=();
      foreach my $m (@mass){
        push (@y1,$avg{$b}{$m}{1}{$nb}{$a});
        push (@err1,$err{$b}{$m}{1}{$nb}{$a});
      }
      my $dataSet1 = Chart::Gnuplot::DataSet->new(
        xdata => \@mass,
        ydata => [\@y1, \@err1],
        title => "Observable: 1",
        style => "yerrorlines",
      );

      my @y2=();
      my @err2=();
      foreach my $m (@mass){
        push (@y2,$avg{$b}{$m}{2}{$nb}{$a});
        push (@err2,$err{$b}{$m}{2}{$nb}{$a});
      }
      my $dataSet2 = Chart::Gnuplot::DataSet->new(
        xdata => \@mass,
        ydata => [\@y2, \@err2],
        title => "Observable: 2",
        style => "yerrorlines",
      );
        
      my @y3=();
      my @err3=();
      foreach my $m (@mass){
        push (@y3,$avg{$b}{$m}{3}{$nb}{$a});
        push (@err3,$err{$b}{$m}{3}{$nb}{$a});
      }
      my $dataSet3 = Chart::Gnuplot::DataSet->new(
        xdata => \@mass,
        ydata => [\@y3, \@err3],
        title => "Observable: 3",
        style => "yerrorlines",
      );
      
      my @y4=();
      my @err4=();
      foreach my $m (@mass){
        push (@y4,$avg{$b}{$m}{4}{$nb}{$a});
        push (@err4,$err{$b}{$m}{4}{$nb}{$a});
      }
      my $dataSet4 = Chart::Gnuplot::DataSet->new(
        xdata => \@mass,
        ydata => [\@y4, \@err4],
        title => "Observable: 4",
        style => "yerrorlines",
      );

      my @y5=();
      my @err5=();
      foreach my $m (@mass){
        push (@y5,$avg{$b}{$m}{5}{$nb}{$a});
        push (@err5,$err{$b}{$m}{5}{$nb}{$a});
      }
      my $dataSet5 = Chart::Gnuplot::DataSet->new(
        xdata => \@mass,
        ydata => [\@y5, \@err5],
        title => "Observable: 5",
        style => "yerrorlines",
      );
      
      my @ypoly=();
      my @errpoly=();
      foreach my $m (@mass){
        push (@ypoly,$avg{$b}{$m}{'POLY'}{$nb}{$a});
        push (@errpoly,$err{$b}{$m}{'POLY'}{$nb}{$a});
      }
      my $dataSetPoly = Chart::Gnuplot::DataSet->new(
        xdata => \@mass,
        ydata => [\@ypoly, \@errpoly],
        title => "Observable: POLY",
        style => "yerrorlines",
      );

      # Plot many data sets on a single chart
      $chart->plot2d($dataSet0, $dataSet1, $dataSet2, $dataSet3, $dataSet4, $dataSet5, $dataSetPoly);
    }
  }
}
