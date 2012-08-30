#!/usr/bin/perl
use warnings;
use strict;

#This program performs MCRG matching for SU(3) with Nf numbers of fermions

#adds packages
use lib '..';
#use lib '/Users/gregpetrop/perl5/lib/perl5/darwin-thread-multi-2level';
use lib '/Users/gregpetrop/perl5/lib/perl5';
use lib::stat_mod qw(:DEFAULT);                                  #tells perl to use the jackknife package
use Chart::Gnuplot;
#static variables
my @Vol= (612,1224,2448);
my %Nb=(
  612  => '1',
  1224 => '2',
  2448 => '3',
);
#user input
my $Alphain=$ARGV[0];
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

fetch();
full();
#-----

#fetches the data and stores it in memory
sub fetch() {
  print"*\t*\t*\nReading in data from files\n";
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
  print"*\t*\t*\nAveraging Data\n";
  foreach my $v (keys %Mass){ 
  my $m=$Mass{$v};         
    foreach my $b (@{$Betavalues{$v}{$m}}){
      foreach my $nb (0..$Nb{$v}){
        foreach my $a (@Alphavalues){
          foreach my $label (qw/0 1 2 3 4 5/){
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
  print"*\t*\t*\nDone\n";
  print"*\t*\t*\nPlotting Data\n";
  foreach my $alpha (@Alphavalues){                       #find delta beta for each alpha and loop
    foreach my $label (qw/0 1 2 3 4 5/){
      foreach my $bb (@{$Betavalues{2448}{$Mass{2448}}}){         #loops over the large volume beta values
        my @x1=@{$Betavalues{612}{$Mass{612}}};
        my @y1=();
        my @e1=();
        my @x2=@{$Betavalues{612}{$Mass{612}}};
        my @y2=();
        my @e2=();
        my @x3=@{$Betavalues{1224}{$Mass{1224}}};
        my @y3=();
        my @e3=();
        my @x4=@{$Betavalues{1224}{$Mass{1224}}};
        my @y4=();
        my @e4=();
        my @x5=@{$Betavalues{1224}{$Mass{1224}}};
        my @y5=();
        my @e5=();
        my @x6=@{$Betavalues{2448}{$Mass{2448}}};
        my @y6=();
        my @e6=();
        my @x7=@{$Betavalues{2448}{$Mass{2448}}};
        my @y7=();
        my @e7=();
        my @x8=@{$Betavalues{2448}{$Mass{2448}}};
        my @y8=();
        my @e8=();
        my @x9=@{$Betavalues{2448}{$Mass{2448}}};
        my @y9=();
        my @e9=();
        foreach my $b (@{$Betavalues{612}{$Mass{612}}}){        #loops over the small volume beta values
          push(@y1,$Avg{612}{$Mass{612}}{$b}{$label}{0}{$alpha});
          push(@e1,$Err{612}{$Mass{612}}{$b}{$label}{0}{$alpha});
        }
        foreach my $b (@{$Betavalues{612}{$Mass{612}}}){        #loops over the small volume beta values
          push(@y2,$Avg{612}{$Mass{612}}{$b}{$label}{1}{$alpha});
          push(@e2,$Err{612}{$Mass{612}}{$b}{$label}{1}{$alpha});
        }
        foreach my $b (@{$Betavalues{1224}{$Mass{1224}}}){        #loops over the small volume beta values
          push(@y3,$Avg{1224}{$Mass{1224}}{$b}{$label}{0}{$alpha});
          push(@e3,$Err{1224}{$Mass{1224}}{$b}{$label}{0}{$alpha});
        }
        foreach my $b (@{$Betavalues{1224}{$Mass{1224}}}){        #loops over the small volume beta values
          push(@y4,$Avg{1224}{$Mass{1224}}{$b}{$label}{1}{$alpha});
          push(@e4,$Err{1224}{$Mass{1224}}{$b}{$label}{1}{$alpha});
        }
        foreach my $b (@{$Betavalues{1224}{$Mass{1224}}}){        #loops over the small volume beta values
          push(@y5,$Avg{1224}{$Mass{1224}}{$b}{$label}{2}{$alpha});
          push(@e5,$Err{1224}{$Mass{1224}}{$b}{$label}{2}{$alpha});
        }
        foreach my $b (@{$Betavalues{2448}{$Mass{2448}}}){        #loops over the small volume beta values
          push(@y6,$Avg{2448}{$Mass{2448}}{$b}{$label}{0}{$alpha});
          push(@e6,$Err{2448}{$Mass{2448}}{$b}{$label}{0}{$alpha});
        }
        foreach my $b (@{$Betavalues{2448}{$Mass{2448}}}){        #loops over the small volume beta values
          push(@y7,$Avg{2448}{$Mass{2448}}{$b}{$label}{1}{$alpha});
          push(@e7,$Err{2448}{$Mass{2448}}{$b}{$label}{1}{$alpha});
        }
        foreach my $b (@{$Betavalues{2448}{$Mass{2448}}}){        #loops over the small volume beta values
          push(@y8,$Avg{2448}{$Mass{2448}}{$b}{$label}{2}{$alpha});
          push(@e8,$Err{2448}{$Mass{2448}}{$b}{$label}{2}{$alpha});
        }
        foreach my $b (@{$Betavalues{2448}{$Mass{2448}}}){        #loops over the small volume beta values
          push(@y9,$Avg{2448}{$Mass{2448}}{$b}{$label}{3}{$alpha});
          push(@e9,$Err{2448}{$Mass{2448}}{$b}{$label}{3}{$alpha});
        }
        my $chart = Chart::Gnuplot->new(                  #Create chart object 
          output => "${Alphain}/overlay/alpha_${alpha}_obs_${label}.png",
          xtics => {
            font => "Times-Roman, 30",
          },
          ytics => {
            font => "Times-Roman, 30",
          },
        );
        $chart->command("unset key");
        $chart->command("set size ratio .714");
        my $dataSet1 = Chart::Gnuplot::DataSet->new(      #Create dataset object for small volumes
          xdata => \@x1,
          ydata => [\@y1, \@e1],
          title => "6x12 blocked 0",
          style => "yerrorlines",
        );
        my $dataSet2 = Chart::Gnuplot::DataSet->new(      #Create dataset object for small volumes
          xdata => \@x2,
          ydata => \@y2,
          title => "6x12 blocked 1",
          style => "lines",
          width => "6",
          color => "blue",
          linetype => "solid",
        );
        my $dataSet3 = Chart::Gnuplot::DataSet->new(      #Create dataset object for small volumes
          xdata => \@x3,
          ydata => [\@y3, \@e3],
          title => "12x24 blocked 0",
          style => "yerrorlines",
        );
        my $dataSet4 = Chart::Gnuplot::DataSet->new(      #Create dataset object for small volumes
          xdata => \@x4,
          ydata => \@y4,
          title => "12x24 blocked 1",
          style => "lines",
          linetype => "dash",
          color => "red",
          linesize => "2",
        );
        my $dataSet5 = Chart::Gnuplot::DataSet->new(      #Create dataset object for small volumes
          xdata => \@x5,
          ydata => \@y5,
          title => "12x24 blocked 2",
          style => "lines",
          width => "6",
          color => "dark-green",
          linetype => "solid",
        );
        my $dataSet6 = Chart::Gnuplot::DataSet->new(      #Create dataset object for small volumes
          xdata => \@x6,
          ydata => [\@y6, \@e6],
          title => "24x48 blocked 0",
          style => "yerrorlines",
        );
        my $dataSet7 = Chart::Gnuplot::DataSet->new(      #Create dataset object for large volume
          xdata => \@x7,
          ydata => \@y7,
          title => "24x48 blocked 1",
          style => "lines",
          linetype => "2dash",
          color => "blue",
        );
        my $dataSet8 = Chart::Gnuplot::DataSet->new(      #Create dataset object for large volume
          xdata => \@x8,
          ydata => \@y8,
          title => "24x48 blocked 2",
          style => "lines",
          linetype => "dash",
          color => "blue",
        );
        my $dataSet9 = Chart::Gnuplot::DataSet->new(      #Create dataset object for large volume
          xdata => \@x9,
          ydata => \@y9,
          title => "24x48 blocked 3",
          style => "lines",
          width => "6",
          color => "red",
          linetype => "solid",
        );
        $chart->plot2d($dataSet2, $dataSet5, $dataSet9);  #plots the chart
      }
    }
  }
  print"*\t*\t*\nDone\n";
}
