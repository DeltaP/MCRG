package stat_mod;
use strict;
use warnings;
#needs serious upgrades in the block function, and jack knife sections.  Standard deviation and average work fine.
sub auto_stdev($\@);
sub block(\$ \@);
sub jack_block_error(\$ \@);
sub jack_error2($\@);

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION    =1.00;
@ISA	    =qw(Exporter);
@EXPORT	    =();
@EXPORT_OK  =qw(block);
%EXPORT_TAGS=(DEFAULT=>[qw(&block &avg &stdev &jack_error &jack_block_error)]);

sub block(\$ \@) {#blocks an array into chuncks given by an argument blocksize forming a new array
  if(@_ != 2){print("WARNING!  The block function expects a reference to a scalar and a reference to an array, you have failed to provide one or both of these references\n");}
  my ($data,$blocksize)=@_;
  my $bs=${$blocksize};
  my @blocked=();
  my $i=0;
  my $total=0;
 
  foreach my $value (@{$data}){
    $total+=$value;
    $i++;
    if ($i==$bs){
      push(@blocked, ($total/$bs));
      $i=0;
      $total=0;
    }
  }
  return @blocked;
}

sub avg {#finds the average of an array of data
  if(@_ < 2){print("WARNING!  The average function expects at least two values, it appears that you are trying to average one number\n");}
  
  my @data=@_;
  my $average=0.0;
  my $n=@data;
  foreach my $i (@data) {$average+=$i;}
  $average/=$n;
  return $average;

}

sub stdev{#finds the standard deviation of an array of data
  if(@_ < 2){print("WARNING!  The stdev function expects at least two values, it appears that you are trying to find the standard deviation of one number\n");}
  
  my @data=@_;
  my $average=&avg(@data);
  my $diff=0;
  foreach my $i (@data) {$diff+=($i-$average)**2;}
  my $error=sqrt($diff/(@data-1));
  return $error;
}

sub auto_stdev($\@){
 if(@_ != 2){print("WARNING!  The auto_stdev function expects at least two values, it appears that you are trying to find the standard deviation of one number\n");}
  
  my $t    = $_[0];
  my @data = @{$_[1]};

  my $average=&avg(@data);
  my $diff=0;
  foreach my $i (@data) {$diff+=($i-$average)**2;}
  my $error=sqrt(2*$t*$diff/(@data));
  return $error;
}

sub jack_error1{#finds the jack knife error of a set of blocked data using the average
  if(@_ < 2){print("WARNING!  The error function requires the following argument (\@data)\n");}
  
  my @data=@_;
  my $average=&avg(@data);
  my $diff=0;
  foreach my $i (@data) {$diff+=($i-$average)**2;}
  my $error=sqrt($diff/(@data*(@data-1)));
  return $error;
}

sub jack_error2($\@){#finds the jack knife error of a set of blocked data using the full data set instead of the average
  if(@_ != 2){print("WARNING!  The jack_error2 function expects a reference to a scalar and a reference to an array, you have failed to provide one or both of these references\n");}
  
  my $full=$_[0];
  my @data=@{$_[1]};
  my $diff=0;
  foreach my $i (@data) {$diff+=($i-$full)**2;}
  my $error=sqrt($diff/(@data*(@data-1)));
  return $error;
}
1;
