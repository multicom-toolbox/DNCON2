# Author: Jesse Eickholt

package MyMathUtils;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(calc_cosine calc_pearson_r round mutual_information calc_mu calc_sigma calc_zscores);
%EXPORT_TAGS = ( DEFAULT => [qw( )]);

my $DEBUG = 1;

#######################################################################
# Name : calc_pearson_r
# Takes: references to two arrays
# Returns: the pearson correlation coefficent for the two arrays 
######################################################################

sub calc_pearson_r { 
  my ($r, $x, $y, $x_avg, $y_avg, $xy_avg, $i, $length);
  my ($x_sigma, $y_sigma);

  ($x, $y) = @_;
  $x_avg = 0; $y_avg = 0; $x_sigma = 0; $y_sigma = 0; $xy_avg = 0; $r = 0;

  if($DEBUG == 1) {
    if(@{$x} != @{$y}) { die "Arrays are of differing length in calc_pearson_r.\n"; }
  }

  $length = @{$x};
  for($i = 0; $i < $length ; $i++) {
    $x_avg += $x->[$i];
    $y_avg += $y->[$i];
  }
  $x_avg = $x_avg/$length;
  $y_avg = $y_avg/$length;
  
  for($i = 0; $i < $length; $i++) {
    $x_sigma += ($x->[$i] - $x_avg) * ($x->[$i] - $x_avg);
    $y_sigma += ($y->[$i] - $y_avg) * ($y->[$i] - $y_avg);
    $xy_avg += ($x->[$i] - $x_avg) * ($y->[$i] - $y_avg);
  }

  $x_sigma = sqrt($x_sigma);
  $y_sigma = sqrt($y_sigma);

  if($x_sigma > 0 && $y_sigma > 0) {
    $r = ($xy_avg)/ ($x_sigma * $y_sigma);
  }

  # Make sure r is between [-1,1]
  $r = (abs($r) < 1) ? $r : abs($r)/$r;
  return $r;

}

#######################################################################
# Name : calc_cosine 
# Takes: Two references to arrays representing two vectors
# Returns: The angle between the two vectors
######################################################################
sub calc_cosine {

  my ($dot_prod, $x, $y, $x_norm, $y_norm, $cos, $i, $length);
  ($x, $y) = @_; 
 
  if($DEBUG == 1) {
    if(@{$x} != @{$y}) { die "Arrays are of differing length in calc_pearson_r.\n"; }
  }

  $length = @{$x}; 

  for($i = 0; $i < $length; $i++) {
    $dot_prod += $x->[$i] * $y->[$i];
    $x_norm += $x->[$i] * $x->[$i];
    $y_norm += $y->[$i] * $y->[$i];
  }
  $x_norm = sqrt($x_norm);
  $y_norm = sqrt($y_norm);

  $cos = $dot_prod / ($x_norm * $y_norm);

  return $cos;
}

#######################################################################
# Name : round
# Takes : a number and the number of decimals to keep
# Returns : the rounded value
######################################################################

sub round {
  my $round_up = 0;
  my($places, $value);
  ($value, $places) = @_;

  if (! defined($places)) { $places = 3; }
  if(($value * (10 ** ($places+1)) % 10) >= 5) {
    $round_up = 1;
  }

  $value = int($value * (10 ** $places));
  if($round_up) {
    $value++;
  }
  $value = $value / (10 ** $places);

  return $value;
}

#########################################################################
# Name : mutual_information
# Takes : two references to probability vectors
# Returns: the calculated mutual information between the vectors
########################################################################

sub mutual_information {
  my($x, $y, @xy, $i, $j, $mutual_information);
  ($x, $y) = @_;
  @xy = ();
  $mutual_information = 0;

  # Calculate joint dist and mutual info in one pass.
  for($i = 0; $i < @{$x}; $i++) {
    for($j = 0; $j < @{$y}; $j++) {
      $xy[$i][$j] = $x->[$i] * $y->[$j];
      if($xy[$i][$j] > 0) {
        $mutual_information += ($xy[$i][$j] * log($xy[$i][$j] /$x->[$i] / $y->[$j]));
      }
    }
  }

  return ($mutual_information * (10**16));
  
}

########################################################################
# Name : calc_mu 
# Takes : an array of values 
# Returns: the calculated mean of the data points 
########################################################################

sub calc_mu {

  my @data = @_;
  my $sum = 0;

  foreach my $value (@data) {
    $sum += $value;
  }
  
  if(@data > 0) {
    return $sum / @data;
  } else {
    return 0;
  }

}

#########################################################################
# Name : calc_sigma
# Takes : an array of values
# Returns : the calculated standard deviation of the data points
########################################################################

sub calc_sigma {
  my @data = @_;
  my $deviations = 0;
  my $mu = calc_mu(@data);

  foreach my $value (@data) {
    $deviations += ($value - $mu) ** 2;
  }
 
  if(@data > 1) {
    sqrt($deviations / (@data -1));
  } else {
    return 0;
  }

}

#########################################################################
# calc_zscores
# Takes : an array of values
# Returns : an array where each value has been replaced by its zscore
#########################################################################

sub calc_zscores {
  my @data = @_;
  my @zscores = ();
  my $mu = calc_mu(@data);
  my $sigma = calc_sigma(@data);

  foreach my $value (@data) {
    if($sigma > 0) {
      push (@zscores, ($value - $mu)/$sigma);
    } else { 
      push (@zscores, 0);
    }
  }

  return @zscores;

}

1;


