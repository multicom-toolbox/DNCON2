# Author: Jesse Eickholt

package ProteinUtils;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(residue_prob_from_msa scaled_log_odds_n_inf calc_entropy hydrophobicity residue_propensity residue_residue_contacts 
                braun_contact_potential levitt_contact_potential lu_contact_potential atchley_factor scaled_ordered_mean pssm_counts 
                aa_hot_encoder aa_hot_decoder strip_and_number_vector flesh_out_vector gen_aa_hot_window estimate_joint_entropy
                scld_residue_residue_contacts scld_lu_contact_potential scld_levitt_contact_potential scld_braun_contact_potential);
%EXPORT_TAGS = ( DEFAULT => [qw( )]);


#####################################################################
# Name : aa_hot_encoder
# Takes: AA 
# Returns: A binary vector representing the AA  
####################################################################

sub aa_hot_encoder {

  my $aa_str = "ACDEFGHIKLMNPQRSTVWY"; #20 std aa
  my $encoding ="0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ";
  my $features = uc($_[0]);

  my $idx = index($aa_str, $features);
  if ($idx < 0) { # unknonw        
    $idx = 20;
  }
  
  substr($encoding, $idx * 2, 1, '1');
  return $encoding;
}

#####################################################################
# Name : gen_aa_hot_window
# Takes: A reference to a lookup table for a residue, length, index, 
# and window size
# Returns: A string representing the feature vector  
####################################################################

sub gen_aa_hot_window {
  my ($seq_ref, $length, $index, $window_size) = @_;
  my $vector = '';

  if($window_size % 2 == 0) { $window_size++; }
  
  for( my $offset = -1 * ($window_size - 1)/2; $offset <= ($window_size -1)/2; $offset++) {
    if( $index+$offset >= 1 && $index+$offset <= $length) {
      $vector .= aa_hot_encoder($seq_ref->[$index+$offset]);
    } else {
      $vector .= aa_hot_encoder(' ');
    }
  }

  #$vector =~ s/\s+$//;
  return $vector;
}

#####################################################################
# Name : aa_hot_decoder
# Takes: An aa_hot_encoded vector
# Returns: The AA or blank  
####################################################################

sub aa_hot_decoder {
  
  #20 std aa
  my @residues = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y', '_');  
  my $encoding = $_[0]; 
  my @values = split(/\s+/, $encoding);

  for(my $i=0; $i<=20; $i++) {
    if($values[$i] eq 1) {
#print "|| " . $values[$i] . ", $i ->  ";
      return $residues[$i];
    }
  }
  
  return "_";

}

######################################################################
# Name : strip_and_number_vector
# Takes: A string of features seperated by spaces
# Returns: Each no zero vector will be proceded by a count starting at
#  1.
######################################################################

sub strip_and_number_vector {
  my @features = split(/\s+/, $_[0]);
  my $count = 1;
  my $vector = '';

  foreach my $features_t (@features) {
    if($features_t != 0) {
      $vector .= "$count:" . $features_t . " ";
    } 
    $count++;
  }
  
  $vector =~ s/\s+$/ /;
  return $vector;
}

######################################################################
# Name : flesh_out_vector
# Takes : A string with feature pairs containing feature number and 
#  value and a max feature number
# Returns: A string with feature values filling in any missing features
#  with zeros
######################################################################

sub flesh_out_vector {
  my @pairs = split(/\s+/, $_[0]);
  my $max_feature_number = $_[1];
  my %features = ();
  my $vector = ''; 

  foreach my $pair_t (@pairs) {
    my ($num, $value) = split(/\:/, $pair_t);
    $features{$num} = $value;
  }  


  for(my $i = 1; $i <= $max_feature_number; $i++) {
    if(exists $features{$i}) {
      $vector .= $features{$i} . " ";
    } else {
      $vector .= "0 ";
    }
  }
  
  $vector =~ s/\s+$/ /;
  return $vector;

}


######################################################################
# Name : residue_prob_from_msa
# Takes: number of sequences in msa and array of sequences
# Returns: a 2d array indexed first by residue and second by probability of 
#  aa at that position.  There are 20 std aa and index 20 is used for gap.
#  The probabilities are initialized to 0.
######################################################################

sub residue_prob_from_msa { 
  my($num_of_sequences, @sequences);
  my ($aa_str, @residue_prob, $i, $j, $aa, $idx, $length);

  ($num_of_sequences, @sequences) = @_;
  chomp($num_of_sequences);
  chomp(@sequences);

  $aa_str = "ACDEFGHIKLMNPQRSTVWY"; #20 std aa
  $length = length($sequences[0]);

  #profile size: 20 + one gap, if all zero, means not exist(for pos outside of window)
  @residue_prob = ();
  for ($i = 0; $i < $length; $i++) {
    for ($j = 0; $j < 21; $j++) {
     $residue_prob[$i][$j] = 0;
    }
  }


  for ($i = 0; $i < $length; $i++) {
    for($j = 0; $j < $num_of_sequences; $j++) {
      $aa = substr($sequences[$j], $i, 1);
      $aa = uc($aa);
      $idx = index($aa_str, $aa);
      if ($idx < 0) { #gap case or unknonw        
        #treated as a gap
        $idx = 20;
      }
      $residue_prob[$i][$idx] +=  (1 / $num_of_sequences);
    }
  }

  return @residue_prob;
}

#######################################################################
# Name : scaled_log_odds_n_inf 
# Takes: An array of lines coming from a pssm
# Returns: A 2d array for scaled versions of the log odds of each 
#  residue at each positions and an array of scaled information  
#  values at each position
######################################################################
sub scaled_log_odds_n_inf {
  my @lines = @_;
  my @pssm =();
  my @inf = ();

  # Get rid of header
  shift @lines; shift @lines; shift @lines;

  while(my $line_t = shift @lines) {
    if($line_t =~ m/^\s*$/) {
      last; 
    }
    my @fields = ();
    for(my $i = 0; $i < 20; $i++) {
       my $val = substr($line_t, 9+$i*3, 3);
       push(@fields, sprintf("%.4f", 1/(1+exp(-1*$val))));
    }
    push(@pssm, [@fields]);

    my $inf_t = substr($line_t, 151, 5);
    if($inf_t > 6) {
      push(@inf, 1.0);
    } elsif( $inf_t < 0) {
      push(@inf, 0);
    } else {
      push(@inf, sprintf("%.4f", ($inf_t/6.0)));
    }
  }

  return(\@pssm, \@inf);

}

#######################################################################
# Name : estimate_joint_entropy 
# Takes: An array of lines from a pssm
# Returns: A 2d array for scaled versions of the joint entropy  
#  here we assume p(x,y) = p(x) * p(y)
######################################################################
sub estimate_joint_entropy {
  my @lines = @_;
  my @counts =();

  # Get rid of header
  shift @lines; shift @lines; shift @lines;

  while(my $line_t = shift @lines) {
    if($line_t =~ m/^\s*$/) {
      last; 
    }
    my @fields = ();
    for(my $i = 0; $i < 20; $i++) {
       # offset is to jump over to the counts section of a pssm
       my $val = substr($line_t, 71+$i*4, 3);
       push(@fields, sprintf("%d", $val));
    }
    push(@counts, [@fields]);

  }

  my @joint_entropy = ();
  for(my $i = 0; $i < @counts; $i++) {
    $joint_entropy[$i] = ();
    for(my $j = 0; $j < @counts; $j++) {
      $joint_entropy[$i][$j] = 0;
    }
  }  

  for(my $i = 0; $i < @counts; $i++) {
    for(my $j = $i+1; $j < @counts; $j++) {
 
      my $sum = 0;     
      for(my $k = 0; $k < 20; $k++) {
        for(my $l = 0; $l < 20; $l++) {
          my $joint = ($counts[$i][$k]*$counts[$j][$l])/10000;
          if($joint > 1) {
            $joint = 1;
          }
          if($joint > 0) {
            $sum -= ($joint * log($joint));
          }
        }
      }
      # 6 is what you get if you have a uniform dist over the joint, add some extract just in case 
      $joint_entropy[$i][$j] = $sum/6;
      $joint_entropy[$j][$i] = $joint_entropy[$i][$j];
    }
  }

  return @joint_entropy;
  
}


#######################################################################
# Name : pssm_counts 
# Takes: Lines read in from a pssm file
# Returns: A 2d matrix containing the residue type counts in the pssm 
######################################################################
sub pssm_counts {
  my @lines = @_;
  my @counts =();

  # Get rid of header
  shift @lines; shift @lines; shift @lines;

  while(my $line_t = shift @lines) {
    if($line_t =~ m/^\s*$/) {
      last; 
    }
    my @fields = ();
    for(my $i = 0; $i < 20; $i++) {
       # offset is to jump over to the counts section of a pssm
       my $val = substr($line_t, 71+$i*4, 3);
       push(@fields, sprintf("%d", $val));
    }
    push(@counts, [@fields]);

  }

  return @counts;

}

#######################################################################
# Name : calc_entropy 
# Takes: An array of probablities
# Returns: The negated sum of the probablilites of each aa times the 
#  log of the probablity.  
######################################################################
sub calc_entropy {

  my @vector = @_;
  my $entropy = 0;
  my $i;

  for ($i = 0; $i < @vector; $i++) {
    if ($vector[$i] > 0) {

      $entropy -= ($vector[$i] * log($vector[$i]));
    }
  }
  return $entropy;
}


########################################################################
# Name : scaled_ordered_mean 
# Takes: An array reference, start, stop, order and character to match
#  for d(k) to be one
# Returns: A scaled value for the nth order weighted mean
#  The value came from Li et al., doi:10.1093/bioinformatics/btr579
#######################################################################

sub scaled_ordered_mean {
  my ($values_ref, $i, $j, $order, $type) = @_;
  my @values = @{$values_ref};

  my $total = 0.0; my $max_val = 0.0;
  for(my $m = $i; $m <= $j; $m++) {
    my $weight = (1/($j-$i))*(($m-$i)/($j-$i))**($order);
    $max_val += $weight;
    if($values[$m] eq $type) {
      $total += $weight;
    } 

  }  

  return (1.0*$total/$max_val);

}

########################################################################
# Name : atchley_factor 
# Takes: one letter identifier for aa and factor number
# Returns: value for factor and aa
#######################################################################

sub atchley_factor {
  my ($aa, $index) = @_;

  my @factor = ({}, 
                {"A" => -0.591, "C" => -1.343, "D" =>1.050, "E" => 1.357, "F" => -1.006, "G" => -0.384, 
                 "H" => 0.336, "I" => -1.239, "K" => 1.831, "L" => -1.019, "M" => -0.663, "N" => 0.945, 
                 "P" => 0.189, "Q" => 0.931, "R" => 1.538, "S" => -0.228, "T" => -0.032, "V" => -1.337, 
                 "W" => -0.595, "Y" => 0.260},
                 {"A" => -1.302, "C" => 0.465, "D" => 0.302, "E" => -1.453, "F" => -0.590, "G" => 1.652, 
                  "H" => -0.417, "I" => -0.547, "K" => -0.561, "L" => -0.987, "M" => -1.524, "N" => 0.828, 
                  "P" => 2.081, "Q" => -0.179, "R" => -0.055, "S" => 1.399, "T" => 0.326, "V" => -0.279, 
                  "W" => 0.009, "Y" => 0.830},
                 {"A" => -0.733, "C" => -0.862, "D" => -3.656, "E" => 1.477, "F" => 1.891, "G" => 1.330, 
                  "H" => -1.673, "I" => 2.131, "K" => 0.533, "L" => -1.505, "M" => 2.219, "N" => 1.299, 
                  "P" => -1.628, "Q" => -3.005, "R" => 1.502, "S" => -4.760, "T" => 2.213, "V" => -0.544, 
                  "W" => 0.672, "Y" => 3.097},
                 {"A" => 1.570, "C" => -1.020, "D" => -0.259, "E" => 0.113, "F" => -0.397, "G" => 1.045, 
                  "H" => -1.474, "I" => 0.393, "K" => -0.277, "L" => 1.266, "M" => -1.005, "N" => -0.169, 
                  "P" => 0.421, "Q" => -0.503, "R" => 0.440, "S" => 0.670, "T" => 0.908, "V" => 1.242, 
                  "W" => -2.128, "Y" => -0.838},
                 {"A" => -0.146, "C" => -0.255, "D" => -3.242, "E" => -0.837, "F" => 0.412, "G" => 2.064, 
                  "H" => -0.078, "I" => 0.816, "K" => 1.648, "L" => -0.912, "M" => 1.212, "N" => 0.933, 
                  "P" => -1.392, "Q" => -1.853, "R" => 2.897, "S" => -2.647, "T" => 1.313, "V" => -1.262, 
                  "W" => -0.184, "Y" => 1.512});

  return (defined($factor[$index]) && defined($factor[$index]->{$aa})) ? $factor[$index]->{$aa} : 0;

}



########################################################################
# Name : hydrophobicity
# Takes: one letter identifier for aa
# Returns: hydrophobiicy of aa
#######################################################################

sub hydrophobicity {

  my %hydrophobicity = ("I", "0.73", "F", "0.61", "V", "0.54", "L", "0.53", "W", "0.37", "M", "0.26",
                      "A", "0.25", "G", "0.16", "C", "0.04", "Y", "0.02", "P", "-0.07", "T", "-0.18",
                      "S", "-0.26", "H", "-0.40", "E", "-0.62", "N", "-0.64", "Q", "-0.69", "D", "-0.72",
                      "K", "-1.10", "R", "-1.76");
  ($_) = @_;
  return (defined($hydrophobicity{$_})) ? $hydrophobicity{$_} : 0;

}

#######################################################################
# Name : residue_propensity
# Takes: one letter identifier for aa
# Returns: residue_propensity
#######################################################################

######################################################################
# NOTA GRANDE ! These values are were calculated using Heterodimers
#  Should I be using something else?  Also, they have been scaled 
#  down (divided by 100)
#####################################################################

sub residue_propensity {

  my %residue_propensity = ("I", "0.28", "F", "0.68", "V", "0.24", "L", "0.22", "W", "0.35", "M", "0.10",
                      "A", "0.08", "G", "0.15", "C", "0.16", "Y", "0.58", "P", "0.07", "T", "0.11",
                      "S", "0.10", "H", "0.08", "E", "0.20", "N", "0.10", "Q", "0.09", "D", "0.24",
                      "K", "0.24", "R", "0.63");

  ($_) = @_;
  return (defined($residue_propensity{$_})) ? $residue_propensity{$_} : 0;  


}

#####################################################################
# Name : residue_residue_contacts
# Takes: two residues
# Returns: a preference score scaled down (divided by 100).  Values 
#  come from Residue Frequencies and Pairing Preferences at Protein Protein Interfaces
#  Fabian Glaser, David M. Steinberg, Ilya A. Vakser, and Nir Ben-Tal
####################################################################

sub residue_residue_contacts {
  my($i, $j, @preference_score, $alphabet, $i_index, $j_index);
  ($i, $j) = @_;

  $alphabet = "IVLFCMAGTSWYPHEQDNKR";
  $i_index = index($alphabet, $i);
  $j_index = index($alphabet, $j);

  # if an unknown letter for an aa is recieved, set indicies to return something close to 0.
  if($i_index < 0 || $j_index < 0) {
    $i_index = 10; $j_index = 10;  

  }

# Create a array of preference scores.   
#                I     V     L    F    C    M     A     G     T    S    W    Y     P    H    E    Q    D    N    K    R
#
  @preference_score
      = ( [ qw( .78  1.52  1.37  .82  .23  .50  1.74  1.55  1.09  .87  .32  .77  1.10  .35  .73  .56  .70  .56  .49  .60 ) ], # I
          [ qw(      1.79  1.93 1.09  .46  .63  2.52  1.82  1.60 1.84  .23  .81  1.56  .52 1.13  .79  .99  .82 1.00 1.01 ) ], # V
          [ qw(            1.80 1.10  .45  .76  2.56  1.78  1.30 1.43  .43  .83  1.38  .74 1.07  .81  .85  .99  .72 1.18 ) ], # L
          [ qw(                  .62  .27  .38  1.36  1.01   .88  .78  .22  .61  1.04  .27  .51  .49  .39  .60  .40  .53 ) ], # F
          [ qw(                       .43  .11   .61   .59   .33  .59  .06  .18   .47  .20  .30  .16  .21  .17  .18  .23 ) ], # C
          [ qw(                            .28   .72   .75   .41  .74  .11  .30   .53  .22  .40  .30  .21  .31  .27  .27 ) ], # M
          [ qw(                                 2.28  2.45  2.03 2.15  .47 1.06  1.95  .83 1.47 1.03 1.52 1.63 1.08 1.10 ) ], # A
          [ qw(                                       1.92  2.31 1.98  .43 1.15  1.88  .84 1.16 1.17 1.65 1.40 1.29 1.47 ) ], # G
          [ qw(                                             1.23 1.82  .42  .42   .74 1.62  .51 1.15  .63 1.71  .92 1.01 ) ], # T
          [ qw(                                                  1.47  .32  .78  1.53  .42 1.38  .84 1.76 1.27  .95 1.04 ) ], # S
          [ qw(                                                        .07  .21   .76  .17  .11  .08  .18  .21  .21  .43 ) ], # W
          [ qw(                                                             .55   .91  .43  .66  .26  .41  .60  .52  .56 ) ], # Y
          [ qw(                                                                   .97  .51 1.18  .89  .94 1.29  .90 1.02 ) ], # P
          [ qw(                                                                        .28  .30  .31  .69  .34  .22  .39 ) ], # H
          [ qw(                                                                             .56  .42  .46  .79  .87 1.03 ) ], # E
          [ qw(                                                                                  .36  .67  .66  .40  .54 ) ], # Q
          [ qw(                                                                                       .55 1.22  .74 1.01 ) ], # D
          [ qw(                                                                                            .93  .59  .74 ) ], # N
          [ qw(                                                                                                 .36  .31 ) ], # K
          [ qw(                                                                                                      .83 ) ] ); # R

  return ($j_index >= $i_index) ? $preference_score[$i_index][$j_index - $i_index] : $preference_score[$j_index][$i_index - $j_index]; 
}

sub scld_residue_residue_contacts {
  my ($i, $j) = @_;

  my $r_val = (residue_residue_contacts($i,$j) - 0.06 ) / 2.50;
  return $r_val;

}

#####################################################################
# Name : lu_contact_potential
# Takes: two residues
# Returns: a contact potential score.  
####################################################################

sub lu_contact_potential {
  my($i, $j, @contact_potential, $alphabet, $i_index, $j_index);
  ($i, $j) = @_;

  $alphabet = "GASCVTIPMDNLKEQRHFYW";
  $i_index = index($alphabet, $i);
  $j_index = index($alphabet, $j);


  # if an unknown letter for an aa is recieved, set indicies to return something close to 0.
  if($i_index < 0 || $j_index < 0) {
    $i_index = 0; $j_index = 16;  
  }


       # GLY   ALA   SER   CYS   VAL   THR   ILE   PRO   MET   ASP   ASN   LEU   LYS   GLU   GLN   ARG   HIS   PHE   TYR   TRP
  @contact_potential = ( 
  [ qw(  0.4   0.4   0.7  -1.0  -0.1   0.6  -0.4   0.6  -0.4   1.3   0.7  -0.4   1.1   1.3   0.6   0.1  -0.0  -0.8  -0.8  -1.0 ) ], # G  
  [ qw(       -0.4   0.4  -1.2  -0.7   0.4  -0.9   0.5  -0.9   1.2   0.5  -1.0   1.3   1.1   0.5  -0.0  -0.1  -1.2  -1.0  -0.9 ) ], # A
  [ qw(              0.1  -1.2  -0.2   0.5  -0.6   0.6  -0.2   0.9   0.6  -0.5   1.1   0.8   0.6   0.1  -0.4  -0.9  -0.7  -0.8 ) ], # S
  [ qw(                   -4.4  -1.9  -0.8  -2.2  -1.0  -2.2   0.0  -0.3  -2.0  -0.2   0.2  -0.8  -1.6  -1.9  -2.4  -2.3  -2.6 ) ], # C
  [ qw(                         -2.0  -0.5  -1.9  -0.3  -1.7   0.4  -0.2  -1.9   0.3   0.4  -0.3  -0.8  -0.9  -2.1  -1.8  -2.0 ) ], # V
  [ qw(                               -0.1  -0.8   0.6  -0.5   1.0   0.5  -0.7   0.8   1.0   0.5   0.1  -0.3  -0.9  -0.7  -1.1 ) ], # T
  [ qw(                                     -2.7  -0.6  -2.0   0.3  -0.3  -2.3  -0.0  -0.0  -0.6  -1.1  -1.3  -2.4  -2.0  -2.2 ) ], # I
  [ qw(                                            0.1  -0.8   1.3   0.8  -0.8   1.1   1.3   0.4   0.1   0.1  -0.8  -1.0  -1.0 ) ], # P
  [ qw(                                                 -2.9   0.3  -0.4  -2.2  -0.0  -0.0  -0.7  -0.7  -1.3  -2.4  -2.1  -2.6 ) ], # M
  [ qw(                                                        1.1   1.0   0.2   0.7   1.9   0.9  -0.2   0.3   0.1  -0.4  -0.2 ) ], # D
  [ qw(                                                              0.3  -0.3   1.0   1.1   0.6   0.1   0.3  -0.7  -0.5  -0.6 ) ], # N
  [ qw(                                                                   -2.7   0.0  -0.0  -0.6  -1.0  -1.1  -2.3  -2.1  -2.4 ) ], # L
  [ qw(                                                                          1.6   0.8   1.0   0.8   0.7  -0.2  -0.2  -0.2 ) ], # K
  [ qw(                                                                                1.2   1.1   0.0   0.2  -0.1  -0.4   0.0 ) ], # E
  [ qw(                                                                                     -0.2   0.1  -0.0  -0.9  -0.6  -1.2 ) ], # Q
  [ qw(                                                                                           -0.6  -0.4  -1.2  -1.3  -1.4 ) ], # R
  [ qw(                                                                                                 -1.3  -1.4  -1.6  -1.9 ) ], # H
  [ qw(                                                                                                       -3.0  -2.3  -2.5 ) ], # F
  [ qw(                                                                                                             -2.3  -2.4 ) ], # Y
  [ qw(                                                                                                                   -3.3 ) ]); # W

  # Matrix is upper triangluar so we need to use smaller index first and then offset the second index
  return ($j_index >= $i_index) ? $contact_potential[$i_index][$j_index - $i_index] : $contact_potential[$j_index][$i_index - $j_index];
}

sub scld_lu_contact_potential {
  my ($i, $j) = @_;

  my $r_val = (lu_contact_potential($i,$j) + 4.4  ) / 6.3;
  return $r_val;

}


#####################################################################
# Name : levitt_contact_potential
# Takes: two residues
# Returns: a contact potential score.  
####################################################################

sub levitt_contact_potential {

  my($i, $j, $alphabet, @contact_potential, $i_index, $j_index);
  ($i, $j) = @_;

  $alphabet = "GAVLIPDENQKRSTMCYWHF";
  $i_index = index($alphabet, $i);
  $j_index = index($alphabet, $j);

  # if an unknown letter for an aa is recieved, set indicies to return something close to 0.
  if($i_index < 0 || $j_index < 0) {
    $i_index = 0; $j_index = 4;  
  }


  # Table from SVNcon
  #     Gly   Ala   Val   Leu   Ile   Pro   Asp   Glu   Asn   Gln   Lys   Arg   Ser   Thr   Met   Cys   Tyr   Trp   His   Phe
  @contact_potential = (
  [ qw(  .1    .7    .1    .1   0      .5    .4    .6    .1   0      .4  -0.1    .4    .2  -0.1  -0.1  -0.4  -0.7   0    -0.3 ) ], 
  [ qw(        .5  -0.3  -0.4  -0.4    .6    .3    .6    .3   0     1.0    .2    .5   0    -0.5    .3  -0.7  -0.8   0    -0.8 ) ],
  [ qw(             1.1  -1.2  -1.2   0      .4   0     0    -0.4   0.1  -0.5   0    -0.3  -1.0  -0.5  -1.2  -1.6  -0.5  -1.5 ) ],
  [ qw(                  -1.4  -1.4 -0.1    0    -0.1  -0.1  -0.6   0.1  -0.6   0    -0.3  -1.3  -0.8  -1.4  -1.7  -0.7  -1.6 ) ],
  [ qw(                        -1.5 -0.1    0    -0.2  -0.1  -0.4   0    -0.7  -0.1  -0.6  -1.4  -0.8  -1.4  -1.8  -0.8  -1.7 ) ],
  [ qw(                               .1     .1    .1  -0.1  -0.3    .6  -0.2    .2   0    -0.5   0    -1.0  -1.3  -0.4  -0.7 ) ],
  [ qw(                                     0     0    -0.6  -0.3  -1.0  -1.4  -0.3  -0.3   0.1   0    -1.0  -0.6  -1.1  -0.3 ) ],
  [ qw(                                           0.1  -0.6  -0.4  -1.1  -1.5  -0.2  -0.3  -0.3   0.1  -1.0  -0.8  -1.0  -0.5 ) ],
  [ qw(                                                -0.7  -0.7  -0.3  -0.8  -0.1  -0.4  -0.3   0    -0.8  -0.8  -0.8  -0.6 ) ],
  [ qw(                                                      -0.5  -0.4  -0.9   0    -0.5  -0.6  -0.2  -1.1  -1.0  -0.5  -0.8 ) ],
  [ qw(                                                              .7    .1    .1   0    -0.1    .5  -1.0  -0.8   0    -0.4 ) ],
  [ qw(                                                                  -0.9  -0.4  -0.6  -0.5   0    -1.4  -1.3  -1.0  -0.9 ) ],
  [ qw(                                                                         0    -0.2  -0.1  -0.1  -0.6  -0.6  -0.6  -0.4 ) ],
  [ qw(                                                                              -0.5  -0.6  -0.3  -0.8  -0.9  -0.7  -0.7 ) ],
  [ qw(                                                                                    -1.5  -0.8  -1.5  -2.0  -0.9  -1.9 ) ],
  [ qw(                                                                                          -2.7  -0.8  -1.3  -0.6  -1.2 ) ],
  [ qw(                                                                                                -1.6  -1.8  -1.5  -1.7 ) ],
  [ qw(                                                                                                      -2.2  -1.5  -2.0 ) ],
  [ qw(                                                                                                            -1.6  -1.2 ) ],
  [ qw(                                                                                                                  -2.0 ) ] );

  # Matrix is upper triangluar so we need to use smaller index first and then offset the second index
  return ($j_index >= $i_index) ? $contact_potential[$i_index][$j_index - $i_index] : $contact_potential[$j_index][$i_index - $j_index];
}


sub scld_levitt_contact_potential {
  my ($i, $j) = @_;

  my $r_val = (levitt_contact_potential($i,$j) + 2.7) / 3.8;
  return $r_val;

}



#####################################################################
# Name : braun_contact_potential
# Takes: two residues
# Returns: a contact potential score.  
####################################################################

sub braun_contact_potential {

  my($i, $j, @contact_potential, $alphabet, $i_index, $j_index);
  ($i, $j) = @_;


  $alphabet = "GAVLIFYWMCPSTNQHKRDE";
  $i_index = index($alphabet, $i);
  $j_index = index($alphabet, $j);

  # if an unknown letter for an aa is recieved, set indicies to return something close to 0.
  if($i_index < 0 || $j_index < 0) {
    $i_index = 11; $j_index = 0;  
  }

  # Table from SVMcon
  @contact_potential = (
  [ qw( -0.29                                                                                                                 ) ],
  [ qw( -0.14 -0.18                                                                                                           ) ],
  [ qw( -0.10 -0.15 -0.48                                                                                                     ) ],
  [ qw( -0.04 -0.24 -0.29 -0.43                                                                                               ) ],
  [ qw(  0.27 -0.25 -0.31 -0.45 -0.48                                                                                         ) ],
  [ qw( -0.09 -0.16 -0.31 -0.28 -0.05 -0.50                                                                                   ) ],
  [ qw( -0.21 -0.18  0.00 -0.10 -0.34 -0.27 -0.11                                                                             ) ],
  [ qw( -0.34 -0.01  0.18 -0.18 -0.28  0.16 -0.30 -0.53                                                                       ) ],
  [ qw(  0.25 -0.02 -0.02 -0.32  0.21 -0.36  0.01 -0.73 -0.75                                                                 ) ],
  [ qw( -0.42  0.08  0.08  0.36 -0.16 -0.28  0.69 -0.74  0.27 -1.77                                                           ) ],
  [ qw(  0.06  0.28  0.76  0.30  0.99  0.65 -0.02  0.70 -0.78  0.31 -0.78                                                     ) ],
  [ qw(  0.04  0.38  0.18  0.30  0.57  0.15 -0.03  0.44  0.00  0.12  0.21 -0.68                                               ) ],
  [ qw(  0.28  0.06  0.19  0.57  0.34  0.25  0.23  0.74  0.43  0.28  0.04 -0.23 -0.58                                         ) ],
  [ qw(  0.49 -0.04  0.48  0.25  1.45  0.12 -0.14  0.46 -0.52  0.07  0.59 -0.21 -0.06 -0.45                                   ) ],
  [ qw(  0.54  0.35  0.41  0.35  0.44 -0.04 -0.06 -0.09  0.07  0.39  0.73  0.19 -0.31  0.20 -0.17                             ) ],
  [ qw( -0.09  0.44  0.37  0.10  0.24  0.25  0.33 -0.34  1.07 -0.45 -0.21 -0.13 -0.22 -0.56  0.28 -0.15                       ) ],
  [ qw(  0.56  0.28  0.53  0.37 -0.00  0.75 -0.00  0.02  0.44  0.68  0.26 -0.05 -0.26 -0.27  0.05  0.57  0.21                 ) ],
  [ qw(  0.40  0.59  0.43  0.37  0.05  0.31  0.03 -0.20  0.53  0.92  0.34  0.24 -0.31 -0.00  0.56 -0.11  0.58 -0.03           ) ],
  [ qw( -0.26  0.24  0.51  0.80  0.26  0.33  0.61  0.74  0.21  0.53  0.87 -0.03  0.32 -0.43 -0.03 -0.61 -0.43 -0.79 0.11      ) ],
  [ qw(  0.21  0.53  0.37  0.51  0.53  0.38  0.25  1.37  0.44  0.17  0.41  0.10  0.27  0.76 -0.20 -0.14 -1.12 -0.85 0.86  0.58) ] );
  #       G     A     V     L     I     F     Y     W     M     C     P     S     T     N     Q     H     K     R    D     E

 # Here the matrix is lower triangluar.  No need for funny offset business/magic
  return ($i_index >= $j_index) ? $contact_potential[$i_index][$j_index] : $contact_potential[$j_index][$i_index];

}

sub scld_braun_contact_potential {
  my ($i, $j) = @_;

  my $r_val = (braun_contact_potential($i,$j) + 1.77  ) / 3.22;
  return $r_val;

}




1;


#### Test Scripts #####
##!/bin/bash

#len=`wc -l ../datasets/skolnick/ss_sa/5HPGB.align | cut -d ' ' -f 1`;
#echo $len;

#for j in `seq 1 84`; do
#echo "$j: -----";
#for i in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z "\."; do
#  count=`cut -b $j  ../datasets/skolnick/ss_sa/5HPGB.align | grep $i | wc -l`
#  if [[ $count -gt 0 ]]; then
#    echo -n "   $i: "
#    echo "scale=4; $count / $len" | bc;
#  fi
#done

#done

#########################

##!/usr/bin/perl -w

#use strict;

#use ProteinUtils qw( residue_prob_from_msa );
#use MyMathUtils qw( round );
#open FILE, "<" . $ARGV[0];
#my @lines = <FILE>;
#close FILE;

#my $aa_str = "ACDEFGHIKLMNPQRSTVWY-";
#my @matrix = residue_prob_from_msa(@lines, \@lines);

#for(my $i=0; $i < @matrix; $i++) {

#  print $i+1 .  ": -------\n";
#  for(my $j=0; $j < @{$matrix[$i]}; $j++) {
#    if(@{$matrix[$i]}[$j] > 0) {
#      print "   " . substr($aa_str, $j, 1) . ": " . round(@{$matrix[$i]}[$j], 4) . " \n";
#    }
#  }
#  print "\n";
#}



