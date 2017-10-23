#!/usr/bin/perl -w
# Jesse Eickholt (original script)
# Badri Adhikari (added coevolution-based features), 5-16-2017

use Cwd 'abs_path';
use File::Basename;
use strict;
use warnings;
use Carp;
use lib ''.abs_path(dirname($0));
use ProteinUtils qw(estimate_joint_entropy calc_entropy atchley_factor scaled_log_odds_n_inf scld_residue_residue_contacts scld_lu_contact_potential scld_levitt_contact_potential scld_braun_contact_potential scaled_ordered_mean pssm_counts);
use MyMathUtils qw(calc_cosine calc_pearson_r);

if (@ARGV ne 10){
  print STDERR "Usage: $0 <fasta> <pssm> <ss_sa> <colstat> <pairstat> <freecontact> <ccmpred> <psicov> <psipred> <psisolv>\n";
  exit;
}

####################################################################################################
my $fasta_fname       = $ARGV[0];
my $pssm_fname        = $ARGV[1];
my $ss_sa_fname       = $ARGV[2];
my $colstat_fname     = $ARGV[3];
my $pairstat_fname    = $ARGV[4];
my $freecontact_fname = $ARGV[5];
my $ccmpred_fname     = $ARGV[6];
my $psicov_fname      = $ARGV[7];
my $psipred_fname     = $ARGV[8];
my $psisolv_fname     = $ARGV[9];

my $id = $fasta_fname;
$id =~ s/\..*//;
$id =~ s/^.*\///;

####################################################################################################
open FASTA, "<" . $fasta_fname or die "Couldn't open fasta file ".$fasta_fname."\n";
my @lines = <FASTA>;
chomp(@lines);
close FASTA;

shift @lines;
my $seq = join('', @lines);
$seq =~ s/ //g;
my @seq = split(//, $seq);

# Calculate scalars and acthley factors
my $aas = "ACDEFGHIKLMNPQRSTVWY";
my @scalar = ();
my @offset = ();
my @res = split(//, $aas);

for(my $i = 1; $i <= 5; $i++) {
  my $min = 15;
  my $max = - 15;
  for(my $j = 0; $j < @res; $j++) {

    if(atchley_factor($res[$j],$i) < $min) {
      $min = atchley_factor($res[$j], $i);
    }
    if(atchley_factor($res[$j],$i) > $max) {
      $max = atchley_factor($res[$j], $i);
    }

  }
  $scalar[$i] = (1/($max-$min));
  $offset[$i] = ($min/($min-$max));
}

my @factors = ();
for (my $i = 0; $i < @seq; $i++) {
  $factors[$i] = [];
  unshift(@{$factors[$i]}, 0);
}

# Get all factors, once
for (my $i = 0; $i < @seq; $i++) {
	for(my $j = 1; $j <= 5; $j++) {
		$factors[$i][$j] = atchley_factor($seq[$i], $j) * $scalar[$j] + $offset[$j];
	}
}

my $seq_len = length($seq);
# Determine composition of protein by amino acid type
my @comp_counts = ();
for(my $i = 0; $i < length($aas); $i++) {
  $comp_counts[$i] = 0;
}
for(my $i = 0; $i < @seq; $i++) {
  $comp_counts[index($aas, $seq[$i])]++;
}
for(my $i = 0; $i < length($aas); $i++) {
  $comp_counts[$i] = $comp_counts[$i]/$seq_len;
}

####################################################################################################
open PSSM, "<" . $pssm_fname or die "Couldn't open pssm file $pssm_fname ".$!."\n";
@lines = <PSSM>;
chomp(@lines);
close PSSM;

my ($odds_ref, $inf_ref) = scaled_log_odds_n_inf(@lines);
my @odds = @{$odds_ref};
my @inf = @{$inf_ref};
my @joint_entropy = estimate_joint_entropy(@lines);
my @pssm_counts = pssm_counts(@lines);

my @pssm_sums = ();
for(my $i = 0; $i < @pssm_counts; $i++) {
  my $sum = 0; 
  for(my $j = 0; $j < @{$pssm_counts[$i]}; $j++) {
    $sum += $pssm_counts[$i][$j];
  }
  $pssm_sums[$i] = $sum;
}

####################################################################################################
open SS_SA, "<" . $ss_sa_fname or die "Couldn't open ss_sa file\n";
my @ss_sa=<SS_SA>;
chomp @ss_sa;
my @ss = split(//, $ss_sa[2]);
my @sa = split(//, $ss_sa[3]);
close SS_SA;

my $beta_count = 0; my $alpha_count = 0;
foreach my $ss_t (@ss) {
  if($ss_t eq 'E') {
    $beta_count++;
  } 
  if($ss_t eq 'H') {
    $alpha_count++;
  }
}

my $exposed_count = 0;
foreach my $sa_t (@sa) {
  if($sa_t eq 'b') {
    $exposed_count++;
  } 
}

####################################################################################################
my @psipredss;
open INPUT, $psipred_fname or confess $!;
while(<INPUT>){
	$_ =~ s/^\s+//;
	next if $_ !~ m/^[0-9]/;
	my @columns = split(/\s+/, $_);
	$psipredss[$columns[0]][0] = $columns[3];
	$psipredss[$columns[0]][1] = $columns[4];
	$psipredss[$columns[0]][2] = $columns[5];
}
close INPUT;

####################################################################################################
my @psisolv;
open INPUT, $psisolv_fname or confess $!;
while(<INPUT>){
	$_ =~ s/^\s+//;
	next if $_ !~ m/^[0-9]/;
	my @columns = split(/\s+/, $_);
	$psisolv[$columns[0]] = $columns[2];
}
close INPUT;

####################################################################################################
# Initialize coevolutionary features to zero
my %con_ccmpre = %{all_zero_2D_features(0)};
my %con_frecon = %{all_zero_2D_features(1)};
my %con_psicov = %{all_zero_2D_features(1)};
my %pstat_pots = %{all_zero_2D_features(1)};
my %pstat_mimt = %{all_zero_2D_features(1)};
my %pstat_mip  = %{all_zero_2D_features(1)};
my $colstatrow = "1 1";
my @colstat    = split /\s+/, $colstatrow;

# Obtain co-evolution features
if (-f $ccmpred_fname){
	%con_ccmpre = %{ccmpred2hash($ccmpred_fname)}; 
	%con_frecon = %{freecontact2hash($freecontact_fname)};
	%con_psicov = %{psicov2hash($psicov_fname)};
	%pstat_pots = %{pairstat2hash($pairstat_fname, "potsum")};
	%pstat_mimt = %{pairstat2hash($pairstat_fname, "mimat")};
	%pstat_mip  = %{pairstat2hash($pairstat_fname, "mip")}; 
	$colstatrow = colstatfeatures($colstat_fname);
	@colstat    = split /\s+/, $colstatrow;
}
else{
	print STDERR 'Coevolutionary features absent for '. $ccmpred_fname;
}

####################################################################################################
print "# Sequence Length (log)\n";
printf "%.4f\n", log($seq_len);
print "# alignment-count (log)\n";
printf "%.6f\n", log($colstat[0]);
print "# effective-alignment-count (log)\n";
printf "%.6f\n", log($colstat[1]);  
print "# Relative 'b' count\n";
printf "%.2f\n", $exposed_count/$seq_len;
print "# Relative 'H' count\n";
printf "%.2f\n", $alpha_count/$seq_len;
print "# Relative 'E' count\n";
printf "%.2f\n", $beta_count/$seq_len;  

####################################################################################################
print "# AA composition\n";
for(my $i = 0; $i < length($aas); $i++) {
	printf " %.3f", $comp_counts[$i];
	print "\n";
}

####################################################################################################
print "# Atchley factors\n";
for(my $factor = 1; $factor <= 5; $factor++) {
	for(my $i = 0; $i < $seq_len; $i++) {
		printf " %.2f", $factors[$i][$factor];
	}
	print "\n";
}

####################################################################################################
print "# Secondary Structure\n";
my @sslist = ('H', 'E', 'C');
my $sum_HCE = 0;
foreach my $sschar(@sslist){
	for(my $i = 0; $i < $seq_len; $i++){
		if (not defined $ss[$i]){
			warn "Secondary structure not defined for position ".$i."\n";
			warn "".(@ss)."\n";
			exit 1;
		}
		if ($ss[$i] eq $sschar){
			print " 1";
			$sum_HCE ++;
		}
		else{
			print " 0";
		}
	}
	print "\n";
}
if ($sum_HCE != $seq_len){
	print "SS H/E/C feature has some issues, all columns must be either H/E/C!";
	exit 1;
}

####################################################################################################
print "# Solvent accessibility\n"; # 'exposed check' flag
for(my $i = 0; $i < $seq_len; $i++) {
	if ($sa[$i] eq 'e'){
		print " 1";
	}
	else{
		print " 0";
	}
}
print "\n";

####################################################################################################
print "# PSSM inf feature\n";
for(my $i = 0; $i < $seq_len; $i++) {
	print " ".$inf[$i];
}
print "\n";

####################################################################################################
print "# PSSM\n";
for(my $l = 0; $l < @{$odds[0]}; $l++) {
	for(my $i = 0; $i < $seq_len; $i++) {
		#need to make sure that this is correct!!
		my $xx = $odds[$i][$l];
		$xx = 0 if not defined $xx;
		print " $xx";
	}
	print "\n";
}

####################################################################################################
print "# PSSM Sums (divided by 100)\n";
for(my $i = 0; $i < $seq_len; $i++) {
	printf " %.4f", $pssm_sums[$i]/100;
}
print "\n";

####################################################################################################
print "# PSSM sum cosines\n";
for(my $i = 0; $i < $seq_len; $i++) {
	for(my $j = 0; $j < $seq_len; $j++) {
		if($pssm_sums[$i] > 0 && $pssm_sums[$j] > 0) {
			printf " %.2f", (calc_cosine($pssm_counts[$i], $pssm_counts[$j]));
		}
		else{
			print " 0";
		}
	}
}
print "\n";

####################################################################################################
print "# Relative sequence separation\n";
for(my $i = 0; $i < $seq_len; $i++) {
	for(my $j = 0; $j < $seq_len; $j++) {
		printf " %.2f", (abs($j - $i)/$seq_len);
	}
}
print "\n";

####################################################################################################
print "# Sequence separation between 23 and 28\n";
for(my $i = 0; $i < $seq_len; $i++) {
	for(my $j = 0; $j < $seq_len; $j++) {
		my $flag = 0;
		if (abs($j - $i) >= 23 and abs($j - $i) < 28){
			$flag = 1;
		}
		printf " $flag";
	}
}
print "\n";
print "# Sequence separation between 28 and 38\n";
for(my $i = 0; $i < $seq_len; $i++) {
	for(my $j = 0; $j < $seq_len; $j++) {
		my $flag = 0;
		if (abs($j - $i) >= 28 and abs($j - $i) < 38){
			$flag = 1;
		}
		printf " $flag";
	}
}
print "\n";
print "# Sequence separation between 38 and 48\n";
for(my $i = 0; $i < $seq_len; $i++) {
	for(my $j = 0; $j < $seq_len; $j++) {
		my $flag = 0;
		if (abs($j - $i) >= 38 and abs($j - $i) < 48){
			$flag = 1;
		}
		printf " $flag";
	}
}
print "\n";
print "# Sequence separation 48+\n";
for(my $i = 0; $i < $seq_len; $i++) {
	for(my $j = 0; $j < $seq_len; $j++) {
		my $flag = 0;
		if (abs($j - $i) >= 48){
			$flag = 1;
		}
		printf " $flag";
	}
}
print "\n";

####################################################################################################
print "# Psipred\n";
for(my $i = 0; $i < 3; $i++){
	for(my $j = 1; $j <= $seq_len; $j++) {
		printf " %.3f", $psipredss[$j][$i];
	}
	print "\n";
}

####################################################################################################
print "# Psisolv\n";
for(my $i = 1; $i <= $seq_len; $i++) {
	printf " %.3f", $psisolv[$i];
}
print "\n";

####################################################################################################
print "# pref score\n";
for(my $i = 0; $i < $seq_len; $i++){ 
	for(my $j = 0; $j < $seq_len; $j++){ 
		my $xx = scld_residue_residue_contacts($seq[$i], $seq[$j]);
		$xx = 0 if not defined $xx;
		printf " %.3f", $xx;
	}
}
print "\n";

####################################################################################################
print "# scld lu con pot\n";
for(my $i = 0; $i < $seq_len; $i++){ 
	for(my $j = 0; $j < $seq_len; $j++){ 
		my $xx = scld_lu_contact_potential($seq[$i], $seq[$j]);
		$xx = 0 if not defined $xx;
		printf " %.3f", $xx;
	}
}
print "\n";

####################################################################################################
print "# levitt con pot\n";
for(my $i = 0; $i < $seq_len; $i++){ 
	for(my $j = 0; $j < $seq_len; $j++){ 
		my $xx = scld_levitt_contact_potential($seq[$i], $seq[$j]);
		$xx = 0 if not defined $xx;
		printf " %.3f", $xx;
	}
}
print "\n";

####################################################################################################
print "# braun con pot\n";
for(my $i = 0; $i < $seq_len; $i++){ 
	for(my $j = 0; $j < $seq_len; $j++){ 
		my $xx = scld_braun_contact_potential($seq[$i], $seq[$j]);
		$xx = 0 if not defined $xx;
		printf " %.3f", $xx;
	}
}
print "\n";

####################################################################################################
print "# joint entro\n";
for(my $i = 0; $i < $seq_len; $i++){ 
	for(my $j = 0; $j < $seq_len; $j++){ 
		printf " %.3f", $joint_entropy[$i][$j];
	}
}
print "\n";

####################################################################################################
print "# pearson r\n";
for(my $i = 0; $i < $seq_len; $i++){ 
	for(my $j = 0; $j < $seq_len; $j++){ 
		printf " %.3f", (calc_pearson_r($pssm_counts[$i],$pssm_counts[$j])+1)/2.0;
	}
}
print "\n";

####################################################################################################
print "# Shannon entropy sum\n";
my @entropy = colstat_entropy($colstat_fname);
for(my $i = 21; $i <= 21; $i++){
	for(my $j = 0; $j < $seq_len; $j++) {
		printf " %.3f", $entropy[$i][$j];
	}
	print "\n";
}

####################################################################################################
print "# ccmpred\n";
for(my $i = 0; $i < $seq_len; $i++){ 
	for(my $j = 0; $j < $seq_len; $j++){ 
		my $xx = $con_ccmpre{$i." ".$j};
		$xx = 0 if not defined $xx;
		printf " %.4f", $xx;
	}
}
print "\n";

####################################################################################################
print "# freecontact\n";
for(my $i = 1; $i <= $seq_len; $i++){ 
	for(my $j = 1; $j <= $seq_len; $j++){ 
		my $xx = $con_frecon{$i." ".$j};
		$xx = 0 if not defined $xx;
		$xx = 0 if $xx < 0;
		printf " %.4f", $xx;
	}
}
print "\n";

####################################################################################################
print "# psicov\n";
for(my $i = 1; $i <= $seq_len; $i++){ 
	for(my $j = 1; $j <= $seq_len; $j++){ 
		my $xx = $con_psicov{$i." ".$j};
		$xx = 0 if not defined $xx;
		$xx = 0 if $xx < 0;
		printf " %.4f", $xx;
	}
}
print "\n";

####################################################################################################
print "# pstat_pots\n";
for(my $i = 1; $i <= $seq_len; $i++){ 
	for(my $j = 1; $j <= $seq_len; $j++){ 
		my $xx = $pstat_pots{$i." ".$j};
		$xx = 0 if not defined $xx;
		printf " %.4f", (1 + exp(-$xx)) ** -1;
	}
}
print "\n";

####################################################################################################
print "# pstat_mimt\n";
for(my $i = 1; $i <= $seq_len; $i++){ 
	for(my $j = 1; $j <= $seq_len; $j++){ 
		my $xx = $pstat_mimt{$i." ".$j};
		$xx = 0 if not defined $xx;
		printf " %.4f", (1 + exp(-$xx)) ** -1;
	}
}
print "\n";

####################################################################################################
print "# pstat_mip\n";
for(my $i = 1; $i <= $seq_len; $i++){ 
	for(my $j = 1; $j <= $seq_len; $j++){ 
		my $xx = $pstat_mip{$i." ".$j};
		$xx = 0 if not defined $xx;
		printf " %.4f", (1 + exp(-$xx)) ** -1;
	}
}
print "\n";

####################################################################################################
sub ccmpred2hash{
	my $file_ccmpred = shift;
	die "ERROR! file_ccmpred not defined!" if !$file_ccmpred;
	my %conf = ();
	open CCM, $file_ccmpred or die $!." $file_ccmpred";
	my $i = 0;
	while(<CCM>){
		$_ =~ s/^\s+//;
		my @C = split /\s+/, $_;
		for(my $j = 0; $j <= $#C; $j++){
			$conf{$i." ".$j} = $C[$j];
			$conf{$j." ".$i} = $C[$j];
		}
		$conf{$i." ".$i} = 1;
		$i++;
	}
	close CCM;
	return \%conf;
}

####################################################################################################
sub psicov2hash{
	my $file_psicov = shift;
	die "ERROR! file_psicov not defined!" if !$file_psicov;
	my %conf = ();
	open PSI, $file_psicov or die $!." $file_psicov";
	while(<PSI>){
		$_ =~ s/\r//g;
		$_ =~ s/^\s+//;
		next unless $_ =~ /^[0-9]/;
		my @C = split /\s+/, $_;
		$conf{$C[0]." ".$C[1]} = $C[4];
		$conf{$C[1]." ".$C[0]} = $conf{$C[0]." ".$C[1]};
		$conf{$C[0]." ".$C[0]} = 1;
		$conf{$C[1]." ".$C[1]} = 1;
	}
	close PSI;
	return \%conf;
}

####################################################################################################
sub freecontact2hash{
	my $file_fc = shift;
	die "ERROR! file_freecontact not defined!" if !$file_fc;
	my %conf = ();
	open FC, $file_fc or die $!." $file_fc";
	while(<FC>){
		$_ =~ s/\r//g;
		$_ =~ s/^\s+//;
		next unless $_ =~ /^[0-9]/;
		my @C = split /\s+/, $_;
		$conf{$C[0]." ".$C[2]} = $C[5];
		$conf{$C[2]." ".$C[0]} = $C[5];
		$conf{$C[0]." ".$C[0]} = 1;
		$conf{$C[2]." ".$C[2]} = 1;
	}
	close FC;
	return \%conf;
}

####################################################################################################
sub pairstat2hash{
	my $file_pairstat = shift;
	my $option = shift;
	die "ERROR! file_freecontact not defined!" if !$file_pairstat;
	my %pairs = ();
	open PS, $file_pairstat or die $!." $file_pairstat";
	while(<PS>){
		$_ =~ s/\r//g;
		$_ =~ s/^\s+//;
		next unless $_ =~ /^[0-9]/;
		my @C = split /\s+/, $_;
		my $value = 0;
		$value = $C[2] if $option eq "potsum"; # mean contact potential
		$value = $C[3] if $option eq "mimat";  # mutual information
		$value = $C[4] if $option eq "mip";    # normalized mutual information
		$pairs{$C[0]." ".$C[1]} = $value;
		$pairs{$C[1]." ".$C[0]} = $value;
	}
	close PS;
	return \%pairs;
}

####################################################################################################
sub colstatfeatures{
	my $file_colstat = shift;
	my $option = shift;
	die "ERROR! file_freecontact not defined!" if !$file_colstat;
	open CS, $file_colstat or die $!." $file_colstat";
	my @lines = <CS>;
	close CS;
	chomp @lines;
	my $seqlen      = $lines[0];
	my $alignlen    = $lines[1];
	my $effalignlen = $lines[2];
	if(not defined $seqlen){
		return "0 0";
	}
	$alignlen = 0 if not defined $alignlen;
	$effalignlen = 0 if not defined $effalignlen;
	if($seqlen ne length($seq)){
		print STDERR "ERROR! Fasta file reports seqlen is ".length($seq)." but colstat feature file reports ".$seqlen." [fasta: $fasta_fname]\n";
	}
	return ($alignlen)." ".($effalignlen);
}

####################################################################################################
sub colstat_entropy{
	my $file_colstat = shift;
	my $option = shift;
	die "ERROR! file_freecontact not defined!" if !$file_colstat;
	open CS, $file_colstat or die $!." $file_colstat";
	my @lines = <CS>;
	close CS;
	my @aaentropy;
	for (my $l = 4; $l < $seq_len + 4; $l++){
		my @aacomp = split /\s+/, $lines[$l];
		for (my $i = 0; $i <= $#aacomp; $i++){
			$aaentropy[$i][$l - 4] = $aacomp[$i];
			if ($aacomp[$i] eq ""){
				print STDERR "space in aacomp.. something is wrong.. $l $i";
				exit 1;
			}
		}
	}
	return @aaentropy;
}

####################################################################################################
sub all_zero_2D_features{
	my $start_at = shift;
	my $start = 0;
	my $end = $seq_len;
	$start = 1 if $start_at == 0;
	$end   = $end + 1 if $start_at == 0;
	my %feat = ();
	for(my $i = $start; $i < $end; $i++){ 
		for(my $j = $start; $j < $end; $j++){ 
			$feat{$i." ".$j} = 0
		}
	}
	return \%feat;
}
