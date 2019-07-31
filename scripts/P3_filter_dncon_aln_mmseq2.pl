 #! /usr/bin/perl -w
 # perl P3_filter_dncon_aln_mmseq2.pl  T0658.aln  T0658_filt.aln
use Carp;

$num = @ARGV;
if($num != 2)
{
  print "Usage: perl P3_filter_dncon_aln_mmseq2.pl  T0658.aln  T0658_filt.aln\n";
  die "The number of parameter is not correct!\n";
}

$in_aln = $ARGV[0];
$out_aln = $ARGV[1];
$in_aln_id  = substr($in_aln,0,index($in_aln,'.aln'));

use constant{
  MMseqs   => '/home/casp13/tools/MMseqs2/mmseqs2/bin/mmseqs'
};

open(IN,$in_aln) || die "Failed to open file $in_aln\n";
open(OUT,">${in_aln_id}_tmp.fa") || die "Failed to open file ${in_aln_id}_tmp.fa\n";
@content = <IN>;
close IN;

$orig_seq= shift @content;
chomp $orig_seq;

$ind =0;
foreach (@content)
{
  $line=$_;
  chomp $line;
  #if(index($line,'X') >= 0)
  #{
  #  die "!!! the sequence should not have X\n\n";
  #}
  #$line =~ s/\-/X/g;
  $ind++;
  print OUT ">aln_$ind\n$line\n";
}
close IN;
close OUT;

#### use mmseq to filter
system_cmd("rm -rf *_tmp_DB*");
system_cmd(MMseqs." createdb ${in_aln_id}_tmp.fa  ${in_aln_id}_tmp_DB");
system_cmd(MMseqs." linclust ${in_aln_id}_tmp_DB ${in_aln_id}_tmp_DB_90 ${in_aln_id}_tmp_DB_tmp90 --min-seq-id 0.9");

system_cmd(MMseqs." result2repseq ${in_aln_id}_tmp_DB ${in_aln_id}_tmp_DB_90 ${in_aln_id}_tmp_DB_90_seq");

system_cmd(MMseqs." result2flat ${in_aln_id}_tmp_DB ${in_aln_id}_tmp_DB ${in_aln_id}_tmp_DB_90_seq ${in_aln_id}_tmp_DB_90_seq.fasta");

open(OUT,">$out_aln") || die "Failed to open file $out_aln\n";
print OUT "$orig_seq\n";
open(IN,"${in_aln_id}_tmp_DB_90_seq.fasta") || die "Failed to open file ${in_aln_id}_tmp_DB_90_seq.fasta\n";
while(<IN>)
{
  $line=$_;
  chomp $line;
  if(substr($line,0,1) eq '>')
  {
    next;
  }

  print OUT "$line\n";
}

close OUT;

system_cmd("rm -rf  ${in_aln_id}_tmp*");

####################################################################################################
sub system_cmd{
	my $command = shift;
	my $log = shift;
	confess "EXECUTE [$command]?\n" if (length($command) < 5  and $command =~ m/^rm/);
	if(defined $log){
		system("$command &> $log");
	}
	else{
		system($command);
	}
	if($? != 0){
		my $exit_code  = $? >> 8;
		confess "ERROR!! Could not execute [$command]! \nError message: [$!]";
	}
}
