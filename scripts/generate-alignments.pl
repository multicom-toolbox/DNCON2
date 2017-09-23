#!/usr/bin/perl -w
# Badri Adhikari, 5-21-2017

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;

use constant{
	JACKHMMER   => '/home/badri/hmmer-3.1b2-linux-intel-x86_64/binaries/jackhmmer',
	REFORMAT    => abs_path(dirname($0)).'/reformat.pl',
	JACKHMMERDB => '/home/badri/databases/uniref/uniref90pfilt',
	HHBLITS     => '/usr/bin/hhblits',
	HHBLITSDB   => '/home/badri/databases/uniprot20_2016_02/uniprot20_2016_02',
	CPU         => 2
};

my %COVERAGE = ();
$COVERAGE{60} = 1;

my %EVALUE = ();
$EVALUE{'1e-20'} = 1;
$EVALUE{'1e-10'} = 2;
$EVALUE{'1e-4'}  = 3;
$EVALUE{'1'}     = 4;

confess 'Oops!! jackhmmer not found!'   if not -f JACKHMMER;
confess 'Oops!! reformat not found!'    if not -f REFORMAT;
confess 'Oops!! jackhmmerdb not found!' if not -f JACKHMMERDB;
confess 'Oops!! hhblits not found!'     if not -f HHBLITS;
confess 'Oops!! hhblitsdb not found!'   if not -f HHBLITSDB.'_a3m_db';

####################################################################################################
my $fasta  = shift;
my $outdir = shift;

if (not $fasta or not -f $fasta){
	print "Fasta file $fasta does not exist!\n" if ($fasta and not -f $fasta);
	print "Usage: $0 <fasta> <output-directory>\n";
	exit(1);
}

if (not $outdir){
	print 'Output directory not defined!';
	print "Usage: $0 <fasta> <output-directory>\n";
	exit(1);
}

####################################################################################################
my $id = basename($fasta, ".fasta");
system_cmd("mkdir -p $outdir") if not -d $outdir;
system_cmd("cp $fasta $outdir/") if not -f $outdir."/$id.fasta";
chdir $outdir or confess $!;
$fasta = basename($fasta);
my $seq = seq_fasta($fasta);

# check and quit, if there are any results already
my $existing = `find . -name "*.aln" | wc -l`;
$existing = 0 if not $existing;
confess 'Oops!! There are already some alignment file in the ouput directory! Consider running in an empty directory!' if int($existing) > 0;

####################################################################################################
print "Started [$0]: ".(localtime)."\n";

my %jobs = ();
foreach my $c (keys %COVERAGE){
	my $hhbid = "hhb-cov".$c;
	open  JOB, ">$hhbid.sh" or confess "ERROR! Could not open $hhbid.sh $!";
	print JOB "#!/bin/bash\n";
	print JOB "touch $hhbid.running\n";
	print JOB "echo \"running hhblits job $hhbid..\"\n";
	print JOB HHBLITS." -i $fasta -d ".HHBLITSDB." -oa3m $hhbid.a3m -cpu ".CPU." -n 3 -maxfilt 500000 -diff inf -e 0.001 -id 99 -cov $c > $hhbid-hhblits.log\n";
	print JOB "if [ ! -f \"${hhbid}.a3m\" ]; then\n";
	print JOB "   mv $hhbid.running $hhbid.failed\n";
	print JOB "   echo \"hhblits job $hhbid failed!\"\n";
	print JOB "   exit\n";
	print JOB "fi\n";
	print JOB "egrep -v \"^>\" $hhbid.a3m | sed 's/[a-z]//g' > $hhbid.aln\n";
	print JOB "if [ -f \"${hhbid}.aln\" ]; then\n";
	print JOB "   mv $hhbid.running $hhbid.done\n";
	print JOB "   echo \"hhblits $hhbid job done.\"\n";
	print JOB "   exit\n";
	print JOB "fi\n";
	print JOB "echo \"Something went wrong! $hhbid.aln file not present!\"\n";
	print JOB "mv $hhbid.running $hhbid.failed\n";
	close JOB;
	system_cmd("chmod 755 $hhbid.sh");
	$jobs{$hhbid.".sh"} = 1;
}

foreach my $job (sort keys %jobs){
	print "Starting job $job ..\n";
	system "./$job &";
	sleep 1;
}

####################################################################################################
print("Wait until all HHblits jobs are done ..\n");
my $running = `find . -name "*.running" | wc -l`;
chomp $running;
confess 'Oops!! Something went wrong! No jobs are running!' if (int($running) < 0);
print "$running jobs running currently\n";
while (int($running) > 0){
	sleep 2;
	my $this_running = `find . -name "*.running" | wc -l`;
	chomp $this_running;
	$this_running = 0 if not $this_running;
	if(int($this_running) != $running){
		print "$this_running jobs running currently\n";
	}
	$running = $this_running;
}

####################################################################################################
print "\nAlignment Summary:\n";
print 'L = '.length($seq)."\n";
system "wc -l *.aln";
print "\n";

# Apply alignment selection rule to select the best alignment file as $id.aln
# Increasing this threshold to 5 (from 2.5) after observing the case of T0855 where e-10 has 321 rows and e-4 has 11K rows
my $T = 5 * length($seq);
my $found_aln = 0;
foreach my $c (sort {$COVERAGE{$a} <=> $COVERAGE{$b}} keys %COVERAGE){
	last if $found_aln;
	my $hhbid = "hhb-cov".$c;
	confess "Oops!! Expected file $hhbid.aln not found!" if not -f "$hhbid.aln";
	if (count_lines("$hhbid.aln") > $T){
		print("Copying $hhbid.aln as $id.aln\n");
		system_cmd("echo \"cp $hhbid.aln $id.aln\" > result.txt");
		system_cmd("cp $hhbid.aln $id.aln");
		$found_aln = 1;
		last;
	}
}

if($found_aln){
	print "HHblits jobs have enough alignments! Not running JackHmmer!\n";
	print "\nFinished [$0]: ".(localtime)."\n";
	exit 0;
}

####################################################################################################
# JackHmmer needs fasta sequence to be in a single line
open FASTA, ">$id.jh.fasta" or confess "ERROR! Could not open $id.jh.fasta $!";
print FASTA ">$id.jh.fasta\n";
print FASTA "$seq\n";
close FASTA;

%jobs = ();
foreach my $e (keys %EVALUE){
	my $jhmid = "jhm-".$e;
	$jhmid = "jhm-e-0" if $e eq '1';
	open  JOB, ">$jhmid.sh" or confess "ERROR! Could not open $jhmid.sh $!";
	print JOB "#!/bin/bash\n";
	print JOB "touch $jhmid.running\n";
	print JOB "echo \"running jackhmmer job $jhmid..\"\n";
	print JOB JACKHMMER.' --cpu '.CPU." -N 5 -E $e -A $jhmid.ali $id.jh.fasta ".JACKHMMERDB." > $jhmid-jackhmmer.log\n";
	print JOB "if [ ! -f \"${jhmid}.ali\" ]; then\n";
	print JOB "   mv $jhmid.running $jhmid.failed\n";
	print JOB "   echo \"jackhmmer job $jhmid failed!\"\n";
	print JOB "   exit\n";
	print JOB "fi\n";
	print JOB REFORMAT." -l 1500 -d 1500 sto a3m $jhmid.ali $jhmid.a3m\n";
	print JOB "egrep -v \"^>\" $jhmid.a3m | sed 's/[a-z]//g' > $jhmid.aln\n";
	print JOB "if [ -f \"${jhmid}.aln\" ]; then\n";
	print JOB "   rm $jhmid.running\n";
	# jackhmmer log files and .ali files use up a lot of space
	print JOB "   rm $jhmid-jackhmmer.log\n";
	print JOB "   rm $jhmid.ali\n";
	print JOB "   echo \"jackhmmer job $jhmid done.\"\n";
	print JOB "   exit\n";
	print JOB "fi\n";
	print JOB "echo \"Something went wrong! $jhmid.aln file not present!\"\n";
	print JOB "mv $jhmid.running $jhmid.failed\n";
	close JOB;
	system_cmd("chmod 755 $jhmid.sh");
	$jobs{$jhmid.".sh"} = 1;
}

foreach my $job (sort keys %jobs){
	print "Starting job $job ..\n";
	system "./$job &";
	sleep 1;
}

####################################################################################################
print("Wait until all JackHmmer jobs are done ..\n");
$running = `find . -name "*.running" | wc -l`;
chomp $running;
confess 'Oops!! Something went wrong! No jobs are running!' if (int($running) < 0);
print "$running jobs running currently\n";
while (int($running) > 0){
	sleep 2;
	my $this_running = `find . -name "*.running" | wc -l`;
	chomp $this_running;
	$this_running = 0 if not $this_running;
	if(int($this_running) != $running){
		print "$this_running jobs running currently\n";
	}
	$running = $this_running;
}

####################################################################################################
print "\nAlignment Summary:\n";
print 'L = '.length($seq)."\n";
system "wc -l *.aln";
print "\n";

####################################################################################################
# Apply alignment selection rule to select the best alignment file as $id.aln
foreach my $e (sort {$EVALUE{$a} <=> $EVALUE{$b}} keys %EVALUE){
	last if $found_aln;
	my $jhmid = "jhm-".$e;
	$jhmid = "jhm-e-0" if $e eq '1';
	confess "Oops!! Expected file $jhmid.aln not found!" if not -f "$jhmid.aln";
	if (count_lines("$jhmid.aln") > $T){
		print("Copying $jhmid.aln as $id.aln\n");
		system_cmd("echo \"cp $jhmid.aln $id.aln\" > result.txt");
		system_cmd("cp $jhmid.aln $id.aln");
		$found_aln = 1;
		last;
	}
}
if(not $found_aln){
	my $jhmid = "jhm-e-0";
	confess "Oops!! Expected file $jhmid.aln not found!" if not -f "$jhmid.aln";
	print("Copying $jhmid.aln as $id.aln\n");
	system_cmd("echo \"cp $jhmid.aln $id.aln\" > result.txt");
	system_cmd("cp $jhmid.aln $id.aln");
}
if (count_lines("$id.aln") > 75000){
	print("More than 75,000 rows in the alignment file.. trimming..\n");
	system_cmd("head -75000 $id.aln > temp.aln");
	system_cmd("rm $id.aln");
	system_cmd("mv temp.aln $id.aln");
}

####################################################################################################
print "Check sequences that are shorter and throw them away..\n";
my $L = length($seq);
open ALN, "$id.aln" or confess $!;
open TEMP, ">temp.aln" or confess $!;
while (<ALN>){
	chomp $_;
	if (length($_) != $L){
		print "Skipping - $_\n";
		next; 
	}
	print TEMP $_."\n";
}
close TEMP;
close ALN;

system_cmd("mv temp.aln $id.aln");

print "\nFinished [$0]: ".(localtime)."\n";

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

####################################################################################################
sub seq_fasta{
	my $file_fasta = shift;
	confess "ERROR! Fasta file $file_fasta does not exist!" if not -f $file_fasta;
	my $seq = "";
	open FASTA, $file_fasta or confess $!;
	while (<FASTA>){
		next if (substr($_,0,1) eq ">"); 
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$seq .= $_;
	}
	close FASTA;
	return $seq;
}

####################################################################################################
sub count_lines{
	my $file = shift;
	my $lines = 0;
	return 0 if not -f $file;
	open FILE, $file or confess "ERROR! Could not open $file! $!";
	while (<FILE>){
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		next if not defined $_;
		next if length($_) < 1;
		$lines ++;
	}
	close FILE;
	return $lines;
}
