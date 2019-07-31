#!/usr/bin/perl -w
# Badri Adhikari, 5-21-2017
# Tianqi Wu modified alignments pipeline, 4-25-2018

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use Cwd qw();

use constant{
	JACKHMMER   => '/storage/htc/bdm/tools/hmmer-3.1b2-linux-intel-x86_64/binaries/jackhmmer',
	REFORMAT    => abs_path(dirname($0)).'/reformat.pl',
	JACKHMMERDB => '/storage/htc/bdm/tools/databases/uniref90_04_2018/uniref90',
	JACK_HH => abs_path(dirname($0)).'/jack_hhblits.pl',
	Metapsicov => '/storage/htc/bdm/tools/metapsicov-2.0.3/bin',
	ALNSTAT      => '/storage/htc/bdm/tools/metapsicov-2.0.3/bin/alnstats',
	HHBLITS     => '/storage/htc/bdm/tools/hhsuite-3.0-beta.1/bin/hhblits',
	HHBLITS2     => '/storage/htc/bdm/tools/hhsuite-2.0.16-linux-x86_64',
	HHBLITSDB   => '/storage/htc/bdm/tools/databases/uniclust30_2017_10/uniclust30_2017_10',
	CPU         => 8
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
confess 'Oops!! jack_hh scripts not found!' if not -f JACK_HH;
confess 'Oops!! alnstat not found!'     if not -f ALNSTAT;
confess 'Oops!! folder metapsicov not found!' if not -d Metapsicov;
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

my $path = Cwd::cwd();
print "current path:$path\n";
$fasta = basename($fasta);
my $seq = seq_fasta($fasta);

# check and quit, if there are any results already
my $existing = `find . -name "*.aln" | wc -l`;
$existing = 0 if not $existing;
confess 'Oops!! There are already some alignment file in the ouput directory! Consider running in an empty directory!' if int($existing) > 0;

####################################################################################################
print "Started [$0]: ".(localtime)."\n";

my ($jhmid,$hhbid);
my $evalue = '1';
my %jobs = ();

foreach my $c (keys %COVERAGE){
	$hhbid = "hhb-cov".$c;
	open  JOB, ">$hhbid.sh" or confess "ERROR! Could not open $hhbid.sh $!";
	print JOB "#!/bin/bash\n";
	print JOB " export HHLIB=".HHBLITS."/build\n";
	print JOB "PATH=\$PATH:\$HHLIB/bin:\$HHLIB/scripts\n";
	print JOB "touch $hhbid.running\n";
	print JOB "echo \"running hhblits job $hhbid..\"\n";
	print JOB HHBLITS." -i $fasta -d ".HHBLITSDB." -oa3m $id.a3m -cpu ".CPU." -n 3 -maxfilt 500000 -diff inf -e 0.001 -id 99 -cov $c > $hhbid-hhblits.log\n";
	print JOB "export HHLIB=".HHBLITS2."/lib/hh\n";
	print JOB "PATH=\$PATH:".HHBLITS2."/bin:\$HHLIB/scripts\n";
	print JOB "cp $hhbid.a3m $id.a3m\n";
	print JOB "if [ ! -f \"${id}.a3m\" ]; then\n";
	print JOB "   mv $hhbid.running $hhbid.failed\n";
	print JOB "   echo \"hhblits job $hhbid failed!\"\n";
	print JOB "   exit\n";
	print JOB "fi\n";
	print JOB "egrep -v \"^>\" $id.a3m | sed 's/[a-z]//g' > $hhbid.aln\n";
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
my $T = 10 * length($seq);
my $found_aln = 0;
foreach my $c (sort {$COVERAGE{$a} <=> $COVERAGE{$b}} keys %COVERAGE){
	last if $found_aln;
	my $hhbid = "hhb-cov".$c;
	confess "Oops!! Expected file $hhbid.aln not found!" if not -f "$hhbid.aln";
	if ((count_lines("$hhbid.aln") > $T) && (count_lines("$hhbid.aln")>2000)){
		print("Copying $hhbid.aln as $id.aln\n");
		system_cmd("echo \"cp $hhbid.aln $id.aln\" > result.txt");
		system_cmd("cp $hhbid.aln $id.aln");
		system_cmd("cp $hhbid.aln ${id}_orig.aln");
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
	$jhmid = "jhm-".$e;
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
my ($jhm_n,$jhm_neff,$jackaln_neff);
foreach my $e (sort {$EVALUE{$a} <=> $EVALUE{$b}} keys %EVALUE){
	last if $found_aln;
	$jhmid = "jhm-".$e;
	$jhmid = "jhm-e-0" if $e eq '1';
	confess "Oops!! Expected file $jhmid.aln not found!" if not -f "$jhmid.aln";
	if ((count_lines("$jhmid.aln") > $T) && (count_lines("$jhmid.aln")>2000)){
		print("Copying $jhmid.aln as $id.aln\n");
		system_cmd("echo \"cp $jhmid.aln $id.aln\" > result.txt");
		$found_aln = 1;
		if($e eq '1e-20'){
			$evalue = '1e-20';
		}
		if($e eq '1e-10'){
			$evalue = '1e-10';
		}
		if($e eq '1e-4'){
			$evalue = '1e-4';
		}
		if($e eq '1'){
			$evalue = '1';
		}

		if (count_lines("$jhmid.aln") > 120000){
			print("More than 120,000 rows in the alignment file.. trimming..\n");
			system_cmd("head -120000 $jhmid.aln > temp.aln");
			system_cmd("rm $jhmid.aln");
			system_cmd("mv temp.aln $jhmid.aln");
		}

    system_cmd("cp $jhmid.aln ${id}_orig.aln");

		print "Check sequences $jhmid.aln that are shorter and throw them away..\n";
		open ALN, "$jhmid.aln" or confess $!;
		open TEMP, ">temp.aln" or confess $!;
		while (<ALN>){
			chomp $_;
			if (length($_) != length($seq)){
				print "Skipping - $_\n";
				next;
			}
			print TEMP $_."\n";
		}
		close TEMP;
		close ALN;
		system_cmd("mv temp.aln $jhmid.aln");

		system_cmd(ALNSTAT." $jhmid.aln $jhmid.colstats $jhmid.pairstats");
	   open COLSTATS, "<$jhmid.colstats" or confess $!;
		while(<COLSTATS>){
			$jhm_n=<COLSTATS>;
			$jhm_neff =<COLSTATS>;
			print "jhm_neff:$jhm_neff\n";
			last;
		}
		last;
	}
}
if(not $found_aln){
	$jhmid = "jhm-e-0";
	$evalue='1';
	confess "Oops!! Expected file $jhmid.aln not found!" if not -f "$jhmid.aln";
	#print("Copying $jhmid.aln as $id.aln\n");
	system_cmd("echo \"cp $jhmid.aln $id.aln\" > result.txt");

	if (count_lines("$jhmid.aln") > 120000){
		print("More than 120,000 rows in the alignment file.. trimming..\n");
		system_cmd("head -120000 $jhmid.aln > temp.aln");
		system_cmd("rm $jhmid.aln");
		system_cmd("mv temp.aln $jhmid.aln");
	}

	system_cmd("cp $jhmid.aln ${id}_orig.aln");
	print "Check sequences $jhmid.aln that are shorter and throw them away..\n";
	open ALN, "$jhmid.aln" or confess $!;
	open TEMP, ">temp.aln" or confess $!;
	while (<ALN>){
		chomp $_;
		if (length($_) != length($seq)){
			print "Skipping - $_\n";
			next;
		}
		print TEMP $_."\n";
	}
	close TEMP;
	close ALN;
	system_cmd("mv temp.aln $jhmid.aln");

	system_cmd(ALNSTAT." $jhmid.aln $jhmid.colstats $jhmid.pairstats");
	open COLSTATS, "<$jhmid.colstats" or confess $!;
	while(<COLSTATS>){
		$jhm_n=<COLSTATS>;
		$jhm_neff =<COLSTATS>;
		print "jhm_neff:$jhm_neff\n";
		last;
	}
	#system_cmd("cp $jhmid.aln $id.aln");
}

print "evalue:$evalue\n";
if($jhm_n > 20000){
      #combine jhmid,aln with hhbidaln
		system_cmd("cp $jhmid.aln $id.jackaln");
		system_cmd("cp $jhmid.colstats $id.jackaln.colstats");
}else{
	print "evalue:$evalue\n";
    system_cmd(JACK_HH." $id ".Metapsicov." $path ".JACKHMMERDB." $evalue > $path/$id.jacklog");
    system_cmd(ALNSTAT." $id.jackaln $id.jackaln.colstats $id.jackaln.pairstats");
}

open COLSTATS, "<$id.jackaln.colstats" or confess $!;
while(<COLSTATS>){
	$jackaln_neff=<COLSTATS>;
	$jackaln_neff =<COLSTATS>;
	print "jackaln_neff:$jackaln_neff\n";
	last;
}

if ($jackaln_neff>$jhm_neff){
	if (count_lines("$hhbid.aln") < count_lines("$id.jackaln")){
			system_cmd("cp $id.jackaln $id.aln");
	}
	else{
			system_cmd("cp $hhbid.aln $id.aln");
	}
}
else{
		system_cmd("cp $jhmid.aln $id.aln");
}


if (count_lines("$id.aln") > 120000){
	print("More than 120,000 rows in the alignment file.. trimming..\n");
	system_cmd("head -120000 $id.aln > temp.aln");
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
