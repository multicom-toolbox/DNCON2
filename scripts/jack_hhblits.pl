#!/usr/bin/perl -w
# Tianqi Wu, 04-15-201

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use Cwd qw();

my $jobid  = shift;
my $bindir = shift;
my $tmpdir = shift;
my $db = shift;
my $e = shift;

use constant{
	JACKHMMER   => '/storage/htc/bdm/tools/hmmer-3.1b2-linux-intel-x86_64/',
	HHDIR    => '/storage/htc/bdm/tools/hhsuite-2.0.16-linux-x86_64'
};

chdir $tmpdir or confess $!;
system_cmd(JACKHMMER."/src/jackhmmer --cpu 8 -N 5 -E $e --noali --tblout $jobid.tbl $jobid.fasta $db");
system_cmd(JACKHMMER."/easel/miniapps/esl-sfetch -f $db $jobid.tbl > $jobid.fseqs");
system_cmd("cat $jobid.fasta >> $jobid.fseqs");
system_cmd("mkdir -p $jobid-mya3m");
system_cmd("mkdir -p $jobid-mydb");
chdir "$jobid-mya3m" or confess $!;
#system_cmd("cd $jobid-mya3m");
system_cmd("$bindir/fasta2a3msplit < $tmpdir/$jobid.fseqs");
#system_cmd("cd ..");
chdir $tmpdir or confess $!;
system_cmd("perl ".HHDIR."/lib/hh/scripts/hhblitsdb.pl -cpu 8 -o $jobid-mydb/mydb -ia3m $jobid-mya3m > /dev/null");
system_cmd(HHDIR."/bin/hhblits -i $jobid.a3m -d $jobid-mydb/mydb -oa3m $jobid.a3m -maxfilt 500000 -e 1e-3 -n 3 -cpu 8 -diff inf -id 99 -cov 50");
system_cmd("rm -rf $jobid-mya3m &");
system_cmd("rm -rf $jobid-mydb &");

system_cmd("egrep -v \"^>\" $jobid.a3m | sed 's/[a-z]//g' > $jobid.jackaln");

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
