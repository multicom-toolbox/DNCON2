#!/usr/bin/perl -w
# Badri Adhikari, 5-21-2017

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use LWP::UserAgent;
use Time::Piece;

####################################################################################################
use constant{
	FREECONTACT=> '/usr/bin/freecontact',
	PSICOV    => '/home/badri/psicov/psicov',
	CCMPRED   => '/home/badri/CCMpred/bin/ccmpred',
	HOURLIMIT => 24,
	NPROC     => 8
};
my %PSICOVOPTIONS = ('d0.03' => '-o -d 0.03', 'r0.001' => '-o -r 0.001', 'r0.01' => '-o -r 0.01');

confess 'Oops!! psicov not found at '.PSICOV.' !' if not -f PSICOV;
confess 'Oops!! ccmpred qq not found!'.CCMPRED if not -f CCMPRED;
confess 'Oops!! freecontact not found!'.FREECONTACT if not -f FREECONTACT;

####################################################################################################
my $aln            = shift;
my $psicovdir      = shift;
my $ccmpreddir     = shift;
my $freecontactdir = shift;

if (not $aln or not -f $aln){
	print "Alignment file $aln does not exist!\n" if ($aln and not -f $aln);
	print "Usage: $0 <aln-file> <psicov-output-directory> <ccmpred-output-directory> <freecontact-output-directory>\n";
	exit(1);
}
$aln = abs_path($aln);
if (not $psicovdir){
	print 'PSICOV Output directory not defined!';
	print "Usage: $0 <aln-file> <psicov-output-directory> <ccmpred-output-directory> <freecontact-output-directory>\n";
	exit(1);
}
system_cmd("mkdir -p $psicovdir");
$psicovdir = abs_path($psicovdir);

if (not $ccmpreddir){
	print 'CCMpred Output directory not defined!';
	print "Usage: $0 <aln-file> <psicov-output-directory> <ccmpred-output-directory> <freecontact-output-directory>\n";
	exit(1);
}
system_cmd("mkdir -p $ccmpreddir");
$ccmpreddir = abs_path($ccmpreddir);

if (not $freecontactdir){
	print 'FreeContact Output directory not defined!';
	print "Usage: $0 <aln-file> <psicov-output-directory> <ccmpred-output-directory> <freecontact-output-directory>\n";
	exit(1);
}
system_cmd("mkdir -p $freecontactdir");
$freecontactdir = abs_path($freecontactdir);

####################################################################################################
my $id = basename($aln, ".aln");
$aln = abs_path($aln);

# check and quit, if there are any results already
my $existing = `find $psicovdir -name "*.rr" | wc -l`;
$existing = 0 if not $existing;
confess 'Oops!! There are already some rr files in the PSICOV ouput directory! Consider running in an empty directory!' if int($existing) > 0;
$existing = `find $ccmpreddir -name "*.rr" | wc -l`;
$existing = 0 if not $existing;
confess 'Oops!! There are already some rr files in the CCMPRED ouput directory! Consider running in an empty directory!' if int($existing) > 0;

####################################################################################################
print "Started [$0]: ".(localtime)."\n";

foreach my $opt (keys %PSICOVOPTIONS){
	my $psicovopt = "$id-".$opt;
	system_cmd("mkdir -p $psicovdir/$opt");
	chdir "$psicovdir/$opt" or confess $!;
	system_cmd("cp $aln ./");
	open  JOB, ">$psicovopt.sh" or confess "ERROR! Could not open $psicovopt.sh $!";
	print JOB "#!/bin/bash\n";
	print JOB "touch $psicovopt.running\n";
	print JOB "echo \"running $psicovopt ..\"\n";
	print JOB "date\n";
	print JOB PSICOV." ".$PSICOVOPTIONS{$opt}." $id.aln > ".$psicovopt.".psicov\n";
	print JOB "if [ -s \"${psicovopt}.psicov\" ]; then\n";
	print JOB "   mv $psicovopt.running $psicovopt.done\n";
	print JOB "   echo \"$psicovopt job done.\"\n";
	print JOB "   date\n";
	print JOB "   exit\n";
	print JOB "fi\n";
	print JOB "mv $psicovopt.running $psicovopt.failed\n";
	print JOB "echo \"psicov job $psicovopt failed!\"\n";
	print JOB "date\n";
	close JOB;
	system_cmd("chmod 755 $psicovopt.sh");
	print "Starting job $psicovopt.sh ..\n";
	system "./$psicovopt.sh > $psicovopt.log &";
	sleep 1;
}

####################################################################################################
chdir $ccmpreddir or confess $!;
system_cmd("cp $aln ./");
open  JOB, ">$id-ccmpred.sh" or confess "ERROR! Could not open $id-ccmpred.sh $!";
print JOB "#!/bin/bash\n";
print JOB "touch ccmpred.running\n";
print JOB "echo \"running ccmpred ..\"\n";
print JOB CCMPRED." -t ".NPROC." $id.aln $id.ccmpred > ccmpred.log\n";
print JOB "if [ -s \"$id.ccmpred\" ]; then\n";
print JOB "   mv ccmpred.running ccmpred.done\n";
print JOB "   echo \"ccmpred job done.\"\n";
print JOB "   exit\n";
print JOB "fi\n";
print JOB "echo \"ccmpred failed!\"\n";
print JOB "mv ccmpred.running ccmpred.failed\n";
close JOB;
system_cmd("chmod 755 $id-ccmpred.sh");
print "Starting job $id-ccmpred.sh ..\n";
system "./$id-ccmpred.sh > $id-ccmpred.log &";
sleep 1;

####################################################################################################
chdir $freecontactdir or confess $!;
system_cmd("cp $aln ./");
open  JOB, ">$id-freecontact.sh" or confess "ERROR! Could not open $id-freecontact.sh $!";
print JOB "#!/bin/bash\n";
print JOB "touch freecontact.running\n";
print JOB "echo \"running freecontact ..\"\n";
print JOB "".FREECONTACT." < $id.aln > $id.freecontact.rr\n";
print JOB "if [ -s \"$id.freecontact.rr\" ]; then\n";
print JOB "   mv freecontact.running freecontact.done\n";
print JOB "   echo \"freecontact job done.\"\n";
print JOB "   exit\n";
print JOB "fi\n";
print JOB "echo \"freecontact failed!\"\n";
print JOB "mv freecontact.running freecontact.failed\n";
close JOB;
system_cmd("chmod 755 $id-freecontact.sh");
print "Starting job $id-freecontact.sh ..\n";
system "./$id-freecontact.sh &";
sleep 1;

####################################################################################################
print("\nWait for max ".HOURLIMIT." hours until all jobs are done ..\n");
my $running = 1;
my $i = 0;
while(int($running) > 0){
	sleep (HOURLIMIT);
	my $psicov_running  = `find $psicovdir/ -name "*.running" | wc -l`;
	my $ccmpred_running = `find $ccmpreddir/ -name "*.running" | wc -l`;
	chomp $psicov_running;
	chomp $ccmpred_running;
	$running = int($psicov_running) + int($ccmpred_running);
	last if $running == 0;
	$i++;
	last if $i == 3600;
}

####################################################################################################
print "\nAttempting to kill psicov processes that aren't finished..\n";
system("ps aux | grep -e ".PSICOV." | grep -e $id | awk '{print \$2}' | grep -v grep | xargs kill -9");


####################################################################################################
print "\nChecking FreeContact prediction..\n";
if(not -f "$freecontactdir/$id.freecontact.rr"){
	confess "Looks like CCMpred did not finish! $freecontactdir/$id.freecontact.rr is absent!\n";
	system_cmd("touch $freecontactdir/$id.freecontact.rr");
}

####################################################################################################
print "\nChecking CCMpred prediction..\n";
if(not -f "$ccmpreddir/$id.ccmpred"){
	confess "Looks like CCMpred did not finish! $ccmpreddir/$id.ccmpred is absent!\n";
	system_cmd("touch $ccmpreddir/$id.ccmpred");
}

####################################################################################################
print "\nChecking PSICOV predictions..\n";
if (-s "$psicovdir/d0.03/$id-d0.03.psicov"){
	print "Looks like PSICOV 'd-0.03' option has already finished!\n";
	system_cmd("cp $psicovdir/d0.03/$id-d0.03.psicov $psicovdir/$id.psicov.rr");
}
elsif (-s "$psicovdir/r0.001/$id-r0.001.psicov"){
	print "Looks like PSICOV 'r-0.001' option has finished!\n";
	system_cmd("cp $psicovdir/r0.001/$id-r0.001.psicov $psicovdir/$id.psicov.rr");
}
elsif (-s "$psicovdir/r0.01/$id-r0.01.psicov"){
	print "Looks like PSICOV 'r-0.01' option has finished!\n";
	system_cmd("cp $psicovdir/r0.01/$id-r0.01.psicov $psicovdir/$id.psicov.rr");
}
else{
	confess "Looks like none of the PSICOV jobs finished!\n";
	system_cmd("touch $psicovdir/$id.psicov.rr");
}

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
