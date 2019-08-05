#!/usr/bin/perl -w
# Badri Adhikari, 9-9-2017
# Tianqi Wu modified alignments pipeline, 7-25-2019

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use LWP::UserAgent;
use Time::Piece;
use Getopt::Long;

# User inputs
my ($aln_method);

GetOptions(
	"aln=s"	=> \$aln_method)
or confess "ERROR! Error in command line arguments!";

# Defaults
$aln_method  = "NOVEL"  if !$aln_method;
####################################################################################################
my $fasta  = shift;
my $outdir = shift;

if (not $fasta or not -f $fasta){
	print "Fasta file $fasta does not exist!\n" if ($fasta and not -f $fasta);
	print "Usage: $0 -aln <NOVEL,CONSTRUCT,CLUSTER> <fasta> <output-directory>\n";
	print "Usage: $0 -aln NOVEL 3e7u.fasta output/3e7u\n";
	exit(1);
}

if (not $outdir){
	print 'Output directory not defined!';
	print "Usage: $0 -aln <NOVEL,CONSTRUCT,CLUSTER> <fasta> <output-directory>\n";
	print "Usage: $0 -aln NOVEL 3e7u.fasta output/3e7u\n";
	exit(1);
}

####################################################################################################
use constant{
	SCRATCH      => '/home/badri/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh',
	BLASTPATH    => '/home/badri/ncbi-blast-2.2.25+/bin', 
	BLASTNRDB    => '/home/badri/databases/nr90-2012',
	PSIPRED      => '/home/badri/metapsicov/runpsipredandsolv',
	ALNSTAT      => '/home/badri/metapsicov/bin/alnstats',
	FEAT2RR      => abs_path(dirname($0)).'/predict-rr-from-features.sh',
	ALNSCRIPT    => abs_path(dirname($0)).'/generate-alignments.pl',
	ALNSCRIPT_NOVEL    => abs_path(dirname($0)).'/generate-alignments_novel.pl',
	ALNSCRIPT_CONSTRUCT    => abs_path(dirname($0)).'/generate-alignments_construct.pl',
	ALNSCRIPT_CLUSTER    => abs_path(dirname($0)).'/generate-alignments_cluster.pl',
	RRPREDSCRIPT => abs_path(dirname($0)).'/run-ccmpred-freecontact-psicov.pl',
	GENFEAT      => abs_path(dirname($0)).'/generate-dncon2-features.pl',
	HOURLIMIT    => 1
};

# confess "Oops!! alignment-script not found at ".ALNSCRIPT   if not -f ALNSCRIPT;
# confess "Oops!! rr-pred-script not found at ".RRPREDSCRIPT  if not -f RRPREDSCRIPT;
# confess "Oops!! psipred program not found at ".PSIPRED      if not -f PSIPRED;
# confess "Oops!! alnstat program not found at ".ALNSTAT      if not -f ALNSTAT;
# confess "Oops!! scratch program not found at ".SCRATCH      if not -f SCRATCH;
# confess "Oops!! GENFEAT program not found at ".GENFEAT      if not -f GENFEAT;
# confess "Oops!! Blast Path not found at ".BLASTPATH         if not -d BLASTPATH;
# confess "Oops!! BlastNRDB Path not found at ".BLASTNRDB     if not -d BLASTNRDB;

####################################################################################################
print "Started [$0]: ".(localtime)."\n";

print "\n";
print "Input: $fasta\n";
print "L    : ".length(seq_fasta($fasta))."\n";
print "Seq  : ".seq_fasta($fasta)."\n\n";

my $id = basename($fasta, ".fasta");
system_cmd("mkdir -p $outdir") if not -d $outdir;
system_cmd("cp $fasta $outdir/");
system_cmd("echo \">$id\" > $outdir/$id.fasta");
system_cmd("echo \"".seq_fasta($fasta)."\" >> $outdir/$id.fasta");
$outdir = abs_path($outdir);
chdir $outdir or confess $!;
$fasta = basename($fasta);

####################################################################################################
print "\n\n";
print "Generating PSSM..\n";
system_cmd("mkdir -p $outdir/pssm");
chdir $outdir."/pssm" or confess $!;
system_cmd("cp ../$fasta ./");
if ( -s "$id.pssm"){
	print "Looks like .pssm file is already here.. skipping..\n";
}
else{
	generate_pssm($fasta, "$id.pssm", "$id.psiblast.output");
	if( ! -e "$id.pssm") {
		print "Running less stringent version of PSSM generation..\n";
		generate_pssm_less_stringent($fasta, "$id.pssm", "$id.psiblast.output");
	}
}
chdir $outdir or confess $!;

####################################################################################################
print "\n\n";
print "Predicting secondary structure and solvent accessibility using SCRATCH..\n";
system_cmd("mkdir -p $outdir/ss_sa");
chdir $outdir."/ss_sa" or confess $!;
system_cmd("cp ../$fasta ./");
if (-s "$id.ss"){
	print "Looks like .aln file is already here.. skipping..\n";
}
else{
	system_cmd(SCRATCH." $fasta $id 4");
}
chdir $outdir or confess $!;
confess "ERROR!! No ss file!! Scratch failed!\n" if not -f "$outdir/ss_sa/$id.ss";
confess "ERROR!! No acc file!! Scratch failed!\n" if not -f "$outdir/ss_sa/$id.acc";
system_cmd("cp $outdir/ss_sa/$fasta $outdir/ss_sa/$id.ss_sa");
system_cmd("echo \">$id\" > $outdir/ss_sa/$id.ss_sa");
system_cmd("echo \"".seq_fasta($fasta)."\" >> $outdir/ss_sa/$id.ss_sa");
system_cmd("tail -n 1 $outdir/ss_sa/$id.ss >> $outdir/ss_sa/$id.ss_sa");
system_cmd("tail -n 1 $outdir/ss_sa/$id.acc >> $outdir/ss_sa/$id.ss_sa");
system_cmd("sed -i 's/-/b/g' $outdir/ss_sa/$id.ss_sa");
print "Predicted SS and SA:\n";
system "cat $outdir/ss_sa/$id.ss_sa";

####################################################################################################
print "\n\n";
print "Predicting secondary structure and solvent accessibility using PSIPRED..\n";
system_cmd("mkdir -p $outdir/psipred");
chdir $outdir."/psipred" or confess $!;
system_cmd("echo \">$id\" > ./$id.fasta");
system_cmd("echo \"".seq_fasta("$outdir/$fasta")."\" >> ./$id.fasta");
if (-s "$id.solv"){
	print "Looks like .solv file is already here.. skipping..\n";
}
else{
	system_cmd(PSIPRED." $fasta");
}
chdir $outdir or confess $!;

####################################################################################################
print "\n\n";
print "Generating alignments..\n";
$aln_method = uc($aln_method);
system_cmd("mkdir -p alignments");
if (-s "alignments/$id.aln"){
	print "Looks like .aln file is already here.. skipping..\n";
}
else{
	if ($aln_method eq "NOVEL") {
		system_cmd(ALNSCRIPT_NOVEL." $id.fasta alignments");
	}elsif ($aln_method eq "CONSTRUCT") {
		system_cmd(ALNSCRIPT_NOVEL." $id.fasta alignments");
	}elsif ($aln_method eq "CLUSTER") {
		system_cmd(ALNSCRIPT_NOVEL." $id.fasta alignments");
	}else{
		system_cmd(ALNSCRIPT." $id.fasta alignments");
	}
}
if (not -s "alignments/$id.aln"){
	print "Warning! Something went wrong! Alignments were not generated!\n";
	system_cmd("touch alignments/$id.aln");
	system_cmd("echo \"".seq_fasta($fasta)."\" > alignments/$id.aln");
}

####################################################################################################
print "\n\n";
print "Generate alignment stats ..\n";
system_cmd("mkdir -p alnstat");
system_cmd("cp alignments/$id.aln ./alnstat/");
system_cmd("cp alignments/$id.aln ./alnstat/");
chdir $outdir."/alnstat" or confess $!;
if (-s "$id.colstats"){
	print "Looks like .colstats file is already here.. skipping..\n";
}
else{
	system_cmd(ALNSTAT." $id.aln $id.colstats $id.pairstats");
}
chdir $outdir or confess $!;

####################################################################################################
print "\n\n";
print "Contact Predictions ..\n";
system_cmd("mkdir -p psicov");
system_cmd("mkdir -p ccmpred");
system_cmd("mkdir -p freecontact");
if (count_lines("alignments/$id.aln") <= 5){
	print "Too few sequences in the alignment! Skipping PSICOV, CCMpred and FreeContact run..\n";
}
elsif (-s "psicov/$id.psicov.rr"){
	print "Looks like contact predictions is already here.. skipping..\n";
}
else{
	print "\nRunning PSICOV, CCmpred, and FreeContact parallely..\n";
	system_cmd(RRPREDSCRIPT." alignments/$id.aln psicov ccmpred freecontact");
}

####################################################################################################
print "\n\n";
print "Verify coevolution-based contact predictions ..\n";
if (not -s "psicov/$id.psicov.rr"){
	warn "Warning! psicov/$id.psicov.rr file is empty!\n";
	system_cmd("touch psicov/$id.psicov.rr");
}
if (not -s "ccmpred/$id.ccmpred"){
	confess "Warning! ccmpred/$id.ccmpred file is empty!\n";
}
if (not -s "freecontact/$id.freecontact.rr"){
	confess "Warning! freecontact/$id.freecontact.rr file is empty!\n";
}

####################################################################################################
print "\n\n";
print "Generating feature file..\n";
system "pwd";
if (-s "feat-$id.txt"){
	print "Looks like feature file is already here.. skipping..\n";
}
else{
	system_cmd(GENFEAT. " ./$fasta ./pssm/$id.pssm ./ss_sa/$id.ss_sa ./alnstat/$id.colstats ./alnstat/$id.pairstats ./freecontact/$id.freecontact.rr ./ccmpred/$id.ccmpred ./psicov/$id.psicov.rr ./psipred/$id.ss2 ./psipred/$id.solv > feat-$id.txt");
}
print "Strip leading spaces from feature files..\n";
system_cmd("sed -i 's/^ *//g' feat-$id.txt");

####################################################################################################
print "\n\n";
print "Predict RR from features\n";
if (not -s "$id.rr.raw"){
	system_cmd(FEAT2RR. " feat-$id.txt $id.rr.raw $id.feat.stage2.txt")
}

####################################################################################################
print "\n\n";
print "Preparing predictions..\n";
system_cmd("mkdir -p predictions");
system_cmd("rm -f predictions/*");
freecontact2rr("$outdir/$fasta", "./freecontact/$id.freecontact.rr", "./predictions/$id.freecontact.rr", 6);
ccmpred2rr("$outdir/$fasta", "./ccmpred/$id.ccmpred", "./predictions/$id.ccmpred.rr", 6);
system_cmd("echo \"".seq_fasta($fasta)."\" > ./predictions/$id.psicov.rr");
system_cmd("cat ./psicov/$id.psicov.rr >> ./predictions/$id.psicov.rr");
system_cmd("echo \"".seq_fasta($fasta)."\" > ./predictions/dncon2.rr");
system_cmd("cat ./$id.rr.raw >> ./predictions/dncon2.rr");

####################################################################################################
print "Add some details to the prediction and prepared CASP RR format..\n";
my $aln  = "N/A";
my $N    = "N/A";
my $Neff = "N/A";
my $aln_line = `cat ./alignments/result.txt`;
chomp $aln_line;
my @C = split /\s+/, $aln_line;
$aln = $C[1] if $C[1];
my @colstats = `head -3 ./alnstat/*.colstats`;
chomp (@colstats);
$N = $colstats[1] if $colstats[1];
$Neff = $colstats[2] if $colstats[2];
add_casp_header_footer("./predictions/dncon2.rr", $id, $N, $Neff, $aln, "./$id.dncon2.rr");
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
		print "[[Executing: $command]]\n";
		system($command);
	}
	if($? != 0){
		my $exit_code  = $? >> 8;
		confess "ERROR!! Could not execute [$command]! \nError message: [$!]";
	}
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
sub freecontact2rr{
	my $file_fasta = shift;
	my $file_freecontact = shift;
	my $file_rr = shift;
	my $seq_sep_th = shift;
	confess "ERROR! file_fasta not defined!" if !$file_fasta;
	confess "ERROR! file_freecontact not defined!" if !$file_freecontact;
	confess "ERROR! file_rr not defined!" if !$file_rr;
	confess "ERROR! seq_sep_th not defined!" if !$seq_sep_th;
	my %conf = ();
	open CCM, $file_freecontact or confess $!;
	while(<CCM>){
		my @C = split /\s+/, $_;
		$conf{$C[0] . " " . $C[2]} = $C[5];
	}
	close CCM;
	open RR, ">$file_rr" or confess $!;
	print RR "".seq_fasta($file_fasta)."\n";
	foreach (sort {$conf{$b} <=> $conf{$a}} keys %conf){
		my @C = split /\s+/, $_;
		next if abs($C[0] - $C[1]) < $seq_sep_th;
		print RR $_." 0 8 ".$conf{$_}."\n";
	}
	close RR;
}

####################################################################################################
sub ccmpred2rr{
	my $file_fasta = shift;
	my $file_ccmpred = shift;
	my $file_rr = shift;
	my $seq_sep_th = shift;
	confess "ERROR! file_fasta not defined!" if !$file_fasta;
	confess "ERROR! file_ccmpred not defined!" if !$file_ccmpred;
	confess "ERROR! file_rr not defined!" if !$file_rr;
	confess "ERROR! seq_sep_th not defined!" if !$seq_sep_th;
	my %conf = ();
	open CCM, $file_ccmpred or confess $!;
	my $i = 1;
	while(<CCM>){
		my @C = split /\s+/, $_;
		for(my $j = 0; $j <= $#C; $j++){
			my $pair = $i." ".($j+1);
			$pair = ($j+1)." ".$i if ($j+1) < $i;
			my $confidence = $C[$j];
			$confidence = $conf{$pair} if (defined $conf{$pair} && $conf{$pair} > $confidence);
			$conf{$pair} = $confidence;
		}
		$i++;
	}
	close CCM;
	open RR, ">$file_rr" or confess $!;
	print RR "".seq_fasta($file_fasta)."\n";
	foreach (sort {$conf{$b} <=> $conf{$a}} keys %conf){
		my @C = split /\s+/, $_;
		next if abs($C[0] - $C[1]) < $seq_sep_th;
		print RR $_." 0 8 ".$conf{$_}."\n";
	}
	close RR;
}

####################################################################################################
sub add_casp_header_footer{
	my $rr     = shift;
	my $target = shift;
	my $N      = shift;
	my $Neff   = shift;
	my $aln    = shift;
	my $out_rr = shift;

	confess "ERROR!" if !$out_rr;
	confess "ERROR!" if !$rr;
	confess "$rr does not exit!"   if not -f $rr;
	confess "target-id not supplied!" if !$target;

	$N    = "??" if ! $N;
	$Neff = "??" if ! $Neff;
	$aln  = "??" if ! $aln;

	open ORR, ">$out_rr" or confess $!;
	print ORR "PFRMAT RR\n";
	print ORR "TARGET $target\n";
	print ORR "AUTHOR DNCON2\n";
	print ORR "METHOD DNCON2\n";
	print ORR "REMARK Number of sequences in the alignment = $N\n";
	print ORR "REMARK Effective number of sequences in the alignment = $Neff\n";
	print ORR "REMARK Alignment generated using $aln\n";
	print ORR "MODEL  1\n";
	print ORR "".wrap_seq(seq_rr($rr))."\n";
	open RR, $rr or confess "ERROR! Could not open $rr $!";
	while(<RR>){
		next if $_ !~ m/^[0-9]/;
		print ORR $_;
	}
	close RR;
	print ORR "END\n";
	close ORR;
}

####################################################################################################
sub seq_rr{
	my $file_rr = shift;
	confess ":(" if not -f $file_rr;
	my $seq;
	open RR, $file_rr or confess "ERROR! Could not open $file_rr! $!";
	while(<RR>){
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$_ =~ s/^\s+//;
		next if ($_ =~ /^>/);
		next if ($_ =~ /^PFRMAT/);
		next if ($_ =~ /^TARGET/);
		next if ($_ =~ /^AUTHOR/);
		next if ($_ =~ /^SCORE/); 
		next if ($_ =~ /^REMARK/);
		next if ($_ =~ /^METHOD/);
		next if ($_ =~ /^MODEL/); 
		next if ($_ =~ /^PARENT/);
		last if ($_ =~ /^TER/);   
		last if ($_ =~ /^END/);
		# Now, I can directly merge to RR files with sequences on top
		last if ($_ =~ /^[0-9]/);
		$seq .= $_;
	}
	close RR;
	$seq =~ s/\s+//;
	confess ":( no sequence header in $file_rr" if not defined $seq;
	return $seq;
}

####################################################################################################
sub wrap_seq{
	my $seq = shift;
	confess ":(" if !$seq;
	my $seq_new = "";
	while($seq){
		if(length($seq) <= 50){
			$seq_new .= $seq;
			last;
		}
		$seq_new .= substr($seq, 0, 50)."\n";
		$seq = substr($seq, 50);
	}
	return $seq_new;
}

####################################################################################################
sub generate_pssm{
	my $FASTA = shift;
	my $PSSM = shift;
	my $REPORT = shift;
	my $DB = "nr90";
	print "Running PSI-Blast with $DB...\n";
	system_cmd(BLASTPATH. "/psiblast -query ".$FASTA." -evalue .001 -inclusion_ethresh .002 -db ".BLASTNRDB. "/$DB -num_iterations 3 -outfmt 0 -out $REPORT -seg yes -out_ascii_pssm $PSSM");
}

####################################################################################################
sub generate_pssm_less_stringent{
	my $FASTA = shift;
	my $PSSM = shift;
	my $REPORT = shift;
	my $DB = "nr90";
	print "Running PSI-Blast with $DB...\n";
	system_cmd(BLASTPATH. "/psiblast -query ".$FASTA." -evalue 10 -inclusion_ethresh 10 -db ".BLASTNRDB. "/$DB -num_iterations 3 -outfmt 0 -out $REPORT -num_alignments 2000 -out_ascii_pssm $PSSM");
}
