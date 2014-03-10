#! /usr/bin/perl

use strict;
use warnings;
use lib qw(/usr/global/blp/perllib);
use Vutil;
use Qsub;
use Cwd;
use Config::General;
use Getopt::Long;

sub usage;

my ($fastaList,$configFileName,$force,$help,$keeplog,$createdatabase);

my $result = GetOptions ( "f=s" => \$fastaList, # list of fasta files (one file name per line) with protein sequences
			  "force" => \$force, # overwrite files
			  "d" => \$createdatabase,
			  "c=s" => \$configFileName,
			  "keeplogs" => \$keeplog, # do not delete log files
			  "h|help" => \$help);

usage() if ($help);

my $path=Vutil::getFileAbsolutePath(__FILE__);
my $inparanoid="${path}inparanoid.pl";
my $multiparanoid="${path}multiparanoid.pl";
my $mp2omcl="${path}mp2omcl.pl";
my $createdb="${path}createdb.pl";
my $cwd=getcwd();

if (not defined $configFileName) {
    Vutil::fileCheck("${path}/config.txt","Default config file required if one is not provided!");
    $configFileName="${path}/config.txt";
}

Vutil::fileCheck($configFileName,"Check the config file!");

my $cfg=new Config::General($configFileName);
my %cfgHash=$cfg->getall;

Vutil::fileCheck($fastaList,"Check the config file!");
open F, "<$fastaList";
my @fastas;
while (my $line=<F>){
	chomp $line;
	Vutil::fileCheck("${cwd}/${line}","Check the $fastaList file!");
	push @fastas,$line;
}
my $nfastas = scalar @fastas;

# make blast database for each fasta file
print "Making BLAST databases...\n";
for my $f (@fastas) {
	system ("$cfgHash{formatdb} -i $f");
}

# Create self-blast files
print "Creating preliminary files...\n";
my $cmd;
my @jobs1;
my $job;
for my $i (0..$nfastas-1) {
	$cmd="perl $inparanoid 1 0 $fastas[$i]";
	$job=new Qsub(name=>"SP1_$i",outFile=>"$fastas[$i].tmp.out",wd=>$cwd,cmd=>$cmd,nproc=>6);
	$job->submit();
	push @jobs1,$job;
}

print "Waiting for jobs to finish\n";
for $job (@jobs1) {
	$job->waitForCompletion();
}

# InParanoid
print "Starting InParanoid...\n";
my @jobs2;
for my $i (0..$nfastas-2) {
    for my $j ($i+1..$nfastas-1) {
        $cmd = "perl $inparanoid 1 1 $fastas[$i] $fastas[$j]";
	$job=new Qsub(name=>"SP2_${i}_${j}",outFile=>"$fastas[$i]-$fastas[$j].tmp.out",wd=>$cwd,cmd=>$cmd,nproc=>6);
	$job->submit();
	push @jobs2,$job;
    }
}
# Wait for InParanoid jobs to run
print "Waiting for InParanoid jobs to complete...\n";
for $job (@jobs2) {
	$job->waitForCompletion();
}
undef(@jobs1);
undef(@jobs2);
print "InParanoid jobs are complete.\n";

# start Multiparanoid
print "Starting MultiParanoid...\n";
$cmd = "perl $multiparanoid -species $fastas[0]";
for my $k (1..$nfastas-1) {
   $cmd .= "+$fastas[$k]";
}
system($cmd);
print "Multiparanoid job is complete.\n";

# create groups.txt
print "Creating groups.txt like orthomcl\n";
my $sqlFile="$fastas[0]";
for my $k (1..$nfastas-1) {
   $sqlFile .= "-$fastas[$k]";
}
$sqlFile .= ".MP.sql";
Vutil::fileCheck($sqlFile,"Multiparanoid did not generated the required .sql file\n");
$cmd = "perl $mp2omcl $sqlFile|sort -k1n > groups.txt";
system($cmd);

#create sqlite database
if (defined $createdatabase) {
	print "Creating sqlite database ...\n";
	$cmd = "perl $createdb groups.txt ".join(" ",@fastas);
	$job=new Qsub(name=>"SP_db",outFile=>"db.tmp.out",wd=>$cwd,cmd=>$cmd,nproc=>1);
	$job->submit();
	print "Waiting for jobs to finish.\n";
	$job->waitForCompletion();
}
print "Done\n";

unless ($keeplog) {
	system("rm -f SP_*");
}

sub usage {
print <<EOF;

Superparanoid is a wrapper script for running inparanoid and multiparanoid in conjuction.

Author: Vishal N. Koparde, Ph. D.
Created: 140220
Modified: 140220

options:
-f file containing list protein-sequences fasta files (one file per line, all files or symbolic links should be in the current folder; fasta filename should atleast one character... all digits may crash)
-force overwrite existing files
-d convert output to sqlite3 database
-c config file to use (optional)
-k keep log files (optional ... logs files are deleted by default)
-h help


EOF
exit 1;
}
