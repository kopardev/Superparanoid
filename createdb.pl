use strict;
use Math::Combinatorics;
use List::Compare;
use Data::Dumper;
use DBI;
use Bio::SeqIO;
use lib qw(/usr/global/blp/perllib);
use Vutil;

# Opem sqlite db

my $dbh = DBI->connect(          
    "dbi:SQLite:dbname=groups.db", 
    "",                          
    "",                          
    { RaiseError => 1 },         
) or die $DBI::errstr;

my @elements;
foreach my $element (@ARGV) {
  push @elements,$element;
}
my $groupsFile=shift @elements;
my $nspecies=(scalar @elements);

# Check that all fastas exists

for my $fasta (@elements) {
  Vutil::fileCheck($fasta);
}

#initialize geneID2groupID
my %sname2ngenes;
my %geneID2groupID;
for my $fasta (@elements) {
  my $f=Bio::SeqIO->new(-file=>"$fasta",-format=>'fasta');
  my $ngenes=0;
  while (my $s=$f->next_seq) {
    $ngenes+=1;
    my $geneid=sprintf"%s##%s",$fasta,$s->id;
    $geneID2groupID{$geneid}=0;
  }
  $sname2ngenes{$fasta}=$ngenes;
}

my %element2speciesnumber;
for (my $i=0;$i<$nspecies;$i++) {
  my $j=$i+1;
  $element2speciesnumber{$elements[$i]}="Species_".$j;
}

# CREATE THE SPECIES TABLE

my %speciesName2ID;
$dbh->do("DROP TABLE IF EXISTS Species");
$dbh->do("CREATE TABLE Species(sid INT PRIMARY KEY,sname TEXT)");
for (my $i=0;$i<$nspecies;$i++) {
  $dbh->do("INSERT INTO Species VALUES($i+1,\'$elements[$i]\')");
  $speciesName2ID{$elements[$i]}=$i+1;
}

# CREATE THE COMBINATIONS.... NO TABLE YET... IDEA IS TO ADD ONLY COMBINATIONS WITH NON-ZERO NUMBER OF CLUSTERS TO
# COMBINATIONS TABLE

my %comboHash;
my %comboComplementHash;
my $ncombinations=0;
for (my $i=1;$i<=$nspecies;$i++) {
my $combinat = Math::Combinatorics->new(count => $i,data => [@elements] );
  while(my @combo = $combinat->next_combination){
    $ncombinations+=1;
    my $lc = List::Compare->new(\@combo, \@elements);
    my @combo_complement=$lc->get_Ronly();
    my $comboName=join("__",@combo);
    $comboHash{$comboName}=();
    foreach my $c (@combo) {
      push @{$comboHash{$comboName}},$element2speciesnumber{$c};
    }
    $comboComplementHash{$comboName}=();
    foreach my $c (@combo_complement) {
      push @{$comboComplementHash{$comboName}},$element2speciesnumber{$c};
    }
  }
}


# READ THE GROUPS.TXT FILE INTO AN ARRAY

open G,"<$groupsFile";
my @groupslines;
while (my $line=<G>) {
  chomp $line;
  push @groupslines,$line;
}
close G;

for my $line (@groupslines) {
  my @cols=split/\t/,$line;
  for (my $i=1;$i<(scalar @cols);$i++) {
    $geneID2groupID{$cols[$i]}=$cols[0];
  }
}


# CREATE THE SINGLETONS TABLE

$dbh->do("DROP TABLE IF EXISTS Singletons");
my $sql_cmd="CREATE TABLE Singletons(singleton_id INT PRIMARY KEY";
$sql_cmd .= ", sid INT, geneid TEXT)";
$dbh->do($sql_cmd);
my $singleton_id=0;
for my $geneid (keys %geneID2groupID) {
  next unless $geneID2groupID{$geneid}==0;
  $singleton_id+=1;
  my @tmp=split/##/,$geneid;
  my $geneID=sprintf "'%s'",$tmp[1];
  my @tmp2=split/_/,$element2speciesnumber{$tmp[0]};
  my $sid=$tmp2[1];
  my @values;
  push @values,$singleton_id;
  push @values,$sid;
  push @values,$geneID;
  my $sql_cmd="INSERT INTO Singletons VALUES(".join(",",@values).")";
  $dbh->do($sql_cmd);
}

open S, ">Singletons.txt";
print S "#Species\tNo_of_singletons\tTotal_no_of_genes\n";
for my $sname (keys %element2speciesnumber){
  my @tmp=split/_/,$element2speciesnumber{$sname};
  my $sid=$tmp[1];
  my $sql_cmd="SELECT COUNT(*) FROM Singletons WHERE sid = $sid";
  my $sth = $dbh->prepare( $sql_cmd );
  my $rv = $sth->execute() or die $DBI::errstr;
  if($rv < 0){
     print $DBI::errstr;
  }
  my @row = $sth->fetchrow_array();
  my $nsingletons=$row[0];
  chomp $nsingletons;
  print S "$sname\t$nsingletons\t$sname2ngenes{$sname}\n";
}
close S;


# CREATE THE CLUSTERS TABLE

$dbh->do("DROP TABLE IF EXISTS Clusters");
my $sql_cmd="CREATE TABLE Clusters(cid INT PRIMARY KEY";
for (my $i=0;$i<$nspecies;$i++) {
  my $y=$i+1;
  $sql_cmd.=",Species_$y INT";
}
for (my $i=0;$i<$nspecies;$i++) {
  my $y=$i+1;
  $sql_cmd.=",Species_${y}_seqid TEXT";
}
$sql_cmd.=")";
$dbh->do($sql_cmd);
#print "$sql_cmd\n"; exit;


my $ngroups=(scalar @groupslines);
my $nexttarget=10;

my $ctr=0;
print "Adding clusters(or groups) to the database: ... \n";
for my $line (@groupslines) {
  $ctr+=1;
  my @values;
  my @values2;
  my @cols=split/\t/,$line;
  push @values,$cols[0];
  for my $x (@elements) {
    my @t=$line=~/(?<!\S)($x\S*)/g;
    #print "line=$line\n";
    #print "t=@t\n";
    push @values,(scalar @t);
    my @y;
    for my $g (@t) {
      #print "g=$g\n";
      my @tmp=split/##/,$g;
      push @y,$tmp[1];
    }
    my $z=sprintf("'%s'",join(":",@y));
    push @values2,$z;
  }
  $sql_cmd="INSERT INTO Clusters VALUES(".join(",",@values).",".join(",",@values2).")";
  #print $sql_cmd,"\n";
  $dbh->do($sql_cmd);
  #print "$cols[0]\n";
  #last if $ctr==10;
  my $percent_complete=$ctr*100.0/$ngroups;
  if ($percent_complete>=$nexttarget) {
    $nexttarget+=10;
    printf "%.0f%\n",$percent_complete;
  }
}

#$dbh->disconnect();
#exit;

# CREATE COMBINATIONS TABLE NOW

print "Trying all combinations ...\n";

$dbh->do("DROP TABLE IF EXISTS Combinations");
my $sql_cmd="CREATE TABLE Combinations(combination_id INT PRIMARY KEY";
$sql_cmd.=",combination_name TEXT,nclusters INT";
for (my $i=0;$i<$nspecies;$i++) {
  my $y=$i+1;
  $sql_cmd.=",Species_$y INT";
}
$sql_cmd.=")";
#print "$sql_cmd\n"; exit;
$dbh->do($sql_cmd);

my $combination_ctr=0;
$ctr=0;
$nexttarget=10;
my $totalClusters=$dbh->selectrow_array("SELECT COUNT(*) FROM Clusters");
my $totalClustersinCombinations=0;
for my $comboName (keys %comboHash) {
  $ctr+=1;
  my $percent_complete=100.0*$ctr/$ncombinations;
  if ($percent_complete>=$nexttarget) {
    $nexttarget+=10;
    printf "%.0f%\n",$percent_complete;
  }
  my @combo=@{$comboHash{$comboName}};
  my @comboComplement=();
  if ( defined $comboComplementHash{$comboName}) {
    @comboComplement=@{$comboComplementHash{$comboName}};
  }
  #print "NAME:$comboName\tCONTAINS:@combo\tDOES_NOT_CONTAIN:@comboComplement\n";
  #print "NAME:$comboName\n";
  my $sql_cmd="SELECT COUNT(*) FROM Clusters WHERE ";
  my $nwheres=0;
  for my $s (@combo) {
    $sql_cmd.=" $s > 0 " if $nwheres==0;
    $sql_cmd.=" AND $s > 0" if $nwheres!=0;
    $nwheres+=1;
  }
  for my $s (@comboComplement) {
    $sql_cmd.=" AND $s = 0"
  }
  my $count = $dbh->selectrow_array($sql_cmd); # COUNT NUMBER OF ROWS RETURNED BY SQL QUERY
  $count=0 if not defined $count;
  next if $count==0;
  $totalClustersinCombinations+=$count;
  #print "$count\t$comboName\t$sql_cmd\n";
  $combination_ctr+=1;
  my @values;
  for (my $i=0;$i<$nspecies;$i++) {
    push @values,0;
  }
  for my $species (@combo) {
    my @tmp=split/_/,$species;
    my $speciesNumber=int($tmp[1]);
    my $j=$speciesNumber-1;
    $values[$j]=1;
  }
  $sql_cmd="INSERT INTO Combinations VALUES($combination_ctr,\'$comboName\',$count,".join(",",@values).")";
  $dbh->do($sql_cmd);
  last if $totalClustersinCombinations == $totalClusters;
  #exit if $ctr==500;
}
#print "$ctr\t$combination_ctr\n"

# WRITE OUT THE COMBINATIONS WITH NONZERO NUMBER OF CLUSTERS

open C, ">Combinations.txt";
print C "#No_of_clusters\tCombination_id\tnSpecies(total=$nspecies)\tSpecies\n";
my $sql_cmd="SELECT * FROM Combinations ORDER BY nclusters DESC";
my $sth = $dbh->prepare( $sql_cmd );
my $rv = $sth->execute() or die $DBI::errstr;
if($rv < 0){
   print $DBI::errstr;
}
while(my @row = $sth->fetchrow_array()) {
  my $combinationNumber=$row[0];
  my @species=split/__/,$row[1];
  my $nclusters=$row[2];
  my $nspeciesinCombination=(scalar @species);
  print C "$nclusters\t$combinationNumber\t$nspeciesinCombination\t",join(",",@species),"\n";
}



$dbh->disconnect();
exit;