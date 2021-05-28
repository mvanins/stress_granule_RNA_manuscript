#!/usr/bin/perl
#
# Michael VanInsberghe 20210207
# removes reads from paired fastqs that match an input hit list
# outputs gzipped fastq with sequences removed
# Inputs: --hits     (gzipped) file containing hit list
#         --reads    (gzipped) fastq containing input reads
#         --prefix   output file prefix
#

use strict;
use warnings;
use Getopt::Long;
use File::Basename;


my @reads;
my $hits;
my $prefix;

GetOptions("reads=s{1,}"=>\@reads,
           "hits=s"=>\$hits,
           "prefix=s"=>\$prefix) or die("Error in command line arguments: $!\n");

die("No fastq reads files supplied") unless @reads;
die("No hits file supplied") unless $hits;


# open read and output read filehandles
my %readhandles;
my %outhandles;

foreach my $read (@reads){
  if($read =~/\.gz$/){
    open $readhandles{$read}, "zcat $read |" or die "Can't open gzipped file $read: $!\n";

    (my $outfn = $read) =~ s{\.[^.]*(?:\.gz)?$}{.depleted.fastq.gz};
    open $outhandles{$read}, "| gzip -c > $outfn" or die "Can't open gzipped output file $outfn: $!\n";
  }else{
    open $readhandles{$read}, "$read" or die "Can't open file $read: $!\n";

    (my $outfn = $read) =~ s{\.[^.]*(?:\.gz)?$}{.depleted.fastq};
    open $outhandles{$read}, "> $outfn" or die "Can't open output file $outfn: $!\n";
  }
}


# read hits
my $hitshandle;
if($hits =~ /\.gz$/){
  open $hitshandle, "zcat $hits |" or die "Can't open gzipped hits file $hits: $!\n";
}else{
  open $hitshandle, "$hits" or die "Can't open hits file $hits: $!\n";
}

my %hitlist;
while(<$hitshandle>){
  chomp;
  next if /^\s*$/; # next if whitespace
  $hitlist{$_} = 1;
}
close $hitshandle;


my %headers;
my %seqs;
my %quals;
my %pluss;
my $endOfReads = 1;

while($endOfReads){
  my @readnames;
  foreach my $read (@reads){
    chomp($headers{$read} = readline($readhandles{$read}));
    chomp($seqs{$read} = readline($readhandles{$read}));
    chomp($pluss{$read} = readline($readhandles{$read}));
    chomp($quals{$read} = readline($readhandles{$read}));

    my $name = (split /\s/, $headers{$read})[0];
    $name =~ s/^@//;
    push @readnames, $name;
  }
  my @unique = do { my %seen; grep { !$seen{$_}++ } @readnames };
  die("Read names not matched!\n") unless scalar @unique < 2;

  if(!exists $hitlist{$unique[0]}){
    foreach my $read (@reads){
      #print $outhandles{$read} "$headers{$read}\n$seqs{$read}\n$pluss{$read}\n$quals{$read}\n";
      print { $outhandles{$read} } "$headers{$read}\n$seqs{$read}\n$pluss{$read}\n$quals{$read}\n";
    }
  }

  foreach my $k (@reads){
    $endOfReads = $endOfReads & !eof($readhandles{$k});
  }
}

close $_ for values %readhandles;
close $_ for values %outhandles;
