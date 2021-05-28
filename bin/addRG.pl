#!/usr/bin/perl
#
# Michael VanInsberghe 20210208
# adds RG header and tags to bam file
# Inputs: --bam      bam file
#

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $bam;
my $prefix;
my $out;
my $pseudobulk = '';

GetOptions("bam=s"=>\$bam,
           "out=s"=>\$out,
           "pseudobulk"=>\$pseudobulk,
           "prefix=s"=>\$prefix) or die("Error in command line arguments: $!\n");

die("No bam file supplied") unless $bam;

my $expt = $prefix || (split "_", $bam)[0];

(my $outfn = $bam) =~ s{\.[^.]*$}{.rg.bam};
$outfn = $out || $outfn;

open my $inbam, "samtools view -h $bam |" or die "Cannot open samtools pipe $bam: $!\n";
open my $outbam, "| samtools view -S -b - > $outfn" or die "Cannot open outfile samtools view -S b - > $outfn: $!\n";

my $header;
my $sam;
my %readgroups;

while (<$inbam>){
  chomp;
  next if /^\s*$/; # next if whitespace
  if (/^@/){
    $header .= "$_\n";
    next;
  }
  my @fields = split "\t";
  my @readfields = split ":", $fields[0];
  # readfields:
  # 0 instrument
  # 1 run number
  # 2 flow cell ID
  # 3 lane
  # 4 tile
  # 5 cluster x location
  # 6 cluster y location

  $fields[0] = join(":",@readfields[0..6]);

  # extract corrected cell barcode sequence (CB), use raw CR if unavailable
  my @CB = grep { /^CB:Z:/ } @fields;
  my $cb = "";
  if(scalar @CB > 0){
    ($cb = $CB[0]) =~ s/^CB:Z://;
  }else{
    @CB = grep { /^CR:Z:/ } @fields;
    if(scalar @CB > 0){
      ($cb = $CB[0]) =~ s/^CR:Z://;
    }
  }

  # extract corrected UMI sequence (UB), use raw UR if unavailable
  my @RX = grep { /^UB:Z:/ } @fields;
  my $rx = "";
  if(scalar @RX > 0){
    ($rx = $RX[0]) =~ s/^UB:Z://;
  }else{
    @RX = grep { /^UR:Z:/ } @fields;
    if(scalar @RX > 0){
      ($rx = $RX[0]) =~ s/^UR:Z://;
    }
  }

  my $rgid;
  my $lbid;
  my $smid;
  if($pseudobulk){
    $rgid = $readfields[2] . "." . $readfields[3] . "." . $expt;
    $lbid = $expt;
    $smid = $expt;
  }else{
    $rgid = $readfields[2] . "." . $readfields[3] . "." . $expt . "." . $cb;
    $lbid = $expt . "." . $cb;
    $smid = $expt . ".". $cb;
  }

  if (!exists($readgroups{$rgid})){
    $readgroups{$rgid} = "\@RG\tID:$rgid\tSM:$smid\tPU:$readfields[2].$readfields[3]\tPL:ILLUMINA\tPM:$readfields[0]\tLB:$lbid";
  }

  push @fields , ("RX:Z:$rx", "RG:Z:$rgid", "LB:Z:$lbid");

  $sam .= join("\t",@fields) . "\n";

}

close $inbam;

for my $key (keys %readgroups){
	$header .= "$readgroups{$key}\n";
}


print $outbam $header;
print $outbam $sam;
close $outbam;
