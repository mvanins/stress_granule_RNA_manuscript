#!/usr/bin/perl
# adds sequence around variant to CONTEXT vcf tag

use warnings;
use strict;
use Getopt::Std;
use POSIX;

use vars qw($opt_r $opt_v $opt_d $opt_o);
getopts("r:v:d:o:");

my $usage = "\nuse me like this:\n$0 -r reference.fa -v variants.vcf.gz -o output.vcf.gz (-d bases_on_either_side)\n";
die $usage unless ($opt_r && $opt_v && $opt_o);

my $distance = $opt_d || 1;

my %chrs;
my $chr;

open REF, $opt_r or die "Cannot open reference $opt_r: $!\n";
while (<REF>){
  chomp;
  if (/^\>/){
    s/^\>//;	# remove leading >
    my @inf = split ' ';
    $chr = $inf[0];
    next;
  }
  $chrs{$chr} .= $_;
}

close REF;

open VAR, "zcat $opt_v |" or die "Cannot open variant $opt_v: $!\n";

open OUT, "| bgzip -c > $opt_o" or die "Cannot open output variant bgzip pipe $opt_o: $!\n";

my $varheader;
my $headerind = 1;
while(<VAR>){
  chomp;
  if (/^##/){
    $varheader .= "$_\n";
    # print "$_\n";
    next;
  }
  if (/^#/){
    $varheader .= '##INFO=<ID=CONTEXT,Number=1,Type=String,Description="Sequence context of mutation, in forward coordinates. Mutation happens to the middle base.">';
    $varheader .= "\n##$0 -r $opt_r -v $opt_v -o $opt_o -d $distance ; Date=" . strftime("%Y-%m-%d %H:%M:%S", localtime);

    print OUT "$varheader\n$_";
    next;
  }

  my @vcfs = split "\t";
  my $chr = $vcfs[0];
  my $location = $vcfs[1];
  # my $annotation = $vcfs[7];

  my $sequence = substr $chrs{$chr}, $location - 1 - $distance, 2*$distance + 1;

  $vcfs[7] .= ";CONTEXT=$sequence";

  print OUT "\n" . join "\t" , @vcfs
}

close VAR;
close OUT;

my $tabixcmd = "tabix -p vcf $opt_o";
system($tabixcmd);
