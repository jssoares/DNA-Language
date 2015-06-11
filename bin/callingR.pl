#!/software/bin/perl-5.14.2

use v5.14.2;
use strict;
use warnings;
use Sanger::CGP::Config::Config qw(/software/CGP/projects/translationDB/perl/config/local_config.ini);
use Sanger::CGP::Database::Conn;
use DBI;
use Carp;
use Data::Dumper;
use Getopt::Long;
use Statistics::R;

#program version
my $VERSION="0.1";

local $| = 1;       #setautoflush, outputs immediately instead of buffering

eval {

  main();

};

warn() if $@;


sub main {

  my $files = getFiles();

  #print Dumper($files);
  #exit;

  #calculateRelativeFreqs($files);
  writeFilesAllK($files);
  exit;

  my @grams = qw(k11);


  for my $file(@$files) {
    for my $gram(@grams) {
      if ($file =~ m/real\/(k\d{0,2})\//) {
	my $k = $1;
	if ($k eq $gram) {
	  #print "$file\n";
	  createGraphs($k,$file);
	}
      }
    }
  }

  #print Dumper(\%hash);
}


sub getFiles {

  my @files;
  push(@files, `find ../dna_language_analysis/real/\* -iname '*.txt' -not -iname '*_all*'`);
  for my $file(@files) {
    $file =~ s/\n//;
  }
  return(\@files);

}


sub writeFilesAllK {

  my ($files) = @_;
  my ($gram,$org,$chr,$df_label);
  my %hash;
  for my $file(@$files) {
    if ($file =~ m/real\/(k\d{0,2})\/(\w*)\/(\w*)\//) {
      $gram = $1;
      $org = $2;
      $chr = $3;
      $df_label = $1 . '_' . $2 . '_' . $3;
      open(my $fh, '<', $file) or croak "can't open .$!";
      my @lines = <$fh>;
      close($fh);
      for (my $i = 0; $i < scalar @lines; $i++) {
	$lines[$i] =~ s/\n//;
	my ($k,$sum,$rel_freq,$label) = split(/\t/, $lines[$i]);
	if ( $i > 0) {
	  push( @{$hash{$gram}{k}}, $k);
	  push( @{$hash{$gram}{sum}}, $sum);
	  push( @{$hash{$gram}{rel_freq}}, $rel_freq);
	  push( @{$hash{$gram}{label}}, $label);
	}
      }
    }
  }
  for $gram(keys %hash) {
    my $file_path = "../dna_language_analysis/real/$gram/$gram" . '_all.txt';
    open(my $fh2,'>', $file_path);
    print $fh2 "$gram\tCount\tRelFreq\tORG_CHR\n";
    for (my $i = 0; $i < scalar @{$hash{$gram}{k}}; $i++) {
      print $fh2 "$hash{$gram}{k}[$i]\t$hash{$gram}{sum}[$i]\t$hash{$gram}{rel_freq}[$i]\t$hash{$gram}{label}[$i]\n";
    }
    close($fh2);
  }
  return;
}

sub calculateRelativeFreqs {

  my ($files) = @_;
  my ($gram,$org,$chr,$df_label);
  my %hash;
  for my $file(@$files) {
    if ($file =~ m/real\/(k\d{0,2})\/(\w*)\/(\w*)\//) {
      $gram = $1;
      $org = $2;
      $chr = $3;
      $df_label = $1 . '_' . $2 . '_' . $3;
      open(my $fh, '<', $file) or croak "can't open .$!";
      my @lines = <$fh>;
      close($fh);
      for (my $i = 0; $i < scalar @lines; $i++) {
	$lines[$i] =~ s/\n//;
	my ($k,$sum,$label) = split(/\t/, $lines[$i]);
	if ( $i > 0) {
	  push(@{ $hash{$gram}{k} }, $k);
	  push(@{ $hash{$gram}{sum} }, $sum);
	  push(@{ $hash{$gram}{label} }, $label);
	}
      }
    }

    my $total_grams_in_file;
    for $gram(keys %hash) {
      for (my $i =0; $i < scalar @{ $hash{$gram}{k} }; $i++) {
	$total_grams_in_file += $hash{$gram}{sum}[$i];
      }
    }
    for $gram(keys %hash) {
      open(my $fh,'>', $file);
      print $fh "$gram\tCount\tRelFreq\tORG_CHR\n";
      for (my $i = 0; $i < scalar @{$hash{$gram}{k}}; $i++) {
	print $fh ("$hash{$gram}{k}[$i]\t$hash{$gram}{sum}[$i]\t", $hash{$gram}{sum}[$i] / $total_grams_in_file, "\t$hash{$gram}{label}[$i]\n");
      }
      close($fh);
    }
  }
  return;
}


sub createGraphs {

  my ($k,$file) = @_;
  # Create a communication bridge with R and start R
  my $R = Statistics::R->new();
  $R->start();

  my $cmds1 = <<EOF;
library(ROracle)
library(proto)
library(grid)
library(reshape)
library(plyr)
library(ggplot2)
drv <- dbDriver("Oracle")
EOF


  my $cmds2 = <<EOF;

$k <- read.table("$file", header=TRUE)
EOF

  my $pdf_name = "/nfs/cancer_translation/Jorge/dna_language_proj/dna_language_analysis/graphs/$k" . "_countAbove500.pdf";

  my $cmds3 = <<EOF;
pdf("$pdf_name")
EOF

  my $r_for_loop = q{};
  $r_for_loop .= 'for (i in unique(' . $k . '$ORG_CHR)) {';
  $r_for_loop .= "if ( $k" . '$count >= 500) {';
  $r_for_loop .= " p <- ggplot($k" . '[' . $k . '$ORG_CHR == i, ], aes(x=';
  $r_for_loop .= $k . ', y=count, group=ORG_CHR));';
  $r_for_loop .= 'p <- p + geom_line(size = 0.5, colour="green3") + facet_wrap(~ORG_CHR) + opts(axis.text.x=theme_text(angle=90, hjust=1, size=3));print(p);};}; dev.off()';

  my $cmds4 = <<EOF;
$r_for_loop
EOF

  #print "$r_for_loop\n";

  my $out2 = $R->run($cmds1,$cmds2,$cmds3,$cmds4);
#my $a = $R->get( 'c' );
#print Dumper($out2);



}
