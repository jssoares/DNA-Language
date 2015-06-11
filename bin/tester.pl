#!/software/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Net::FTP;
use Carp;
use lib 'Utils';
use EnsemblFTP;

local $| = 1;       #setautoflush, outputs immediately instead of buffering

main();

sub main {

  my %opts;

  GetOptions( \%opts,
	      'help=s' => \my $help,
	      'h=s' => \my $host,    #Options: 'ensembl' or 'genomes'
	      'b=s' => \my $branch,     #Options: 'plants' or 'fungi' or 'protist' or 'metazoa' or 'bacteria'
	      'r=s' => \my $rmode,      #Options: 'ro', 'rw'
	      't=i' => \my $test,   #Options: 1 (run the test sub_routine) or 0 (don't run it)
	    ) or pod2usage(2);

  pod2usage(2) if ($help);

  my %args;
  $args{host} = $host; 
  $args{branch} = $branch;
  $args{rmode} = $rmode;
  $args{test} = $test;

  my $ens_ftp = EnsemblFTP->new(%args);
  $ens_ftp->getSpeciesEnsemblFTP;
  $ens_ftp->testFTP;
  print Dumper($ens_ftp);

}
