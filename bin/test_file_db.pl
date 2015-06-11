#!/software/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Storable;
use List::MoreUtils qw(uniq);
use List::Util qw(shuffle);
use Switch;
use Carp;


local $| = 1;       #setautoflush, outputs immediately instead of buffering

main();

sub main {

  my $file11 = '../dna_language_db/ngram11';

  my $array_ref = retrieve($file11);
  print Dumper($array_ref);

#   my $hashref = retrieve($file);
#   for my $b1(sort keys %{$hashref}) {
#     for my $ngram2(sort keys %{$hashref->{$b1}}) {
#       for my $ngram3(sort keys %{$hashref->{$b1}->{$ngram2}}) {
# 	print "$b1\t$ngram2\t$ngram3\n";
#       }
#     }
#   }
}
