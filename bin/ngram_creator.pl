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

#  my $file3 = '../dna_language_db/ngram3';
#   my $file4 = '../dna_language_db/ngram4';
#   my $file5 = '../dna_language_db/ngram5';
#   my $file6 = '../dna_language_db/ngram6';
#   my $file7 = '../dna_language_db/ngram7';
#   my $file8 = '../dna_language_db/ngram8';
#   my $file9 = '../dna_language_db/ngram9';
#   my $file10 = '../dna_language_db/ngram10';
#   my $file11 = '../dna_language_db/ngram11';
#   my $file12 = '../dna_language_db/ngram12';
#   my $file13 = '../dna_language_db/ngram13';
#   my $file14 = '../dna_language_db/ngram14';
   my $file15 = '../dna_language_db/ngram15';
#   my $file16 = '../dna_language_db/ngram16';

  my @ngrams3;
  my @ngrams4;
  my @ngrams5;
  my @ngrams6;
  my @ngrams7;
  my @ngrams8;
  my @ngrams9;
  my @ngrams10;
  my @ngrams11;
  my @ngrams12;
  my @ngrams13;
  my @ngrams14;
  my @ngrams15;
  my @ngrams16;

  my @bases = qw{a c g t};
  for my $b1(@bases) {
    for my $b2(@bases) {
      for my $b3(@bases) {
#	my $g3 = $b1 . $b2 . $b3;
#	push(@ngrams3,$g3);
   	for my $b4(@bases) {
#	  my $g4 = $b1 . $b2 . $b3 . $b4;
# 	  push(@ngrams4,$g4);
   	  for my $b5(@bases) {
# 	    my $g5 = $b1 . $b2 . $b3 . $b4 . $b5;
# 	    push(@ngrams5,$g5);
   	    for my $b6(@bases) {
# 	      my $g6 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6;
# 	      push(@ngrams6,$g6);
   	      for my $b7(@bases) {
# 		my $g7 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7;
# 		push(@ngrams7,$g7);
   		for my $b8(@bases) {
# 		  my $g8 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7 . $b8;
# 		  push(@ngrams8,$g8);
   		  for my $b9 (@bases) {
# 		    my $g9 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7 . $b8 . $b9;
# 		    push(@ngrams9,$g9);
   		    for my $b10 (@bases) {
# 		      my $g10 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7 . $b8 . $b9 . $b10;
# 		      push(@ngrams10,$g10);
 		      for my $b11 (@bases) {
# 			my $g11 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7 . $b8 . $b9 . $b10 . $b11;
# 			push(@ngrams11,$g11);
 			for my $b12 (@bases) {
# 			  my $g12 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7 . $b8 . $b9 . $b10 . $b11 . $b12;
# 			  push(@ngrams12,$g12);
 			  for my $b13 (@bases) {
# 			    my $g13 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7 . $b8 . $b9 . $b10 . $b11 . $b12 . $b13;
# 			    push(@ngrams13,$g13);
 			    for my $b14 (@bases) {
# 			      my $g14 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7 . $b8 . $b9 . $b10 . $b11 . $b12 . $b13 . $b14;
# 			      push(@ngrams14,$g14);
 			      for my $b15 (@bases) {
 				my $g15 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7 . $b8 . $b9 . $b10 . $b11 . $b12 . $b13 . $b14 . $b15;
 				push(@ngrams15,$g15);
# 				for my $b16 (@bases) {
# 				  my $g16 = $b1 . $b2 . $b3 . $b4 . $b5 . $b6 . $b7 . $b8 . $b9 . $b10 . $b11 . $b12 . $b13 . $b14 . $b15 . $b16;
# 				  push(@ngrams16,$g16);

# 				}
 			      }
 			    }
 			  }
 			}
 		      }
   		    }
   		  }
   		}
    	      }
    	    }
 	  }
   	}
      }
    }
  }
#  store \@ngrams3, $file3;
#  store \@ngrams4, $file4;
#   store \@ngrams5, $file5;
#   store \@ngrams6, $file6;
#   store \@ngrams7, $file7;
#   store \@ngrams8, $file8;
#   store \@ngrams9, $file9;
#   store \@ngrams10, $file10;
#   store \@ngrams11, $file11;
#   store \@ngrams12, $file12;
#   store \@ngrams13, $file13;
#   store \@ngrams14, $file14;
   store \@ngrams15, $file15;
#   store \@ngrams16, $file16;
}
