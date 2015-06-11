#!/usr/bin/env perl
use Moose;
use Data::Dumper;

BEGIN { unshift( @INC, './lib' ) }
BEGIN { unshift( @INC, './t/lib' ) }
BEGIN {
  use Test::Most;
  use Test::Exception;
  use_ok('DNA::Language');
}

=head
  
ok ( my $lanalysis =  DNA::Language->new( fasta_file => 't/data/fasta_seqs/very_small.fa', mode => 'real', k => 4, rand_range => 4, tnb => 200, host => '', branch => '', rmode => '', test => 1 ), 'Creating language analysis object' );
$lanalysis->run;
is ( $lanalysis->fasta_file, 't/data/fasta_seqs/very_small.fa', 'Fasta file' );
is_deeply ( $lanalysis->genetic_code, ['a','t','c','g'], 'Genetic code array' );
isa_ok($lanalysis->bio_seq,'Bio::SeqIO');
while (my $seq_obj = $lanalysis->bio_seq->next_seq){
  is( $seq_obj->seq,'AAACTGTAAACTG', 'DNA sequence');
}
#is ( $lanalysis->basic_stats->{A_small_fasta_sequence_for_testing}->{A}, 0.461538461538462, 'Number of As' );
#is ( $lanalysis->basic_stats->{A_small_fasta_sequence_for_testing}->{T}, 0.230769230769231, 'Number of Ts' );
#is ( $lanalysis->basic_stats->{A_small_fasta_sequence_for_testing}->{C}, 0.153846153846154, 'Number of Cs' );
#is ( $lanalysis->basic_stats->{A_small_fasta_sequence_for_testing}->{G}, 0.153846153846154, 'Number of Gs' );



=cut


ok ( my $lanalysis =  DNA::Language->new( fasta_file => 't/data/fasta_seqs/slightly_bigger.fa', mode => 'test', k => 16, rand_range => 4, tnb => 200, host => '', branch => '', rmode => '', test => 1 ), 'Creating language analysis object' );
$lanalysis->run;
print Dumper($lanalysis->dna_cryptic->digit_to_dna);

=head

ok ( $lanalysis =  DNA::Language->new( fasta_file => 't/data/fasta_seqs/A_nuc.fasta', mode => 'test', k => 3, rand_range => 4, tnb => 200, host => '', branch => '', rmode => '', test => 1 ), 'Creating language analysis object' );
$lanalysis->run;
print Dumper($lanalysis);

=cut

done_testing();
