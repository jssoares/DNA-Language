package DNA::Language::Utils::SimpleBaseCounter;

use Moose;
use Bio::Tools::SeqStats;

has 'seq' => ( is => 'rw', isa => 'Bio::Seq', required => 1 );
has 'monomers_freq' => ( is => 'rw', isa => 'HashRef', lazy => 1, builder => '_build_monomers_freq' );

sub _build_monomers_freq {

  my ($self) = @_;
  my $seq_stats = Bio::Tools::SeqStats->new(-seq => $self->seq);
  my $monomers = $seq_stats->count_monomers;
  for my $base(%{$monomers}) {
    if (defined $monomers->{$base}) {
      $monomers->{$base} =  $monomers->{$base}/$self->seq->length;
    }
  }
  return($monomers);
}

no Moose;
__PACKAGE__->meta->make_immutable;
1;
