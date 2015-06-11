package DNA::Language::Utils::FastaParser;

use Moose;
use Bio::SeqIO;

has 'fasta_file' => ( is => 'rw', isa => 'Str', required => 1 );
has 'bio_seq' => ( is => 'rw', isa => 'Bio::SeqIO', lazy => 1, builder => '_build_bio_seq' ); 

sub _build_bio_seq {

  my ($self) = @_;
  my $bio_seq = Bio::SeqIO->new(-file => $self->fasta_file, -format => "fasta" );
  return($bio_seq);
}

no Moose;
__PACKAGE__->meta->make_immutable;
1;

