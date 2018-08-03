=head1 NAME

RosenLab - Base module for all Rosen Lab Perl-related stuff.

=head1 SYNOPSIS

	use RosenLab;

=head1 DESCRIPTION

The module serves mainly as a common root for various other RosenLab modules.

=cut

# module stuff
package RosenLab;
BEGIN
{
	our $VERSION = '0.001'; # use string for MakeMaker (and similar)
	$VERSION = eval "$VERSION"; # convert to number for `use` (and similar)
	use Exporter qw< import >;
	our @ISA = qw< Exporter >;
	our @EXPORT = qw<  >;
	our @EXPORT_OK = qw< projects >;
	our %EXPORT_TAGS =
		(
			'STD' => \@EXPORT,
			'ALL' => [@EXPORT,@EXPORT_OK],
		);
}

# pragmas
use 5.026;
use strict; use warnings;
use experimental qw< signatures postderef >;

# core modules
use File::Spec::Functions qw< catdir catfile catpath rootdir >;
use List::Util qw< reduce pairgrep pairmap >;

=head1 CONSTANTS

=over

=item ROOT

The absolute path to the Rosen Lab shared space. This location serves as the de
facto root for the vast majority of any task(s) performed within it, hence its
name.

=item RAW_DATA_REPOSITORY

The absolute path to the Rosen Lab raw data repository.

=item SOFTWARE_REPOSITORY

The absolute path to the Rosen Lab software repository.

=back

=cut

use constant
	{
		ROOT => catdir(rootdir(),'broad','rosenlab_archive'),
	};

use constant
	{
		RAW_DATA_REPOSITORY => catdir(ROOT,'Data'),
		PROJECT_REPOSITORY => catdir(ROOT,'Projects'),
		SOFTWARE_REPOSITORY => catdir(ROOT,'Software'),
	};

# data sub-repositories
use constant
	{
		GENOME_DATA_REPOSITORY => catdir(RAW_DATA_REPOSITORY,'Genomes'),
		MOTIF_DATA_REPOSITORY => catdir(RAW_DATA_REPOSITORY,'Motifs'),
		SEQUENCING_DATA_REPOSITORY => catdir(RAW_DATA_REPOSITORY,'Sequencing'),
	};

=head1 SUBROUTINES

=cut



=head1 SEE ALSO

L<RosenLab::Pipeline>

=head1 COPYRIGHT & LICENSE

=cut

1;

__END__
