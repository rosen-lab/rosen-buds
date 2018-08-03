
=head1 NAME

RosenLab::Log::Simple - Straightforward logging facilities for Rosen Lab jobs.

=head1 SYNOPSIS

	# using a subset of the logging facilities for reporting when things go bad
	use RosenLab::Log::Simple qw< fatal error warning >;
	warn('report a warning message if RosenLab::Log::Simple's warnings are on");
	warn('report a warning message if ${category}'s warnings are on",$category);
	error('report an error message');
	fatal('report an error message and die (as per Carp::croak)');
	
	# using a subset of the logging facilities for reporting only what's wanted
	use RosenLab::Log::Simple qw< error message debug >;
	local $RosenLab::Log::Simple::verbosity = 1; # default = 0
	error('error messages are always reported');
	message('this message will be reported when verbosity >= 1');
	debug('this message will _not_ be reported when verbosity < 2');
	
	# use the Carp-like logging functions
	use RosenLab::Log::Simple qw< :CARP >;
	carp('report error');
	croak('report error and die');
	cluck('report error with stack trace');
	confess('report error and die with stack trace');
	
	# additional Carp-like logging fucntions
	use RosenLab::Log::Simple qw< :CARP >;
	cavil('carp with verbosity control',$verbosity>=1);
	complain('cluck with verbosity control',$verbosity>=1);
	concede('confess with verbosity control',$verbosity>=1);
	
	# redirect all output to a file
	my $fh = IO::File->new('log.txt','w');
	local $RosenLab::Log::Simple::error_fh = $fh;
	local $RosenLab::Log::Simple::output_fh = $fh;

=head1 DESCRIPTION

This module provides logging as per the core module: L<Log::Message::Simple>. It
has some extended functionality to control its verbosity across multiple scripts
that C<use> it.

This module should not be used for any serious logging activities that require
advanced features such as asynchronous or fail-safe logging. Speed was also not
a consideration in the design of this module.

=cut

# module stuff
package RosenLab::Log::Simple;
BEGIN
{
	our $VERSION = '1.001'; # use string for MakeMaker (and similar)
	$VERSION = eval "$VERSION"; # convert to number for `use` (and similar)
	use Exporter qw< import >;
	our @ISA = qw< Exporter >;
	our @EXPORT = qw< fatal error warning message debug >;
	our @EXPORT_OK = qw< carp croak cluck confess cavil complain concede >;
	our %EXPORT_TAGS =
		(
			'STD' => \@EXPORT,
			'CARP' => \@EXPORT_OK,
			'ALL' => [@EXPORT,@EXPORT_OK],
		);
}

# pragmas
use 5.020;
use strict; use warnings; use warnings::register;
use experimental qw< signatures >;

# core modules
use File::Basename qw< basename >;
use POSIX qw< strftime >;

# cpan modules
use Log::Message private => 0;

=head1 CONSTANTS

=over

=item DEFAULT_ERR_FH

The default file handle to write error messages to. This includes all messages
logged by the functions L<C<fatal>>, L<C<error>>, and L<C<warning>>.

=item DEFAULT_OUT_FH

The default file handle for non-error messages. This includes messages from its
namesake function L<C<message>>, as well as any logging level not explicitly set
to write to L<C<DEFAULT_ERR_FH>>.

=back

=cut

use constant
	{
		DEFAULT_ERR_FH => \*STDERR,
		DEFAULT_OUT_FH => \*STDOUT,
	};

=head1 GLOBAL VARIABLES

=cut

=head2 Log handles

=over

=item $error_fh

This is the file handle that L<C<error>>, L<C<fatal>>, and L<C<warning>> write
to. It defaults to L<C<DEFAULT_ERR_FH>>.

=item $output_fh

This is the file handle that L<C<message>> and L<C<debug>> write to. It defaults
to L<C<DEFAULT_OUT_FH>>.

=back

=cut

local $| = 1;
our $error_fh = DEFAULT_ERR_FH;
our $output_fh = DEFAULT_OUT_FH;

=head2 Log level control

=over

=item $verbosity

The verbosity level. This variable controls which messages are printed to logs
and which are (effectively) ignored.

=back

=cut

our $verbosity = 0;

=begin comment LOCAL VARIABLES

These variables aren't really relevant to usage. They just maintain internal
state. They are documented in this comment block for the benefit of maintainers.

=over

=item log

The logger object.

=back

=end comment

=cut

my $log = new Log::Message;

=head1 CLASS METHODS

=cut

=head2 Methods for interacting with the logging stack.

These functions manipulate the logging stack. The logging stack is implemented
using L<F<Log::Message>>.

=cut

=head3 flush

Removes all of the messages on the stack and returns them (in reverse order) as
L<F<Log::Message::Item>>s. Consult that module's manpage for more detail on how
to use these items.

=cut

sub flush
{
	return reverse $log->flush();
}

=head3 peek

Retrieves the top (last) message on the logging stack (without removing it) as
a L<F<Log::Message::Item>>. Consult that module's manpage for more detail on how
to use these items.

=cut

sub peek
{
	return reverse $log->last();
}

=head3 spy

Retrieves the bottom (first) message on the logging stack (without removing it)
as a L<F<Log::Message::Item>>. Consult that module's manpage for more detail on
how to use these items.

=cut

sub spy
{
	return reverse $log->first();
}

=head3 stack

Retrieves all of the messages on the stack as L<F<Log::Message::Item>>s. Consult
that module's manpage for more detail on how to use these items.

=cut

sub stack
{
	return $log->retrieve( chrono => 1 );
}

=head1 FUNCTIONS

This module provides logging functions as per L<RosenLab::Log::Simple>. It has
some extended functionality to control its verbosity across multiple scripts
that C<use> it.

=cut

# TODO : expose the prefix command for tailoring?
sub _prefix ( $tag )
{
	my $datetime = strftime("%F\t%T",localtime());
	my ($file,$line) = (caller(4))[1,2];
	my $subroutine = (caller(5))[3] // 'main';
	return sprintf("%s\t%s::%s::%s\t[%s]\t",$datetime,basename($file),$subroutine,$line,$tag);
}


=head2 Simple logging functions

These functions are exported by default. They use the L<Log::Message> as per the
L<Log::Message::Simple> module.

=cut

=head3 debug

	debug("MESSAGE",$VERBOSITY);

Prints a message to L<C<$output_fh>> if the verbosity level is greater than or
equal to two (C<2>). The C<$VERBOSITY> argument is optional and defaults to the
package's verbosity level: L<C<$RosenLab::Log::Simple::verbosity>>.

=cut

sub debug ( $message , $prolixity = $verbosity )
{
	$log->store
		(
			message => $message,
			tag     => 'DEBUG',
			level   => 'debug',
			extra   => [$prolixity],
		);
}

{
	package Log::Message::Handlers;
	
	sub debug ( $self , $prolixity = 0 )
	{
		return unless($prolixity >= 2);
		my $tmp_fh = select($RosenLab::Log::Simple::output_fh);
		say RosenLab::Log::Simple::_prefix($self->tag) . $self->message;
		select($tmp_fh);
		return;
	}
}

=head3 error

	error("MESSAGE");

Prints a message to L<C<$error_fh>>. Error messages are always logged.

=cut

sub error ( $message )
{
	$log->store
		(
			message => $message,
			tag     => 'ERROR',
			level   => 'error',
			extra   => [],
		);
}

{
	package Log::Message::Handlers;
	
	sub error ( $self )
	{
		my $tmp_fh = select($RosenLab::Log::Simple::error_fh);
		say RosenLab::Log::Simple::_prefix($self->tag) . $self->message;
		select($tmp_fh);
		return;
	}
}

=head3 fatal

	fatal("MESSAGE");

Prints a message to L<C<$error_fh>> and then dies. Fatal error messages
are always logged.

Special note for this function: Log messages are passed to L<C<CORE::die>> as
well as logged to the stack (which should flush as the program dies). This may
result in near-duplicated error messages being printed in the same location
(depending on what C<$error_fh> points to and where program output is being
directed).

=cut

sub fatal ( $message )
{
	$log->store
		(
			message => $message,
			tag     => 'FATAL',
			level   => 'fatal',
			extra   => [],
		);
}

{
	package Log::Message::Handlers;
	
	sub fatal ( $self , $verbose = 0 )
	{
		my $tmp_fh = select($RosenLab::Log::Simple::error_fh);
		say RosenLab::Log::Simple::_prefix($self->tag) . $self->message;
		select($tmp_fh);
		# `die` like normal (but suppress further output)
		my $code = $!;
		$code = ($? >> 8) unless ($code);
		$code = 255 unless ($code);
		exit($code);
	}
}

=head3 message

	message("MESSAGE",$VERBOSITY);

Prints a message to L<C<$output_fh>> if the verbosity level is greater than or
equal to one (C<1>). The C<$VERBOSITY> argument is optional and defaults to the
package's verbosity level: L<C<$RosenLab::Log::Simple::verbosity>>.

=cut

sub message ( $message , $prolixity = $verbosity )
{
	$log->store
		(
			message => $message,
			tag     => 'MESSAGE',
			level   => 'msg',
			extra   => [$prolixity],
		);
}

{
	package Log::Message::Handlers;
	
	sub msg ( $self , $prolixity = 0 )
	{
		return unless($prolixity >= 1);
		my $tmp_fh = select($RosenLab::Log::Simple::output_fh);
		say RosenLab::Log::Simple::_prefix($self->tag) . $self->message;
		select($tmp_fh);
		return;
	}
}

=head3 warning

	warning("MESSAGE",$CATEGORY);

Prints a message to L<C<$error_fh>> if warnings are enabled for the appropriate
category. The C<$CATEGORY> argument is optional and defaults to the package name
as given by C<caller>.

=cut

sub warning ( $message , $category = caller )
{
	return unless (warnings::enabled($category));
	$log->store(
			message => $message,
			tag     => 'WARNING',
			level   => 'warns',
			extra   => [],
		);
}

{
	package Log::Message::Handlers;
	
	sub warns ( $self )
	{
		my $tmp_fh = select($RosenLab::Log::Simple::error_fh);
		say RosenLab::Log::Simple::_prefix($self->tag) . $self->message;
		select($tmp_fh);
		return;
	}
}

=head2 Carp functions

This module can also optionally expose logging functions for the alternative
warning (and C<die>) functions from Perl's L<F<Carp>> module. These functions
allow log messages to be issued to the stack, but no additional functionality is
provided.

To expose all of these functions at once, use the tag C<:CARP> (or C<:ALL>).

	use RosenLab::Log::Simple qw( :CARP );

Although the function L<C<cluck>> is not exported by the L<F<Carp>> module by
default, it is included in the export tag C<:CARP> (and C<:ALL>).

=head3 carp

	carp("MESSAGE");

=head3 cluck

	cluck("MESSAGE");

=head3 confess

	confess("MESSAGE");

=head3 croak

	croak("MESSAGE");

=cut

foreach my $subroutine ( qw( carp cluck croak confess ) )
{
	no strict 'refs';
	*$subroutine = sub ( $message , @extra )
		{
			local $Carp::CarpLevel += 2;
			$log->store
				(
					message => $message,
					tag     => uc($subroutine),
					level   => $subroutine,
					extra   => [@extra],
				);
		};
}

=head2 Carp-like functions

Three supplementary L<Carp>-like functions are provided. In addition to logging
to the stack, these functions allow messages to be logged conditionally (without
the need to wrap each log statement in an C<if>/C<else> block). Otherwise these
functions are identical to the L<Carp>-like functions they claim to emulate.

=cut

=head3 cavil

	cavil("MESSAGE",$VERBOSE);

L<C<carp>>s, but only if asked to be C<$VERBOSE>.

=cut

sub cavil ( $message , $verbose = 0 , @extra )
{
	return unless($verbose);
	local $Carp::CarpLevel += 2;
	$log->store
		(
			message => $message,
			tag     => 'CARP',
			level   => 'carp',
			extra   => [@extra],
		);
}

=head3 complain

	complain("MESSAGE",$VERBOSE);

L<C<clucks>>s, but only if asked to be C<$VERBOSE>.

=cut

sub complain ( $message , $verbose = 0 , @extra )
{
	return unless($verbose);
	local $Carp::CarpLevel += 2;
	$log->store
		(
			message => $message,
			tag     => 'CLUCK',
			level   => 'cluck',
			extra   => [@extra],
		);
}

=head3 concede

	concede("MESSAGE",$VERBOSE);

L<C<croak>>s, but only if asked to be C<$VERBOSE>.

=cut

sub concede ( $message , $verbose = 0 , @extra )
{
	return unless($verbose);
	local $Carp::CarpLevel += 2;
	$log->store
		(
			message => $message,
			tag     => 'CROAK',
			level   => 'croak',
			extra   => [@extra],
		);
}

=head3 connip

	connip("MESSAGE",$VERBOSE);

L<C<confess>>es, but only if asked to be C<$VERBOSE>.

=cut

sub connip ( $message , $verbose = 0 , @extra )
{
	return unless($verbose);
	local $Carp::CarpLevel += 2;
	$log->store
		(
			message => $message,
			tag     => 'CONFESS',
			level   => 'confess',
			extra   => [@extra],
		);
}

=head1 SEE ALSO

L<Log::Message>

=head1 COPYRIGHT & LICENSE

=cut

1;

__END__
