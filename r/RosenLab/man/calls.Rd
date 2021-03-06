% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calls.r
\name{calls}
\alias{calls}
\alias{caller}
\alias{script}
\title{Functions to interrogate the call stack.}
\usage{
caller(which = 0)

script()
}
\arguments{
\item{which}{The frame number on the call stack. If negative, the number of frames less
than the current frame.}
}
\description{
These functions provide limited access to a function's environment. They are
convenience functions. For more direct access see the help pages linked at
the end of this page.
}
\details{
\code{\link[=caller]{caller()}} gets the character "name" of the calling function. This function
is very similar to the \code{\link[base:sys.call]{base::sys.call()}} function, except that it returns
a function's name (converted to a string with \code{\link[=as.character]{as.character()}}) rather than
the call itself. If the call is \code{NULL}, as when called from an interative R
shell, the string "script" is returned.

As with the \code{\link[base:sys.call]{base::sys.call()}} function, this function will never return
"itself". That is, \code{caller(0)} will return the name of the calling function.

\code{\link[=script]{script()}} gets the character "name" of the executing script. If called from
the interactive command line interface, this function returns \code{NA} (that is,
\code{\link{NA_character_}}), otherwise the name is parsed from the command line
arguments and returned as a string. The name, including the file path, is
returned as-is, without (e.g.) canonicalization.
}
\seealso{
sys.parent
}
