#' Create simple log message "object"s.
#'
#' @description
#' The following functions serve as log message generators. They will create
#' S3 objects with, at least, the classes "message" and "condition" (in that
#' order). These may be used by any of the log messaging functions provided by
#' this package (e.g. [message()] or [print.message()]).
#'
#' @export
#' @usage simpleMessage(message, call = NULL, domain = NULL)
#' @param message
#'   An R object, which is coerced to character (via [as.character()]).
#' @param call
#'   The calling expression. This option exists for compatibility with the
#    built-in function [base::simpleMessage()] (and similar) and should not be
#'   used in most circumstances. If `NULL`, the name of the calling function is
#'   is deduced from the stack frame.
simpleMessage <- function(message , call = NULL , domain = NULL)
{
	# already a message, nothing to do here
	if (inherits(message,"message"))
	{
		return(message)
	}
	# get the call, unless provided
	call <- ifelse(is.null(call) , caller(which = -2) , call)
	# construct the message
	message <- paste0(
		.makeMessagePrefix(level = "MESSAGE" , call = call),
		.makeMessage(message , domain = domain , appendLF = TRUE)
		)
	# build and return the message object
	class <- c("message","condition")
	return(structure(list(message = message , call = call) , class = class))
}

#' @export simpleWarning
#' @export simpleWarningMessage
#' @rdname simpleMessage
#' @usage simpleWarning(message, call = NULL, domain = NULL)
#' @usage simpleWarningMessage(message, call = NULL, domain = NULL)
simpleWarning <- simpleWarningMessage <- function(message , call = NULL , domain = NULL)
{
	# already a message, nothing to do here
	if (inherits(message,"message"))
	{
		return(message)
	}
	# get the call, unless provided
	call <- ifelse(is.null(call) , caller(which = -2) , call)
	# construct the message
	message <- paste0(
		.makeMessagePrefix(level = "WARNING" , call = call),
		.makeMessage(message,domain = domain,appendLF = TRUE)
		)
	# build and return the message object
	class <- c("warning","message","condition")
	return(structure(list(message = message, call = call),class = class))
}

#' @export simpleError
#' @export simpleErrorMessage
#' @rdname simpleMessage
#' @usage simpleError(message, call = NULL, domain = NULL)
#' @usage simpleErrorMessage(message, call = NULL, domain = NULL)
simpleError <- simpleErrorMessage <- function(message , call = NULL , domain = NULL)
{
	# already a message, nothing to do here
	if (inherits(message,"message"))
	{
		return(message)
	}
	# get the call, unless provided
	call <- ifelse(is.null(call) , caller(which = -2) , call)
	message <- paste0(
		.makeMessagePrefix(level = "ERROR"),
		.makeMessage(message,domain = domain,appendLF = TRUE)
		)
	# construct the message
	class <- c("error","message","condition")
	return(structure(list(message = message, call = call),class = class))
}

#' @export simpleDebug
#' @export simpleDebugMessage
#' @rdname simpleMessage
#' @usage simpleDebug(message, call = NULL, domain = NULL)
#' @usage simpleDebugMessage(message, call = NULL, domain = NULL)
simpleDebug <- simpleDebugMessage <- function(message , call = NULL , domain = NULL)
{
	# already a message, nothing to do here
	if (inherits(message,"message"))
	{
		return(message)
	}
	# get the call, unless provided
	call <- ifelse(is.null(call) , caller(which = -2) , call)
	# construct the message
	message <- paste0(
		.makeMessagePrefix(level = "DEBUG" , call = call),
		.makeMessage(message,domain = domain,appendLF = TRUE)
		)
	# build and return the message object
	class <- c("debug","message","condition")
	return(structure(list(message = message, call = call),class = class))
}

# @usage .makeMessagePrefix(...,
#   level = toupper(sys.call(-2)[[1]]),
#   domain = NULL , appendLF = N...,
#   call = NULL)
# @param file
#   The file (connection) to which the log message is written.
# @param verbosity
#   The verbosity level to apply this message. Defaults to grabbing this value
#   from the parent environment (via [`dynGet()`]).
# @param domain
#   Translate log message via Native Language Support. See [base::gettext()]
#   more details about this. If `NA`, messages will not be translated; `NULL`
#   (default) implies domain "R-pkg".
# @param appendLF
#   IGNORED. This option exists for compatibility with built-in functions with
#   the same name as any of these functions. See the details section for more
#   information about line breaks in messages.
.makeMessagePrefix <- function(...,
	level = toupper(sys.call(-2)[[1]]),
	domain = NULL , appendLF = FALSE,
	call = NULL
	)
{
	# get a timestamp for this message
	datetime <- strftime(Sys.time(),"%F%t%T")
	# find the calling function
	call <- ifelse(is.null(call) , caller(which = -3) , call)
	# find the calling file
	script <- RosenLab::script()
	script <- ifelse(is.na(script),"R",script)
	return(sprintf("%s\t%s::%s\t[%s]\t",datetime,script,call,level))
	# TODO find the calling line
	line <- NaN
	return(sprintf("%s\t%s::%s::%s\t[%s]\t",datetime,script,call,line,level))
}
