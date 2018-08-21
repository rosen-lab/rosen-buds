#' Straightforward logging facilities for Rosen Lab jobs.
#'
#' @description
#' These functions aim to provide a base logging interface for users using a
#' familiar "API". They are not intended to be used for any serious logging
#' activities that require advanced features, such as asynchronous or fail-safe
#' logging. Speed was also not a consideration in the design of these functions.
#'
#' @details
#' There are five (5) log levels, with corresponding function names, which in
#' order of increasing verbosity are: [fatal()], [error()], [warning()],
#' [message()], and [debug()]. Errors, fatal errors, and warnings "always" write
#' their messages to the appropriate connection (see below). Messages, debugging
#' or otherwise, will only write their messages to this connection if value of
#' their [`verbosity`] argument is high enough.
#'
#' By default, all messaging functions will "inherit" their message connection
#' and the verbosity level from the calling environment (via [base::dynGet()]).
#' The message connection will inherit from the variable `message_connection`,
#' while the verbosity level inherits from the `verbosity` variable. Whenever
#' either of these variables are not defined (or otherwise uninheritable), the
#' values [stderr()] and `0` are substituted, respectively. These values may
#' be changed ona per-call basis, via the `file` and `verbosity` arguments.
#'
#' Some of these functions mask built-in functions in the base R package. Both
#' [warning()] and [message()] *should* be seamlessly compatible with these
#' built-ins ([base::warning()] and [base::message()], respectively). The
#' [debug()] function is not backwards compatible, although this function
#' should probably(?) not be used in any released code; developers can still
#' use [base::debug()] for this purpose.
#'
#' The built-in functions [suppressWarnings()] and [suppressMessages()],
#' likewise, *should* work as before with these functions. Their use is
#' discouraged, however.
#'
#' @export
#' @usage message(...,
#'   file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION),
#'   verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
#'   domain = NULL , appendLF = NULL)
#' @param file
#'   The file (connection) to which the log message is written.
#' @param verbosity
#'   The verbosity level to apply this message. Defaults to grabbing this value
#'   from the parent environment (via [`dynGet()`]).
#' @param domain
#'   Translate log message via Native Language Support. See [base::gettext()]
#'   more details about this. If `NA`, messages will not be translated; `NULL`
#'   (default) implies domain "R-pkg".
#' @param appendLF
#'   IGNORED. This option exists for compatibility with built-in functions with
#'   the same name as any of these functions. See the details section for more
#'   information about line breaks in messages.
message <- function(...,
	file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION),
	verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
	domain = NULL , appendLF = NA
	)
{
	# translate objects into a message (condition) via simpleMessage
	arguments <- list(...)
	conditions <- lapply(arguments , RosenLab::simpleMessage , call = caller(which = -4))
	# specify the default handler
	withRestarts(
		lapply(conditions,
			function(condition)
				{
					signalCondition(condition)
					print.message(condition)
				}
			),
		muffleMessage = function() { return(NULL) }
		)
	return(invisible(sapply(conditions,`[[`,"message")))
}

#' @export
#' @rdname message
#' @usage warning(...,
#'   file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION),
#'   verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
#'   call. = TRUE , immediate. = FALSE , noBreaks. = TRUE
#'   domain = NULL , appendLF = NULL)
#' @param call.
#'   IGNORED. This named argument exists only for "compatibility" with existing
#'   calls to [base::warning()], [base::stop()], etc. The calling context is
#'   always logged.
#' @param immediate.
#'   IGNORED. This named argument exists only for "compatibility" with existing
#'   calls to [base::warning()], [base::stop()], etc. Warnings are always logged
#'   immediately.
#' @param noBreaks.
#'   IGNORED. This named argument exists only for "compatibility" with existing
#'   calls to [base::warning()], [base::stop()], etc. See the details section
#'   for more information about line breaks in messages.
warning <- warningMessage <- function(...,
	file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION),
	verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
	call. = TRUE , immediate. = FALSE , noBreaks. = FALSE,
	domain = NULL , appendLF = NA
	)
{
	# translate objects into a warning (condition) via simpleWarning
	arguments <- list(...)
	conditions <- lapply(arguments , RosenLab::simpleWarning , call = caller(which = -4))
	# specify the default handler
	withRestarts(
		lapply(conditions,
			function(condition)
				{
					message <- conditionMessage(condition)
					call <- conditionCall(condition)
					.Internal(.signalCondition(condition,message,call))
					print.message(condition)
					#.Internal(.dfltWarn(condition$message,condition$call))
				}
			),
		muffleMessage = function() { return(NULL) }
		)
	return(invisible(sapply(conditions,`[[`,"message")))
}

#' @export
#' @rdname message
#' @usage error(...,
#'   file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION),
#'   verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
#'   call. = TRUE , immediate. = FALSE , noBreaks. = TRUE
#'   domain = NULL , appendLF = NULL)
error <- errorMessage <- function(...,
	file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION),
	verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
	call. = TRUE , immediate. = FALSE , noBreaks. = FALSE,
	domain = NULL , appendLF = NA
	)
{
	# translate objects into a warning (condition) via simpleError
	arguments <- list(...)
	conditions <- lapply(arguments , RosenLab::simpleError , call = caller(which = -4))
	# specify the default handler
	withRestarts(
		lapply(conditions,
			function(condition)
				{
					message <- conditionMessage(condition)
					call <- conditionCall(condition)
					.Internal(.signalCondition(condition,message,call))
					print.message(condition)
					#.Internal(.dfltWarn(condition$message,condition$call))
				}
			),
		muffleMessage = function() { return(NULL) }
		)
	return(invisible(sapply(conditions,`[[`,"message")))
}

#' @export
#' @rdname message
#' @usage fatal(... , status = ...,
#'   file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION),
#'   verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
#'   call. = TRUE , domain = NULL , appendLF = NULL)
#' @param status
#'   The exit status. See [base::quit()] for more details about appropriate
#'   values for this argument.
fatal <- fatalMessage <- fatalErrorMessage <- function(... , status = 255,
	file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION),
	verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
	call. = TRUE , domain = NULL , appendLF = NA
	)
{
	# translate objects into a warning (condition) via simpleWarning
	arguments <- list(...)
	conditions <- lapply(arguments , RosenLab::simpleWarning , call = caller(which = -4))
	# specify the default handler
	withRestarts(
		lapply(conditions,
			function(condition)
				{
					message <- conditionMessage(condition)
					call <- conditionCall(condition)
					.Internal(.signalCondition(condition,message,call))
					print.message(condition)
					#.Internal(.dfltStop(condition$message,condition$call))
				}
			),
		muffleMessage = function() { return(NULL) }
		)
	quit(save = "no" , status = status)
}

#' @export
#' @rdname message
#' @usage debug(...,
#'   file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION),
#'   verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
#'   domain = NULL , appendLF = NULL)
debug <- debugMessage <- function(...,
	file = dynGet("logfile" , inherits = TRUE , RosenLab:::DEFAULT_MESSAGE_CONNECTION),
	verbosity = dynGet("verbosity" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_VERBOSITY),
	domain = NULL , appendLF = NA
	)
{
	# translate objects into a message (condition) via simpleDebug
	arguments <- list(...)
	conditions <- lapply(arguments , RosenLab::simpleDebug , call = caller(which = -4))
	# specify the default handler
	withRestarts(
		lapply(conditions,
			function(condition)
				{
					signalCondition(condition)
					print.message(condition)
				}
			),
		muffleMessage = function() { return(NULL) }
		)
	return(invisible(sapply(conditions,`[[`,"message")))
}

# The default message connection
DEFAULT_MESSAGE_CONNECTION <- stderr()

# The default verbosity level
DEFAULT_VERBOSITY <- 0
