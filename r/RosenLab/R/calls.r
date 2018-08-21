#' Functions to interrogate the call stack.
#'
#' These functions provide limited access to a function's environment. They are
#' convenience functions. For more direct access see the help pages linked at
#' the end of this page.
#'
#' @seealso sys.parent

#' @details
#' [caller()] gets the character "name" of the calling function. This function
#' is very similar to the [base::sys.call()] function, except that it returns
#' a function's name (converted to a string with [as.character()]) rather than
#' the call itself. If the call is `NULL`, as when called from an interative R
#' shell, the string "script" is returned.
#'
#' As with the [base::sys.call()] function, this function will never return
#' "itself". That is, `caller(0)` will return the name of the calling function.
#'
#' @export
#' @usage caller(which = 0)
#' @param which
#'   The frame number on the call stack. If negative, the number of frames less
#'   than the current frame.
caller <- function(which = 0)
{
	# get the call (but looking back one)
	nframes <- sys.nframe()
	which <- ifelse(nframes == which,
		which + 1, # skip this function call when looking "forward"
		ifelse(which <= 0,
			which - 1, # remove this function stack when looking "backward"
			which
			)
		)
	call <- sys.call(which = which)
	# extract the call name as a string
	caller <- ifelse(is.null(call) , "script" , as.character(call[[1]]))
	return(caller)
}

#' @details
#' [script()] gets the character "name" of the executing script. If called from
#' the interactive command line interface, this function returns `NA` (that is,
#' [`NA_character_`]), otherwise the name is parsed from the command line
#' arguments and returned as a string. The name, including the file path, is
#' returned as-is, without (e.g.) canonicalization.
#'
#' @export
#' @rdname caller
#' @usage script()
script <- function()
{
	# grab the argument (options) provided to the script when called
	options <- commandArgs(trailingOnly = FALSE)
	# comb through the options and filter out `--file`, then remove that prefix
	file.pattern <- '--file='
	name <- sub(file.pattern,'',options[grep(file.pattern,options)])
	# if there was no matching pattern, name will be empty but *not* NA or NULL
	name <- ifelse(length(name) == 0 , NA_character_ , name)
}
