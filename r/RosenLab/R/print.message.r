#' @inherit base::print
#' @export
#' @param file
#'   The file (connection) to which the log message is written.
print.message <- function(message,...,
	file = dynGet("logfile" , inherits = TRUE , ifnotfound = RosenLab:::DEFAULT_MESSAGE_CONNECTION)
	)
{
	cat(message$message, file = file, sep = "")
	return(invisible(message))
}
