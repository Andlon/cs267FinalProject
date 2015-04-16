"gstat.debug" <-
function(level = 0)
{
#ifdef R
	invisible(.Call("gstat_debug_level", as.integer(level), PACKAGE = "gstat"))
#else
	invisible(.Call("gstat_debug_level", as.integer(level)))
#endif
}
