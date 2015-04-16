".First.lib" <-
function(lib, pkg) {
#ifdef R
	require(lattice)
	library.dynam("gstat", pkg, lib)
	.Call("gstat_init", as.integer(1), PACKAGE = "gstat")
#endif
#ifdef SP5
	setOldClass(c("variogram", "data.frame"))
	setOldClass(c("variogram.cloud", "data.frame"))
	setOldClass(c("variogram.model", "data.frame"))
	setOldClass(c("point.pairs", "data.frame"))
	setOldClass(c("gstat", "list"))
	.Call("gstat_init", as.integer(1))
#endif
}
variogram <- function(object, ...) UseMethod("variogram")
