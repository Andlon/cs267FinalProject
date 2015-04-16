"vgm.panel.xyplot" <-
function (x, y, type = "p", pch = plot.symbol$pch, col, 
	col.line = plot.line$col, col.symbol = plot.symbol$col, lty = 
	plot.line$lty, cex = plot.symbol$cex, lwd = plot.line$lwd, 
	model = model, labels, ...) 
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (length(x) > 0) {
        if (!missing(col)) {
            if (missing(col.line)) 
                col.line <- col
            if (missing(col.symbol)) 
                col.symbol <- col
        }
        plot.symbol <- trellis.par.get("plot.symbol")
        plot.line <- trellis.par.get("plot.line")
#ifdef R
        lpoints(x = x, y = y, cex = cex, col = col.symbol, pch = pch)
#else
        points(x = x, y = y, cex = cex, col = col.symbol, pch = pch)
#endif
        if (!missing(model)) {
            ret <- variogram.line(model, max(x))
#ifdef R
            llines(x = ret$dist, y = ret$gamma, lty = lty, col = col.line, 
                lwd = lwd)
#else
            lines(x = ret$dist, y = ret$gamma, lty = lty, col = col.line, 
                lwd = lwd)
#endif
        }
		if (!is.null(labels))
#ifdef R
			ltext(x = x + 0.03 * max(x), y = y, labels = labels)
#else
			text(x = x + 0.03 * max(x), y = y, labels = labels)
#endif
    }
}
