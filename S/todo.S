two levelplots on one "page"

add na.action and subset arguments???

DO FIRST:
	residuals.gstat:
		cross validation
		w.r.t. trend

TEST:
	verify heavy simulation -- memory leaks?

add: 
  ?abline(h=0) for cross variogram plotting
  verification to gstat objects, should never contain nonsense when
    called through gstat()
  ? gstat object a; a[c(2:3)]<- and so on.

Later:
	lhs?
	user-supplied variogram models -- how about fitting them?
	passing through as function or as table?
	fitting a Matern model -- derivative?
	fitting using R functions (writing R extensions, 5.8)
