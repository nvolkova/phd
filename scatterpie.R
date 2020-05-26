corr_scatterpie <- function (x, y, p, r, xlab = "", ylab = "", circles = FALSE, 
          lwd.circle = rep(1, length(x)), lty.circle = rep(1, length(x)), 
          add = FALSE, col.circle='black', ...) 
{
  if (!add) 
    plot(x, y, xlab = xlab, ylab = ylab, pch = NA)
  for (i in seq_along(x)) {
    mg14:::.pie(p[i, ], x0 = x[i], y0 = y[i], radius = r[i], add = TRUE, 
         ...)
    if (circles) {
      u <- par("usr")
      pr <- (u[2] - u[1])/(u[4] - u[3])
      fr <- par("pin")[1]/par("pin")[2]
      polygon(x[i] + cos(seq(0, 2 * pi, l = 100)) * r[i], 
              y[i] + sin(seq(0, 2 * pi, l = 100)) * r[i]/pr * 
                fr, col = NA, lty = lty.circle[i], lwd = lwd.circle[i], border=col.circle)
    }
  }
}