rm(list=ls())
library(Matrix)

# Number of points and dimensions
I <- 25
J <- 3

# Create a random 25x3 matrix
A <- matrix(runif(I * J, min = 0.3, max = 0.7), nrow = I, ncol = J)

# Compute the rank
rankMatrix(A)

# Normalized data such that each row sums to 1
A_norm <- diag(1/apply(A, 1, sum)) %*% A


ternary_to_xy_A <- t(apply(A_norm, 
                           1, 
                           function(x) {
                             x1 <- 0.5 * (2*x[2] + x[3])
                             x2 <- (sqrt(3)/2) * x[3]
                             return(c(x1, x2))
                           }))

basis <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, byrow = TRUE)
rownames(basis) <- c("e1", "e2", "e3")

ternary_to_xy_basis <- t(apply(basis, 
                               1, 
                               function(x) {
                                 x1 <- 0.5 * (2*x[2] + x[3])
                                 x2 <- (sqrt(3)/2) * x[3]
                                 return(c(x1, x2))
                               }))

# draw the inscribed circle
v1 <- ternary_to_xy_basis[1, ]
v2 <- ternary_to_xy_basis[2, ]
v3 <- ternary_to_xy_basis[3, ]
a <- sqrt(sum((v2 - v3)^2))
b <- sqrt(sum((v1 - v3)^2))
c <- sqrt(sum((v1 - v2)^2))
px <- (a * v1[1] + b * v2[1] + c * v3[1]) / (a + b + c)
py <- (a * v1[2] + b * v2[2] + c * v3[2]) / (a + b + c)
s <- (a + b + c) / 2
area <- sqrt(s * (s - a) * (s - b) * (s - c))
r <- area / s

incenter <- (a*v1 + b*v2 + c*v3) / (a + b + c)
theta <- seq(0, 2*pi, length.out = 300)
incircle_x <- incenter[1] + r * cos(theta)
incircle_y <- incenter[2] + r * sin(theta)

png(file=paste0("plots//halfpurity.png"),
    width=500, height=500)
par(mar = rep(0.2, 4))
plot(ternary_to_xy_basis, asp = 1, type="n", axes = F, xlim=c(-0.3,1.2),ylim=c(-0.3,1.2),xlab="",ylab="", cex.lab = 2.75, cex = 2.75)

points(ternary_to_xy_basis, cex = 2.75, lwd = 3, col = "red", pch = 4)
lines(rbind(ternary_to_xy_basis, ternary_to_xy_basis[1, ]), cex = 2.75, lwd = 1.5, col = "red", lty = 2)
points(ternary_to_xy_A, cex = 2.75, lwd = 3, col = "blue")

lines(incircle_x, incircle_y, col = "green", cex = 2.75, lty = 3, lwd = 2)

text(ternary_to_xy_basis[1,1]-0.1, ternary_to_xy_basis[1,2], expression(bold(e)[1]), cex = 2.75, lwd = 3)
text(ternary_to_xy_basis[2,1]+0.12, ternary_to_xy_basis[2,2], expression(bold(e)[2]), cex = 2.75, lwd = 3)
text(ternary_to_xy_basis[3,1], ternary_to_xy_basis[3,2]+0.1, expression(bold(e)[3]), cex = 2.75, lwd = 3)

points_on_line <- function(p1, p2, t_vals = c(1/3, 2/3)) {
  sapply(t_vals, function(t) t*p1 + (1-t)*p2)
}

# Compute points for each edge
edge12 <- t(points_on_line(ternary_to_xy_basis[1, ], ternary_to_xy_basis[2, ], t_vals = c(7/10)))
edge21 <- t(points_on_line(ternary_to_xy_basis[1, ], ternary_to_xy_basis[2, ], t_vals = c(1-8/10)))

edge23 <- t(points_on_line(ternary_to_xy_basis[2, ], ternary_to_xy_basis[3, ], t_vals = c(7/10)))
edge32 <- t(points_on_line(ternary_to_xy_basis[2, ], ternary_to_xy_basis[3, ], t_vals = c(1-8/10)))


edge31 <- t(points_on_line(ternary_to_xy_basis[3, ], ternary_to_xy_basis[1, ], t_vals = c(9/12)))
edge13 <- t(points_on_line(ternary_to_xy_basis[3, ], ternary_to_xy_basis[1, ], t_vals = c(1-10/12)))

# Get x and y coordinates
extra_points <- rbind(edge12, edge21, edge23, edge32, edge31, edge13)

# Plot the extra points
points(extra_points, cex = 2.75, lwd = 3, col="blue")
lines(rbind(extra_points, extra_points[1, ]), cex = 2.75, lwd = 1)

dev.off()

