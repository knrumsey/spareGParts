# Test functions

X <- lhs::randomLHS(300, 3)
y <- apply(X, 1, duqling::ishigami)

gp1 <- mpgp(X, y)
plot(gp1)

gp2 <- svecgp(X, y)
plot(gp2)
