# Constrained Gomez and Levy function (modified) in two variables
#
#
# min  f(X) = 4x^2 - 2.1x^4 + (1/3)x^6 + xy - 4y^2 + 4y^4
#  X
#
# X = [x, y]^T
#
# s.t: -sin(4 pi x) + 2 sin(2 pi y)^2 - 1.5 <= 0
#
# Bounds:
# -1.0 <= x <= 0.75
# -1.0 <= y <= 1.0
#
# Global minimun:
# X*     = [0.08984201, -0.7126564]^T
# f(X*)  = -1.031628453
#
# =============================================================================
# DATE:    November 2021
# WHO:     Steven Vanegas Giraldo
# EMAIL:   stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# * https://en.wikipedia.org/wiki/Test_functions_for_optimization
#
# -----------------------------------------------------------------------------

# ===================== data of optimization problem ==========================

# objective function (Gomez and levy function)
function f(X, opti)

    x = X[1];
    y = X[2];
    return 4*x^2 - 2.1*x^4 + (1/3)*x^6 + x*y - 4*y^2 + 4*y^4
end

# constraints function
function f_c(X, opti)

    x = X[1];
    y = X[2];
    return [-sin(4*π*x) + 2*(sin(2*π*y))^2 - 1.5]
end

# variable bounds
bounds = [-1.0   0.75;    # ----> x
          -1.0   1.0 ];   # ----> y

# X* (optimal)
#               x           y
X_optimal = [0.08984201 -0.7126564];
