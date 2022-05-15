# Constrained Rosembrock function in two variables, case 1
#
#
# min  f(X) = (1 - x)^2 + 100*(y - x^2)^2
#  X
#
# X = [x, y]^T
#
# s.t: --> (x - 1)^3 - y + 1 <= 0
#      --> x + y - 2 <= 0
#
# Bounds (or search domain):
# -1.5  <= x <= 1.5
# -0.5  <= y <= 2.5
#
# Global minimun:
# X*      = [1.0, 1.0]^T
# f(X*)   = 0
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

# objective function (Rosenbrock function)
function f(X, opti)

    x = X[1];   y = X[2];
    return (1 - x)^2 + 100*(y - x^2)^2
end;

# constraints function
function f_c(X, opti)

    x = X[1];   y = X[2];
    return [(x - 1)^3 - y + 1, x + y - 2]
end;

# variable bounds
#          min       max
bounds = [-1.5       1.5     # ----> x
          -0.5       2.5];   # ----> y
