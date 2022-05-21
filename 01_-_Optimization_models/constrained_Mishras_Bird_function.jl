# Constrained Mishra's Bird function in two variables
#
#
# min  f(X) = sin(y)e^((1 - cos(x))^2) + cos(x)e^((1 - sin(y))^2) + (x - y)^2
#  X
#
# X = [x, y]^T
#
# s.t: (x + 5)^2 + (y + 5)^2 - 25 < 0
#
# Bounds (or search domain):
# -10.0  <= x <= 0.0
#  -6.5  <= y <= 0.0
#
# Global minimun:
# X*      = [-3.1302468, -1.5821422]^T
# f(X*)   = -106.7645367
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

# objective function (Mishra's Bird function)
function f(X, opti)

    x = X[1];   y = X[2];
    part_1 = sin(y)*exp((1 - cos(x))^2);
    part_2 = cos(x)*exp((1 - sin(y))^2);
    return part_1 + part_2 + (x - y)^2
end

# constraints function
function f_c(X, opti)

    x = X[1];   y = X[2];
    return [(x + 5)^2 + (y + 5)^2 - 25]
end

# variable bounds
#          min       max
bounds = [-10.0      0.0;    # ----> x
           -6.5      0.0];   # ----> y

# X* (optimal)
#               x           y
X_optimal = [-3.1302468 -1.5821422];
