# Constrained Townsend function (modified) in two variables
#
#
# min  f(X) = -(cos((x - 0.1) * y))^2 - x*sin(3x + y)
#  X
#
# X = [x, y]^T
#
# s.t: x^2 + y^2 -(2cos(t) - 0.5cos(2t) - 0.25cos(3t) - 0.125cos(4t))^2 - (2sin(t))^2 < 0
# t = atan2(x,y)
#
# Bounds:
# -2.25 <= x <= 2.25
# -2.5  <= y <= 1.75
#
# Global minimun:
# X*     = [2.0052938, 1.1944509]^T
# f(X*)  = -2.0239884
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

# objective function (Townsend function)
function f(X, opti)

    x = X[1];   y = X[2];
    return -(cos((x - 0.1) * y))^2 - x*sin(3*x + y)
end

# constraints function
function f_c(X, opti)

    x = X[1];   y = X[2];
    t = atan(y, x);
    return [x^2 + y^2 - (2*cos(t) - 0.5*cos(2*t) - 0.25*cos(3*t) - 0.125*cos(4*t))^2 - (2*sin(t))^2]
end

# variable bounds
bounds = [-2.25   2.25;    # ----> x
          -2.5    1.75];   # ----> y
