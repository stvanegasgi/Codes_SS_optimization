# Structural model unknowns for 3D truss, load case 1, 2 and 3 for truss with
# 200 bars with the PSO flyback method
#
#
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 200
#  A
#
# A = [A1, A2, ..., A200]^T
#
# s.t:
# Stress limits: (compresive -) and (tensile +)
# Tensile limits:     all bars 10 [ksi]
# Compressive limits: all bar -10 [ksi]
#
# Limits displacement:
# No restricted
#
# Bounds:
# --> 0.1  [in^2] <= A_i <= 20.0 [in^2] for i=1, 2, ..., 200
#
# Number of restrictions: 1200
# * Displacement: = 0
# * Stress: 2 (Compressive and Tensile) * 200 Bars * 3 load case = 1200
#
# keep in mind the units of unknows
# Areas:     [in^2]
# Forces:    [kips = 1_000 lb_f]
# Distances: [in]
# Stress:    [ksi = 1_000psi = 1_000 lb_f/in^2]
# Weight:    [lbf]
# Density:   [lb/in^3]
#
#
# =============================================================================
# DATE:    December 2021
# WHO:     Steven Vanegas Giraldo
# EMAIL:   stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# *(1)  Li LJ, Huang ZB, Liu F, Wu QH. A heuristic particle swarm optimizer for
#       optimization of pin connected structures.
#       Comput Struct 2007;85(7–8):340–9
#
# -----------------------------------------------------------------------------


# ============================ packages =======================================

using Plots # for plots
ENV["GKSwstype"] = "100";

# ======================= PSO functions and model =============================

include("../../../01_-_Optimization_models/truss_200_bars_2D.jl"); # model
include("./PSO_flyback_functions.jl");                             # load the PSO flyback functions

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by W and f_c respectivily

# variable bounds
bounds = truss_model_200_bars.variables_bounds;

num_var, _ = size(bounds); # number of variables

vmax = abs(20.0 - 0.1);
vbounds = [-vmax*ones(num_var)  vmax*ones(num_var)]; # velocity limits

opt_arg = truss_model_200_bars;  # optional arguments
num_particles = 50;              # number of particles
k_max = 3000;                    # maximunm iteration
ω   =  [0.9, 0.4];               # inertia weight
c1  =  [0.8, 0.8];               # acceleration constants
c2  =  [0.8, 0.8];

values_f, X_g_value, population_data, Swarm_data, fun_evals, const_evals = PSO_flyback(W, f_c, bounds,
                                                                                       opt_arg,
                                                                                       vbounds,
                                                                                       num_particles,
                                                                                       k_max, ω, c1, c2);

# solution
solution    = Swarm_data.X_global;
f_solution  = Swarm_data.f_X_global;

println("\n\n=================================================================")
println("Truss 200 bars --> PSO flyback")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, values_f, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
