# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 25
# bars with the PSOPC flyback method
#
#
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 25
#  A
#
# A = [A1, A2, ..., A25]^T
#        
# s.t:
# Stress limits:
#-------------------------------------------------------------------------
#    Variables         Compressive stress [ksi]      Tensile stress [ksi]
#-------------------------------------------------------------------------
#   1       A1                 35.092                      40.0
#   2    A2-A5                 11.590                      40.0
#   3    A6-A9                 17.305                      40.0
#   4  A10-A11                 35.092                      40.0
#   5  A12-A13                 35.092                      40.0
#   6  A14-A17                  6.759                      40.0
#   7  A18-A21                  6.959                      40.0
#   8  A22-A25                 11.082                      40.0
#
# All limits displacement in three direction: +- 0.35 [in]
#
# Bounds:
# 0.01 [in^2] <= Ai <= 3.5 [in], for i=1, 2, ..., 25
#
# Number of restrictions: 172
# * Displacement: 2 (+-) * 3 Direction * 6 Nodes * 2 Cases  = 72
# * Stress: 2 (Compressive and Tensile) * 25 Bars * 2 Cases = 100
#
# keep in mind the units of unknows
# Areas:     [in^2]
# Forces:    [kips = 1_000 lb_f]
# Distances: [in]
# Stress:    [ksi = 1_000psi = 1_000 lb_f/in^2]
# Weight:    [lbf]
# Density:   [lb/in^3]
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
# *(1)  Li LJ, Huang ZB, Liu F, Wu QH. A heuristic particle swarm optimizer for
#       optimization of pin connected structures.
#       Comput Struct 2007;85(7–8):340–9
#
# -----------------------------------------------------------------------------

# ============================ packages =======================================

using Plots # for plots
ENV["GKSwstype"] = "100";

# ======================= PSOPC functions and model ===========================

include("../../../01_-_Optimization_models/truss_25_bars_3D.jl"); # model truss
include("./PSOPC_flyback_functions.jl");                          # load the PSOPC flyback functions

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by W and f_c respectivily

# variable bounds
bounds = truss_model_25_bars.variables_bounds;

num_var, _ = size(bounds); # number of variables

vmax = abs(3.5-0.01);
vbounds = [-vmax*ones(num_var)  vmax*ones(num_var)]; # velocity limits

opt_arg = truss_model_25_bars;  # optional arguments
num_particles = 50;             # number of particles
k_max = 3000;                   # maximum iteration
ω   =  [0.9, 0.4];              # inertia weight
c1  =  [0.8, 0.8];              # acceleration constants
c2  =  [0.8, 0.8];
c3  =  [0.2, 0.2];

values_f, X_g_value, population_data, Swarm_data, fun_evals, const_evals = PSOPC_flyback(W, f_c, bounds,
                                                                                         opt_arg,
                                                                                         vbounds,
                                                                                         num_particles,
                                                                                         k_max, ω, c1, c2, c3);

# solution
solution    = Swarm_data.X_global;
f_solution  = Swarm_data.f_X_global;

println("\n\n=================================================================")
println("Truss 25 bars --> PSOPC flyback")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, values_f, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
