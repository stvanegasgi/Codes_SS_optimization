# Test solver 3D truss, load case 1 for truss with 200 bars
#
#
# keep in mind the units of unknows
# Areas:     [in^2]
# Forces:    [kips=1000lb_f]
# Distances: [in]
# Stress:    [ksi=1000psi=1000lb_f/in^2]
#
# =============================================================================
# ------- December 2021
# -----------------------------------------------------------------------------
# by 
# ------- Steven Vanegas Giraldo -------> stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# *(1) Awad, R. (2021, October). Sizing optimization of truss structures using
#      the political optimizer (PO) algorithm. In Structures
#      (Vol. 33, pp. 4871-4894). Elsevier.
#
# *(2) Degertekin, S. O., & Hayalioglu, M. S. (2013). Sizing truss structures 
#      using teaching-learning-based optimization.
#      Computers & Structures, 119, 177-188.
#
# -----------------------------------------------------------------------------

# ================================ load packages ==============================

using Plots # load the package to plot

# ================== load the unknowns model and functions ====================

include("./load_case_1_2D_truss_200_bars.jl"); # load the model
include("./structure_solver_functions.jl");    # load solver functions

# ============================ solver structure ===============================

a, q, N_e, sigma_e, epsilon_e, a_loc = solver_truss(xcor,
                                                    elements,
                                                    material_properties,
                                                    forces,
                                                    restricted_dof);

# ============================== the results ==================================

println("\nDisplacements\n")
show(stdout, "text/plain", a)

println("\n\n\nVector of equilibrium nodal forces\n")
show(stdout, "text/plain", q)

println("\n\n\nElements normal forces\n")
show(stdout, "text/plain", N_e)

println("\n\n\nElements stress\n")
show(stdout, "text/plain", sigma_e)

println("\n\n\nElements strains\n")
show(stdout, "text/plain", epsilon_e)

println("\n\n\nLocal displacements\n")
show(stdout, "text/plain", a_loc)
println("\n")

# ============================== plot the structue ============================

# indexes for better reading and writing in the rest of the code
# -----------------------------------------------------------------------------

# for coordinates
X = 1; # x coordinate
Y = 2; # y coordinate
Z = 3; # z coordinate

# for nodes
N_start = 1; # element start node
N_end   = 2; # element end node

# ========================== information of structure =========================

L2G = elements[:, 1:2];             # local to global matrix
num_elem, _     = size(L2G);        # number of elements
num_nodes, _    = size(xcor);       # number of nodes

# =============================== plot structure ==============================

# plot the elements
for element in 1:num_elem

    p = plot!([xcor[L2G[element, N_start], X], xcor[L2G[element, N_end], X]],
              [xcor[L2G[element, N_start], Y], xcor[L2G[element, N_end], Y]],
              label="",
              color=:blue,
              xlabel="x-coordinate",
              ylabel="y-coordinate",
              aspect_ratio=:equal,
              markersize = 2,
              show=true)

#   the number of element
    annotate!((xcor[L2G[element, N_start], X] + xcor[L2G[element, N_end], X])/2,
              (xcor[L2G[element, N_start], Y] + xcor[L2G[element, N_end], Y])/2,
              text("$element", :red, :center, 3), label="")
end

# plot the nodes of the structure
scatter!(xcor[:, X], xcor[:, Y],     # data
         marker=:c,                  # form marker (cicle)
         color=:white,               # color marker
         markerstrokecolor=:black,   # marker stroker color
         markerstrokewidth=1,        # marker stroker width
         label="",                   # no label
         show=true)

# nuumber nodes
for node in 1: num_nodes

    annotate!(xcor[node, X], xcor[node, Y], text("$node", :black, :center, 3), label="")
end

savefig("./structure_200_bars_case_1.svg")
