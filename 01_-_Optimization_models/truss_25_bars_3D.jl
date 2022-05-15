# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 25
# bars
#
#             n
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 25
#  A         i=1
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
# *(1) Li, H.-S., Au, S.-K (2010). Desing optimization using subset simulation
#      algorithm. Structural Safety, 32(6), 384-392.
#
# *(2) Haftka RT, Gurdal Z. Elements of structural optimization. 3rd
#      ed. Dordrecht: Kluwer Academic Publishers; 1992. (pag-245)
#
# -----------------------------------------------------------------------------

# ============== load the model data and truss solver =========================

include("./load_case_1_2_3D_truss_25_bars.jl"); # load the variables of model
include("../00_-_Truss_structure_solver/structure_solver_functions.jl"); # load the structural solver

# ===================== data of optimization problem ==========================

# objective function (weight of structure)
function W(A, model)

    L_i           = model.length_bars; # the length of each bar
    Type_material = model.info_elements[:, 3]; # type of material
    ρ_i           = model.info_material_propieties[Type_material, 4]; # density
    A_i           = A[Type_material]; # areas for each element

    return sum(ρ_i .* L_i .* A_i);
end;

# constraints function
function f_c(A, model)

#   material propieties with A to test
    material_properties_load_1_2 = deepcopy(model.info_material_propieties[:, 1:3]);
    material_properties_load_1_2[:, 2] = A;

#   solution FEM load case 1
    δ1, _, _, σ1, _, _ = solver_truss(model.coordinates,
                                      model.info_elements,
                                      material_properties_load_1_2,
                                      model.load_cases[1],
                                      model.restricted_dof);

#   stress constraints for the load case 1
    σ1_T = σ1 - model.axial_sigma_restricted[:, 2];   # σ - σmax <= 0
    σ1_C = model.axial_sigma_restricted[:, 3]  - σ1;  # σmin - σ <= 0

#   displacement constraint
    nodes           = model.constraint_displacement[:, 1];
    direction       = model.constraint_displacement[:, 2];
    constraint_dof  = Int64.(3*nodes .- 3 .+ direction);

#   positive and negative displacement
#   δ - δmax <= 0
    δ1_p = δ1[constraint_dof] - model.constraint_displacement[:, 3];
#   δmin - δ <= 0
    δ1_n = model.constraint_displacement[:, 4] - δ1[constraint_dof];

#   solution FEM load case 2
    δ2, _, _, σ2, _, _ = solver_truss(model.coordinates,
                                      model.info_elements,
                                      material_properties_load_1_2,
                                      model.load_cases[2],
                                      model.restricted_dof);

#   stress constraints for the load case 2
    σ2_T = σ2 - model.axial_sigma_restricted[:, 2];   # σ - σmax <= 0
    σ2_C = model.axial_sigma_restricted[:, 3]  - σ2;  # σmin - σ <= 0

#   positive and negative displacement
#   δ - δmax <= 0
    δ2_p = δ2[constraint_dof] - model.constraint_displacement[:, 3];
#   δmin - δ <= 0
    δ2_n = model.constraint_displacement[:, 4] - δ2[constraint_dof];

    return [σ1_T; σ1_C; δ1_p; δ1_n; σ2_T; σ2_C; δ2_p; δ2_n]
end;

# variable bounds
bounds = truss_model_25_bars.variables_bounds;
