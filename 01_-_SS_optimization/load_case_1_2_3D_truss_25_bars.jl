# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 25
# bars
#
#
# Load conditions:
# * Load case 1 for truss with 25 bars
# * Load case 2 for truss with 25 bars
#-----------------------------------------------------
#  Node      Case(1) [kips]         Case(2) [kips]
#            Px     Py     Pz       Px     Py     Pz
#-----------------------------------------------------
#   1        0.0   20.0   -5.0      1.0   10.0   -5.0
#   2        0.0  -20.0   -5.0      0.0   10.0   -5.0
#   3        0.0    0.0    0.0      0.5    0.0    0.0
#   6        0.0    0.0    0.0      0.5    0.0    0.0
#-----------------------------------------------------
#
# Density of bar material: 0.1 [lb/in^3]
#
# Modulus of elasticity: 10_000 [ksi] (all bars)
#
# All areas bounds: 0.01 [in^2] <= A_i <= 3.5 [in^2]
#
# All limits displacement in three direction: +- 0.35 [in]
#
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
#
# The reference [1] has some error in stress limits, see the reference [2]
#
# =============================================================================
# ------- November 2021
# -----------------------------------------------------------------------------
# by 
# ------- Steven Vanegas Giraldo -------> stvanegasgi@unal.edu.co
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

# node coordinates
#                     x         y         z
xcor =   Float64[   -37.5      0.0       200.0;  #  1
                     37.5      0.0       200.0;  #  2
                    -37.5     37.5       100.0;  #  3
                     37.5     37.5       100.0;  #  4
                     37.5    -37.5       100.0;  #  5
                    -37.5    -37.5       100.0;  #  6
                   -100.0    100.0         0.0;  #  7
                    100.0    100.0         0.0;  #  8
                    100.0   -100.0         0.0;  #  9
                   -100.0   -100.0         0.0]; # 10

# elements information
#                N_start   N_end   Type_material
elements = Int64[    1       2          1;       #  1
                     1       4          2;       #  2
                     2       3          2;       #  3
                     1       5          2;       #  4
                     2       6          2;       #  5
                     2       4          3;       #  6
                     2       5          3;       #  7
                     1       3          3;       #  8
                     1       6          3;       #  9
                     3       6          4;       # 10
                     4       5          4;       # 11
                     3       4          5;       # 12
                     5       6          5;       # 13
                     3      10          6;       # 14
                     6       7          6;       # 15
                     4       9          6;       # 16
                     5       8          6;       # 17
                     4       7          7;       # 18
                     3       8          7;       # 19
                     5      10          7;       # 20
                     6       9          7;       # 21
                     6      10          8;       # 22
                     3       7          8;       # 23
                     4       8          8;       # 24
                     5       9          8];      # 25

# materials properties
#                              Type   Areas     E        ρ
material_properties = Float64[   1      1     10_000    0.1;    # A1
                                 2      1     10_000    0.1;    # A2-A5
                                 3      1     10_000    0.1;    # A6-A9
                                 4      1     10_000    0.1;    # A10-A11
                                 5      1     10_000    0.1;    # A12-A13
                                 6      1     10_000    0.1;    # A14-A17
                                 7      1     10_000    0.1;    # A18-A21
                                 8      1     10_000    0.1];   # A22-A25

# load case 1
# forces information (Direction X=1, Y=2, Z=3)
#                      Node     Direction   Magnitude
forces_case_1 = Float64[ 1           2          20;
                         1           3          -5;
                         2           2         -20;
                         2           3          -5];

# load case 2
# forces information (Direction X=1, Y=2, Z=3)
#                       Node     Direction   Magnitude
forces_case_2 = Float64[  1          1          1.0;
                          1          2         10.0;
                          1          3         -5.0;
                          2          2         10.0;
                          2          3         -5.0;
                          3          1          0.5;
                          6          1          0.5];

# restricted degree of freedom (Direction X=1, Y=2, Z=3)
#                      Node   Direction
restricted_dof = Int64[  7        1;
                         7        2;
                         7        3;
                         8        1;
                         8        2;
                         8        3;
                         9        1;
                         9        2;
                         9        3;
                        10        1;
                        10        2;
                        10        3];

# minimun and maximum bounds for variables (areas of structure)
#                           Min       Max
bounds_variables = Float64[ 0.01      3.5;  # A1        [in^2]
                            0.01      3.5;  # A2-A5     [in^2]
                            0.01      3.5;  # A6-A9     [in^2]
                            0.01      3.5;  # A10-A11   [in^2]
                            0.01      3.5;  # A12-A13   [in^2]
                            0.01      3.5;  # A14-A17   [in^2]
                            0.01      3.5;  # A18-A21   [in^2]
                            0.01      3.5]; # A22-A25   [in^2]

# compressive and tensile limits in axial forces for each element
# (compresive -) and (tensile +)
#                                  Element   Tensile   Compresive
axial_stress_restricted = Float64[   1       40.0      -35.092;   # A1
                                     2       40.0      -11.590;   # A2
                                     3       40.0      -11.590;   # A3
                                     4       40.0      -11.590;   # A4
                                     5       40.0      -11.590;   # A5
                                     6       40.0      -17.305;   # A6
                                     7       40.0      -17.305;   # A7
                                     8       40.0      -17.305;   # A8
                                     9       40.0      -17.305;   # A9
                                    10       40.0      -35.092;   # A10
                                    11       40.0      -35.092;   # A11
                                    12       40.0      -35.092;   # A12
                                    13       40.0      -35.092;   # A13
                                    14       40.0       -6.759;   # A14
                                    15       40.0       -6.759;   # A15
                                    16       40.0       -6.759;   # A16
                                    17       40.0       -6.759;   # A17
                                    18       40.0       -6.959;   # A18
                                    19       40.0       -6.959;   # A19
                                    20       40.0       -6.959;   # A20
                                    21       40.0       -6.959;   # A21
                                    22       40.0      -11.082;   # A22
                                    23       40.0      -11.082;   # A23
                                    24       40.0      -11.082;   # A24
                                    25       40.0      -11.082];  # A25

# restricted displacement in the estructure (Direction X=1, Y=2, Z=3)
# all the directions except the dof restricted 
#                                 Node  Direction u_max     u_min
restricted_displacement = Float64[  1      1       0.35     -0.35;
                                    1      2       0.35     -0.35;
                                    1      3       0.35     -0.35;
                                    2      1       0.35     -0.35;
                                    2      2       0.35     -0.35;
                                    2      3       0.35     -0.35;
                                    3      1       0.35     -0.35;
                                    3      2       0.35     -0.35;
                                    3      3       0.35     -0.35;
                                    4      1       0.35     -0.35;
                                    4      2       0.35     -0.35;
                                    4      3       0.35     -0.35;
                                    5      1       0.35     -0.35;
                                    5      2       0.35     -0.35;
                                    5      3       0.35     -0.35;
                                    6      1       0.35     -0.35;
                                    6      2       0.35     -0.35;
                                    6      3       0.35     -0.35];

# connection matrix (local to global nodes), column 1 start node, column 2
# end node, each row indicates a elements
L2G_truss_25 = elements[:, 1:2];

# difference coordinates end and start node for each element
dif_x = xcor[L2G_truss_25[:, 2], 1] - xcor[L2G_truss_25[:, 1], 1];
dif_y = xcor[L2G_truss_25[:, 2], 2] - xcor[L2G_truss_25[:, 1], 2];
dif_z = xcor[L2G_truss_25[:, 2], 3] - xcor[L2G_truss_25[:, 1], 3];

# vector with the length of each element
L_e = sqrt.(dif_x.^2 + dif_y.^2 + dif_z.^2);

# model struct
mutable struct structural_model
    coordinates             ::Array{Float64,2} # the nodes coordinates
    info_elements           ::Array{Int64,2}   # elements information
    info_material_propieties::Array{Float64,2} # material properties information
    restricted_dof          ::Array{Int64,2}   # degree of freedom restricted
    number_load_cases       ::Int64            # number of loads cases
    load_cases              ::Array{Array{Float64,2},1} # the information of loads cases in an array
    variables_bounds        ::Array{Float64,2} # the maximum and minimun values of variables
    axial_sigma_restricted  ::Array{Float64,2} # the maximum and minimun values of axial stress
    constraint_displacement ::Array{Float64,2} # the maximum and minimun values of displacement and degree of freedom
    length_bars             ::Vector{Float64}  # length of each bar
end

# struct with the model information
truss_model_25_bars = structural_model(xcor,
                                       elements,
                                       material_properties,
                                       restricted_dof,
                                       2,
                                       [forces_case_1, forces_case_2],
                                       bounds_variables,
                                       axial_stress_restricted,
                                       restricted_displacement,
                                       L_e);
