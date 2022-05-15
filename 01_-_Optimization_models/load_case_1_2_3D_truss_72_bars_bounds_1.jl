# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 72
# bars, for the areas bound (case 1) 0.1 [in^2] <= A_i <= 5.0 [in^2]
#
#
# Load conditions:
# * Load case 1 for truss with 72 bars
# * Load case 2 for truss with 72 bars
#-----------------------------------------------------
#  Node      Case(1) [kips]          Case(2) [kips]
#            Px     Py     Pz       Px     Py     Pz
#-----------------------------------------------------
#   17       5.0    5.0   -5.0      0.0    0.0   -5.0
#   18       0.0    0.0    0.0      0.0    0.0   -5.0
#   19       0.0    0.0    0.0      0.0    0.0   -5.0
#   20       0.0    0.0    0.0      0.0    0.0   -5.0
#-----------------------------------------------------
#
# Density of bar material: 0.1 [lb/in^3]
#
# Modulus of elasticity: 10_000 [ksi] (all bars)
#
# All areas bounds:
# --> (case 1) 0.1  [in^2] <= A_i <= 5.0 [in^2]
#
# Limits displacement:
#-----------------------------------------
#           Displacement limits [in]
#   Node     x [in]      y [in]     z [in]
#-----------------------------------------
#    17       +-0.25     +-0.25       -
#    18       +-0.25     +-0.25       -
#    19       +-0.25     +-0.25       -
#    20       +-0.25     +-0.25       -
#
# Stress limits: (compresive -) and (tensile +)
# Tensile limits: all bars 25 [ksi]
# Compressive limits: all bar -25 [ksi]
#
# Number of restrictions: 320
# * Displacement: 2 (+-) * 2 Direction * 4 Nodes * 2 Cases  = 32
# * Stress: 2 (Compressive and Tensile) * 72 Bars * 2 Cases = 288
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
# The reference [1] has some error in direction displacement constraints,
# see the reference [2]
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

# node coordinates
#                     x         y          z
xcor =   Float64[     0.0       0.0       0.0;  #  1
                    120.0       0.0       0.0;  #  2
                    120.0     120.0       0.0;  #  3
                      0.0     120.0       0.0;  #  4
                      0.0       0.0      60.0;  #  5
                    120.0       0.0      60.0;  #  6
                    120.0     120.0      60.0;  #  7
                      0.0     120.0      60.0;  #  8
                      0.0       0.0     120.0;  #  9
                    120.0       0.0     120.0;  # 10
                    120.0     120.0     120.0;  # 11
                      0.0     120.0     120.0;  # 12
                      0.0       0.0     180.0;  # 13
                    120.0       0.0     180.0;  # 14
                    120.0     120.0     180.0;  # 15
                      0.0     120.0     180.0;  # 16
                      0.0       0.0     240.0;  # 17
                    120.0       0.0     240.0;  # 18
                    120.0     120.0     240.0;  # 19
                      0.0     120.0     240.0]; # 20

# elements information
#                N_start   N_end   Type_material
elements = Int64[   1        5          1;    #    1
                    2        6          1;    #    2
                    3        7          1;    #    3
                    4        8          1;    #    4
                    1        6          2;    #    5
                    2        5          2;    #    6
                    2        7          2;    #    7
                    3        6          2;    #    8
                    3        8          2;    #    9
                    4        7          2;    #   10
                    4        5          2;    #   11
                    1        8          2;    #   12
                    5        6          3;    #   13
                    6        7          3;    #   14
                    7        8          3;    #   15
                    8        5          3;    #   16
                    5        7          4;    #   17
                    6        8          4;    #   18
                    5        9          5;    #   19
                    6       10          5;    #   20
                    7       11          5;    #   21
                    8       12          5;    #   22
                    5       10          6;    #   23
                    6        9          6;    #   24
                    6       11          6;    #   25
                    7       10          6;    #   26
                    7       12          6;    #   27
                    8       11          6;    #   28
                    8        9          6;    #   29
                    5       12          6;    #   30
                    9       10          7;    #   31
                   10       11          7;    #   32
                   11       12          7;    #   33
                   12        9          7;    #   34
                    9       11          8;    #   35
                   10       12          8;    #   36
                    9       13          9;    #   37
                   10       14          9;    #   38
                   11       15          9;    #   39
                   12       16          9;    #   40
                    9       14         10;    #   41
                   10       13         10;    #   42
                   10       15         10;    #   43
                   11       14         10;    #   44
                   11       16         10;    #   45
                   12       15         10;    #   46
                   12       13         10;    #   47
                    9       16         10;    #   48
                   13       14         11;    #   49
                   14       15         11;    #   50
                   15       16         11;    #   51
                   16       13         11;    #   52
                   13       15         12;    #   53
                   14       16         12;    #   54
                   13       17         13;    #   55
                   14       18         13;    #   56
                   15       19         13;    #   57
                   16       20         13;    #   58
                   13       18         14;    #   59
                   14       17         14;    #   60
                   14       19         14;    #   61
                   15       18         14;    #   62
                   15       20         14;    #   63
                   16       19         14;    #   64
                   16       17         14;    #   65
                   13       20         14;    #   66
                   17       18         15;    #   67
                   18       19         15;    #   68
                   19       20         15;    #   69
                   20       17         15;    #   70
                   17       19         16;    #   71
                   18       20         16];   #   72

# materials properties
#                               Type   Areas    E        ρ
material_properties = Float64[    1     2     10_000    0.1;    # A1-A4
                                  2     2     10_000    0.1;    # A5-A12
                                  3     2     10_000    0.1;    # A13-A16
                                  4     2     10_000    0.1;    # A17-A18
                                  5     2     10_000    0.1;    # A19-A22
                                  6     2     10_000    0.1;    # A23-A30
                                  7     2     10_000    0.1;    # A31-A34
                                  8     2     10_000    0.1;    # A35-A36
                                  9     2     10_000    0.1;    # A37-A40
                                 10     2     10_000    0.1;    # A41-A48
                                 11     2     10_000    0.1;    # A49-A52
                                 12     2     10_000    0.1;    # A53-A54
                                 13     2     10_000    0.1;    # A55-A58
                                 14     2     10_000    0.1;    # A59-A66
                                 15     2     10_000    0.1;    # A67-A70
                                 16     2     10_000    0.1];   # A71-A72

# load case 1
# forces information (Direction X=1, Y=2, Z=3)
#                       Node     Direction   Magnitude
forces_case_1 = Float64[ 17          1          5;
                         17          2          5;
                         17          3         -5];

# load case 2
# forces information (Direction X=1, Y=2, Z=3)
#                      Node     Direction   Magnitude
forces_case_2 = Float64[ 17          3         -5;
                         18          3         -5;
                         19          3         -5;
                         20          3         -5];

# restricted degree of freedom (Direction X=1, Y=2, Z=3)
#                      Node   Direction
restricted_dof = Int64[  1       1;
                         1       2;
                         1       3;
                         2       1;
                         2       2;
                         2       3;
                         3       1;
                         3       2;
                         3       3;
                         4       1;
                         4       2;
                         4       3];

# minimun and maximum bounds for variables (areas of structure)
#                           Min      Max
bounds_variables = Float64[ 0.1      5.0;    # A1-A4       [in^2]
                            0.1      5.0;    # A5-A12      [in^2]
                            0.1      5.0;    # A13-A16     [in^2]
                            0.1      5.0;    # A17-A18     [in^2]
                            0.1      5.0;    # A19-A22     [in^2]
                            0.1      5.0;    # A23-A30     [in^2]
                            0.1      5.0;    # A31-A34     [in^2]
                            0.1      5.0;    # A35-A36     [in^2]
                            0.1      5.0;    # A37-A40     [in^2]
                            0.1      5.0;    # A41-A48     [in^2]
                            0.1      5.0;    # A49-A52     [in^2]
                            0.1      5.0;    # A53-A54     [in^2]
                            0.1      5.0;    # A55-A58     [in^2]
                            0.1      5.0;    # A59-A66     [in^2]
                            0.1      5.0;    # A67-A70     [in^2]
                            0.1      5.0];   # A71-A72     [in^2]

# compressive and tensile limits in axial forces for each element
# (compresive -) and (tensile +)
#                                 Element   Tensile   Compresive
axial_stress_restricted = Float64[   1       25.0      -25.0;    #    A1
                                     2       25.0      -25.0;    #    A2
                                     3       25.0      -25.0;    #    A3
                                     4       25.0      -25.0;    #    A4
                                     5       25.0      -25.0;    #    A5
                                     6       25.0      -25.0;    #    A6
                                     7       25.0      -25.0;    #    A7
                                     8       25.0      -25.0;    #    A8
                                     9       25.0      -25.0;    #    A9
                                    10       25.0      -25.0;    #   A10
                                    11       25.0      -25.0;    #   A11
                                    12       25.0      -25.0;    #   A12
                                    13       25.0      -25.0;    #   A13
                                    14       25.0      -25.0;    #   A14
                                    15       25.0      -25.0;    #   A15
                                    16       25.0      -25.0;    #   A16
                                    17       25.0      -25.0;    #   A17
                                    18       25.0      -25.0;    #   A18
                                    19       25.0      -25.0;    #   A19
                                    20       25.0      -25.0;    #   A20
                                    21       25.0      -25.0;    #   A21
                                    22       25.0      -25.0;    #   A22
                                    23       25.0      -25.0;    #   A23
                                    24       25.0      -25.0;    #   A24
                                    25       25.0      -25.0;    #   A25
                                    26       25.0      -25.0;    #   A26
                                    27       25.0      -25.0;    #   A27
                                    28       25.0      -25.0;    #   A28
                                    29       25.0      -25.0;    #   A29
                                    30       25.0      -25.0;    #   A30
                                    31       25.0      -25.0;    #   A31
                                    32       25.0      -25.0;    #   A32
                                    33       25.0      -25.0;    #   A33
                                    34       25.0      -25.0;    #   A34
                                    35       25.0      -25.0;    #   A35
                                    36       25.0      -25.0;    #   A36
                                    37       25.0      -25.0;    #   A37
                                    38       25.0      -25.0;    #   A38
                                    39       25.0      -25.0;    #   A39
                                    40       25.0      -25.0;    #   A40
                                    41       25.0      -25.0;    #   A41
                                    42       25.0      -25.0;    #   A42
                                    43       25.0      -25.0;    #   A43
                                    44       25.0      -25.0;    #   A44
                                    45       25.0      -25.0;    #   A45
                                    46       25.0      -25.0;    #   A46
                                    47       25.0      -25.0;    #   A47
                                    48       25.0      -25.0;    #   A48
                                    49       25.0      -25.0;    #   A49
                                    50       25.0      -25.0;    #   A50
                                    51       25.0      -25.0;    #   A51
                                    52       25.0      -25.0;    #   A52
                                    53       25.0      -25.0;    #   A53
                                    54       25.0      -25.0;    #   A54
                                    55       25.0      -25.0;    #   A55
                                    56       25.0      -25.0;    #   A56
                                    57       25.0      -25.0;    #   A57
                                    58       25.0      -25.0;    #   A58
                                    59       25.0      -25.0;    #   A59
                                    60       25.0      -25.0;    #   A60
                                    61       25.0      -25.0;    #   A61
                                    62       25.0      -25.0;    #   A62
                                    63       25.0      -25.0;    #   A63
                                    64       25.0      -25.0;    #   A64
                                    65       25.0      -25.0;    #   A65
                                    66       25.0      -25.0;    #   A66
                                    67       25.0      -25.0;    #   A67
                                    68       25.0      -25.0;    #   A68
                                    69       25.0      -25.0;    #   A69
                                    70       25.0      -25.0;    #   A70
                                    71       25.0      -25.0;    #   A71
                                    72       25.0      -25.0];   #   A72

# restricted displacement in the estructure (Direction X=1, Y=2, Z=3)
# all the directions except the dof restricted 
#                                 Node  Direction   u_max     u_min
restricted_displacement = Float64[  17      1       0.25     -0.25;
                                    17      2       0.25     -0.25;
                                    18      1       0.25     -0.25;
                                    18      2       0.25     -0.25;
                                    19      1       0.25     -0.25;
                                    19      2       0.25     -0.25;
                                    20      1       0.25     -0.25;
                                    20      2       0.25     -0.25];

# connection matrix (local to global nodes), column 1 strat node, column 2
# end node, each row indicates a elements
L2G_truss_72_c1 = elements[:, 1:2];

# difference coordinates end and start node for each element
dif_x = xcor[L2G_truss_72_c1[:, 2], 1] - xcor[L2G_truss_72_c1[:, 1], 1];
dif_y = xcor[L2G_truss_72_c1[:, 2], 2] - xcor[L2G_truss_72_c1[:, 1], 2];
dif_z = xcor[L2G_truss_72_c1[:, 2], 3] - xcor[L2G_truss_72_c1[:, 1], 3];

# vector with the length of each element
L_e = hypot.(dif_x, dif_y, dif_z);

# model struct
mutable struct structural_model
    coordinates             ::Array{Float64,2} # the nodes coordinates
    info_elements           ::Array{Int64,2}   # elements information
    info_material_propieties::Array{Float64,2} # material properties information
    restricted_dof          ::Array{Int64,2}   # degree of freedom restricted
    number_load_cases       ::Int64            # number of loads cases
    load_cases              ::Array{Array{Float64,2},1} # the information of loads cases in a array
    variables_bounds        ::Array{Float64,2} # the maximum and minimun values of variables
    axial_sigma_restricted  ::Array{Float64,2} # the maximum and minimun values of axial stress
    constraint_displacement ::Array{Float64,2} # the maximum and minimun values of displacement and degree of freedom
    length_bars             ::Vector{Float64}  # length of each bar
end

# struct with the model information
truss_model_72_bars_bounds_1 = structural_model(xcor,
                                                elements,
                                                material_properties,
                                                restricted_dof,
                                                2,
                                                [forces_case_1, forces_case_2],
                                                bounds_variables,
                                                axial_stress_restricted,
                                                restricted_displacement,
                                                L_e);
