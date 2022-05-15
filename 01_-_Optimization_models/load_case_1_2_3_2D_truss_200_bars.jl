# Structural model unknowns for 2D truss, load case 1, 2 and 3 for truss with
# 200 bars
#
#
# Load conditions:
# * Load case 1 for truss with 200 bars
# * Load case 2 for truss with 200 bars
# * Load case 3 for truss with 200 bars
#-----------------------------------------------------
# Load case 1
#-----------------------------------------------------
# Load = 1.0 [kip]
# Direction: positive x-direction
# Nodes: 1, 6, 15, 20, 29, 34, 43, 48, 57, 62 and 71
#-----------------------------------------------------
# Load case 2
#-----------------------------------------------------
# Load = 10.0 [kips]
# Direction: negative y-direction
# Nodes: 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26,
# 28, 29, 30, 31, 32, 33, 34, 36, 38, 40, 42, 43, 44, 45, 46, 47, 48, 50, 52,
# 54, 56, 57, 58, 59, 60, 61, 62, 64, 66, 68, 70, 71, 72, 73, 74 and 75
#-----------------------------------------------------
# Load case 3
#-----------------------------------------------------
# load case (1) and (2) acting together
#-----------------------------------------------------
#
# Density of bar material: 0.283 [lb/in^3]
#
# Modulus of elasticity: 30_000 [ksi] (all bars)
#
# All areas bounds: 0.1 [in^2] <= A_i <= 20.0 [in^2]
#
# No limits displacement
#
# Stress limits: (compresive -) and (tensile +)
# Tensile limits:     all bars 10 [ksi]
# Compressive limits: all bars -10 [ksi]
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
# *(1) Awad, R. (2021, October). Sizing optimization of truss structures using
#      the political optimizer (PO) algorithm. In Structures
#      (Vol. 33, pp. 4871-4894). Elsevier.
#
# *(2) Degertekin, S. O., & Hayalioglu, M. S. (2013). Sizing truss structures 
#      using teaching-learning-based optimization.
#      Computers & Structures, 119, 177-188.
#
# -----------------------------------------------------------------------------

# generate the coordinates
#=
xcor_p = zeros((75, 3));
change_flod = 10;
x, y = [], [];
for i in collect(10.0:-1.0:0)*144.0 .+ 360.0
    global xcor_p, change_flod, x, y

    if Int32(change_flod%2) == 0
        
        X = [0.0, 1.0, 2.0, 3.0, 4.0]*240.0;
        Y = ones(5)*i;
        x = [x; X];
        y = [y; Y];
    end

    if Int32(change_flod%2) != 0
        X = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]*(240.0/2);
        Y = ones(9)*i;
        x = [x; X];
        y = [y; Y];
    end
    change_flod -= 1;
end

xcor_p[:, 1] = Float64.(x);
xcor_p[:, 2] = Float64.(y);
xcor_p[:, 3] = collect(1:75);

nodes_76_77 = [240.0   0.0   76;
               720.0   0.0   77];

xcor_p = [     xcor_p;
          nodes_76_77];

println("\n\n\nLocal displacements\n")
show(stdout, "text/plain", xcor_p)
println("\n")
=#

# node coordinates
#                   x      y
xcor =   Float64[  0.0  1800.0;  #   1
                 240.0  1800.0;  #   2
                 480.0  1800.0;  #   3
                 720.0  1800.0;  #   4
                 960.0  1800.0;  #   5
                   0.0  1656.0;  #   6
                 120.0  1656.0;  #   7
                 240.0  1656.0;  #   8
                 360.0  1656.0;  #   9
                 480.0  1656.0;  #  10
                 600.0  1656.0;  #  11
                 720.0  1656.0;  #  12
                 840.0  1656.0;  #  13
                 960.0  1656.0;  #  14
                   0.0  1512.0;  #  15
                 240.0  1512.0;  #  16
                 480.0  1512.0;  #  17
                 720.0  1512.0;  #  18
                 960.0  1512.0;  #  19
                   0.0  1368.0;  #  20
                 120.0  1368.0;  #  21
                 240.0  1368.0;  #  22
                 360.0  1368.0;  #  23
                 480.0  1368.0;  #  24
                 600.0  1368.0;  #  25
                 720.0  1368.0;  #  26
                 840.0  1368.0;  #  27
                 960.0  1368.0;  #  28
                   0.0  1224.0;  #  29
                 240.0  1224.0;  #  30
                 480.0  1224.0;  #  31
                 720.0  1224.0;  #  32
                 960.0  1224.0;  #  33
                   0.0  1080.0;  #  34
                 120.0  1080.0;  #  35
                 240.0  1080.0;  #  36
                 360.0  1080.0;  #  37
                 480.0  1080.0;  #  38
                 600.0  1080.0;  #  39
                 720.0  1080.0;  #  40
                 840.0  1080.0;  #  41
                 960.0  1080.0;  #  42
                   0.0   936.0;  #  43
                 240.0   936.0;  #  44
                 480.0   936.0;  #  45
                 720.0   936.0;  #  46
                 960.0   936.0;  #  47
                   0.0   792.0;  #  48
                 120.0   792.0;  #  49
                 240.0   792.0;  #  50
                 360.0   792.0;  #  51
                 480.0   792.0;  #  52
                 600.0   792.0;  #  53
                 720.0   792.0;  #  54
                 840.0   792.0;  #  55
                 960.0   792.0;  #  56
                   0.0   648.0;  #  57
                 240.0   648.0;  #  58
                 480.0   648.0;  #  59
                 720.0   648.0;  #  60
                 960.0   648.0;  #  61
                   0.0   504.0;  #  62
                 120.0   504.0;  #  63
                 240.0   504.0;  #  64
                 360.0   504.0;  #  65
                 480.0   504.0;  #  66
                 600.0   504.0;  #  67
                 720.0   504.0;  #  68
                 840.0   504.0;  #  69
                 960.0   504.0;  #  70
                   0.0   360.0;  #  71
                 240.0   360.0;  #  72
                 480.0   360.0;  #  73
                 720.0   360.0;  #  74
                 960.0   360.0;  #  75
                 240.0     0.0;  #  76
                 720.0     0.0]; #  77

# elements information
#                N_start   N_end   Type_material
elements = Int64[    1        2          1;       #   1
                     2        3          1;       #   2
                     3        4          1;       #   3
                     4        5          1;       #   4
                     1        6          2;       #   5
                     1        7          6;       #   6
                     2        7          6;       #   7
                     2        8          2;       #   8
                     2        9          6;       #   9
                     3        9          6;       #  10
                     3       10          2;       #  11
                     3       11          6;       #  12
                     4       11          6;       #  13
                     4       12          2;       #  14
                     4       13          6;       #  15
                     5       13          6;       #  16
                     5       14          2;       #  17
                     6        7          4;       #  18
                     7        8          3;       #  19
                     8        9          3;       #  20
                     9       10          3;       #  21
                    10       11          3;       #  22
                    11       12          3;       #  23
                    12       13          3;       #  24
                    13       14          4;       #  25
                     6       15          5;       #  26
                     7       15          6;       #  27
                     7       16          6;       #  28
                     8       16          5;       #  29
                     9       16          6;       #  30
                     9       17          6;       #  31
                    10       17          5;       #  32
                    11       17          6;       #  33
                    11       18          6;       #  34
                    12       18          5;       #  35
                    13       18          6;       #  36
                    13       19          6;       #  37
                    14       19          5;       #  38
                    15       16          7;       #  39
                    16       17          7;       #  40
                    17       18          7;       #  41
                    18       19          7;       #  42
                    15       20          8;       #  43
                    15       21         11;       #  44
                    16       21         11;       #  45
                    16       22          8;       #  46
                    16       23         11;       #  47
                    17       23         11;       #  48
                    17       24          8;       #  49
                    17       25         11;       #  50
                    18       25         11;       #  51
                    18       26          8;       #  52
                    18       27         11;       #  53
                    19       27         11;       #  54
                    19       28          8;       #  55
                    20       21          4;       #  56
                    21       22          9;       #  57
                    22       23          9;       #  58
                    23       24          9;       #  59
                    24       25          9;       #  60
                    25       26          9;       #  61
                    26       27          9;       #  62
                    27       28          4;       #  63
                    20       29         10;       #  64
                    21       29         11;       #  65
                    21       30         11;       #  66
                    22       30         10;       #  67
                    23       30         11;       #  68
                    23       31         11;       #  69
                    24       31         10;       #  70
                    25       31         11;       #  71
                    25       32         11;       #  72
                    26       32         10;       #  73
                    27       32         11;       #  74
                    27       33         11;       #  75
                    28       33         10;       #  76
                    29       30         12;       #  77
                    30       31         12;       #  78
                    31       32         12;       #  79
                    32       33         12;       #  80
                    29       34         13;       #  81
                    29       35         16;       #  82
                    30       35         16;       #  83
                    30       36         13;       #  84
                    30       37         16;       #  85
                    31       37         16;       #  86
                    31       38         13;       #  87
                    31       39         16;       #  88
                    32       39         16;       #  89
                    32       40         13;       #  90
                    32       41         16;       #  91
                    33       41         16;       #  92
                    33       42         13;       #  93
                    34       35          4;       #  94
                    35       36         14;       #  95
                    36       37         14;       #  96
                    37       38         14;       #  97
                    38       39         14;       #  98
                    39       40         14;       #  99
                    40       41         14;       # 100
                    41       42          4;       # 101
                    34       43         15;       # 102
                    35       43         16;       # 103
                    35       44         16;       # 104
                    36       44         15;       # 105
                    37       44         16;       # 106
                    37       45         16;       # 107
                    38       45         15;       # 108
                    39       45         16;       # 109
                    39       46         16;       # 110
                    40       46         15;       # 111
                    41       46         16;       # 112
                    41       47         16;       # 113
                    42       47         15;       # 114
                    43       44         17;       # 115
                    44       45         17;       # 116
                    45       46         17;       # 117
                    46       47         17;       # 118
                    43       48         18;       # 119
                    43       49         21;       # 120
                    44       49         21;       # 121
                    44       50         18;       # 122
                    44       51         21;       # 123
                    45       51         21;       # 124
                    45       52         18;       # 125
                    45       53         21;       # 126
                    46       53         21;       # 127
                    46       54         18;       # 128
                    46       55         21;       # 129
                    47       55         21;       # 130
                    47       56         18;       # 131
                    48       49          4;       # 132
                    49       50         19;       # 133
                    50       51         19;       # 134
                    51       52         19;       # 135
                    52       53         19;       # 136
                    53       54         19;       # 137
                    54       55         19;       # 138
                    55       56          4;       # 139
                    48       57         20;       # 140
                    49       57         21;       # 141
                    49       58         21;       # 142
                    50       58         20;       # 143
                    51       58         21;       # 144
                    51       59         21;       # 145
                    52       59         20;       # 146
                    53       59         21;       # 147
                    53       60         21;       # 148
                    54       60         20;       # 149
                    55       60         21;       # 150
                    55       61         21;       # 151
                    56       61         20;       # 152
                    57       58         22;       # 153
                    58       59         22;       # 154
                    59       60         22;       # 155
                    60       61         22;       # 156
                    57       62         23;       # 157
                    57       63         26;       # 158
                    58       63         26;       # 159
                    58       64         23;       # 160
                    58       65         26;       # 161
                    59       65         26;       # 162
                    59       66         23;       # 163
                    59       67         26;       # 164
                    60       67         26;       # 165
                    60       68         23;       # 166
                    60       69         26;       # 167
                    61       69         26;       # 168
                    61       70         23;       # 169
                    62       63          4;       # 170
                    63       64         24;       # 171
                    64       65         24;       # 172
                    65       66         24;       # 173
                    66       67         24;       # 174
                    67       68         24;       # 175
                    68       69         24;       # 176
                    69       70          4;       # 177
                    62       71         25;       # 178
                    63       71         26;       # 179
                    63       72         26;       # 180
                    64       72         25;       # 181
                    65       72         26;       # 182
                    65       73         26;       # 183
                    66       73         25;       # 184
                    67       73         26;       # 185
                    67       74         26;       # 186
                    68       74         25;       # 187
                    69       74         26;       # 188
                    69       75         26;       # 189
                    70       75         25;       # 190
                    71       72         27;       # 191
                    72       73         27;       # 192
                    73       74         27;       # 193
                    74       75         27;       # 194
                    71       76         28;       # 195
                    72       76         29;       # 196
                    73       76         28;       # 197
                    73       77         28;       # 198
                    74       77         29;       # 199
                    75       77         28];      # 200

# materials properties
#                              Type   Areas     E         ρ
material_properties = Float64[   1      1     30_000    0.283;
                                 2      1     30_000    0.283;
                                 3      1     30_000    0.283;
                                 4      1     30_000    0.283;
                                 5      1     30_000    0.283;
                                 6      1     30_000    0.283;
                                 7      1     30_000    0.283;
                                 8      1     30_000    0.283;
                                 9      1     30_000    0.283;
                                10      1     30_000    0.283;
                                11      1     30_000    0.283;
                                12      1     30_000    0.283;
                                13      1     30_000    0.283;
                                14      1     30_000    0.283;
                                15      1     30_000    0.283;
                                16      1     30_000    0.283;
                                17      1     30_000    0.283;
                                18      1     30_000    0.283;
                                19      1     30_000    0.283;
                                20      1     30_000    0.283;
                                21      1     30_000    0.283;
                                22      1     30_000    0.283;
                                23      1     30_000    0.283;
                                24      1     30_000    0.283;
                                25      1     30_000    0.283;
                                26      1     30_000    0.283;
                                27      1     30_000    0.283;
                                28      1     30_000    0.283;
                                29      1     30_000    0.283];

# load case 1
# forces information (Direction X=1, Y=2, Z=3)
#                       Node     Direction   Magnitude
forces_case_1 = Float64[  1           1         1.0;
                          6           1         1.0;
                         15           1         1.0;
                         20           1         1.0;
                         29           1         1.0;
                         34           1         1.0;
                         43           1         1.0;
                         48           1         1.0;
                         57           1         1.0;
                         62           1         1.0;
                         71           1         1.0];

# load case 2
# forces information (Direction X=1, Y=2, Z=3)
#                       Node     Direction   Magnitude
forces_case_2 = Float64[  1          2         -10.0;
                          2          2         -10.0;
                          3          2         -10.0;
                          4          2         -10.0;
                          5          2         -10.0;
                          6          2         -10.0;
                          8          2         -10.0;
                         10          2         -10.0;
                         12          2         -10.0;
                         14          2         -10.0;
                         15          2         -10.0;
                         16          2         -10.0;
                         17          2         -10.0;
                         18          2         -10.0;
                         19          2         -10.0;
                         20          2         -10.0;
                         22          2         -10.0;
                         24          2         -10.0;
                         26          2         -10.0;
                         28          2         -10.0;
                         29          2         -10.0;
                         30          2         -10.0;
                         31          2         -10.0;
                         32          2         -10.0;
                         33          2         -10.0;
                         34          2         -10.0;
                         36          2         -10.0;
                         38          2         -10.0;
                         40          2         -10.0;
                         42          2         -10.0;
                         43          2         -10.0;
                         44          2         -10.0;
                         45          2         -10.0;
                         46          2         -10.0;
                         47          2         -10.0;
                         48          2         -10.0;
                         50          2         -10.0;
                         52          2         -10.0;
                         54          2         -10.0;
                         56          2         -10.0;
                         57          2         -10.0;
                         58          2         -10.0;
                         59          2         -10.0;
                         60          2         -10.0;
                         61          2         -10.0;
                         62          2         -10.0;
                         64          2         -10.0;
                         66          2         -10.0;
                         68          2         -10.0;
                         70          2         -10.0;
                         71          2         -10.0;
                         72          2         -10.0;
                         73          2         -10.0;
                         74          2         -10.0;
                         75          2         -10.0];

# load case 3
# forces information
# the load case 1 and 2 acting togethers
forces_case_3 = [copy(forces_case_1);
                 copy(forces_case_2)];

# restricted degree of freedom (Direction X=1, Y=2, Z=3)
#                       Node   Direction
restricted_dof = Int64[  76        1;
                         76        2;
                         77        1;
                         77        2];

# minimun and maximum bounds for variables (areas of structure)
#                           Min       Max
bounds_variables = Float64[ 0.1      20.0;  # A1  [in^2]
                            0.1      20.0;  # A2  [in^2]
                            0.1      20.0;  # A3  [in^2]
                            0.1      20.0;  # A4  [in^2]
                            0.1      20.0;  # A5  [in^2]
                            0.1      20.0;  # A6  [in^2]
                            0.1      20.0;  # A7  [in^2]
                            0.1      20.0;  # A8  [in^2]
                            0.1      20.0;  # A9  [in^2]
                            0.1      20.0;  # A10  [in^2]
                            0.1      20.0;  # A11  [in^2]
                            0.1      20.0;  # A12  [in^2]
                            0.1      20.0;  # A13  [in^2]
                            0.1      20.0;  # A14  [in^2]
                            0.1      20.0;  # A15  [in^2]
                            0.1      20.0;  # A16  [in^2]
                            0.1      20.0;  # A17  [in^2]
                            0.1      20.0;  # A18  [in^2]
                            0.1      20.0;  # A19  [in^2]
                            0.1      20.0;  # A20  [in^2]
                            0.1      20.0;  # A21  [in^2]
                            0.1      20.0;  # A22  [in^2]
                            0.1      20.0;  # A23  [in^2]
                            0.1      20.0;  # A24  [in^2]
                            0.1      20.0;  # A25  [in^2]
                            0.1      20.0;  # A26  [in^2]
                            0.1      20.0;  # A27  [in^2]
                            0.1      20.0;  # A28  [in^2]
                            0.1      20.0]; # A29  [in^2]

# compressive and tensile limits in axial forces for each element
# (compresive -) and (tensile +)
elements_i =  collect(1:200);
Tensile    =  10.0*ones(200);
Compresive = -10.0*ones(200);

#                           Element   Tensile  Compresive
axial_stress_restricted = [elements_i Tensile Compresive];

# restricted displacement in the estructure (Direction X=1, Y=2, Z=3)
# no restricted
restricted_displacement = Array{Float64, 2}(undef, 2, 2);

# connection matrix (local to global nodes), column 1 strat node, column 2
# end node, each row indicates a elements
L2G_truss_200 = elements[:, 1:2];

# difference coordinates end and start node for each element
dif_x = xcor[L2G_truss_200[:, 2], 1] - xcor[L2G_truss_200[:, 1], 1];
dif_y = xcor[L2G_truss_200[:, 2], 2] - xcor[L2G_truss_200[:, 1], 2];

# vector with the length of each element
L_e = hypot.(dif_x, dif_y);

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
truss_model_200_bars = structural_model(xcor,
                                        elements,
                                        material_properties,
                                        restricted_dof,
                                        3,
                                        [forces_case_1, forces_case_2, forces_case_3],
                                        bounds_variables,
                                        axial_stress_restricted,
                                        restricted_displacement,
                                        L_e);
