# Structural model unknowns for 3D truss, load case 2 for truss with 25 bars
#
#
# keep in mind the units of unknows
# Areas:     [in^2]
# Forces:    [kips=1000lb_f]
# Distances: [in]
# Stress:    [ksi=1000psi=1000lb_f/in^2]
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
# * Li, H.-S., Au, S.-K (2010). Desing optimization using subset simulation
#   algorithm. Structural Safety, 32(6), 384-392.
#
# -----------------------------------------------------------------------------

# node coordinates
#                     x         y         z
xcor =   Float64[   -37.5      0.0       200.0;     #  1
                     37.5      0.0       200.0;     #  2
                    -37.5     37.5       100.0;     #  3
                     37.5     37.5       100.0;     #  4
                     37.5    -37.5       100.0;     #  5
                    -37.5    -37.5       100.0;     #  6
                   -100.0    100.0         0.0;     #  7
                    100.0    100.0         0.0;     #  8
                    100.0   -100.0         0.0;     #  9
                   -100.0   -100.0         0.0];    # 10

# elements information
#                N_start   N_end   Type_material
elements = Int64[    1       2          1;      #  1
                     1       4          2;      #  2
                     2       3          2;      #  3
                     1       5          2;      #  4
                     2       6          2;      #  5
                     2       4          3;      #  6
                     2       5          3;      #  7
                     1       3          3;      #  8
                     1       6          3;      #  9
                     3       6          4;      # 10
                     4       5          4;      # 11
                     3       4          5;      # 12
                     5       6          5;      # 13
                     3      10          6;      # 14
                     6       7          6;      # 15
                     4       9          6;      # 16
                     5       8          6;      # 17
                     4       7          7;      # 18
                     3       8          7;      # 19
                     5      10          7;      # 20
                     6       9          7;      # 21
                     6      10          8;      # 22
                     3       7          8;      # 23
                     4       8          8;      # 24
                     5       9          8];     # 25

# materials properties
#                              Type   Areas     E
material_properties = Float64[   1      1     10000;
                                 2      1     10000;
                                 3      1     10000;
                                 4      1     10000;
                                 5      1     10000;
                                 6      1     10000;
                                 7      1     10000;
                                 8      1     10000];

# forces information (Direction X=1, Y=2, Z=3)
#                Node     Direction   Magnitude
forces = Float64[  1          1          1.0;
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
