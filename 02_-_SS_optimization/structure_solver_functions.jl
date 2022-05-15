# Functions for solving trusses structures with FEM
# Truss in 2D and 3D
#
#
# keep in mind the units of unknows
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
# *(1) Análisis estructural - Jairo Uribe Escamilla - 2 edición
#
# *(2) Análisis matricial de estructuras (Curso con MATLAB) - Jorge Eduardo
#      Hurtado
#
# -----------------------------------------------------------------------------

# =========================== main function ===================================
"""
Trusses structure solver

a, q, N_e, sigma_e, epsilon_e, a_loc = solver_truss(xcor,
                                                    elements,
                                                    elemt_prop,
                                                    forces,
                                                    dof_restricted,
                                                    plots)

Parameters:

    xcor           (Array): array with the coordinates of the truss nodes
                            for 2D trusses only put the coordinates x and y,
                            for 3D trusses put the coordinates x, y and z
                            for 2D trusses --> for each node [x_i, y_i]
                            for 3D trusses --> for each node [x_i, y_i, z_i]
                            -- array type -- Array{Float64, 2}

    elements       (Array): array that contains the interconection nodes for
                            and the material type of each element
                            ---------------------------------------------------
                            start_node    end_node    material_type
                            [ N_start_1    N_end_1       Type_1;
                              N_start_2    N_end_2       Type_2;
                                 ...        ...           ...  ;
                              N_start_i    N_end_i       Type_i;
                                 ...        ...           ...  ]
                            -- array type -- Array{Int64, 2}

    element_prop   (Array): array that contains the materials properties for
                            each kind of material
                            ---------------------------------------------------
                             type        area      E_mod
                             [ Type_1    Area_1     E_1;
                               Type_2    Area_2     E_2;
                                 ...        ...     ...;
                               Type_i    Area_i     E_i;
                                 ...        ...     ...]
                            -- array type -- Array{Float64, 2}

    forces         (Array): array that contains the external forces acting
                            in the structure.
                            The direction (X = 1, Y = 2, Z = 3)
                            ---------------------------------------------------
                              node      direction      magnitude
                            [ Node_1     Dir_1           Mag_1;
                              Node_2     Dir_2           Mag_2;
                                ...         ...           ... ;
                              Node_i     Dir_i           Mag_i;
                                ...         ...           ... ]
                            -- array type -- Array{Float64, 2}

    dof_restricted (Array): array that contains the degrees of freedom
                            restricted.
                            The direction (X = 1, Y = 2, Z = 3)
                            ---------------------------------------------------
                              node      direction
                            [ Node_1     Dir_1 ;
                              Node_2     Dir_2 ;
                                ...       ...  ;
                              Node_i     Dir_i ;
                                ...       ...  ]
                            -- array type -- Array{Int64, 2}

    plots           (Bool): flag to indicate plots the structure

Returns:

    a           (Array):array with the global displacements of each degree of
                        freedom of the structure (sorted by dof)
                        -- array type -- Array{Float64, 1}

    N_e         (Array):array with the normal forces of each element of the
                        structure (sorted by element)
                        -- array type -- Array{Float64, 1}
                    
    sigma_e     (Array):array with the stress of each element of the
                        structure (sorted by element)
                        -- array type -- Array{Float64, 1}

    epsilon_e   (Array):array with the strains of each element of the
                        structure (sorted by element)
                        -- array type -- Array{Float64, 1}

    a_loc       (Array):array with the local displacements for each element of
                        the structure (sorted by element), each row is a
                        element, the columns mean:
                        for 2D trusses [u_i, v_i, u_j, v_j]
                        for 3D trusses [u_i, v_i, w_i, u_j, v_j, w_j]
                        -- array type -- Array{Float64, 2}
"""
function solver_truss(xcor              ::Array{Float64, 2},
                      elements          ::Array{Int64,   2},
                      element_prop      ::Array{Float64, 2},
                      forces            ::Array{Float64, 2},
                      dof_restricted    ::Array{Int64,   2})

# =============================== constanst ===================================

#   indexes for better reading and writing in the rest of the code
# -----------------------------------------------------------------------------

#   for coordinates
    X = 1; # x coordinate
    Y = 2; # y coordinate
    Z = 3; # z coordinate

#   for nodes
    N_start = 1; # element start node
    N_end   = 2; # element end node

#   for materials
    area_ind = 2; # area index
    E_ind    = 3; # elasticity modulus index
    type_ind = 3; # type of materials index

#   for forces
    node       = 1; # node index
    direction  = 2; # direction index
    magnitude  = 3; # magnitude index

#   for dof_restricted
    nodes_r     = 1; # for nodes restricted
    direction_r = 2; # for directions restricted

# ========================== information of structure =========================

    _, dimension    = size(xcor);       # dimension of structure
    num_elem, _     = size(elements);   # number of elements
    num_nodes, _    = size(xcor);       # number of nodes

#   reshape matrix xcor in 2D case
    if dimension == 2

        xcor = [xcor zeros(num_nodes)]; # z = 0
    end

# =========================== degrees of freedom ==============================

#   finds the degrees of freedom of each element according to dimension
#   (gdl array mapped to each element)
    dof     = dof_structure(dimension, num_nodes);

    num_dof = length(dof); # number of degrees of freedom

# =========================== unknowns of elements ============================

#   connection matrix (local to global nodes), column 1 start node, column 2
#   end node, each row indicates a elements
    L2G = elements[:, 1:2];

#   difference coordinates end and start node for each element
    dif_x = xcor[L2G[:, N_end], X] - xcor[L2G[:, N_start], X];
    dif_y = xcor[L2G[:, N_end], Y] - xcor[L2G[:, N_start], Y];
    dif_z = xcor[L2G[:, N_end], Z] - xcor[L2G[:, N_start], Z];

#   vector with the length of each element
    length_e = sqrt.(dif_x.^2 + dif_y.^2 + dif_z.^2);

# =========================== element properties ==============================

#   type of materials of elements
    type_e_materials = elements[:, type_ind];

#   elements areas
    areas_e          = element_prop[type_e_materials, area_ind];

#   elasticity modulus
    E_e              = element_prop[type_e_materials, E_ind];

# =============== the directing cosines of each element =======================

    Cx            = dif_x./length_e; # in X
    Cy            = dif_y./length_e; # in Y
    Cz            = dif_z./length_e; # in Z
    directing_cos = [Cx Cy Cz]; # array with directing cosines for each element

# ============================ forces vector ==================================

#   memory is separated for the vector of nodal forces
    f = zeros(num_dof);

#   nodes where forces act
    nodes_force      = Int64.(forces[:, node]);

#   direction of forces, for X = 1, for Y = 2, for Z = 3
    direction_forces = Int64.(forces[:, direction]);

#   magnitudes of forces
    magnitude_forces = forces[:, magnitude];

#   add the forces in the vector
    f[dof[CartesianIndex.(nodes_force, direction_forces)]] += magnitude_forces;

# =================== assemble the global stiffness matrix ====================

#   memory is separated for the global stiffness matrix
    K = zeros((num_dof, num_dof));

#   the matrix T and kloc save the transformation matrix and the local stiffness
#   matrix of each element respectively
    T    = zeros(Float64, (2 * dimension, 2 * dimension, num_elem));
    kloc = zeros(Float64, (2 * dimension, 2 * dimension, num_elem));

#   idx save the dof of each element
    idx = zeros(Int64, (num_elem, 2 * dimension));

#   the global stiffness matrix is assembled
    for e in 1:num_elem

#       the local and global stiffness matrix of the element e and its
#       transformation matrix is compute
        kloc_e, kglo_e, T_e = klocal_T_element(dimension, areas_e, E_e, length_e, directing_cos, e);

#       save the unknowns of elements
        T[:, :, e]    = T_e;        # transformation matrix
        kloc[:, :, e] = kloc_e;     # local stiffness matrix

#       index of dof
#                       dof star_node         dof end_node
        idx[e, :] = [dof[L2G[e, N_start], :]' dof[L2G[e, N_end], :]'];

#       sum the global stiffness matrix (total global matrix)
        K[idx[e, :], idx[e, :]] += kglo_e;
    end

# ========= known and unknown degrees of freedom are defined ==================

#   c : knowns
#   d : unknowns

    restricted_nodes     = dof_restricted[:, nodes_r];

    restricted_direction = dof_restricted[:, direction_r];

#   dof knowns
    c = dof[CartesianIndex.(restricted_nodes, restricted_direction)];

#   dof unknows
    d = setdiff(dof, c);

# ========================= solution of system of equations ===================

#   f = vector of equivalent nodal forces
#   q = vector of equilibrium nodal forces of the element
#   a = global displacements of all the dof

#   | qd |   | Kcc Kcd || ac |   | fd |  qc = 0 (always)
#   |    | = |         ||    | - |    |
#   | qc |   | Kdc Kdd || ad |   | fc |

#   memory for the arrays a and q
    a = zeros(num_dof);
    q = zeros(num_dof);

#   build the submatrices
    Kcc = K[c, c];
    Kcd = K[c, d];
    Kdc = K[d, c];
    Kdd = K[d, d];

#   build the subvectors
    fc = f[d];
    fd = f[c];
    ac = a[c];

#   solution of system of equations
    ad = Kdd\(fc - Kdc * ac); # (qc = 0)
    qd = Kcc * ac + Kcd * ad - fd;

#   integrate solution
    a[d] = ad;
    q[c] = qd;

# ====== the local internal forces of each element are calculated =============

#   internal forces and local displacements 
    int_forces, a_loc = local_unknowns(dimension, num_elem, T, idx, a, kloc);

#   for each element compute the internal forces and unknowns
    N_e         = int_forces[:, dimension + 1]; # normal forces, (+) tensile and
#                                                 (-) compression
    sigma_e     = N_e./areas_e; # stress
    epsilon_e   = sigma_e./E_e; # strain

#   returns the results
    return a, q, N_e, sigma_e, epsilon_e, a_loc
end


"""
Determinate the matrix of degrees of freedom of a structure

dof = dof_structure(dimension, num_nodes)

Parameters:

    dimension   (Int64):    integer with the dimension of structure
                            dimension = 2 ("2D"), dimension = 3 ("3D")

    num_nodes   (Int64):    integer with the number of nodes of structure

Returns:

    dof         (Array):    array with de degree of freedom for each node
                            in 2D structures dof: (num_nodes x 2)
                            displacement for X and Y direction for each node
                            in 3D structures dof: (num_nodes x 3)
                            displacement for X, Y and Z direction for each node
                            -- array type -- Array{Int64, 2}
"""
function dof_structure(dimension::Int64, num_nodes::Int64)

    if dimension == 2 # 2D structure

        num_dof = dimension * num_nodes; # degrees of freedom in X and Y

#       each row is a node, column 1 dof in X, column 2 dof in Y
        dof = reshape(collect(1:num_dof), (2, num_nodes))';
#       a simple calculation is followed to determine the degrees of freedom
#       for each node in 2D in trusses. i is the node number (start in 1).
#       --> in X direction dof: 2*i - 1
#       --> in Y direction dof: 2*i

    elseif dimension == 3 # 3D structure

        num_dof = dimension * num_nodes; # degrees of freedom in X, Y and Z

#       each row is a node, column 1 dof in X, column 2 dof in Y,
#       column 3 dof in Z
        dof = reshape(collect(1:num_dof), (3, num_nodes))';
#       a simple calculation is followed to determine the degrees of freedom
#       for each node in 3D in trusses. i is the node number (start in 1).
#       --> in X direction dof: 3*i - 2
#       --> in Y direction dof: 3*i - 1
#       --> in Z direction dof: 3*i
    else

#       invalid dimension
        println("Invalid option.")
        dof = nothing;
    end

    return dof # matrix of degree of freedom
end


# function that calculates the local and global stiffness matrix of each element
# and its corresponding transformation matrix
"""
Function that compute the local and global stiffness matrix of each element and
the matrix transformation

kloc_e, kglo_e, T_e = klocal_T_element(dimension, A, E, L, cosines_d, e)

Parameters:

    dimension   (Int64):    integer with the dimension of structure
                            dimension = 2 ("2D"), dimension = 3 ("3D")

    A           (Array):    array with the areas of elements
                            -- array type -- Array{Float64, 1}

    E           (Array):    array with the elasticity modulus of elements
                            -- array type -- Array{Float64, 1}

    L           (Array):    array with the lengths of elements
                            -- array type -- Array{Float64, 1}

    cosines_d   (Array):    array with the directing cosines of elements
                            cosines_d = [Cx_i,  Cy_i,  Cz_i] each row is a
                            element
                            -- array type -- Array{Float64, 2}

    e           (Int64):    current item number

Returns:

    kloc_e      (Array):    element local stiffness matrix
                            -- array type -- Array{Float64, 2}

    kglo_e      (Array):    element global stiffness matrix
                            -- array type -- Array{Float64, 2}

    T_e         (Array):    transformation matrix of element
                            -- array type -- Array{Float64, 2}
"""
function klocal_T_element(dimension     ::Int64,
                          A             ::Array{Float64, 1},
                          E             ::Array{Float64, 1},
                          L             ::Array{Float64, 1},
                          cosines_d     ::Array{Float64, 2},
                          e             ::Int64)

#   define director cosines
    Cx = cosines_d[e, 1]; # X direction
    Cy = cosines_d[e, 2]; # Y direction
    Cz = cosines_d[e, 3]; # Z direction

#   define element properties
    area_e  = A[e]; # element area
    L_e     = L[e]; # element length
    E_e     = E[e]; # element elasticity modulus

#   element stiffness
    k_e = (area_e * E_e)/L_e;

    if dimension == 2 # 2D structure

#       |N_i|   | Cx, Cy,   0,  0| |X_i| Cx, Cy: directing cosines
#       |V_i| = |-Cy, Cx,   0,  0| |Y_i|
#       |N_j|   |  0,  0,  Cx, Cy| |X_j|
#       |V_j|   |  0,  0, -Cy, Cx| |Y_j|

#       |u_i|   | Cx, Cy,   0,  0| |ug_i|
#       |v_i| = |-Cy, Cx,   0,  0| |vg_i|
#       |u_j|   |  0,  0,  Cx, Cy| |ug_j|
#       |v_j|   |  0,  0, -Cy, Cx| |vg_j|

#       N_i, N_j   : axial forces on the element at nodes i, j
#                    respectively (local axes)
#       V_i, V_j   : shear forces on the element at nodes i, j
#                    respectively (local axes)
#       X_i, X_j   : forces in x on the element at nodes i, j
#                    respectively (global axes)
#       Y_i, Y_j   : forces in y on the element at nodes i, j
#                    respectively (global axes)
#       u_i, u_j   : displacements in x in the element at nodes i, j
#                    respectively (local axes)
#       v_i, v_j   : displacements in y in the element at nodes i, j
#                    respectively (local axes)
#       ug_i, ug_j : displacements in x in the element at nodes i, j
#                    respectively (global axes)
#       vg_i, vg_j : displacements in y in the element at nodes i, j
#                    respectively (global axes)

#       T_e: transformation matrix of element

        T_e = [ Cx Cy   0  0;
               -Cy Cx   0  0;
                 0  0  Cx Cy;
                 0  0 -Cy Cx];

#                      ^ y
#                      |
#                      |_________________________
#       N_i, u_i ----> |_________________________|---->-----> x
#                                                 N_j, u_j

#       |N_i|         | 1,  0, -1,  0| |u_i|
#       |V_i| = K_e * | 0,  0,  0,  0| |v_i|  K_e = (E_e * A_e)/L_e
#       |N_j|         |-1,  0,  1,  0| |u_j|  V_i, V_j don't exist in truss
#       |V_j|         | 0,  0,  0,  0| |v_j|

#       N_i, N_j   : axial forces on the element at nodes i, j
#                    respectively (local axes)
#       V_i, V_j   : shear forces on the element at nodes i, j
#                    respectively (local axes)
#       u_i, u_j   : displacements in x in the element at nodes i, j
#                    respectively (local axes)
#       v_i, v_j   : displacements in y in the element at nodes i, j
#                    respectively (local axes)

#       kloc_e: stiffness matrix of element (local)

        kloc_e = k_e * [ 1  0 -1  0;
                         0  0  0  0;
                        -1  0  1  0;
                         0  0  0  0];

#       kglo_e: stiffness matrix of element (global)
#       kglo_e = T' * kloc_e * T

        kglo_e = k_e * [  Cx^2    Cx*Cy    -Cx^2   -Cx*Cy;
                         Cx*Cy     Cy^2   -Cx*Cy    -Cy^2;
                         -Cx^2   -Cx*Cy     Cx^2    Cx*Cy;
                        -Cx*Cy    -Cy^2    Cx*Cy     Cy^2];

    elseif dimension == 3 # 3D structure

#       only are necessary Cx, Cy and Cz
#       |N_i|   | Cx, Cy, Cz,  0,  0,  0| |X_i|
#       |V_i|   |  0,  0,  0,  0,  0,  0| |Y_i| Cx, Cy, Cz: directing cosines
#       |G_i| = |  0,  0,  0,  0,  0,  0| |Z_i|
#       |N_j|   |  0,  0,  0, Cx, Cy, Cz| |X_j|
#       |V_j|   |  0,  0,  0,  0,  0,  0| |Y_j|
#       |G_j|   |  0,  0,  0,  0,  0,  0| |Z_j|

#       |u_i|   | Cx, Cy, Cz,  0,  0,  0| |ug_i|
#       |v_i|   |  0,  0,  0,  0,  0,  0| |vg_i|
#       |w_i| = |  0,  0,  0,  0,  0,  0| |wg_i|
#       |u_j|   |  0,  0,  0, Cx, Cy, Cz| |ug_j|
#       |v_j|   |  0,  0,  0,  0,  0,  0| |vg_j|
#       |w_j|   |  0,  0,  0,  0,  0,  0| |wg_j|

#       N_i, N_j   : axial forces on the element at nodes i, j
#                    respectively (local axes)
#       V_i, V_j   : shear forces on the element at nodes i, j
#                    respectively (local axes)
#       G_i, G_j   : shear forces in the direction orthogonal to the element
#                    at nodes i, j respectively (local axis)
#       X_i, X_j   : forces in x on the element at nodes i, j
#                    respectively (global axes)
#       Y_i, Y_j   : forces in y on the element at nodes i, j
#                    respectively (global axes)
#       Z_i, Z_j   : forces in z on the element at nodes i, j
#                    respectively (global axes)
#       u_i, u_j   : displacements in x in the element at nodes i, j
#                    respectively (local axes)
#       v_i, v_j   : displacements in y in the element at nodes i, j
#                    respectively (local axes)
#       w_i, w_j   : displacements in z in the element at nodes i, j
#                    respectively (local axes)
#       ug_i, ug_j : displacements in x in the element at nodes i, j
#                    respectively (global axes)
#       vg_i, vg_j : displacements in y in the element at nodes i, j
#                    respectively (global axes)
#       wg_i, wg_j : displacements in x in the element at nodes i, j
#                    respectively (global axes)

#       T_e: transformation matrix of element

        T_e = [Cx  Cy  Cz    0   0   0;
                0   0   0    0   0   0;
                0   0   0    0   0   0;
                0   0   0   Cx  Cy  Cz;
                0   0   0    0   0   0;
                0   0   0    0   0   0];

#              ^ y
#              |
#            +-|-------------------------------------------+
#            |\|                                           |\
#            | \                                           | \
#            |  +------------------------------------------+--+
#            |  |                                          |  |
#      --->  |  |                                          | -|--->---> x
#  N_i, u_i  |  |                                          |  |  N_j, u_j
#            +--+------------------------------------------+  |
#             \ |\                                          \ |
#              \| \                                          \|
#               +--\------------------------------------------+
#                   \z

#       |N_i|         | 1,  0,  0, -1,  0,  0| |u_i| K_e = (E_e * A_e)/L_e
#       |V_i|         | 0,  0,  0,  0,  0,  0| |v_i| V_i, V_j don't exist in
#       |G_i| = K_e * | 0,  0,  0,  0,  0,  0| |w_i| truss
#       |N_j|         |-1,  0,  0,  1,  0,  0| |u_j| G_i, G_j don't exist in
#       |V_j|         | 0,  0,  0,  0,  0,  0| |v_j| truss
#       |G_j|         | 0,  0,  0,  0,  0,  0| |w_j|

#       N_i, N_j   : axial forces on the element at nodes i, j
#                    respectively (local axes)
#       V_i, V_j   : shear forces on the element at nodes i, j
#                    respectively (local axes)
#       G_i, G_j   : shear forces in the direction orthogonal to the element
#                    at nodes i, j respectively (local axis)
#       u_i, u_j   : displacements in x in the element at nodes i, j
#                    respectively (local axes)
#       v_i, v_j   : displacements in y in the element at nodes i, j
#                    respectively (local axes)
#       w_i, w_j   : displacements in z in the element at nodes i, j
#                    respectively (local axes)

#       kloc_e: stiffness matrix of element (local)

        kloc_e = k_e * [ 1   0   0  -1   0   0;
                         0   0   0   0   0   0;
                         0   0   0   0   0   0;
                        -1   0   0   1   0   0;
                         0   0   0   0   0   0;
                         0   0   0   0   0   0];

#       kglo_e: stiffness matrix of element (global)
#       kglo_e = T' * kloc_e * T

        kglo_e = k_e * [  Cx^2   Cx*Cy   Cx*Cz   -Cx^2  -Cx*Cy  -Cx*Cz;
                         Cx*Cy    Cy^2   Cy*Cz  -Cx*Cy   -Cy^2  -Cy*Cz;
                         Cx*Cz   Cy*Cz    Cz^2  -Cx*Cz  -Cy*Cz   -Cz^2;
                         -Cx^2  -Cx*Cy  -Cx*Cz    Cx^2   Cx*Cy   Cx*Cz;
                        -Cx*Cy   -Cy^2  -Cy*Cz   Cx*Cy    Cy^2   Cy*Cz;
                        -Cx*Cz  -Cy*Cz   -Cz^2   Cx*Cz   Cy*Cz    Cz^2];

    else # invalid input

        println("invalid input.")
        kloc_e = nothing;
        kglo_e = nothing;
        T_e    = nothing;

    end

    return kloc_e, kglo_e, T_e # local and global stiffness matrix of the
#                                element and element transformation matrix

end


"""
Function that compute the internal forces and the local displacements for each
element

int_forces, a_loc = local_unknowns(dimension, num_elements, T_e, idx, a, kloc_e)

Parameters:

    dimension   (Int64):    integer with the dimension of structure
                            dimension = 2 ("2D"), dimension = 3 ("3D")

    num_elements(Int64):    number of elements

    T_e         (Array):    set of transformation matrices for each element
                            -- array type -- Array{Float64, 3}

    idx         (Array):    array with the dof index for each element
                            -- array type -- Array{Int64, 2}

    a           (Array):    array with the global displacements for each dof
                            -- array type -- Array{Float64, 1}

    kloc_e      (Array):    array with the local stiffness matrices for each
                            element
                            -- array type -- Array{Float64, 3}

Returns:

    int_force   (Array):    array with the local forces for each element
                            -- array type -- Array{Float64, 1}

    a_loc       (Array):    local displacements for each element
                            for each element returns
                            [u_i, v_i, u_j, v_j]           --> 2D truss
                            [u_i, v_i, w_i, u_j, v_j, w_j] --> 3D truss
                            -- array type -- Array{Int64, 2}
"""
function local_unknowns(dimension       ::Int64,
                        num_elements    ::Int64,
                        T_e             ::Array{Float64, 3},
                        idx             ::Array{Int64,   2},
                        a               ::Array{Float64, 1},
                        kloc_e          ::Array{Float64, 3})

#   separate memory for internal forces and local displacements
    a_loc       = zeros((num_elements, 2*dimension));
    int_forces  = zeros((num_elements, 2*dimension));

#   for each element compute the interna forces and local displacements
    for e in 1:num_elements

#       local displacements for each element
        a_loc[e, :]      = T_e[:, :, e] * a[idx[e, :]];

#       internal forces for each element
        int_forces[e, :] = kloc_e[:, :, e] * a_loc[e, :];
    end

#   return the internal forces and the local displacements for each element
    return int_forces, a_loc

end
