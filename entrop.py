import numpy as np

def RB_Calc(r, min_dist, max_dist, rb_size):
    scaling = 1.0
    if (r < min_dist):
        exit("RadialBasis: r<min_dist. r = " + str(r) + ", min_dist = " + str(min_dist) + '\n')
    if (r > max_dist):
        exit("RadialBasis: r>MaxDist !!!. r = " + str(r) + ", min_dist = " + str(min_dist) + '\n')

    mult = 2.0 / (max_dist - min_dist)
    ksi = (2 * r - (min_dist + max_dist)) / (max_dist - min_dist)

    rb_vals = np.zeros((rb_size))
    rb_ders = np.zeros((rb_size))
    rb_vals[0] = scaling * (1 * (r - max_dist)*(r - max_dist))
    rb_ders[0] = scaling * (0 * (r - max_dist)*(r - max_dist) + 2 * (r - max_dist))
    rb_vals[1] = scaling * (ksi*(r - max_dist)*(r - max_dist))
    rb_ders[1] = scaling * (mult * (r - max_dist)*(r - max_dist) + 2 * ksi*(r - max_dist))
    i = 2
    while (i < rb_size):
        rb_vals[i] = 2 * ksi*rb_vals[i - 1] - rb_vals[i - 2]
        rb_ders[i] = 2 * (mult * rb_vals[i - 1] + ksi * rb_ders[i - 1]) - rb_ders[i - 2]
        i += 1
    return rb_vals, rb_ders


alpha_index_basic_8 = np.array([[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 2, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1], [0, 0, 2, 0], [0, 0, 1, 1], [0, 0, 0, 2], [1, 0, 0, 0]])
alpha_index_times_8 = np.array([[0, 0, 1, 11], [1, 1, 1, 12], [2, 2, 1, 12], [3, 3, 1, 12], [4, 4, 1, 13], [5, 5, 2, 13], [6, 6, 2, 13], [7, 7, 1, 13], [8, 8, 2, 13], [9, 9, 1, 13], [0, 10, 1, 14], [0, 11, 1, 15], [0, 12, 1, 16], [0, 15, 1, 17]])
alpha_moment_mapping_8 = np.array([0, 10, 11, 12, 13, 14, 15, 16, 17])

# Get Radial Basis vals and calculate
def calc_basis(js,  
    species_count=2, 
    radial_func_count=2,
    alpha_moments_count=18, 
    alpha_index_basic=alpha_index_basic_8, 
    alpha_index_basic_count=11, 
    alpha_index_times=alpha_index_times_8,
    alpha_index_times_count=14,
    alpha_moment_mapping=alpha_moment_mapping_8,
    alpha_scalar_moments=9,
    rb_size=2, 
    min_dist=2, 
    max_dist=5.5):

    moment_vals = np.zeros(alpha_moments_count)
    max_alpha_index_basic = 0
    for i in range(alpha_index_basic_count):
        max_alpha_index_basic = max(max_alpha_index_basic, alpha_index_basic[i][1] + alpha_index_basic[i][2] + alpha_index_basic[i][3])
    max_alpha_index_basic += 1
	# dist_powers_.resize(max_alpha_index_basic)
	# coords_powers_.resize(max_alpha_index_basic)

    # type_central = Neighborhood.my_type

    # if (type_central>=species_count):
        # exit("Too few species count in the MTP potential!")

    central_atom = np.array([0,0,0])
    for j in js:
        rijx = j[0] - central_atom[0]
        rijy = j[1] - central_atom[1]
        rijz = j[2] - central_atom[2]
        rij2 = rijx*rijx + rijy*rijy + rijz*rijz
        rt_rij2 = np.sqrt(rij2)

        rb_vals, rb_ders = RB_Calc(rt_rij2, min_dist, max_dist, rb_size)

        dist_powers_ = np.zeros((max_alpha_index_basic))
        dist_powers_[0] = 1
        coords_powers_ = np.zeros((max_alpha_index_basic,3))
        for a in range(3):
            coords_powers_[0][a] = 1
        for k in range(1, max_alpha_index_basic):
            dist_powers_[k] = dist_powers_[k - 1] *rt_rij2
            for a, vec in enumerate([rijx, rijy, rijz]):
                coords_powers_[k][a] = coords_powers_[k - 1][a] * vec

        for i in range(alpha_index_basic_count):
            # p_rb_vals = rb_vals[0]
            # How to include radial coeffs??
            val = 0.0
            for rb_ind in range(rb_size):
                val += rb_vals[rb_ind]
                # val += radial_coeffs(alpha_index_basic[i][0], rb_ind) * p_rb_vals[rb_ind]

            k = alpha_index_basic[i][1] + alpha_index_basic[i][2] + alpha_index_basic[i][3]
            powk = 1.0 / dist_powers_[k]
            val *= powk

            pow0 = coords_powers_[alpha_index_basic[i][1]][0]
            pow1 = coords_powers_[alpha_index_basic[i][2]][1]
            pow2 = coords_powers_[alpha_index_basic[i][3]][2]

            mult0 = pow0*pow1*pow2

            moment_vals[i] += val * mult0

    for i in range(alpha_index_times_count):
        val0 = moment_vals[alpha_index_times[i][0]]
        val1 = moment_vals[alpha_index_times[i][1]]
        val2 = alpha_index_times[i][2]
        moment_vals[alpha_index_times[i][3]] += val2 * val0 * val1

    basis_vals = np.zeros(alpha_scalar_moments + 1)
    basis_vals[0] = 1.0
    for i in range(alpha_scalar_moments):
        basis_vals[i+1] = moment_vals[alpha_moment_mapping[i]]

    return moment_vals, basis_vals


def calc_basis_ders(js,  
    species_count=2, 
    radial_func_count=2,
    alpha_moments_count=18, 
    alpha_index_basic=alpha_index_basic_8, 
    alpha_index_basic_count=11, 
    alpha_index_times=alpha_index_times_8,
    alpha_index_times_count=14,
    alpha_moment_mapping=alpha_moment_mapping_8,
    alpha_scalar_moments=9,
    rb_size=8, 
    min_dist=2, 
    max_dist=5.5):

    natoms = len(js) + 1    
    moment_ders = np.zeros((alpha_moments_count, natoms, 3))
    moment_vals = np.zeros(alpha_moments_count)
  
    # max_alpha_index_basic calculation
    max_alpha_index_basic = 0
    for i in range(alpha_index_basic_count):
        max_alpha_index_basic = max(max_alpha_index_basic, alpha_index_basic[i][1] + alpha_index_basic[i][2] + alpha_index_basic[i][3])
    max_alpha_index_basic += 1
    dist_powers_ = np.zeros(max_alpha_index_basic)
    coords_powers_ = np.zeros((max_alpha_index_basic, 3))
    
    central_atom = np.array([0,0,0])
    for j_count, j in enumerate(js):
        rijx = j[0] - central_atom[0]
        rijy = j[1] - central_atom[1]
        rijz = j[2] - central_atom[2]
        NeighbVect_j = [rijx, rijy, rijz]
        rij2 = rijx*rijx + rijy*rijy + rijz*rijz
        rt_rij2 = np.sqrt(rij2)

        # calculates vals and ders for j-th atom in the nbh
        rb_vals, rb_ders = RB_Calc(rt_rij2, min_dist, max_dist, rb_size)

        dist_powers_[0] = 1
        coords_powers_[0][0] = 1
        coords_powers_[0][1] = 1
        coords_powers_[0][2] = 1
        for k in range(1, max_alpha_index_basic):
            dist_powers_[k] = dist_powers_[k - 1] * rt_rij2
            for a in range(3):
                coords_powers_[k][a] = coords_powers_[k - 1][a] * NeighbVect_j[a]

        for i in range(alpha_index_basic_count):
            val = 0.0
            der = 0.0
            for rb_ind in range(rb_size):
                # print(alpha_index_basic[i][0], rb_ind)
                val += rb_vals[rb_ind]
                der += rb_ders[rb_ind]
            	# val += radial_coeffs(alpha_index_basic[i][0], rb_ind) * p_rb_vals[rb_ind]
                # der += radial_coeffs(alpha_index_basic[i][0], rb_ind) * p_rb_ders[rb_ind]

            k = alpha_index_basic[i][1] + alpha_index_basic[i][2] + alpha_index_basic[i][3]
            powk = 1.0 / dist_powers_[k]
            val *= powk
            der = der * powk - k * val / rt_rij2

            pow0 = coords_powers_[alpha_index_basic[i][1]][0]
            pow1 = coords_powers_[alpha_index_basic[i][2]][1]
            pow2 = coords_powers_[alpha_index_basic[i][3]][2]

            mult0 = pow0*pow1*pow2

            moment_vals[i] += val * mult0
            mult0 *= der / rt_rij2
            moment_ders[i][j_count][0] += mult0 * NeighbVect_j[0]
            moment_ders[i][j_count][1] += mult0 * NeighbVect_j[1]
            moment_ders[i][j_count][2] += mult0 * NeighbVect_j[2]

            if (alpha_index_basic[i][1] != 0):
                moment_ders[i][j_count][0] += val * alpha_index_basic[i][1] * coords_powers_[alpha_index_basic[i][1] - 1][0] * pow1 * pow2
            if (alpha_index_basic[i][2] != 0):
                moment_ders[i][j_count][1] += val * alpha_index_basic[i][2] * pow0 * coords_powers_[alpha_index_basic[i][2] - 1][1] * pow2
            if (alpha_index_basic[i][3] != 0):
                moment_ders[i][j_count][2] += val * alpha_index_basic[i][3] * pow0 * pow1 * coords_powers_[alpha_index_basic[i][3] - 1][2]

    for i in range(alpha_index_times_count):
        val0 = moment_vals[alpha_index_times[i][0]]
        val1 = moment_vals[alpha_index_times[i][1]]
        val2 = alpha_index_times[i][2]
        moment_vals[alpha_index_times[i][3]] += val2 * val0 * val1
        for j in range(natoms):
            for a in range(3):
                moment_ders[alpha_index_times[i][3]][j][a] += val2 * (moment_ders[alpha_index_times[i][0]][j][a] * val1 + val0 * moment_ders[alpha_index_times[i][1]][j][a])

#Next: copying all b_i corresponding to scalars into separate arrays,
#basis_vals and basis_ders
    basis_vals = np.zeros(alpha_scalar_moments + 1)
    basis_vals[0] = 1.0 # setting the constant basis function

    #memset(&basis_ders(0, 0, 0), 0, 3*nbh.count*sizeof(double))
    basis_ders = np.zeros((alpha_scalar_moments + 1, natoms, 3))

    for i in range(alpha_scalar_moments):
        # memcpy(&basis_ders(1+i, 0, 0), &moment_ders(alpha_moment_mapping[i], 0, 0), 3*nbh.count*sizeof(double))
        basis_vals[1 + i] = moment_vals[alpha_moment_mapping[i]]
        for j in range(natoms):
            for a in range(3):
                basis_ders[i+1][j][a] = moment_ders[i][j][a]
    
    return moment_vals, moment_ders, basis_vals, basis_ders
        


# Tests
# neigh1 = [[2,2,2], [1,2,-1.5], [3,-2,0]]
# neigh2 = [[-2,-2,-2], [-1,-2,1.5], [-3,2,0]]

# m1, b1 = calc_basis(neigh1)
# m2, b2 = calc_basis(neigh2)

# mv1, md1, bv1, bd1 = calc_basis_ders(neigh1)
# mv2, md2, bv2, bd2 = calc_basis_ders(neigh2)

# assert m1.all() == m2.all()
# assert m1.all() == mv1.all()
# assert m1.all() == mv2.all()
# assert b1.all() == b2.all()
# assert b1.all() == bv1.all()
# assert b1.all() == bv2.all()
# assert bd1.all() == bd2.all()


def plt_md():
    import matplotlib.pyplot as plt
    x = 2
    xs = []
    ys = []
    ds = []
    while x < 5:
        j = [[0, x, 0]]
        mv1, md1, bv1, bd1 = calc_basis_ders(j, rb_size=2)
        xs.append(x)
        # print(md1[:,:,1].sum())
        ys.append(md1[:,:,1].sum())
        ds.append(mv1.sum())
        x += 0.05

    print((ys[11]-ys[10])/(xs[11]-xs[10]))
    print(xs[11])
    plt.plot(xs, ys)
    plt.plot(xs, ds)
    plt.show()

# plt_md()

# netx = 0
# nety = 0
# netz = 0
# print(bd1.shape)
# for b in bd1:
#     for x, y, z in b:
#         netx += x
#         nety += y
#         netz += z

# print(netx, nety, netz)

# Solve for required basis sets for level 16
# Calculate and compress necessary Moments
# Tensor = (rij, angl1, angl2, angl3)

def plt_basis(min, max, size):
    import matplotlib.pyplot as plt
    c = 0 
    r = min
    while r < max:
        c += 1
        r += 0.01

    plt_vals = np.zeros((size, c))
    r = min
    c = 0
    rs = []
    tots = []
    while r < max:
        vals, ders = RB_Calc(r, min, max, size)
        tmp = 0
        for i in range(size):
            plt_vals[i,c] = vals[i]
            tmp += vals[i]
        rs.append(r)
        tots.append(tmp)
        r += 0.01
        c += 1
    

    # for i in range(size):
    #     plt.plot(rs, plt_vals[i])
    
    plt.plot(rs, tots)
    
    plt.show()

# plt_basis(2, 5, 2)