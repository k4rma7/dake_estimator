load("../leaky_estimator/proba_utils.sage")
load("../leaky_estimator/DBDD_predict_diag.sage")
load("../leaky_estimator/DBDD_predict.sage")
load("../leaky_estimator/DBDD.sage")
load("../leaky_estimator/DBDD_optimized.sage")
load("../leaky_estimator/ntru.sage")


n = 256

# Maul Parameters

Maul_512 = {
"q"     : 7681,
"k"     : 2,
"eta"     : 4,
}

Maul_768 = {
"q"     : 7681,
"k"     : 3,
"eta"     : 4,
}

Maul_1024 = {
"q"     : 9473,
"k"     : 4,
"eta"     : 4,
}

def negacyclic(r):
    a=reversed(list(r[1:]))
    c=[r[0]]+list(-vector(a))
    return matrix.toeplitz(c,r[1:])

def LWE_to_DBDD(dbdd_class, n, q, m, D_e, D_s, A=None, diag=False, verbosity=1):
    """
    constructor that builds a DBDD instance from a LWE instance
    :n: (integer) size of the secret s
    :q: (integer) modulus
    :m: (integer) size of the error e
    :D_e: distribution of the error e (dictionnary form)
    :D_s: distribution of the secret s (dictionnary form)
    return dbdd instance
    """
    if verbosity:
        logging("     Build DBDD from LWE     ", style="HEADER")
        logging("n=%3d \t m=%3d \t q=%d" % (n, m, q), style="VALUE")
    # define the mean and sigma of the instance
    if (A == None):
        A = matrix([[randint(0, q) for _ in range(n)] for _ in range(m)])
    mu_e, s_e = average_variance(D_e)
    mu_s, s_s = average_variance(D_s)
    mu = vec(m * [mu_e] + n * [mu_s] + [1])
    S = diagonal_matrix(m * [s_e] + n * [s_s] + [0])
    # draw matrix A and define the lattice
    B = build_LWE_lattice(-A, q) # primal
    D = build_LWE_lattice(A/q, 1/q) # dual
    # draw the secrets
    s = vec([draw_from_distribution(D_s) for _ in range(n)])
    e = vec([draw_from_distribution(D_e) for _ in range(m)])
    # compute the public value t and build a target
    b = (s * A.T + e) % q
    tar = concatenate([b, [0] * n])
    B = kannan_embedding(B, tar)
    D = kannan_embedding(D, concatenate([-b/q, [0] * n])).T
    u = concatenate([e, s, [1]])
    return  b, dbdd_class(B, S, mu, u, verbosity=verbosity, D=D, Bvol=m*log(q))

def maul_security_estimation(Maul, n, verbo=0):
    """
    Estimate the security level of Maul parameter sets before reduction, with hint, and after reduction
    with the estimation of the associate DBDD problem
    :Maul: (dictionary) Maul parameter set
    :n: (integer) dimension parameter
    returns: (beta_ori, beta_hint, beta_reduc) security levels in bikz
    """
    sigma1 = sqrt(Maul["eta"] / 2)

    # Build H
    D_ct = build_centered_binomial_law(Maul["eta"])
    H = matrix(ZZ, negacyclic([draw_from_distribution(D_ct) for _ in range(n)]))
    for _ in range(2 * Maul["k"] - 1):
        H = H.augment(negacyclic([draw_from_distribution(D_ct) for _ in range(n)]))

    # ---- Post-Heuristic-Reduction LWE ----
    sigma = sigma1 / sqrt(2)
    D = build_Gaussian_law(sigma, int(ceil(32. * sigma)))

    _, dbdd = LWE_to_DBDD( DBDD_predict, Maul["k"] * n, Maul["q"], Maul["k"] * n, D, D, diag=False, verbosity=verbo)
    beta_reduc, _ = dbdd.estimate_attack()

    # ---- Standard LWE ----
    D = build_Gaussian_law(sigma1, int(ceil(32. * sigma1)))

    _, dbdd = LWE_to_DBDD( DBDD_predict, Maul["k"] * n, Maul["q"], Maul["k"] * n, D, D, diag=False, verbosity=verbo)
    beta_ori, _ = dbdd.estimate_attack()

    # ---- LWE with Side Information ----
    sigma2 = sigma1 ** 2 * sqrt(n)
    D_y = build_Gaussian_law(sigma2, int(ceil(32. * sigma2)))

    y = vec([draw_from_distribution(D_y) for _ in range(n)])
    z = vec(dbdd.u[0][:-1]) * H.T + y

    covariance = sigma2 ** 2 * identity_matrix(ZZ, n)
    dbdd.integrate_approx_hint_fulldim(H, z, True, covariance, aposteriori=False)
    beta_hint, _ = dbdd.estimate_attack()

    return beta_ori, beta_hint, beta_reduc

if __name__ == "__main__":

    print("Starting security estimation for Maul512.")
    beta512 = maul_security_estimation(Maul_512,  n)
    print(f"Results (Maul512): " f"before reduction = {beta512[0]} bikz, "f"with hints = {beta512[1]} bikz, "f"after reduction = {beta512[2]} bikz")
    print("Starting security estimation for Maul768.")
    beta768 = maul_security_estimation(Maul_768, n)
    print(f"Results (Maul768): " f"before reduction = {beta768[0]} bikz, "f"with hints = {beta768[1]} bikz, "f"after reduction = {beta768[2]} bikz")
    print("Starting security estimation for Maul1024.")
    beta1024 = maul_security_estimation(Maul_1024, n)
    print(f"Results (Maul1024): " f"before reduction = {beta1024[0]} bikz, "f"with hints = {beta1024[1]} bikz, "f"after reduction = {beta1024[2]} bikz")