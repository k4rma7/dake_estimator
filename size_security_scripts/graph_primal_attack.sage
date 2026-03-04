load("../leaky_estimator/proba_utils.sage")
load("../leaky_estimator/DBDD_predict_diag.sage")
load("../leaky_estimator/DBDD_predict.sage")
load("../leaky_estimator/DBDD.sage")
load("../leaky_estimator/DBDD_optimized.sage")
load("../leaky_estimator/ntru.sage")
import matplotlib.pyplot as plt


def negacyclic(r):
    a=reversed(list(r[1:]))
    c=[r[0]]+list(-vector(a))
    return matrix.toeplitz(c,r[1:])


def LWE(dbdd_class, n, q, m, D_e, D_s, A=None, diag=False, verbosity=1):
    """
    constructor that builds a DBDD instance from a LWE instance
    :n: (integer) size of the secret s
    :q: (integer) modulus
    :m: (integer) size of the error e
    :D_e: distribution of the error e (dictionnary form)
    :D_s: distribution of the secret s (dictionnary form)
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

n = 256
q = 7681
k = 2

sigma1 = 0.5
SIGMA = []
L_after_reduc = []
L_before = []
L_after_hints = []

D_ct = build_centered_binomial_law(4)
H = matrix(ZZ, negacyclic([draw_from_distribution(D_ct) for _ in range(n)]))
for i in range(2*k - 1):
    H = H.augment(negacyclic([draw_from_distribution(D_ct) for _ in range(n)]))

while(sigma1 < 6):
    sigma = sigma1/sqrt(2)
    D = build_Gaussian_law(sigma, int(ceil(32.*sigma)))
    b, dbdd = LWE(DBDD_predict, k*n, q, k*n, D, D, diag=False, verbosity=1)
    L_after_reduc.append(dbdd.estimate_attack()[0])

    D = build_Gaussian_law(sigma1, int(ceil(32.*sigma1)))
    b, dbdd = LWE(DBDD_predict, k*n, q, k*n, D, D, diag=False, verbosity=1)
    L_before.append(dbdd.estimate_attack()[0])

    sigma2 = sigma1 ** 2 * sqrt(n)
    D_y = build_Gaussian_law(sigma2, int(ceil(32.*sigma2)))
    y = vec([draw_from_distribution(D_y) for _ in range(n)])
    z = vec(dbdd.u[0][:-1]) * H.T + y

    covariance = sigma2 ** 2 * identity_matrix(ZZ, n)
    dbdd.integrate_approx_hint_fulldim(H, z, True, covariance, aposteriori=False)
    L_after_hints.append(dbdd.estimate_attack()[0])

    SIGMA.append(sigma1)
    sigma1 += 0.1


with open("list.csv", "w") as f:
    f.write("SIGMA,L_before,L_after_reduc,L_after_hints\n")
    for i in range(len(SIGMA)):
        f.write(f"{SIGMA[i]},{L_before[i]},{L_after_reduc[i]},{L_after_hints[i]}\n")


plt.plot(SIGMA, L_before, label=r'$\sigma_1$', color='blue')
plt.plot(SIGMA, L_after_reduc, label=r'$\sigma_1/\sqrt{2}$', color='red')
plt.plot(SIGMA, L_after_hints, label=r'$\Sigma_Hints$', color='green')

plt.xlabel(r'$\sigma_1$')
plt.ylabel('Bikz')
plt.legend()
plt.grid(True)

plt.show()