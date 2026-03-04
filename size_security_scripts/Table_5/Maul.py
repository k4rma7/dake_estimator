from math import log, sqrt, ceil
from Maul_failure import p2_cyclotomic_error_probability
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet
from proba_util import build_mod_switching_error_law

class MaulParameterSet:
    def __init__(self, n, m, ks, ke,  q, rqk, rqc, rq2, ke_ct=None):
        if ke_ct is None:
            ke_ct = 2 * n * ks**2
        self.n = n
        self.m = m
        self.ks = ks     # binomial parameter for the secret key
        self.ke = ke    # binomial parameter for the public key errors
        self.ke_ct = ke_ct    # binomial parameter for the ciphertext errors
        self.q = q
        self.rqk = rqk  # 2^(bits in the public key)
        self.rqc = rqc  # 2^(bits in the first ciphertext)
        self.rq2 = rq2  # 2^(bits in the second ciphertext)


def Maul_to_MLWE(kps):
    if kps.ks != kps.ke:
        raise "The security script does not handle different error parameter in secrets and errors (ks != ke) "

    # Check whether ciphertext error variance after rounding is larger than secret key error variance
    Rc = build_mod_switching_error_law(kps.q, kps.rqc)
    var_rounding = sum([i*i*Rc[i] for i in Rc.keys()])

    if kps.ke_ct/2. + var_rounding < kps.ke/2.:
        raise "The security of the ciphertext MLWE may not be weaker than the one of the public key MLWE"   

    return MLWEParameterSet(kps.n, kps.m, kps.m + 1, sqrt(kps.ks/4.), kps.q, "gaussian")


def communication_costs(ps):
    """ Compute the communication cost of a parameter set
    :param ps: Parameter set (ParameterSet)
    :returns: (cost_Alice, cost_Bob) (in Bytes)
    """
    A_space = ps.n * ps.m * log(ps.rqk)/log(2)
    B_space = ps.n * ps.m * log(ps.rqc)/log(2)
    C_space = ps.n * log(ps.rq2)/log(2)
    return (int(round(A_space))/8., int(round(B_space))/8., int(round(C_space))/8.)


def summarize(ps):
    print ("params: ", ps.__dict__)
    print ("com costs: ", communication_costs(ps))
    F, f = p2_cyclotomic_error_probability(ps)
    print ("failure: %.1f = 2^%.1f"%(f, log(f + 2.**(-300))/log(2)))


if __name__ == "__main__":
    # Maul parameter sets
    ps_light = MaulParameterSet(256, 2, 4, 4, 7681, 7681, 2**10, 2**4)
    ps_recommended = MaulParameterSet(256, 3, 4, 4, 7681, 7681, 2**11, 2**6)
    ps_paranoid = MaulParameterSet(256, 4, 4, 4, 9473, 9473, 2**12, 2**5)


    # Analyses
    print ("Maul512 (light):")
    print ("--------------------")
    print ("security:")
    MLWE_summarize_attacks(Maul_to_MLWE(ps_light))
    summarize(ps_light)
    print ()

    print ("Maul768 (recommended):")
    print ("--------------------")
    print ("security:")
    MLWE_summarize_attacks(Maul_to_MLWE(ps_recommended))
    summarize(ps_recommended)
    print ()

    print ("Maul1024 (paranoid):")
    print ("--------------------")
    print ("security:")
    MLWE_summarize_attacks(Maul_to_MLWE(ps_paranoid))
    summarize(ps_paranoid)
    print ()
