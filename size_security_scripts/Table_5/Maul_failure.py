import operator as op
#from math import factorial as fac
from math import sqrt, log, exp
import sys
from proba_util import *
from decimal import Decimal, getcontext

getcontext().prec = 300

def fac(n):
    r = 1.
    for i in range(1, n+1):
        r *= i
    return r

def p2_cyclotomic_final_error_distribution(ps):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    chis = build_centered_binomial_law(ps.ks)          # LWE secret law for sC, sL, and sR
    chie = build_centered_binomial_law(ps.ke)          # LWE error law for eC, eL, and eR
    #chif = build_centered_binomial_law(ps.ke_ct)        # LWE error law for fL and fR
    chif = build_gaussian_law(sqrt(ps.ke_ct/2))        # LWE error law for fL and fR

    Ru = build_mod_switching_error_law(ps.q, ps.rqc)    # Rounding error for cC
    Rv = build_mod_switching_error_law(ps.q, ps.rq2)    # rounding error for cL and cR

    chiRc = law_convolution(chie, Ru)                   # eC + ecC

    B1 = law_product(chis, chiRc)                       # (LWE+Rounding error) * LWE (as in a E*S product)
    B2 = law_product(chis, chie)

    C1 = iter_law_convolution(B1, ps.m * ps.n) # to check
    C2 = iter_law_convolution(B2, ps.m * ps.n) # to check

    C=law_convolution(C1, C2)

    F = law_convolution(chif, Rv)                   # LWE + rounding error cL and cR

    D = law_convolution(C, F)                           # Final error
    return D
    
    """
    chis = build_centered_binomial_law(ps.ks)           # LWE error law for the key
    chie = build_centered_binomial_law(ps.ke_ct)        # LWE error law for the ciphertext
    chie_pk = build_centered_binomial_law(ps.ke)
    Rk = build_mod_switching_error_law(ps.q, ps.rqk)    # Rounding error public key
    Rc = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext
    chiRs = law_convolution(chis, Rk)                   # LWE+Rounding error key
    chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext

    B1 = law_product(chie_pk, chiRs)                       # (LWE+Rounding error) * LWE (as in a E*S product)
    B2 = law_product(chis, chiRe)

    C1 = iter_law_convolution(B1, ps.m * ps.n)
    C2 = iter_law_convolution(B2, ps.m * ps.n)

    C=law_convolution(C1, C2)

    R2 = build_mod_switching_error_law(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error
    D = law_convolution(C, F)                           # Final error
    return D
    """

def p2_cyclotomic_error_probability(ps):
    F = p2_cyclotomic_final_error_distribution(ps)
    proba = tail_probability(F, ps.q/4)
    return F, 2.*ps.n*proba
