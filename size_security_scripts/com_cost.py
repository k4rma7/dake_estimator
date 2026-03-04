from math import log, sqrt
import numpy as np

cc_dillithium1 = (1312, 2420)
cc_dillithium3 = (1952, 3309)
cc_dillithium5 = (2592, 4627)

cc_kyber1 = (800, 768)
cc_kyber3 = (1184, 1088)
cc_kyber5 = (1568, 1568)


class ParameterSet:
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

def print_cost(cc, name):
    print(f"{name} Communication costs: pk: {cc[0]} Bytes, ct: {cc[1]} Bytes")

def maul_communication_costs(ps):
    """ Compute the communication cost of a parameter set
    :param ps: Parameter set (ParameterSet)
    :returns: (cost_Alice, cost_Bob) (in Bytes)
    """
    A_space = ps.n * ps.m * log(ps.rqk)/log(2)
    B_space = ps.n * ps.m * log(ps.rqc)/log(2)
    C_space = ps.n * log(ps.rq2)/log(2)
    return ( int(round(A_space))/8., int(round(B_space))/8. + 2 * int(round(C_space))/8.)

def kyber_communication_costs(ps):
    """ Compute the communication cost of a parameter set
    :param ps: Parameter set (ParameterSet)
    :returns: (cost_Alice, cost_Bob) (in Bytes)
    """
    A_space = ps.n * ps.m * log(ps.rqk)/log(2)
    B_space = ps.n * ps.m * log(ps.rqc)/log(2)
    C_space = ps.n * log(ps.rq2)/log(2)
    return ( int(round(A_space))/8., 2 * int(round(B_space))/8. + 2 * int(round(C_space))/8.)


if __name__ == "__main__":

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # Maul parameter sets
    maul_ps_NIST1 = ParameterSet(256, 2, 4, 4, 7681, 7681, 2**10, 2**4)
    maul_ps_NIST3 = ParameterSet(256, 3, 4, 4, 7681, 7681, 2**11, 2**6)
    maul_ps_NIST5 = ParameterSet(256, 4, 4, 4, 9473, 9473, 2**12, 2**5)

    maul_cc_NIST1 = maul_communication_costs(maul_ps_NIST1)
    maul_cc_NIST3 = maul_communication_costs(maul_ps_NIST3)
    maul_cc_NIST5 = maul_communication_costs(maul_ps_NIST5)

    print_cost(maul_cc_NIST1, "Maul (NIST1)")
    print_cost(maul_cc_NIST3, "Maul (NIST3)")
    print_cost(maul_cc_NIST5, "Maul (NIST5)")

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # Double Kyber parameter sets for comparison
    kyber_ps_NIST1 = ParameterSet(256, 2, 3, 3, 3329, 2**12, 2**10, 2**4, ke_ct=2)
    kyber_ps_NIST3 = ParameterSet(256, 3, 2, 2, 3329, 2**12, 2**10, 2**4, ke_ct=2)
    kyber_ps_NIST5 = ParameterSet(256, 4, 2, 2, 3329, 2**12, 2**11, 2**5, ke_ct=2)
    kyber_cc_NIST1 = kyber_communication_costs(kyber_ps_NIST1)
    kyber_cc_NIST3 = kyber_communication_costs(kyber_ps_NIST3)
    kyber_cc_NIST5 = kyber_communication_costs(kyber_ps_NIST5)

    print_cost(kyber_cc_NIST1, "Double Kyber (NIST1)")
    print_cost(kyber_cc_NIST3, "Double Kyber (NIST3)")
    print_cost(kyber_cc_NIST5, "Double Kyber (NIST5)")

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    twinkyber_cc = (1568, 2464) # From TwinKyber paper
    print_cost(twinkyber_cc, "Twinkyber (NIST3)")

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # DAKE communication cost

    tag = (32, 48, 64) # 2* lambda / 8 correspond also to identity and Hash

    #DAKE

    #message 1 : id + pk_Maul + ct_Kyber

    DAKE_message1_size = [0, 0, 0]
    DAKE_message1_size[0] = tag[0] + maul_cc_NIST1[0] + cc_kyber1[1]
    DAKE_message1_size[1] = tag[1] + maul_cc_NIST3[0] + cc_kyber3[1]
    DAKE_message1_size[2] = tag[2] + maul_cc_NIST5[0] + cc_kyber5[1]
    print(f"DAKE Message 1 sizes (NIST1, NIST3, NIST5): {DAKE_message1_size} Bytes")

    #message 2 : ct_Maul + tag

    DAKE_message2_size = [0, 0, 0]
    DAKE_message2_size[0] = maul_cc_NIST1[1] + tag[0]
    DAKE_message2_size[1] = maul_cc_NIST3[1] + tag[1]
    DAKE_message2_size[2] = maul_cc_NIST5[1] + tag[2]

    print(f"DAKE Message 2 sizes (NIST1, NIST3, NIST5): {DAKE_message2_size} Bytes")

    #message 3 : tag
    DAKE_message3_size = tag

    #Total communication cost
    DAKE_total_size =  [DAKE_message1_size[i] + DAKE_message2_size[i] + DAKE_message3_size[i] for i in range(3)]
    print(f"DAKE total sizes (NIST1, NIST3, NIST5): {DAKE_total_size} Bytes")


    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    #DAKESign
    # message 1 : id + pk_Maul
    DAKESign_message1_size = [0, 0, 0]
    DAKESign_message1_size[0] = tag[0] + maul_cc_NIST1[0]
    DAKESign_message1_size[1] = tag[1] + maul_cc_NIST3[0]
    DAKESign_message1_size[2] = tag[2] + maul_cc_NIST5[0]
    print(f"DAKESign Message 1 sizes (NIST1, NIST3, NIST5): {DAKESign_message1_size} Bytes")

    #message 2 : ct_Maul + signature Dillithium
    DAKESign_message2_size = [0, 0, 0]
    DAKESign_message2_size[0] = maul_cc_NIST1[1] + cc_dillithium1[1]
    DAKESign_message2_size[1] = maul_cc_NIST3[1] + cc_dillithium3[1]
    DAKESign_message2_size[2] = maul_cc_NIST5[1] + cc_dillithium5[1] 

    print(f"DAKESign Message 2 sizes (NIST1, NIST3, NIST5): {DAKESign_message2_size} Bytes")

    #message 3 tag
    DAKESign_message3_size = tag
    #Total communication cost
    DAKESign_total_size =  [DAKESign_message1_size[i] + DAKESign_message2_size[i] + DAKESign_message3_size[i] for i in range(3)]
    print(f"DAKESign total sizes (NIST1, NIST3, NIST5): {DAKESign_total_size} Bytes")


    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    #DUAKE

    # message 1 = id + pk_Maul
    DUAKE_message1_size = [0, 0, 0]
    DUAKE_message1_size[0] = tag[0] + maul_cc_NIST1[0]
    DUAKE_message1_size[1] = tag[1] + maul_cc_NIST3[0]
    DUAKE_message1_size[2] = tag[2] + maul_cc_NIST5[0]
    print(f"DUAKE Message 1 sizes (NIST1, NIST3, NIST5): {DUAKE_message1_size} Bytes")

    # message 2 = ct_maul + tag(hash)

    DUAKE_message2_size = [0, 0, 0]
    DUAKE_message2_size[0] = maul_cc_NIST1[1] + tag[0]
    DUAKE_message2_size[1] = maul_cc_NIST3[1] + tag[1]
    DUAKE_message2_size[2] = maul_cc_NIST5[1] + tag[2]
    print(f"DUAKE Message 2 sizes (NIST1, NIST3, NIST5): {DUAKE_message2_size} Bytes")

    # message 3 = tag
    DUAKE_message3_size = tag
    #Total communication cost
    DUAKE_total_size =  [DUAKE_message1_size[i] + DUAKE_message2_size[i] + DUAKE_message3_size[i] for i in range(3)]
    print(f"DUAKE total sizes (NIST1, NIST3, NIST5): {DUAKE_total_size} Bytes")


    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # ML-KEM based AKE

    # total : 1 pk_Kyber + 3 ct_Kyber
    MLKEM_total_size = [0, 0, 0]
    MLKEM_total_size[0] = cc_kyber1[0] + 3 * cc_kyber1[1]
    MLKEM_total_size[1] = cc_kyber3[0] + 3 * cc_kyber3[1]
    MLKEM_total_size[2] = cc_kyber5[0] + 3 * cc_kyber5[1]
    print(f"ML-KEM based AKE total sizes (NIST1, NIST3, NIST5): {MLKEM_total_size} Bytes")
    print(f"Comparaison DAKE and ML-KEM AKE (NIST1, NIST3, NIST5): {[(MLKEM_total_size[i]-DAKE_total_size[i])/MLKEM_total_size[i] for i in range(3)]}")


    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # ML-KEM based UAKE

    # total : 1 pk_Kyber + 2 ct_Kyber
    MLKEM_U_total_size = [0, 0, 0]
    MLKEM_U_total_size[0] = cc_kyber1[0] + 2 * cc_kyber1[1]
    MLKEM_U_total_size[1] = cc_kyber3[0] + 2 * cc_kyber3[1]
    MLKEM_U_total_size[2] = cc_kyber5[0] + 2 * cc_kyber5[1]
    print(f"ML-KEM based UAKE total sizes (NIST1, NIST3, NIST5): {MLKEM_U_total_size} Bytes")
    print(f"Comparaison DUAKE and ML-KEM UAKE (NIST1, NIST3, NIST5): {[(MLKEM_U_total_size[i]-DUAKE_total_size[i])/MLKEM_U_total_size[i] for i in range(3)]}")

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # AKE of Xu et al.
    # total : 1 pk_twinkyber + 2 ct_twinkyber
    Twinkyber_AKE_total_size = twinkyber_cc[0] + 2 * twinkyber_cc[1]

    print(f"Twinkyber AKE total sizes (NIST5): {Twinkyber_AKE_total_size} Bytes")
    print(f"Comparaison DAKE and Twinkyber AKE (NIST5): {(Twinkyber_AKE_total_size-DAKE_total_size[2])/DAKE_total_size[2]}")

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"Gain in average DAKE: {np.mean([(MLKEM_total_size[i]-DAKE_total_size[i])/MLKEM_total_size[i] for i in range(3)] + [(Twinkyber_AKE_total_size-DAKE_total_size[2])/DAKE_total_size[2]])*100:.2f} %")
    print(f"Gain in average DUAKE: {np.mean([(MLKEM_U_total_size[i]-DUAKE_total_size[i])/MLKEM_U_total_size[i] for i in range(3)])*100:.2f} %")