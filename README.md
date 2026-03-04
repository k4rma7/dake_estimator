# Security and size Estimates

This repository contains scripts used to evaluate the **communication costs** and **security estimates** for several **Double-KEM-based AKE constructions**. It also includes tools to analyze the impact of hints and heuristics on **Hint-(M)LWE security**.

The four scripts used to generate the main tables and figures in the associated paper are located in the size_security_scripts folder and are described below.


## Summary

| Script                     | Purpose                                              | Output        |
|----------------------------|------------------------------------------------------|---------------|
| `com_cost.py`              | AKE communication cost evaluation                    | Tables 2,3    |
| `Table_5/Maul.py`          | Kyber-like security estimation for Maul              | Table 5       |
| `maul_estimation.sage`     | Best-known attack analysis against Maul              | Table 6       |
| `graph_primal_attack.sage` | Validation of Heuristics 2--3 via best-known attacks | Figure 9      |


## size_security_scripts contents

### `com_cost.py`

To run it, execute:
```bash
python3 com_cost.py
```

This script **recomputes the communication costs** of different authenticated key exchange (AKE) constructions, based on:

- Maul  
- Double-Kyber  
- Twin-Kyber  

It provides a detailed comparison of the message sizes exchanged by each scheme and protocol.

**Used to generate:**
- Tables 2 and 3 of the submission.


---

### `Table_5/Maul.py`

To run it, execute:
```bash
python3 Maul.py
```
These scripts are modified versions of those used for Kyber security estimation in the repository https://github.com/pq-crystals/security-estimates. They provide the **decryption failure probability for Maul** and its **security estimates** (for both Primal and Dual attacks) via the MLWE instance resulting from the reduction described in Theorem 4 of our submission under Heuristics 1--3.

**Used to generate:**
- Table 5 of the submission.

---

### `maul_estimation.sage`

To run it, execute:
```bash
sage maul_estimation.sage
```

This script provides **security estimates based on the DDGR20 framework** for Maul using three distinct approaches. The goal is to verify that the reduction under Heuristics 2--3 is consistent with the best-known attack against Hint-LWE---specifically, the LWE with Side Information framework:

1. __Standard LWE (Baseline)__: We analyze the Hint-LWE as a standard LWE instance, disregarding the hints entirely. This naive approach serves as a baseline to demonstrate how incorporating hints reduces the security level.
2. __LWE with Side Information__: We study the LWE instance by treating the hints as approximate hints within the LWE with Side Information framework [DDGR20], encoded by the covariance matrix Σ_Hint.
3. __Post-Heuristic-Reduction LWE__: We analyze the LWE instance resulting from the reduction described in Theorem 4 of our submission under Heuristics 2--3.

This script relies on the **LWE with Side Information (leaky-LWE)** estimator from [DDGR20], available at:
https://github.com/lducas/leaky-LWE-Estimator

The corresponding files are located in the `leaky_estimator` folder.

We introduce a minor modification to the original framework to update the integration of approximate hints. This enables the use of **full-dimensional hint matrices** and improves **integration efficiency**.

The modified function is:
- `integrate_approx_hint_fulldim`

**Used to generate:**
- Table 6 of the submission.

---

### `graph_primal_attack.sage`

To run it, execute:
```bash
sage 
load("graph_primal_attack.sage")
```

This script provides **security estimates based on the DDGR20 framework** (similar to maul_estimation) across a wide **range of standard deviations** under the smoothing parameter.

It is used to:
- Track the **security of Hint-(M)LWE** instances.
- Check whether the reduction under **Heuristics 2--3 degrades security** compared to the **best-known attacks** within the [DDGR20] framework.

**Used to generate:**
- Figure 9 of the submission.

---


## Bibliography

[DDGR20] Dana Dachman-Soled, Léo Ducas, Huijing Gong, Mélissa Rossi: LWE with Side Information: Attacks and
Concrete Security Estimation. CRYPTO 2020.
https://eprint.iacr.org/2020/292
