# expintschurkrylov
**The source code of "An Exponential Integrator with Schur-Krylov Approximation to Accelerate Combustion Chemistry Computation"**

Developed by 

    Zaigang LIU
        Institute of Engineering Thermophysics, Chinese Academy of Sciences, Beijing 100190, China
        
    Jean-Louis CONSALVI     
        IUSTI UMR 7343, CNRS, Aix-Marseille Université, Marseille 13453, France
        
    Wenjun KONG             
        Institute of Engineering Thermophysics, Chinese Academy of Sciences, Beijing 100190, China

Please cite this reference:

    Zaigang Liu, Jean-L. Consalvi, Wenjun Kong,
    An Exponential Integrator with Schur–Krylov Approximation to accelerate combustion chemistry 
    computation,
    Combustion and Flame,
    Volume 203,
    2019,
    Pages 180-189,
    ISSN 0010-2180,
    https://doi.org/10.1016/j.combustflame.2019.01.031.
    (http://www.sciencedirect.com/science/article/pii/S0010218019300495)
    Abstract: The Exponential Integrator with Schur–Krylov Approximation (EISKA) algorithm was 
    developed for combustion applications. This algorithm combines the advantages of the explicit 
    large step advancement of the exponential schemes and the dimension reduction effect of the 
    Krylov subspace approximation, and was improved by introducing the Schur decomposition to control 
    the rounding error. The EISKA based on the SpeedCHEM (SC) package was implemented to simulate a 
    methane partially stirred reactor (PaSR) with pair-wise mixing model by considering the 
    mechanisms of Li et  al., GRI-Mech 3.0 and USC Mech II. Accuracy and computational efficiency 
    of EISKA are systematically compared with those of DVODE. In the case of the Li mechanism which 
    is a priori sufficiently small to be handled directly in combustion simulations, the 
    computations were accelerated by a factor of 1.99 without losing accuracy. In the cases of 
    GRI-Mech 3.0 and USC Mech II which are significantly larger than the Li mechanism, chemical 
    reduction methods, namely the Correlated Dynamic Adaptive Chemistry (CoDAC) and the 
    Multi-timescale (MTS) method were coupled with either DVODE or EISKA. The results show that 
    the EISKA is faster than DVODE either with or without chemical reduction methods. Model results 
    show that the best strategy is to use EISKA without any reduction method which leads to the 
    same accuracy as compared to DVODE and acceleration factors of 2.61 and 2.19 for GRI-Mech 3.0 
    and USC Mech II, respectively.
    Keywords: Chemistry computation acceleration; Exponential integrator; Krylov subspace; Mechanism reduction
