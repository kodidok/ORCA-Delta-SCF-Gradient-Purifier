# ORCA Delta-SCF Gradient Purifier

The deltascf_wrapper.py enables geometry optimization of minima, TS and S0/S1 Minimum Energy Crossing Points (MECI) using the ΔSCF method with gradients corrected for spin-contamination. 

The script acts as an external optimizer interface for ORCA 6.1.0, calculating corrected gradients according to the formula by Ziegler et al. (https://doi.org/10.1007/BF00551551):

```math
\nabla E_{pure} =2 \nabla E_{BS} -\nabla E_{T}
```

The script reads the files basename_EXT.extinp.tmp and basename.xyz provided by ORCA called with the `ExtOpt` option, invokes S1 and T1 calculations, computes purified S1 gradients and feeds them back to the ORCA internal optimizer.

As described in the ORCA 6.1.0 documentation (see ORCA manual), you must specify the location of the wrapper script using one of the following methods:
1. As a file or link named otool_external in the same directory as the ORCA executables.
2. By assigning the EXTOPTEXE environment variable to the full path of the external program.
3. Via the ORCA input:

```
%method
  ProgExt "/full/path/to/script"
  Ext_Params "optional command line arguments"
end
```

The deltascf_wrapper.py runs gradient calculations in the directory ./0 at the location of the script. Make sure you created this directory prior to running the calculation.

##  Acknowledgments

This project is a fork and modification of the [ExtOptORCA-OpenQP](https://github.com/CrespiLab/ExtOptORCA-OpenQP) interface developed by **CrespiLab**. 

The original codebase was designed for MRSF-TDDFT gradients using OpenQP. It has been extensively modified to compute spin-contamination corrected gradients for ΔSCF calculations within ORCA 6.1.0. 


---

##  Features
* **Minimum and TS Optimizations:** Find excited state minima and transition states using spin-corrected ΔSCF.
* **MECI Optimization:** Optimize S0/S1 Minimum Energy Crossing Points using corrected gradients.

Since the ORCA optimizer takes the whole task of changing coordinates (based on gradient), it can also be used to run such jobs as constrained optimizations and potential energy surface scans.

---

## Arguments & Parameters



### Original Parameters (Inherited from ExtOptORCA-OpenQP)
* `--basis_set` / `-b` : The basis set in ORCA format
* `--charge` / `-c` : System charge, default = 0
* `--conical` : Calculates the so-called "penalty function" for the S0/S1 MECI search. Read more at https://sharc-md.org/?page_id=1454#tth_sEc8.20
or for the full paper: https://doi.org/10.1021/jp0761618. 
Associated arguments for the penalty function construction are `--alpha` and `--sigma`



### New Parameters (Added for ΔSCF Correction)
* `--functional` / `-f` : Exchange-correlation functional
* `--aux_basis_set` : Auxilary basis set for RI-MP2
* `--ncpus` / `-cpu` : Number of cpus allocated for SCF and ΔSCF computations, default = 12
* `--gs_orbs_at_each_step` : Enables ground-state UKS orbitals recalculation at each step of geometry optimization to form the input non-Aufbau configuration for ΔSCF, default = False
* `--converger` : Converging algorothm for ΔSCF. See available algorithms in ORCA manual, default = MOM
* `--do_not_purify` : Turns off the spin purification procedure. For MECI optimization using spin-contaminated gradients, default = False

---

## Examples

1. *HBDI_S1_Min_Opt* - Optimization of the planar S1 minimum of HBDI anion using gradients corrected for spin-contamination.

2. *Azobenzene_S0/S1_MECI_Opt* - Optimization of the S0/S1 MECI of azobenzene using penalty function method with parameters sigma=3.5 and alpha=0.02. 

