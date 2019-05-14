# DMFT Lanczos Exact Diagonaization

A Lanczos based solver for the Dynamical Mean-Field Theory using the N_up:N_dw implementation.  
Git repo: [https://github.com/aamaricci/dmft-lanc-ed](https://github.com/aamaricci/dmft-lanc-ed)

**This code solves the normal (N_up, N_dw) case only.**
(*For a generic implemenation which allows for superconductivity (s-wave) or SOC see here:*
[https://github.com/aamaricci/dmft-ed](https://github.com/aamaricci/dmft-ed)

The code depends on:  
* SciFortran ( available at [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran) )  

* DMFT_Tools ( available at [https://github.com/aamaricci/DMFTtools](https://github.com/aamaricci/DMFTtools) ) [For the driver part, see below]

* MPI [although not mandatory]


The code structure is as follow:

* The set of modules compile into a top layer named `DMFT_ED.f90`  
* The actual implementation of the DMFT equations is case by case performed in a driver program, usually placed in the directory `drivers`. 
* In the driver code the user must includes the `DMFT_ED` module and call the necessary procedures to solve the DMFT equations.

An example, solving the Hubbard model on the Bethe lattice, is contained in the file `edn_hm_bethe.f90`.

To install the code

- `git clone https://github.com/aamaricci/dmft-lanc-ed DMFT_LANC_ED`  

- `cd DMFT_LANC_ED` 

* `mkdir build` 
* `cd build` 
* `cmake .. -DEXE=edn_hm_bethe`
