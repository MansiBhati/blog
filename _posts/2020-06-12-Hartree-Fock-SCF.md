---
published: true
title: Hartree-Fock SCF
---
<div style="text-align: justify"> I believe the best way to understand the Hartree-Fock Self-Consistent Field Procedure is by actually applying it to a system. Therefore, I have taken an example of the simplest closed-shell heteronuclear molecule, the H—He+ molecule, and solved for its energies using the STO-1G basis set. Here, I have written python code to approach this problem. </div>
* * *

The bond length is taken to be 0.800 Å i.e. 1.5117 a.u. (bohr). Here I will use the simplest Gaussian basis set, i.e. 1s atomic orbital on each of the atoms is approximated by one s type Gaussian function. An s-type Gaussian function: 

$$\phi=a.\exp(-br^2)$$

This is called an STO-1G basis set, which stands for Slater-type orbitals-one Gaussian since we are approximating a slater-type 1s orbital with one Gaussian function. The STO-1G approximations used here are: 

$$\phi(H)=\phi_1=0.3696 \exp(-0.4166 |r-R_1|^2)$$

$$\phi(He)=\phi_2=0.5881 \exp(-0.7739 |r-R_2|^2)$$

Here **|r-R<sub>i</sub>|** is the distance of the electron in $$\phi_i$$ from i<sup>th</sup> nucleus.
The b constant in the helium exponent is larger when compared to hydrogen (0.7739 vs. 0.4166) since an electron in $$\phi_2$$ is more tightly bound to the doubly charged helium nucleus when compared to singly charged hydrogen nucleus. Thus, the electron density of the Helium nucleus falls off more quickly.   