---
published: true
title: Hartree-Fock SCF
---
<div style="text-align: justify"> I believe the best way to understand the Hartree-Fock Self-Consistent Field Procedure is by actually applying it to a system. Therefore, I have taken an example of the simplest closed-shell heteronuclear molecule, the H—He+ molecule, and solved for its energies using the STO-1G basis set. Here, I have written python code to approach this problem.</div>
* * *

### Step 1: Specify the geometry and basis set for the molecule.


The bond length is taken to be 0.800 Å i.e. 1.5117 a.u. (bohr). Here I will use the simplest Gaussian basis set, i.e. 1s atomic orbital on each of the atoms is approximated by one s type Gaussian function. An s-type Gaussian function: 

$$\phi=a.\exp(-br^2)$$

This is called an STO-1G basis set, which stands for Slater-type orbitals-one Gaussian since we are approximating a slater-type 1s orbital with one Gaussian function. The STO-1G approximations used here are: 

$$\phi(H)=\phi_1=0.3696 \exp(-0.4166 |r-R_1|^2)$$

$$\phi(He)=\phi_2=0.5881 \exp(-0.7739 |r-R_2|^2)$$

Here **(r-Ri)** is the distance of the electron in $$\phi_i$$ from i<sup>th</sup> nucleus.The b constant in the helium exponent is larger when compared to hydrogen (0.7739 vs. 0.4166) since an electron in $$\phi_2$$ is more tightly bound to the doubly charged helium nucleus when compared to singly charged hydrogen nucleus. Thus, the electron density of the Helium nucleus falls off more quickly.

### Step 2: Calculation of one electron integrals.

Here we extract all the 1-electron integrals needed to form the Fock matrix from their respective files and then store them as matrices. To form core Hamiltonian H<sup>core</sup> we add the kinetic energy integral to the two potential energy intervals, V(H) and V(He). 
They are called one-electron integrals since they involve coordinates of only one electron at a time. 

$$T_{rs} = \int \phi_r(r)\left(-\frac{1}{2}\nabla_r^2\right)\phi_s(r)dr$$

$$V_{rs} = \int \phi_r(r)\left(-\sum_A^N\frac{Z}{r_A}\right)\phi_s(r)dr$$

$$H^{core}_{rs} =T_{rs} + V_{rs}(H) +V_{rs}(He)$$

```python
import numpy as np
import scipy.linalg as la


Vf= np.genfromtxt('v.txt',dtype=None)
Tf= np.genfromtxt('t.txt',dtype=None)
Sf= np.genfromtxt('s.txt',dtype=None)
Enuc = np.genfromtxt('enuc.txt',dtype=None)
erif= np.genfromtxt('eri.txt',dtype=None)


for i in Vf:
    m= np.amax(i[1])

#print(m)
#m gives the number of basis functions.
    
#Number of electrons in the given molecule.
Ne = 2

#The extracted data from files is stored as numpy arrays.
V = np.zeros((m,m))
T = np.zeros((m,m))
S = np.zeros((m,m))

#The matrices are symmetric, therefore we convert the given lower triangular matrix
#into full matrix for convenience.

for i in Vf: 
    V[i[0]-1,i[1]-1] = i[2]
    V[i[1]-1,i[0]-1] = i[2]

#print('Nuclear Attraction Intergrals:\n\n', V,'\n\n')

for i in Tf: 
    T[i[0]-1,i[1]-1] = i[2]
    T[i[1]-1,i[0]-1] = i[2]

#print('Kinetic-Energy Integrals:\n\n', T,'\n\n')

H=np.zeros((m,m))

for i in range(m):    
    for j in range(m): 
        H[i,j] = V[i,j] + T[i,j] 

#print('Core Hamiltoninan:\n\n', H,'\n\n')  
```

```output
Nuclear Attraction Intergrals:
-2.2855 -1.5555
-1.5555 -3.4639


Kinetic-Energy Integrals:
0.6249 0.2395
0.2395 1.1609


Core Hamiltoninan:
-1.6606 -1.3160
-1.3160 -2.3030
```

Here we observe that T<sub>11</sub> is smaller than T<sub>22</sub> i.e. the kinetic energy of an electron in $$\phi_1 \:(\phi(H))$$ is smaller than an electron in $$\phi_2 \:(\phi(He))$$ which is expected as the electron around the stronger-pulling He nucleus needs to move faster to stay in its orbit around it. While T<sub>11</sub> denotes the kinetic energy of the electron in the H(1s)-He(1s) overlap region. 

$$V(H)= \left( \begin{array}{cc} -1.0300 & -0.4445 \\ -0.4445 & -0.6563 \end{array} \right) $$
$$V(He)= \left( \begin{array}{cc} -1.2555 & -1.1110 \\ -1.1110 & -2.8076 \end{array} \right) $$

The total potential energy matrix can be disintegrated into the magnitude of attraction of the electron to the H and He nucleus. Here V<sub>11</sub>(H) represents the attraction of an electron in $$\phi_1$$ to the hydrogen nucleus which is understandably larger than the V<sub>22</sub>(H) representing the attraction of an electron in $$\phi_2$$ to the hydrogen nucleus since the electron in $$\phi_1 \:(\phi(H))$$ is attracted more strongly to the H nucleus than the electron in $$\phi_2 \:(\phi(He))$$. While V<sub>12</sub>(H) represents the attraction of an electron present in the H(1s)-He(1s) overlap region to the H nucleus. Similar arguments are true for V(He) also. 


