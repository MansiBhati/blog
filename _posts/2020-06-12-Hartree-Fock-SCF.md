---
published: true
title: Hartree-Fock SCF
---
<div style="text-align: justify"> I believe the best way to understand the Hartree-Fock Self-Consistent Field Procedure is by actually applying it to a system. Therefore, I have taken an example of the simplest closed-shell heteronuclear molecule, the H—He+ molecule, and solved for its energies using the STO-1G basis set. Here, I have written python code to approach this problem.</div>
* * *

### Step1: Specify the geometry and basis set for the molecule.


The bond length is taken to be 0.800 Å i.e. 1.5117 a.u. (bohr). I will use the simplest Gaussian basis set, i.e. 1s atomic orbital on each of the atoms is approximated by one s type Gaussian function. An s-type Gaussian function: 

$$\phi=a.\exp(-br^2)$$

This is called an STO-1G basis set, which stands for Slater-type orbitals-one Gaussian since we are approximating a slater-type 1s orbital with one Gaussian function. The STO-1G approximations used here are: 

$$\phi(H)=\phi_1=0.3696 \exp(-0.4166 |r-R_1|^2)$$

$$\phi(He)=\phi_2=0.5881 \exp(-0.7739 |r-R_2|^2)$$

Here **(r-Ri)** is the distance of the electron in $$\phi_i$$ from i<sup>th</sup> nucleus as shown in the schematic. The b constant in the helium exponent is larger when compared to hydrogen (0.7739 vs. 0.4166) since an electron in $$\phi_2$$ is more tightly bound to the doubly charged helium nucleus when compared to a singly charged hydrogen nucleus. Thus, the electron density of the Helium nucleus falls off more quickly.

![Image]({{site.baseurl}}/assets/images/1.jpg "Coordinate system")

### Step2: Calculation of one electron integrals.

The first step in the pythppn program is to extract all the one-electron integrals needed to form the Fock matrix from their respective files and then store them as matrices. To form core Hamiltonian H<sup>core</sup> we add the kinetic energy integral of the electron to the two potential energy integrals, V(H) and V(He). They are called one-electron integrals since they involve coordinates of only one electron at a moment. They are calculated using:

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

```python
#Output:
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

Here we observe that T<sub>11</sub> is smaller than T<sub>22</sub> i.e. the kinetic energy of an electron in $$\phi_1 \:(\phi(H))$$ is smaller than an electron in $$\phi_2 \:(\phi(He))$$ which is expected as the electron around the stronger-pulling He nucleus needs to move faster to stay in its orbit. While T<sub>12</sub> denotes the kinetic energy of the electron in the H(1s)-He(1s) overlap region. 

Also, the total potential energy matrix can be disintegrated into the magnitude of attraction of the electron to the H and He nucleus.

$$V(H)= \left( \begin{array}{cc} V_{11} & V_{12} \\ V_{21} & V_{22} \end{array} \right)=\left( \begin{array}{cc} -1.0300 & -0.4445 \\ -0.4445 & -0.6563 \end{array} \right) $$

$$V(He)= \left( \begin{array}{cc} -1.2555 & -1.1110 \\ -1.1110 & -2.8076 \end{array} \right) $$

Here V<sub>11</sub>(H) represents the attraction of an electron in $$\phi_1$$ to the hydrogen nucleus which is understandably larger than the V<sub>22</sub>(H) representing the attraction of an electron in $$\phi_2$$ to the hydrogen nucleus since the electron in $$\phi_1 \:(\phi(H))$$ is attracted more strongly to the H nucleus than the electron in $$\phi_2 \:(\phi(He))$$. While V<sub>12</sub>(H) represents the attraction of an electron present in the H(1s)-He(1s) overlap region to the H nucleus. Similar arguments are true for V(He) also. 


### Step3: Calculation of orthogonalising matrix. The S-1/2 matrix.

1. Overlap intergrals are calculated using the given formula, and are already stored in the .txt file. Read the file and store it as a NumPy array.

	$$S_{rs} = \int \phi_r(r)\phi_s(r)dr$$
2. Diagonalise this overlap matrix S.

 	$$S = P.D.P^{-1}$$
3. Calculate D<sup>-1/2</sup>

4. Build the symmetric orthogonalisation matrix S<sup>-1/2</sup>

	$$S^{-1/2} = PD^{-1/2}P^{-1}$$

```python
for i in Sf: 
    S[i[0]-1,i[1]-1] = i[2]
    S[i[1]-1,i[0]-1] = i[2]

#print('The overlap matrix is:\n\n', S,'\n\n')
    
eigval = la.eig(S)
Pi = eigval[1]
Di = np.diag(eigval[0].real**(-0.5))
Sh = Pi @ Di @ la.inv(Pi)
#print('S^-1/2 Matrix:\n\n', Sh,'\n\n') 
```
```python
The overlap matrix is:
1.0 0.5071
0.5071 1.0


S^-1/2 Matrix:
1.11946687 -0.30489583
-0.30489583 1.11946687
```

### Step4: Calculation of initial guess density matrix
 
1. Form an initial Fock matrix F, in the orthonormal basis using the core Hamiltonian H<sub>core</sub>. $$F'_0 = \tilde{S}^{-1/2} H^{core} S^{-1/2}$$

2. Diagonalise F

$$F'_0C'_0 = C'_0\epsilon$$

3. Transform the resulting eigenvectors into the original (non-orthonormal) basis

$$C_0= S^{-1/2} C'_0$$

4. Construct the initial-guess density matrix. 

$$D_{rs}^0 = 2 \sum_m^{occ}(C_0)_r^m(C_0)_s^m$$

### Step5: Calculation of two-electron integrals

The two-electron integrals giving the repulsion between the electrons are calculated by the given formula:

$$P_{rs}= \sum_{t=1}^m\sum_{u=1}^m P_{tu}[(rs|tu)-\frac{1}{2}(ru|ts)]$$

Therefore, each element of the electron repulsion matrix G has eight 2-electron repulsion integrals, making up a total of 32 integrals. However, many of these are the same, and there exists only six unique 2-electron repulsion integrals which are given in the file. 

$$(rs|tu)=(rs|ut)=(sr|tu)=(sr|ut)=(tu|rs)=(tu|sr)=(ut|rs)=(ut|sr)$$

![Image]({{site.baseurl}}/assets/images/overlap.jpg "2-electron repulsion")

1. Read these unique two-integrals from the file and store them in a one-dimensional array using a compound index.

1. Form the new Fock matrix including the two-electron contribution.

$$F_{rs}=H^{core}_{rs} + \sum_{t=1}^m\sum_{u=1}^m P_{tu}[(rs|tu)-\frac{1}{2}(ru|ts)]
=T_{rs} + V_{rs}(H) +V_{rs}(He)+G_{rs}$$

### Step6: The Self-Consistent Field Iteration

1. Transform the new Fock matrix to the orthonormal basis.

$$F' = \tilde{S}^{-1/2} H^{core} S^{-1/2}$$

2. Transform the resulting eigenvectors into the original (non-orthonormal) basis

$$C= S^{-1/2} C'$$

3. Construct the new-guess density matrix.

$$D_{rs}^0 = 2 \sum_m^{occ}(C)_r^m(C)_s^m$$

4. Test for convergence using delta and keep repeating this procedure till you get a converged value of delta and energy. 

$$\left[\frac{\sum_{rs}^{AO}(D_{rs}^i-D_{rs}^{i-1})}{4}\right]^{1/2} <\delta_1$$

$$E_{elec}^i = \sum_{rs}^{AO}\frac{1}{2}P_{rs}^i(H_{rs}^{core}+F_{rs})$$

5. After the convergence is achieved, display the final energy by computing total electronic energy and adding to it the nuclear energy which is given in the text file.

$$E_{total}^0 = E_{elec}^0 +E_{nuc}$$

```python
#This function is to form initial Fock matrix in the orthonormal basis
def Fini(F):
    Fo = np.transpose(Sh) @ F @ Sh
    return Fo

#This is to tranform the resulting eigenvectors into the original
#(non-orthonormal) basis    
def Cmat(M):
    Foeig = la.eigh(M)
    C = Sh @ Foeig[1]
    return C,Foeig[0]

#This is to construct density matrix and store the older values for convergence
#test
def Dmat(C,P):
    OldD = np.zeros((m,m))
    for i in range(m):    
        for j in range(m): 
            OldD[i,j]=P[i,j]
            P[i,j] = 0.0
            for s in range(0,Ne//2):
                P[i,j] = P[i,j] + 2*C[i,s]*C[j,s]
                #print(C[i,s]*C[j,s])
    return P,OldD

#This maps the permutaionally unique two-electron integrals into a compound 
#index.
def eint(a,b,c,d):  
    if a > b: ab = a*(a+1)/2 + b  
    else: ab = b*(b+1)/2 + a  
    if c > d: cd = c*(c+1)/2 + d  
    else: cd = d*(d+1)/2 + c  
    if ab > cd: abcd = ab*(ab+1)/2 + cd  
    else: abcd = cd*(cd+1)/2 + ab  
    return abcd

#It is done to assign value to the two electron integral.
def tei(a,b,c,d):  
    return twoe.get(eint(a,b,c,d),0.0)

#This function gives final Fock matrix.
def Fmat(H, P):
    F = np.zeros((m,m))  
    for i in range(0,m):  
        for j in range(0,m):  
            F[i,j] = H[i,j]  
            for k in range(0,m):  
                for l in range(0,m):  
                    F[i,j] = F[i,j] + P[k,l]*(tei(i+1,j+1,k+1,l+1)-0.5e0*tei(i+1,k+1,j+1,l+1))
    return F

#It calculates delta value used for the test of convergence
def Delta(P,OldD):
   delta = 0.0e0  
   for i in range(0,m):  
       for j in range(0,m):  
           delta = delta+((P[i,j]-OldD[i,j])**2)  
   delta = (delta/4)**(0.5)  
   return delta

#It calculates the electronic energy
def Eelec(P, H, F):
    E = 0.0  
    for i in range(0,m):  
        for j in range(0,m):  
            E = E + 0.5*P[i,j]*(H[i,j] + F[i,j])  
    return E  

twoe = {eint(i[0],i[1],i[2],i[3]) : i[4] for i in erif}
    
delta =1.0
convergence=0.00000001 


i=0
#to define maximum number of iterations
maxit=20
P = np.zeros((m,m))

while(i<maxit):
    i +=1
    
    if (delta>convergence):
        F=Fmat(H, P)
        #print('F',i,'=',F)
        Fi=Fini(F)
        #print('Fi',i,'=',Fi)
        C,Foeig=Cmat(Fi)
        #print('C',i,'=',C)
        #print('C',i,'=',Foeig)
        P,OldD=Dmat(C,P)
        #print('Do',i,'=',P)
        #print('Dold',i,'=',OldD)
        delta=Delta(P,OldD)
        #print('delta=',delta)
        E=Eelec(P, H, F)
        #print('E',i,'=',E+Enuc)
        continue
    
    else:
        print('Convergence achieved!')
        break
```
```python
#Output:
E 1 = -3.3370117315590293
E 2 = -2.40903883091145
E 3 = -2.4313259673325143
E 4 = -2.435384003832221
E 5 = -2.436097100960776
E 6 = -2.4362218038169208
E 7 = -2.4362435940809237
E 8 = -2.4362474011412343
E 9 = -2.436248066271592
E 10 = -2.4362481824758366
E 11 = -2.436248202777747
E 12 = -2.43624820632467
Convergence achieved!
```
***

One should realise that modern ab-initio software doesn’t rigidly use the above SCF procedure. To speed up the calculations, they incorporate a variety of other techniques such as:

1. Use the symmetry of the molecules to avoid duplicate calculations.

1. Extremely small two-electron integrals can be neglected, this happens in case for the functions with far separated nuclei, this decreases the calculation time from n<sup>4</sup> to about n<sup>2.3</sup> dependence.

1. Use of pseudo-spectral method which represents the molecular orbitals as a set of grid points, in addition to a basis set expansion is done. This surpasses the need to explicitly calculate the two-electron integrals, which increases the computational speed by a factor of three-four. 

1. Use of Fast-multipole method is employed, which considers the Coulomb repulsion for distant regions in very large systems as the center of the regions.
