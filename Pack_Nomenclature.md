Brief summary on Packing and Rearrangement Nomenclature.

#################   Static Packings   #################

Static packings are organized in sets defined by the particle number N and the target pressure P0.
Each packing have a position file, but the whole set have a log file.

- Packing name : the name of the packing contains the N-P set informations and a number to identify it. 
It should be like :

name = N(Particle number)~P(Pressure)~(ID number)

Example : N16~P1e-1~0001     = 16 particles, target pressure of 1e-1, packing number 0001 of the set

- Packing contents : 
N  number of particles
L1  first lattice vector
L2  second lattice vector
P0  target pressure
P  real pressure
X	Y	R

- Log name : same as packings, but with "log" at the beginning and no identification number at the end.
Example :  log~N16~P1e-1

- Log contents : 

Packing number / seed
N  number of particles
P0  target pressure
P  real pressure
alpha  
delta  
L  box size
phi
z
#ratt  number of rattlers
sigmaxx
sigmayy
sigmaxy
U
dU
H
dH
t  run time
#FIRE  number of FIRE iterations
#CG  number of CG iterations
gg  max gradient component (in absolute value)
hh/mm/ss jj/mm/yyyy  creation date

#################   Shear / Compression   #################

The shearing or compression occurs in different steps, according to different protocols.
The positions at each step are recorded into a single position file.
The relevant quantities are recorded for each step into a log file and a data file.



- position file name : 
The name should contain the particle number, the original target pressure, the protocol, the number of steps and the ID number of the original packing :

name = N(Particle number)~P(Pressure)~(protocol)~step(number of steps)~(ID number)

the protocol names are :
SRn = Shear, stop at n rearrangements
SSx = Shear, stop at strain value x
CRn = compression, stop at n rearrangements
CSx = compression, stop at strain value x
DRn = Decompression, stop at n rearrangements
DSx = Decompress, stop at strain value x

note : x should be a dimensionless value, always defined with respect to the original state.

Example : N16~P1e-1~CS1e-4~step58~0001    
 = from the packing above (16 particles, P0=1e-1), you do a compression strain 1e-4 in 48 steps.

- position file contents : 
several static position files concatenated.
Example : 


N	L1	L2	P0	P
X	Y	R

N	L1	L2	P0	P
X	Y	R

N	L1	L2	P0	P
X	Y	R

...


- Log file name : "log" + name of the position file.

Log file contents : same as the static log file, except some minor differences : 
.Instead of ID number / seed, we have step number 
.We add G at the end

- data file name : 
- data file contents : for each step we need 

epsilon  the strain difference with the original packing 
sigma  the actual stress (corresopnding to the shear direction)
G  the shear modulus
Nc  number of contacts
N+  number of contact created (with respect to the previous step)
N-  number of contact broken (with respect to the previous step)



















