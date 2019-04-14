import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import block_diag

        
def lorentzian(gamma,E):
    return 0.5*gamma/(E**2+(0.5*gamma)**2)*(1/np.pi)


#The number of sites 
N = 10
#The number of particles
N2 = (N*N)
#The size of the hamiltonian
H = np.zeros((N2,N2))
# the bunding energy
t = -1.0


#The position lattice
pos = np.zeros((N, N, 2))


#filling the hamiltonian
k = 1
c1 = 0
for i in range(1, N2 + 1 - N):
    if(i > k*N + 1 and i <(k+1)*N):
        H[i-1, (i-1)-1] = t 
        H[i-1, (i-1)+1] = t 
        H[i-1, (i-1)+N] = t 
        H[i-1, (i-1)-N] = t
        c1 = c1 + 1
        if(c1 == N-2):
            k = k + 1
            c1 = 0



H[0,1] = t
H[0,N-1] = t 
H[0,N*(N-1)] = t 
H[0,N] = t

H[N-1, 0] = t 
H[N-1, N-2] = t
H[N-1, N*N - 1] = t     
H[N-1, N*2 - 1] = t

H[N*(N-1), N*N - 1] = t 
H[N*(N-1), N*(N-1) + 1] = t 
H[N*(N-1), N*(N-1) - N] = t 
H[N*(N-1), 0] = t

H[N*N - 1, N*N - 2] = t 
H[N*N - 1, N*N - N - 1] = t 
H[N*N - 1, N - 1] = t 
H[N*N - 1, N*(N-1)] = t 


for i in range(2, N):
    H[i-1, i-2] = t 
    H[i-1, i] = t 
    H[i-1, i + N*(N-1) - 1] = t 
    H[i-1, i-1 + N] = t

c1 = N*(N-1) + 2

for i in range(c1, N*N):
    H[i-1, i-2] = t 
    H[i-1, i] = t 
    H[i-1, i-1-N] = t 
    H[i-1, i-N*(N-1) - 1] = t 

for i in range (1,N-1):
    H[i*N, i*N + 1] = t
    H[i*N, (i+1)*N - 1] = t
    H[i*N, (i-1)*N] = t
    H[i*N, (i+1)*N] = t

    H[N*(i+1)-1,N*(i+1)-2 ] = t 
    H[N*(i+1)-1, i*N] = t 
    H[N*(i+1)-1, N*(i+1)-1 + N] = t 
    H[N*(i+1)-1, N*(i+1)-1 - N] = t 

H_full = block_diag(H,H)

E,psi = np.linalg.eigh(H_full)

E_arr = np.linspace(min(E),max(E),500)


amp_arr = np.zeros(len(E_arr))

for i in range(len(E)):
    if(E[i] > 1.0):
        min_E = E[i]


for i in range(len(E)):
    for j in range(len(E)):
        new_min = abs(E[i]-E[j])
        if(abs(new_min) > 0.5 and new_min < min_E):
            min_E = new_min




g = min_E/(N)



for i in range(len(E_arr)):
    amp_arr[i]=sum(lorentzian(g,E_arr[i]-Eval) for Eval in E)
    
plt.title("2-D TB :"+str(N*N)+" sites")
plt.xlabel('Energy')
plt.ylabel('DOS')
plt.plot(E_arr,amp_arr)
plt.show()
