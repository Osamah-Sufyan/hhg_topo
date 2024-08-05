import networkx as nx
import matplotlib.pyplot as plt
from ase.units import Bohr, Ha
from ase.build import graphene
import numpy as np
import math
from math import e
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import eigh
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from numba import njit
from scipy.ndimage import gaussian_filter1d



#######section1: cereating a hexagonal lattice flake with arbitrary side size m and filling the time-indpendent hmailtonian in real space, this part you can skip as it has been tested already to work#################


def node_dist(x,y, cx, cy):
    """Distance of each node from the center of the innermost layer"""
    return abs(cx-x) + abs(cy-y)


def remove_unwanted_nodes(G, m):
    """Remove all the nodes that don't belong to an m-layer hexagonal ring."""
    
    cx, cy = m-0.5, 2*m -(m%2) #odd is 2m-1, even is 2m
    
    
    unwanted = []
    for n in G.nodes:    
        x,y = n
        
        if node_dist(x,y, cx, cy) > 2*m:
            unwanted.append(n)

    
    for n in unwanted:
        G.remove_node(n)
        
    return G



scale = 2.68

m = 2 # flake with side size m=2

G = nx.hexagonal_lattice_graph(2*m-1,2*m-1, periodic=False, 
                               with_positions=True, 
                               create_using=None)
G = remove_unwanted_nodes(G, m)
pos = nx.get_node_attributes(G, 'pos')

for node, position in pos.items():
    pos[node] = tuple(round(i * scale,4) for i in position)

dict0= pos


q= list(dict0.values())



r1=[]
for i in np.arange(0,len(dict0)):
    s=q[i][0]
    r1.append(s)
r2=[]
for i in np.arange(0,len(dict0)):
    s=q[i][1]
    r2.append(s)
    
    
A=np.zeros((len(dict0),3))
for i in np.arange(0, len(dict0)):
    A[i][1]= r1[i]
    A[i][2]= r2[i]
    A[i] =i
for i in np.arange(0, len(dict0)):
    A[i][1]= r1[i]
    A[i][2]= r2[i]


Asx= np.arange(0.5,102.5,1.5)
Asx2 = np.arange(1.34, 1.34*1000, 4.0200000000000005)

Acells=[]
Bcells =[]
for i in np.arange(0,len(dict0)):
    if A[i][1] in np.round(Asx2,4):
    
        Acells.append(A[i])
    else:
        Bcells.append(A[i])
        
Acellsin=[]
for i in np.arange(0, len(Acells)):
    x= Acells[i][0]
    Acellsin.append(x)
    
Bcellsin=[]
for i in np.arange(0, len(Bcells)):
    x= Bcells[i][0]
    Bcellsin.append(x)

a_1= (3/2 * scale , np.sqrt(3)/2*scale)
a_2= (-3/2* scale, np.sqrt(3)/2* scale)
a_3= (0,-1*scale*np.sqrt(3))
pAcells=[]
for i in np.arange(0,len(Acells)):
    for j in np.arange(0,len(Acells)):
        dif =[Acells[i][1]-Acells[j][1], Acells[i][2]-Acells[j][2]]
        if round(dif[0]-a_1[0],2)==0 and round(dif[1]-a_1[1])==0:
           x = Acells[i][0]
           y = Acells[j][0]
           pAcells.append((x,y))
        
        
for i in np.arange(0,len(Acells)):
    for j in np.arange(0,len(Acells)):
        dif =[Acells[i][1]-Acells[j][1], Acells[i][2]-Acells[j][2]]
        if round(dif[0]-a_2[0],2)==0 and round(dif[1]-a_2[1])==0:
           x = Acells[i][0]
           y = Acells[j][0]
           pAcells.append((x,y))
for i in np.arange(0,len(Acells)):
    for j in np.arange(0,len(Acells)):
        dif =[Acells[i][1]-Acells[j][1], Acells[i][2]-Acells[j][2]]
        if round(dif[0]-a_3[0],2)==0 and round(dif[1]-a_3[1])==0:
           x = Acells[i][0]
           y = Acells[j][0]
           pAcells.append((x,y))
        
        

pBcells=[]
a_1b= (-3/2*scale, -np.sqrt(3)/2*scale)
a_2b= (3/2*scale, -1*scale*np.sqrt(3)/2)
a_3b= (0,scale*np.sqrt(3))

for i in np.arange(0,len(Bcells)):
    for j in np.arange(0,len(Bcells)):
        dif =[Bcells[i][1]-Bcells[j][1], Bcells[i][2]-Bcells[j][2]]
        if round(dif[0]-a_1b[0],2)==0 and round(dif[1]-a_1b[1])==0:
           x = Bcells[i][0]
           y = Bcells[j][0]
           pBcells.append((x,y))
        
        
for i in np.arange(0,len(Bcells)):
    for j in np.arange(0,len(Bcells)):
        dif =[Bcells[i][1]-Bcells[j][1], Bcells[i][2]-Bcells[j][2]]
        if round(dif[0]-a_2b[0],2)==0 and round(dif[1]-a_2b[1])==0:
           x = Bcells[i][0]
           y = Bcells[j][0]
           pBcells.append((x,y))
for i in np.arange(0,len(Bcells)):
    for j in np.arange(0,len(Bcells)):
        dif =[Bcells[i][1]-Bcells[j][1], Bcells[i][2]-Bcells[j][2]]
        if round(dif[0]-a_3b[0],2)==0 and round(dif[1]-a_3b[1])==0:
           x = Bcells[i][0]
           y = Bcells[j][0]
           pBcells.append((x,y))
        
        
def NN():
    near= []
    for i in np.arange(0,len(dict0)):
        for j in np.arange(0,len(dict0)):
            if round(np.sqrt((A[i][1]-A[j][1])**2+(A[i][2]-A[j][2])**2), 4)==2.68:
                near.append((i,j))
    return near


def NNN():
    nnear= []
    for i in np.arange(0,len(dict0)):
        for j in np.arange(0,len(dict0)):
            if round(np.sqrt((A[i][1]-A[j][1])**2+(A[i][2]-A[j][2])**2), 4)==4.636:
                nnear.append((i,j))
    return nnear

def hamiltonian(M,t_1,t_2):
    hamiltonian = np.zeros((len(dict0),len(dict0)),dtype=np.complex_)
    j=complex(0,1)
    for k in np.arange(0,len(dict0)):
        for l in np.arange(0,len(dict0)):
            if k in Acellsin:
                hamiltonian[k,k] = -M
            
                
            elif k in Bcellsin:
               
                hamiltonian[k,k] = M


    for k in np.arange(0,len(dict0)):
        for l in np.arange(0,len(dict0)):        
            if (k,l) in pAcells:
                hamiltonian[k,l] = t_2*j
                hamiltonian[l,k] = -t_2*j
                
            if (k,l) in pBcells:
                hamiltonian[k,l] = t_2*j
                hamiltonian[l,k] = -t_2*j
            

    for (k,l) in NN():
        hamiltonian[k,l] = t_1
        hamiltonian[l,k] = t_1
    return hamiltonian



#############section 2: solving the time-indepednent Schrödinger equation in real space and plotting the eigenvalues vs space index and spatial distribution of the eigenvectors, this part has been tested so you can skip it as well#################


H=hamiltonian(0,-0.1,0)
ew, ev = LA.eigh(H)
# s= np.argsort(ew[:]) 
# ew[:]=ew[s]
# ev[:]=ev[s]

# print(ew)
# x= np.linspace(0,len(A),len(A))
# start = 295
# end = 309
# plt.plot(x[:start], ew[:start], color='blue')

# # Plot the data in the interval in red
# plt.plot(x[start-5:end+10], ew[start-5:end+10], color='red')

# # Plot the data after the interval in blue
# plt.plot(x[end:], ew[end:], color='blue')

# # Plot the data after the interval in blue
# plt.plot(x[end:], ew[end:], color='blue')
# plt.xlabel('i', fontsize=20)
# plt.ylabel(r'$E_i$', fontsize=20)
# plt.savefig('energy spectrum.pdf')
# plt.show()
# figure size



# plt.figure(figsize=(9,9))
# bond_length = 2.68
# for i in np.arange(0,len(dict0)):
#     for j in np.arange(0,len(dict0)):
#         # Calculate the distance between the two atoms
#         distance = round(np.sqrt((A[i][1]-A[j][1])**2+(A[i][2]-A[j][2])**2),2)
        
#         # If the distance is less than the bond length, draw a line between the atoms
#         if distance == bond_length:
#             plt.plot([A[i][1], A[j][1]], [A[i][2], A[j][2]], 'k-')


# x_coords = [A[i][1] for i in range(len(A))]
# y_coords = [A[i][2] for i in range(len(A))]
# plt.scatter(x_coords, y_coords, s=50, c='blue', marker='o', zorder = 2)      
# # scale the plot
# plt.axis('scaled')
# plt.savefig('lattice.pdf')
# plt.show()

# eigenvectors_squared = [abs(ev[i])**2 for i in range(len(ev))]

# proj =[]
# for i in np.arange(0,len(dict0)):
#     x= np.abs(ev[i,298])**2
#     proj.append(x)
    

# plt.figure(figsize=(9,9))
# bond_length = 2.68

#         # Calculate the distance between the two atoms
        
#         #plt.scatter(x_coords, y_coords, s=eigenvectors_squared[i], c='blue', marker='o', zorder = 2)
        
#         # If the distance is less than the bond length, draw a line between the atoms
        


# # Plot the dots with size proportional to the absolute value squared of the eigenvectors
# x_coords = [A[i][1] for i in range(len(A))]
# y_coords = [A[i][2] for i in range(len(A))]

# for i in np.arange(0,len(dict0)):
#     plt.scatter(x_coords[i], y_coords[i], s=10000*proj[i], c='blue', marker='o', zorder = 2)   
   

# # scale the plot
# plt.axis('scaled')
# plt.savefig('occupation.pdf')
# plt.show()








############SECTION 3: time-propagtion using the time-depednent Hamiltonian given in Equation 5 in the paper#################




###initializing the laser field (vector potential A)

j=complex(0,1)
omega_0 = 7.5e-3  ### angualar frequency in atomic units of the 
A_0 = 20
dt=1000
A_t= []
tint= len(np.arange(0,2*np.pi*10/omega_0,2*np.pi*10/omega_0/dt))
time_steps = np.arange(0, 2 * np.pi * 10 / omega_0, 2 * np.pi * 10 / omega_0 / dt)



t = np.arange(0, 2*np.pi*10/omega_0, 2*np.pi*10/omega_0/dt)
Att = A_0 * (np.sin(omega_0 * t / 20) ** 2) * np.sin(omega_0 * t)





############building the time-dependnet hamiltonan, here a problem can occur####################


Ht=np.zeros((len(t),len(dict0),len(dict0)),dtype=np.complex_)

for t_i, ts in enumerate(t):
    for k in range(len(H)):
        for j in range(len(H[0])):
            Ht[t_i, k, j] = H[k, j] * np.exp(-1j * (A[k][2] - A[j][2]) * Att[t_i])




# plt.figure(figsize=(10, 6))
# plt.plot(t, Att, label='Att vs. Time')
# plt.xlabel('Time')
# plt.ylabel('Att')
# plt.title('Att in Time')
# plt.legend()
# plt.grid(True)
# plt.savefig('laser_pulse.pdf')
# plt.show()





## converting angualr freequency in atomic units to wavelength#######
wavelength_um = 0.04564 / omega_0
print(f"Angular frequency: {omega_0} a.u. -> Wavelength: {wavelength_um:.2f} µm")






def Ha(t):
    # Find the index of the closest time step
    t_index = np.argmin(np.abs(time_steps - t))
    return Ht[t_index]




############# solving the schrödinger equation using the scipy pacage solve_ivp#################

# def schrodinger(t, psi):
#     return -1j * Ha(t) @ psi




# t_start = time_steps[0]
# t_end = time_steps[-1]
# solutions = []

# dim = len(ev)
# for i in np.arange(0,(int(dim))):
    
#         psi0 = ev[:,i]
#         sol = solve_ivp(schrodinger, [t_start, t_end], psi0, t_eval=time_steps, vectorized=True)
#         solutions.append(sol.y)
   





################solving the schrödinger equation using own RK4 code ####################

def schrodinger(t, psi):
    return -1j * Ha(t) @ psi


def RK4_step(f, t, y, dt):
    k1 = dt * f(t, y)
    k2 = dt * f(t + 0.5 * dt, y + 0.5 * k1)
    k3 = dt * f(t + 0.5 * dt, y + 0.5 * k2)
    k4 = dt * f(t + dt, y + k3)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

# Time steps and initial conditions
t_start = time_steps[0]
t_end = time_steps[-1]
dt = time_steps[1] - time_steps[0]
solutions = []

dim = len(ev)
for i in range(dim):
    psi0 = ev[:, i]
    psi_t = psi0
    psi_t_list = [psi0]
    for t in time_steps[:-1]:
        psi_t = RK4_step(schrodinger, t, psi_t, dt)
        psi_t_list.append(psi_t)
    solutions.append(np.array(psi_t_list).T)


def calculate_current(solutions, t_index):
    J = np.zeros(2, dtype=complex)

    H = Ha(time_steps[t_index])
    for l in range(len(solutions)-1):
        for i in range(len(dict0)):
            for j in range(len(dict0)):
                J[0] += (A[i][1] - A[j][1]) * -1j * np.conj(solutions[l][i][t_index]) * H[i][j] * solutions[l][j][t_index]
                J[1] += (A[i][2] - A[j][2]) * -1j * np.conj(solutions[l][i][t_index]) * H[i][j] * solutions[l][j][t_index]
                
                
    return J


J_t = [calculate_current(solutions, t_index) for t_index in range(len(time_steps))]

J_y = [J[1] for J in J_t]

dJ_y_dt = np.gradient(J_y, time_steps)



# Apply a Hann window
hann_window = np.hanning(len(dJ_y_dt))
dJ_y_dt_windowed = dJ_y_dt * hann_window

# Perform FFT on the windowed data
fft_dJ_y_dt = np.fft.fft(dJ_y_dt_windowed)
abs_fft_dJ_y_dt = np.abs(fft_dJ_y_dt)**2

# Calculate frequencies
frequencies = np.fft.fftfreq(len(dJ_y_dt), d=time_steps[1]-time_steps[0])
positive_frequencies =  2*np.pi*frequencies[frequencies >= 0]
abs_fft_dJ_y_dt_positive = abs_fft_dJ_y_dt[frequencies >= 0]

max_frequency = positive_frequencies[-1]


x_ticks = np.arange(0, max_frequency, omega_0)

# Plot the results
plt.figure(figsize=(16, 8))
plt.plot(positive_frequencies, 1e12 * abs_fft_dJ_y_dt_positive)
plt.yscale('log')
plt.xlabel('harmonic order')
plt.ylabel(r'$|FFT(\dot{J})|^2$')
plt.xticks(x_ticks, [f'{i}' for i in range(len(x_ticks))])
plt.title(f'linearly polaralized pulse, wavelength = {wavelength_um} nm, size of flake = {m}')
plt.savefig(f'fft_{wavelength_um}.pdf')
plt.grid(True)
plt.show()











