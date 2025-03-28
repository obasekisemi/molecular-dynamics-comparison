#!/usr/bin/bin/python3
import numpy as np
from math import *
import random

# System parameters
__d__ = 3           # dimensions
NS = 1000          # number of steps
N = 100            # number of particles
dt = 1.E-4         # timestep
T_0 = 300.         # target temperature
sig = 1.           # LJ sigma
eps = 1.           # LJ epsilon
r_ctf = 2.5*sig    # cutoff distance
bs = 10.*sig       # box size
vol = bs**3        # volume
rho = N/vol        # density
ign = 20          # equilibration steps to ignore

# Calculate cutoff potential and force
u_at_ctf = 4.*eps*((sig/r_ctf)**12 - (sig/r_ctf)**6)
du_at_ctf = 24.*eps*(-2.*sig**12/r_ctf**13 + sig**6/r_ctf**7)

# Initialize arrays
pos = np.zeros([__d__, N])
vel = np.zeros([__d__, N])
acc = np.zeros([__d__, N])
R = np.zeros([__d__, N, N])
r = np.zeros([__d__, N, N])
r2 = np.zeros([N, N])

# Initialize positions randomly
for j in range(N):
    for i in range(__d__):
        pos[i,j] = random.uniform(0, bs)
        vel[i,j] = random.uniform(-0.5, 0.5)  # Changed velocity range

# Scale positions to box size and center system
pos = pos/bs

# Remove center of mass motion
com = np.zeros(__d__)
for i in range(__d__):
    com[i] = np.mean(pos[i,:])
    pos[i,:] -= com[i]

# Remove center of mass velocity
for i in range(__d__):
    vel[i,:] -= np.mean(vel[i,:])

# Scale initial velocities to target temperature
v2 = np.sum(vel**2)
T_current = v2/(3*N - 3)
scale_factor = np.sqrt(T_0/T_current)
vel *= scale_factor

# Initialize statistics
P_S = 0.
k_S = 0.
p_S = 0.
T_S = 0.

# Open output files
f_tpz = open('tpz.out', 'w')
f_kpe = open('kpe.out', 'w')
f_xyz = open('pos.xyz', 'w')

def sign(a, b):
    return a if b >= 0 else -a

# Main MD loop
print("Starting simulation...")
for step in range(1, NS + 1):
    # Periodic boundary conditions
    pos = pos % 1.0
    pos[pos > 0.5] -= 1.0
    
    # Reset forces
    acc.fill(0)
    pot = np.zeros(N)
    vrl = 0.
    
    # Force calculation
    for i in range(N):
        for j in range(i + 1, N):
            r2[i, j] = 0
            rij = np.zeros(__d__)
            
            for k in range(__d__):
                rij[k] = pos[k,i] - pos[k,j]
                # Minimum image convention
                if abs(rij[k]) > 0.5:
                    rij[k] -= sign(1.0, rij[k])
                R[k,i,j] = rij[k]
                r[k,i,j] = bs * rij[k]
                r2[i,j] += r[k,i,j]**2
            
            if r2[i,j] < r_ctf*r_ctf:
                r1 = np.sqrt(r2[i,j])
                ri2 = 1.0/r2[i,j]
                ri6 = ri2**3
                ri12 = ri6**2
                
                # LJ force calculation
                force_mag = 24.0*eps*ri2*(2.0*ri12 - ri6)
                
                for k in range(__d__):
                    f = force_mag * R[k,i,j]
                    acc[k,i] += f
                    acc[k,j] -= f
                
                # Potential energy
                u = 4.0*eps*(ri12 - ri6) - u_at_ctf - r1*du_at_ctf
                pot[i] += 0.5*u
                pot[j] += 0.5*u
                vrl += force_mag * r2[i,j]
    
    vrl = -vrl/__d__
    
    # Velocity Verlet integration
    vel += 0.5*dt*acc
    pos += dt*vel
    
    # Calculate kinetic energy and temperature
    v2 = np.sum(vel**2)
    kin_energy = 0.5*v2*bs*bs
    T_i = 2.0*kin_energy/(3*N)
    
    # Temperature rescaling (simple thermostat)
    scale = np.sqrt(T_0/T_i)
    vel *= scale
    
    # Complete velocity Verlet step
    vel += 0.5*dt*acc
    
    # Calculate averages
    k_AVG = kin_energy/N
    p_AVG = np.sum(pot)/N
    etot_AVG = k_AVG + p_AVG
    P = rho*T_i + vrl/vol
    Z = P*vol/(N*T_i)
    
    # Write output
    if step % 10 == 0:  # Reduced output frequency
        f_tpz.write(f"{step*dt:e} {T_i:e} {P:e} {Z:e}\n")
        f_kpe.write(f"{step*dt:e} {k_AVG:e} {p_AVG:e} {etot_AVG:e}\n")
        f_xyz.write(f"{N}\nStep {step}\n")
        for i in range(N):
            f_xyz.write(f"{pos[0,i]*bs:e} {pos[1,i]*bs:e} {pos[2,i]*bs:e}\n")
    
    # Accumulate statistics
    if step > ign:
        P_S += P
        k_S += k_AVG
        p_S += p_AVG
        T_S += T_i
    
    # Print progress
    if step % 100 == 0:
        print(f"Step {step}/{NS}")
        print(f"Temperature: {T_i:.2f}")
        print(f"Total Energy: {etot_AVG:.2f}")
        print("-" * 20)

# Print final averages
n_avg = NS - ign
print("\nFinal Averages:")
print(f"Temperature: {T_S/n_avg:.3f}")
print(f"Kinetic Energy: {k_S/n_avg:.3f}")
print(f"Potential Energy: {p_S/n_avg:.3f}")
print(f"Total Energy: {(k_S + p_S)/n_avg:.3f}")
print(f"Pressure: {P_S/n_avg:.3f}")

# Close files
f_tpz.close()
f_kpe.close()
f_xyz.close()

print("\nSimulation completed!")