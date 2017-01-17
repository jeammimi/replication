# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:36:10 2017

@author: jarbona
"""

import hoomd
from hoomd import data,init,md,group, dump , deprecated, analyze
import numpy as np
import scipy.linalg as linalg

hoomd.context.initialize("--mode=cpu")

Np = 16   #Number of polymer chains
npp = 150 * 2 * 2  # size of each polymer chain
R = 16 * 2
data_folder = "./data/"
N_diffu = 40 #Number of diffusing elements x2



#########################################
#Define polymer bonding and positions

snapshot = data.make_snapshot(N=Np * npp + N_diffu*2, box=data.boxdim(L=2*R), bond_types=['polymer'])

snapshot.bonds.resize(Np * (npp - 1) + N_diffu)

Npp = [ npp for i in range(Np)]

snapshot.bonds.types = ['polymer_A']
plist = ['A','Ori',"Diffu"]
snapshot.particles.types = plist


offset_bond = 0
offset_particle = 0

#Polymer chains
for i in range(Np):
    npp = Npp[i] # Number of particles
    pos_origin = [np.random.randint(npp)] #Position of origin of replication
    
    initp = 2*np.random.rand(3)-1
    for p in range(npp-1):
        snapshot.bonds.group[offset_bond + p] = [offset_particle + p, offset_particle + p + 1]
        snapshot.bonds.typeid[offset_bond + p] = 0  # polymer_A
    offset_bond += npp - 1
     
    for p in range(npp):
    #print(offset_bond, offset_bond + p)
        new = 2*(2*np.random.rand(3)-1)
        while linalg.norm(initp + new) > R-1:
            new = 2*(2*np.random.rand(3)-1)
        
        initp +=new
        snapshot.particles.position[offset_particle + p ] = initp
        
        if p in pos_origin:
            snapshot.particles.typeid[offset_particle + p ] = 1  #Ori
        else:
            snapshot.particles.typeid[offset_particle + p ] = 0  #A
            

        
    offset_particle += npp 
    
#Diffusing element  
for i in range(N_diffu):
    npp = 2 # Number of particles
    
    initp = 2*np.random.rand(3)-1
    for p in range(npp-1):
        snapshot.bonds.group[offset_bond + p] = [offset_particle + p, offset_particle + p + 1]
        snapshot.bonds.typeid[offset_bond + p] = 0  # polymer_A
    offset_bond += npp - 1
     
    for p in range(npp):
    #print(offset_bond, offset_bond + p)
        new = 2*(2*np.random.rand(3)-1)
        while linalg.norm(initp + new) > R-1:
            new = 2*(2*np.random.rand(3)-1)
        
        initp +=new
        snapshot.particles.position[offset_particle + p ] = initp
        snapshot.particles.typeid[offset_particle + p ] = 2  #Diffu
            

        
    offset_particle += npp 

#Load the configuration

system = init.read_snapshot(snapshot)
###############################################


###############################################
#Defining force field:
harmonic=md.bond.harmonic()
harmonic.bond_coeff.set('polymer_A', k=330.0, r0=1)


nl = md.nlist.tree(check_period=1)

#Potential for warmup
gauss = md.pair.gauss(r_cut=3.0, nlist=nl)

gauss.pair_coeff.set(plist, plist, epsilon=1.0, sigma=1.0)
#gauss.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
#gauss.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)


#Spherical confinement
sphere = md.wall.group()
sphere.add_sphere(r=R, origin=(0.0, 0.0, 0.0), inside=True)

wall_force_slj=md.wall.slj(sphere, r_cut=3.0)
wall_force_slj.force_coeff.set(plist, epsilon=1.0, sigma=1.0,r_cut=1.12)

#Group;
all_beads = group.all()


#Log
logger = analyze.log(filename=data_folder + 'mylog.log', period=1000, quantities=['temperature','potential_energy','kinetic_energy','volume','pressure'],overwrite=True)


#Warmup
method=md.integrate.mode_minimize_fire(group=all_beads,dt=0.05)
while not(method.has_converged()):
   hoomd.run(100)
   
   
gauss.disable()

slj=md.pair.slj(r_cut=2, nlist=nl)
slj.pair_coeff.set(plist,plist,sigma=1,epsilon=1,r_cut=1.12)
print("Second minimizing")
method=md.integrate.mode_minimize_fire(group=all_beads,dt=0.05)
while not(method.has_converged()):
   hoomd.run(100)
#hoomd.run(1000000)
#method.disable()


#Dumping
xml = deprecated.dump.xml(filename=data_folder + "atoms.xml",period=None,group=all_beads,vis=True)
gsd = dump.gsd(filename=data_folder + "atoms.gsd",period=None,group=all_beads)
dcd = dump.dcd(filename=data_folder + 'poly.dcd', period=100,overwrite=True)


#Dynamics

import time
t0 = time.time()
md.integrate.mode_standard(dt=0.005)
method=md.integrate.langevin(group=all_beads,kT=1,seed=np.random.randint(0,100))
snp = system.take_snapshot()

for i in range(10):
    #system.restore_snapshot(snp)
    hoomd.run(1000)
    #snp = system.take_snapshot()

print(gauss.get_energy(all_beads),wall_force_slj.get_energy(all_beads))
print(time.time() -t0)