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

Np = 16
npp = 150 * 2 * 2
R = 16 * 2



#########################################
#Define polymer bonding and positions

snapshot = data.make_snapshot(N=Np * npp, box=data.boxdim(L=2*R), bond_types=['polymer'])

snapshot.bonds.resize(Np * (npp - 1))

Npp = [ npp for i in range(Np)]

snapshot.bonds.types = ['polymer_A']
snapshot.particles.types = ["A"]


offset_bond = 0
offset_particle = 0

for i in range(Np):
    npp = Npp[i] 
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
        snapshot.particles.typeid[offset_particle + p ] = 0  #A

        
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
gauss.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)


#Spherical confinement
sphere = md.wall.group()
sphere.add_sphere(r=R, origin=(0.0, 0.0, 0.0), inside=True)

wall_force_slj=md.wall.slj(sphere, r_cut=3.0)
wall_force_slj.force_coeff.set('A', epsilon=1.0, sigma=1.0,r_cut=1.12)

#Group;
all_beads = group.all()


#Log
logger = analyze.log(filename='mylog.log', period=1000, quantities=['temperature','potential_energy','kinetic_energy','volume','pressure'],overwrite=True)


#Warmup
method=md.integrate.mode_minimize_fire(group=all_beads,dt=0.05)
while not(method.has_converged()):
   hoomd.run(100)
   
   
gauss.disable()

slj=md.pair.slj(r_cut=2, nlist=nl)
slj.pair_coeff.set("A","A",sigma=1,epsilon=1,r_cut=1.12)
print("Second minimizing")
method=md.integrate.mode_minimize_fire(group=all_beads,dt=0.05)
while not(method.has_converged()):
   hoomd.run(100)
#hoomd.run(1000000)
#method.disable()


#Dumping
xml = deprecated.dump.xml(filename="atoms.xml",period=None,group=all_beads,vis=True)
gsd = dump.gsd(filename="atoms.gsd",period=None,group=all_beads)
dcd = dump.dcd(filename='poly.dcd', period=100,overwrite=True)


#Dynamics

import time
t0 = time.time()
md.integrate.mode_standard(dt=0.005)
method=md.integrate.langevin(group=all_beads,kT=1,seed=np.random.randint(0,100))
snp = system.take_snapshot()

for i in range(10):
    system.restore_snapshot(snp)
    hoomd.run(1000)
    snp = system.take_snapshot()

print(gauss.get_energy(all_beads),wall_force_slj.get_energy(all_beads))
print(time.time() -t0)