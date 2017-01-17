# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:36:10 2017

@author: jarbona
"""

import hoomd
from hoomd import data,init,md,group, dump , deprecated, analyze
import numpy as np
import scipy.linalg as linalg
from scipy.spatial.distance import cdist

hoomd.context.initialize("--mode=cpu")

Np = 16   #Number of polymer chains
npp = int(150 * 2 * 2   )# size of each polymer chain
R = 16 * 2
data_folder = "./data/"
N_diffu = 100 #Number of diffusing elements x2
cut_off_inte = 10 #1.5
p_inte = 0. # 1/5
dt = 0.1

#########################################
#Define polymer bonding and positions

snapshot = data.make_snapshot(N=Np * npp + N_diffu*2, box=data.boxdim(L=2*R), bond_types=['polymer'])

snapshot.bonds.resize(Np * (npp - 1) + N_diffu)

Npp = [ npp for i in range(Np)]

bond_list = ['Mono_Mono','Diff_Diff','Mono_Diff']
snapshot.bonds.types = bond_list
plist = ['Mono','Ori',"Diff","A_Ori","S_Diff",'F_Diff']

snapshot.particles.types = plist

class Origin:
    def __init__(self,tag):
        self.tag = tag
        self.move=False

class Fork:
    def __init__(self,tag,position,bond_tag):
        self.tag = tag  #particle tagg
        self.position = position #Position on the chromosome
        self.bond_tag = bond_tag #Bond between diffusive elements and monomec
        self.move=True
        self.update_bond=False

    def update_position(self,dt):
        oldp = int(self.position)
        self.position += self.d * dt    
        if oldp != int(self.position):
            self.update_bond =  True
        else:
            self.update_bond =  False
        
class RFork(Fork):
    def __init__(self,tag,position,bond_tag):
        Fork.__init__(self,tag,position,bond_tag)
        self.d = 1
        
class LFork(Fork):
    def __init__(self,tag,position,bond_tag):
        Fork.__init__(self,tag,position,bond_tag)
        self.d = -1

class Polymer():
    def __init__(self,start,end,origins):
        self.start = start
        self.end = end
        self.origins = origins
        self.modules = [Origin(tag) for tag in origins]
    def has_origin(self,ptag):
        if ptag in self.origins:
            return True
        
    def add_fork(self,ptags,otag,new_btag):
        for i,mod in enumerate(self.modules):
            if mod.tag == otag:
                break
        self.modules.insert(i,RFork(ptags[1],otag,new_btag[1]))
        self.modules.insert(i,LFork(ptags[0],otag,new_btag[0]))
        
    def increment_time(self,dt):
        update_bond = []
        for m in self.modules:
            if self.move:
                m.update_position(dt)
                if m.update_bond:
                    update_bond.append([m.bond_tag,int(self.position)])
                    
        #chek for colisions:
        #for m in  
              
        return [],update_bond,[],[],[]
        #bind_diff,shifted_bonds,passivated_origin,to_release,alone = P.increment_time(dt)

        

offset_bond = 0
offset_particle = 0
lPolymers = []
#Polymer chains
for i in range(Np):
    npp = Npp[i] # Number of particles
    pos_origins = [np.random.randint(npp)] #Position of origin of replication
    
    initp = 2*np.random.rand(3)-1
    for p in range(npp-1):
        snapshot.bonds.group[offset_bond + p] = [offset_particle + p, offset_particle + p + 1]
        snapshot.bonds.typeid[offset_bond + p] = bond_list.index('Mono_Mono')  # polymer_A
    offset_bond += npp - 1
     
    for p in range(npp):
    #print(offset_bond, offset_bond + p)
        new = 2*(2*np.random.rand(3)-1)
        while linalg.norm(initp + new) > R-1:
            new = 2*(2*np.random.rand(3)-1)
        
        initp +=new
        snapshot.particles.position[offset_particle + p ] = initp
        
        if p in pos_origins:
            snapshot.particles.typeid[offset_particle + p ] = 1  #Ori
        else:
            snapshot.particles.typeid[offset_particle + p ] = 0  #A
            

    lPolymers.append(Polymer(offset_particle,
                             offset_particle + npp,
                             [po + offset_particle for po in pos_origins]))
    offset_particle += npp 

############################################################
#Diffusing elements
#Defining useful classes



       
    
#Defining particles and bonds for the simulation
    
for i in range(N_diffu):
    npp = 2 # Number of particles
    
    initp = 2*np.random.rand(3)-1
    for p in range(npp-1):
        snapshot.bonds.group[offset_bond + p] = [offset_particle + p, offset_particle + p + 1]
        snapshot.bonds.typeid[offset_bond + p] = bond_list.index('Diff_Diff') # Diff_Diff
    offset_bond += npp - 1
     
    for p in range(npp):
    #print(offset_bond, offset_bond + p)
        new = 2*(2*np.random.rand(3)-1)
        while linalg.norm(initp + new) > R-1:
            new = 2*(2*np.random.rand(3)-1)
        
        initp +=new
        snapshot.particles.position[offset_particle + p ] = initp
        snapshot.particles.typeid[offset_particle + p ] = plist.index("Diff")  #Diffu
            

    
    offset_particle += npp 
    
#Load the configuration

system = init.read_snapshot(snapshot)


for i,p in enumerate(system.particles):
    #print(p)
    #exit()
    assert p.tag == i
    
for i,b in enumerate(system.bonds):
    #print(p)
    #exit()
    assert b.tag == i
###############################################




###############################################
#Defining force field:
harmonic=md.bond.harmonic()
harmonic.bond_coeff.set(bond_list, k=330.0, r0=1)

harmonic.bond_coeff.set('Mono_Diff', k=10.0, r0=1)

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
snp = system#.take_snapshot()



    
 
def Change_type(typep,particle_list,snp):
    #print(particle_list)
    for p in particle_list:
        snp.particles[p].type = typep
        

def Bind(typeb,bondlist,snp):
    btags = []
    for b1,b2 in bondlist:
        btags.append(snp.bonds.add(typeb,b1,b2))
    return btags
    
def Release(btags,snp):
    for bt in btags:
        snp.bonds.remove(bt)
        
def Shift(bonds):
    for tag,new in bonds:
        b = snp.bond.get(tag)
        b.a = new
        

for i in range(5):
    #system.restore_snapshot(snp)
    hoomd.run(1000)
    
    #snp = system.take_snapshot()
    
    
    # update the position of the monomer by updating bonds
    
    group_diffu = group.type(name="Diff", type='Diff')
    group_origin = group.type(name="Ori", type='Ori')
    
    

    for P in lPolymers:
        bind_diff,shifted_bonds,passivated_origin,to_release,alone = P.increment_time(1)
        
        Bind("Diff_Diff",bind_diff,snp) # Pair of diffu to attach
        Shift(shifted_bonds)
        Release.append(to_release)  #Bond tags to release (Alone particle)
        Change_type.append("A_ori",passivated_origin,snp)
        Change_type.append("Mono_Diffu",alone,snp)
        
        
    #Update Type because of (Ori to passivated)

    #Update group


    # First check if Dimer are close from one origin  

    p_diffu = np.array([p.position for p in group_diffu])
    tag_diffu = [p.tag for p in group_diffu]
    p_origin = np.array([p.position for p in group_origin])
    tag_origin = [p.tag for p in group_origin]
    
    distances = cdist(p_diffu,p_origin) 
    print(distances.shape)
    #Reorder the distances with the dimer tags
    Indexes = []
    PTags = []
    #t0 = time.time()
    Btags = []
    for b in system.bonds:
        if b.type == 'Diff_Diff':
            Indexes.append(tag_diffu.index(b.a))
            Indexes.append(tag_diffu.index(b.b))
            Btags.append(b.tag)
            PTags.append([b.a,b.b])
        
    #print(time.time() -t0)
    
    d2 = distances[Indexes][::2] / 2 + distances[Indexes][1::2] /2
    activated = []
    for iD,(btag,ptags) in enumerate(zip(Btags,PTags)):
        #print(d2.shape)
        #print(d2[iD])
        for iorigin,di in enumerate(d2[iD]):
            if iorigin in activated:
                #Needed because we don't want an origin to be activated twice
                continue
            if di < cut_off_inte:
                if np.random.rand() > p_inte:
                    
                    Release([btag],snp) #Break the dimer
                    
                    for P in lPolymers:
                        if P.has_origin(tag_origin[iorigin]):
                            
                            Change_type('F_Diff',ptags,snp)
                            particular_origin = tag_origin[iorigin]
                            Change_type("A_Ori", [particular_origin], snp)
                            new_btags = Bind("Mono_Diff", [[particular_origin, ptags[0]],
                                              [particular_origin ,ptags[1]]], snp)
                            activated.append(iorigin)
                            
                            P.add_fork(ptags,tag_origin,new_btags)

                            break
                    break
    
    # Then if it is the case attach them according to p law to the origin
    



print(gauss.get_energy(all_beads),wall_force_slj.get_energy(all_beads))
print(time.time() -t0)