# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 12:42:19 2017

@author: jarbona
"""


                 
from mayavi import mlab
from tvtk.tools import visual
from colour import Color
import mdtraj as md
import numpy as np
# Create a figure
f = mlab.figure(size=(500,500))
# Tell visual to use this as the viewer.
visual.set_viewer(f)

# A silly visualization.

# Even sillier animation.


data_folder = "../data/"
t = md.load(data_folder + 'poly.dcd', top=data_folder + "atoms.hoomdxml")

class Representation:
    def __init__(self,traj,what,by_resid=True,tube=True,t=0,update=False,time_color=None):
        self.first = True
        self.traj = traj
        self.selections = []
        self.time_color = time_color
        self.tube = tube 
        if by_resid:
            n_chains = len(list(self.traj.topology.chains))
            
            for i in range(n_chains):
                add = ""
                if what != "":
                    add = " and %s"%what
                sele = "resid %i"%i+add
                print(sele)
                sel = self.traj.topology.select(sele)
                
                
                if len(sel) != 0:
                    self.selections.append(sel)
                    
                    
    def draw(self,time):
        if self.first:
            
            self.rep = []
            
        for isel,sel in enumerate(self.selections):
            x,y,z = self.traj.xyz[time, sel].T
            if self.time_color != None:
                colors = self.time_color(time)
            else:
                colors = [1 for i in range(len(x))]
                
            #print(len(x),len(colors))
            
            if self.first:  
                if self.tube:
                    self.rep.append(mlab.plot3d(x, y, z,colors ,vmin=1,vmax=2.1,tube_radius=0.03,opacity=1))
                else:
                    self.rep.append(mlab.points3d(x, y, z,colors ,scale_factor=0.1,color=Color("red").rgb))
            else:
                
                self.rep[isel].mlab_source.set(x=x,y=y,z=z,scalars=colors)#,scalars=replic)
        self.first = False
#print(x.shape,y.shape)
import _pickle as cPickle
import copy
data_folder = "../data/"
with open(data_folder+"polymer_timing.dat","rb") as f:
        lPolymers = cPickle.load(f)
P = np.array(lPolymers[0].get_replication_profile())

def time_color(time):
    replic = copy.deepcopy(P)
    replic[replic > time / 10. ] = 0
    replic[replic != 0 ] = 2
    replic[replic == 0 ] = 1
    #print(len(replic))
    return replic
    
Reprs = [Representation(traj=t,what="name Diff",tube=False),Representation(traj=t,what="not name Diff",time_color=time_color)]
for r in Reprs:
    r.draw(0)

#Load the replication information:



print(P)
i = 10
replic = copy.deepcopy(P)
replic[replic > i / 10. ] = 0
replic[replic != 0 ] = 2
replic[replic == 0 ] = 1
print(replic)
print(len(t))
i = 0
@mlab.show
@mlab.animate(delay=100)
def anim():
    """Animate the b1 box."""
    i = 20
    while 1:

        for r in Reprs:
            r.draw(i)
        i = i+ 1
        i = i % len(t)
        yield
        
# Run the animation.
anim()