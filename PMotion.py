# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 09:25:03 2017

@author: jarbona
"""

class Origin:
    def __init__(self,tag):
        self.tag = tag
        self.position = self.tag
        self.move=False
        self.origin=True
        self.passivated=False
        self.activated=False
    def __repr__(self):
        add = ""
        if self.passivated:
            add="Passi"
        if self.activated:
            add="Activ"
        return "Origin %i %s"%(self.position,add)

class Fork:
    def __init__(self,tag,position,bond_tag):
        self.tag = tag  #particle tagg
        self.position = position #Position on the chromosome
        self.bond_tag = bond_tag #Bond between diffusive elements and monomec
        self.move=True
        self.update_bond=False
        self.origin = False

    def update_position(self,dt):
        #print(self.position)
        oldp = int(self.position)
        self.position += self.d * dt    
        if oldp != int(self.position):
            self.update_bond =  True
        else:
            self.update_bond =  False
            
    def __repr__(self):
        add = "L"
        if self.d == 1:
            add="R"
        return "%sFork %i"%(add,self.position)
        
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
        self.modules.insert(i+1,RFork(ptags[1],otag,new_btag[1]))
        self.modules.insert(i,LFork(ptags[0],otag,new_btag[0]))
        self.modules[i+1].passivated=True
        self.modules[i+1].activated=True

        #print(self.modules)
        
    def increment_time(self,dt,verbose=False):
        update_bond = []
        alone = []
        to_release = []
        diff_diff = []
        bind_diff = []
        passivated_origin = []
        if verbose and self.modules != []:
            print(self.start,self.end)
            print(self.modules)
        for m in self.modules:
            if m.move:
                m.update_position(dt)
        N_mod = len(self.modules)
        im = 0
        to_remove = []
        while im < N_mod:
            #Take care of leaving origins
            m = self.modules[im]
            if m.move:
                if m.position < self.start or m.position > self.end:
                    alone.append(m.tag)
                    to_release.append(m.bond_tag)
                    to_remove.append(im)
                    m.update_bond = False # because we will delete it
            im += 1
        im = 0 
        while im < N_mod:   
        #Take care of passivated Origin Left to right
            m = self.modules[im]
            if m.move:
                if im != N_mod -1 and m.position > self.modules[im + 1].position:
                    if self.modules[im + 1].origin:
                        passivated_origin.append(self.modules[im + 1].tag)
                        self.modules[im + 1].passivated = True
                        a,b = self.modules[im:im+2]
                        self.modules[im:im+2] = b,a
                        im += 1
                            
            im += 1
                
        im = N_mod-1
        while im > 0:   
        #Take care of passivated Origin Right to left
            m = self.modules[im]
            if m.move:
                if im != 0 and m.position < self.modules[im - 1].position:
                    if self.modules[im - 1].origin:
                        passivated_origin.append(self.modules[im - 1].tag)
                        self.modules[im - 1].passivated = True
                        a,b = self.modules[im-1:im+1]
                        self.modules[im-1:im+1] = b,a
                        im -= 1
                            
            im -= 1
                
                
        im = 0 
        while im < N_mod: 
            #Take care of fork collision and bond motion
            m = self.modules[im]
            if m.move:
                if im != N_mod -1 and m.position > self.modules[im + 1].position:
                    if self.modules[im + 1].move:
                        #Collision
                        to_release.append(m.bond_tag)
                        to_release.append(self.modules[im + 1].bond_tag)
                        diff_diff.append(m.tag)
                        diff_diff.append(self.modules[im + 1].tag)
                        bind_diff.append([m.tag,self.modules[im + 1].tag])
                        to_remove.append(im)
                        to_remove.append(im + 1)
                        im += 1
                elif m.update_bond:
                    update_bond.append([m.bond_tag,int(m.position)])
            im += 1
        for im in to_remove[::-1]:
            m = self.modules.pop(im)
            assert(m.move)
        #chek for colisions:
        #for m in  
        if verbose and self.modules != []:
            print(self.modules)   
        return bind_diff, diff_diff,update_bond,passivated_origin,to_release,alone
if __name__ == "__main__":
    P = Polymer(0,30,[5,10,20])
    P.add_fork([0,1],10,["a","b"])
    for i in range(11):
        print(P.increment_time(1,verbose=True))
    