# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 10:56:46 2017
This script calculates the dynamic network from trajectory 
@author: ashutosh
"""

import os
import subprocess
import shutil
from distutils import spawn
from Bio.PDB.PDBParser import PDBParser # for parsing PDB file
#import command_args
def find_contact(fname):
    if os.path.isfile(fname)==True:                                    #checks if the name saved in file_list is a file
        print fname
        parser=PDBParser(PERMISSIVE=1)                                     #parser for PDB file
        temp_struct=parser.get_structure(fname[0:4], fname)                # parsing of PDB file
        contact_count=0
        contacts=[]
        for atom1 in temp_struct.get_atoms():
            for atom2 in temp_struct.get_atoms():
                resid1=atom1.get_parent()
                resid2=atom2.get_parent()
                if (2.5 <atom1-atom2 < 4.5) and (abs(int(resid1.get_id()[1])-resid2.get_id()[1])>2):
                    contacts.append((atom1.get_id(),resid1.get_id()[1],atom2.get_id(),resid2.get_id()[1])) 
                    contact_count=contact_count+1
    #print(contact_count)
    return contacts 
    
def extract_frame(xtc,tpr,ndx,i):
    print("Trajectory file read: %s",xtc)
    print("Structure file read: %s",tpr)
    print("Index file read: %s",ndx)    
    choice=open('choice.txt','r')                            # File with option to write the frames
    fr=str(i)                                                #Frame number to be dumped converted to string for use within subprocess
    fname='frame'+fr+'.pdb'                                  #PDB file in which frame is dumped
    sout=open(os.devnull,'w')                                #dumping standard output from gromacs
    if spawn.find_executable("gmx"):
        subprocess.call(['gmx','trjconv', '-f', xtc, '-s', tpr, '-dump', fr,'-o', fname, '-n', ndx],stdin=choice, stdout=sout,stderr=subprocess.STDOUT)
    choice.close()
    print("Frame %s dumped in file %s",fr,fname)
    return fname