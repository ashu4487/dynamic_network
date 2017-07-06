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
        print(fname)
        parser=PDBParser(PERMISSIVE=1)                                     #parser for PDB file
        temp_struct=parser.get_structure(fname[0:4], fname)                # parsing of PDB file
        model=temp_struct[0]
        contact_count=0
        contacts=[]
        for residue1 in model.get_residues():                      #loop over all the residues in model    
            chain1=residue1.get_parent()                           #chain id for first loop
            if chain1.get_id()=='A' or chain1.get_id()=='B':       # residues belonging only to chain A or B processed
                #count1=count1+1
                #print count1
                #count2=0
                for residue2 in model.get_residues():
                    chain2=residue2.get_parent()  
                    if chain2.get_id()=='A' or chain2.get_id()=='B':
                        #count2=count2+1
                        
                        #print count2
                        temp_resid1=residue1.get_id()
                        temp_resid2=residue2.get_id()
                        if temp_resid1[0]==" " and temp_resid2[0]==" ": # only aa residues, excludes Water and other het atoms
                            for atom1 in residue1:
                                for atom2 in residue2:
                                    if (atom1-atom2 < 4.5) and (abs(int(residue1.get_id()[1])-residue2.get_id()[1])>2):
                                        contact=(residue1.get_id()[1],chain1.get_id(),residue2.get_id()[1],chain2.get_id())
                                        contact_rev=(residue2.get_id()[1],chain2.get_id(),residue1.get_id()[1],chain1.get_id())
                                        if contact in contacts or contact_rev in contacts :
                                            pass
                                        else:
                                            #contacts.append((atom1.get_id(),resid1.get_id()[1],atom2.get_id(),resid2.get_id()[1]))
                                            contacts.append(contact)
                                            contact_count=contact_count+1
                                    else:
                                       pass 
                        else:
                            pass
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
    else:
        subprocess.call(['trjconv', '-f', xtc, '-s', tpr, '-dump', fr,'-o', fname, '-n', ndx],stdin=choice,stderr=subprocess.STDOUT)
    choice.close()
    print("Frame %s dumped in file %s",fr,fname)
    return fname