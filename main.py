# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 11:08:15 2017

@author: ashutosh
"""

import os
import sys
import pickle
import command_args
import create_network
import numpy as np
import matplotlib.pyplot as plt


#Parse the command line arguments
xtc,tpr,ndx,skip,out=command_args.parseargs()

if os.path.isfile("all_contacts.txt"):
    print("Contacts already existing in all_contacts.txt")
    with open('all_contacts.txt','r') as all_con:
        all_contacts=pickle.load(all_con)
else:
    all_contacts=set()                                      #create empty set to keep all contacts
    for i in range(0,10000,int(skip)):
        fname=create_network.extract_frame(xtc,tpr,ndx,i) #Extracts frame and saves the PDB file with fname
        contacts=create_network.find_contact(fname) # Find the contacts ib protein structure
        outf=open(fname[:-4]+'contacts.txt','w')
        pickle.dump(contacts,outf)
        outf.close()
        scon=set()
        for con in contacts:
            scon.add(con)
        all_contacts=set(all_contacts).union(scon) #Add new contacts to all contact list
        #print(list(all_contacts)[0:5])
        #print(con_count)
        #print(len(all_contacts))
        #print(sys.getsizeof(all_contacts))
        os.remove(fname)                            #removes the PDB file to save space    
    with open('all_contacts.txt','w') as all_out:
        pickle.dump(all_contacts,all_out)

    
if os.path.isfile("contact_matrix.txt"):
    print("Contact matrix already exixts in folder. Now reading it")
    contact_mat=np.loadtxt('contact_matrix.txt')
    
else:
    contact_mat=np.zeros([len(all_contacts),5])
    count=0
    for i1 in range(0,1000,int(skip)):
        with open('frame'+str(i1)+'contacts.txt','r') as rcon:
            con_read=pickle.load(rcon)
        for item in con_read:
            if item in list(all_contacts):
                index=list(all_contacts).index(item)
                contact_mat[index,count]=1
        count=count+1
    np.savetxt('contact_matrix.txt',contact_mat,fmt="%d")
    
#plt.imshow(contact_mat,aspect="auto",interpolation='nearest')
[num_cont,num_frames]=np.shape(contact_mat)
sum_contact=np.sum(contact_mat,axis=1)
cont_prob=sum_contact/num_frames    #Calculate contact probability
plt.plot(range(len(cont_prob)),sorted(cont_prob,reverse=True))
plt.show()         
