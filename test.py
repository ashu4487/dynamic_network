# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 12:59:53 2017

@author: ashutosh
"""
import numpy as np
import matplotlib.pyplot as plt

contact_mat=np.loadtxt('contact_matrix.txt')
sum_contact=np.sum(contact_mat,axis=1)
#print(np.shape(sum_contact))
#print
plt.bar(range(len(sum_contact)),sorted(sum_contact,reverse=True))

#plt.imshow(contact_mat,aspect="auto",interpolation='nearest')
plt.show()