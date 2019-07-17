#!/usr/bin/env python
# coding: utf-8

# In[1]:


import psi4
import numpy as np


# In[2]:


# Initial setup
psi4.set_memory('2 GB')
psi4.set_num_threads(2)

file_prefix = 'methane_HF-DZ'

ch4 = psi4.geometry("""
symmetry c1
0 1
   C       -0.85972        2.41258        0.00000
   H        0.21028        2.41258        0.00000
   H       -1.21638        2.69390       -0.96879
   H       -1.21639        3.11091        0.72802
   H       -1.21639        1.43293        0.24076
""")

# Geometry optimization
psi4.set_output_file(file_prefix + '_geomopt.dat', False)
psi4.set_options({'g_convergence': 'gau_tight'})
psi4.optimize('scf/cc-pVDZ', molecule=ch4)


# In[3]:


# Run vibrational frequency analysis
psi4.set_output_file(file_prefix + '_vibfreq.dat', False)
scf_energy, scf_wfn = psi4.frequency('scf/cc-pVDZ', molecule=ch4, return_wfn=True, dertype='gradient')

# Save "raw" frequencies into a variable
Freq = scf_wfn.frequency_analysis
print(Freq) # this command is just to get you started!


# In[13]:


# Eliminate imaginary parts of frequencies,
# round the frequencies (to the nearest whole number),
# and extract only the *non-zero* frequencies

nonZero = np.round(np.real(Freq["omega"][2]))[6:]
print(nonZero)


# In[26]:


# Determine the unique non-zero frequencies and 
# the number of times each such frequency occurs;
# store these in a NumPy array in the format: 
# {frequency, count} (i.e, one line per freq.)
freqCount = np.array(np.unique(nonZero, return_counts=True))
freqCount = freqCount.T    


# In[ ]:


## Save the NumPy array with frequency and count data
# to a text file

np.savetxt('freqCount.txt', freqCount, fmt='%d %d', header='Freq Degen')

