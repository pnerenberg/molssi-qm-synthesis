import psi4
import numpy as np

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


# Run vibrational frequency analysis
psi4.set_output_file(file_prefix + '_vibfreq.dat', False)
scf_energy, scf_wfn = psi4.frequency('scf/cc-pVDZ', molecule=ch4, return_wfn=True, dertype='gradient')

# Save "raw" frequencies into a variable
omegaEntry = scf_wfn.frequency_analysis["omega"].data
omegaArrayReal = np.real(omegaEntry)

#only keep real entries with a real value of greater than 10
clipped_omegaArray = omegaArrayReal[omegaArrayReal>10]

#round each entry to an integer
rounded_omegaArray = np.rint(clipped_omegaArray)

#create an array with only unique frequencies
unique_omegaArray = np.unique(rounded_omegaArray)

#define an empty two-column array
freqValArray = np.empty((0,2), int)
for i in range(len(unique_omegaArray)):
    freq = np.count_nonzero(rounded_omegaArray == unique_omegaArray[i])
    freqValArray = np.append(freqValArray, np.array([[freq, unique_omegaArray[i]]]), axis=0)

#Print data into text file called frequency_list.txt
print("List of Frequencies for %s \n\n" %(file_prefix), file=open("frequency_list.txt", "w"))

for i in range(len(unique_omegaArray)):
    print("%s:   Degeneracy: %s   Frequency: %s" %(i, freqValArray[i,0],freqValArray[i,1]), file=open("frequency_list.txt", "a"))
