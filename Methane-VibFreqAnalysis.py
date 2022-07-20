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
frequencies = scf_wfn.frequency_analysis['omega'].data # this command is just to get you started!

# Eliminate imaginary parts of frequencies,
frequencies2 = np.real(frequencies)
# round the frequencies (to the nearest whole number),
frequencies_rounded = np.round(frequencies2)
# and extract only the *non-zero* frequencies
frequencies_abovezero = frequencies_rounded[frequencies_rounded>0]
# Determine the unique non-zero frequencies and
# the number of times each such frequency occurs;
# store these in a NumPy array in the format: 
# {frequency, count} (i.e, one line per freq.)
(frequency, count)= np.unique(frequencies_abovezero, return_counts = True)
frequencies_final = np.asarray((frequency, count)).T
# Save the NumPy array with frequency and count data
# to a text file with the header line: 'freq degen'
np.savetxt('frequency degeneracies.txt', frequencies_final, header = 'freq degen')

