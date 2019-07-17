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
scf_energy, scf_wfn = psi4.frequency('scf/cc-pVDZ', molecule=ch4, return_wfn=True, dertype='gradient')

# Save "raw" frequencies into a variable
ifrequency=scf_wfn.frequency_analysis.get("omega").data

# Eliminate non-zero frequency.  NOTE, All must be in order
frequency=np.round(ifrequency.real)[6:]

# Condense the frequency lists to elimintat, but count duplicates
tabfreq=np.array(np.unique(frequency,return_counts=True)).T

#TODO one should and could try to pull more information, such as intensities or motions.

#Print to file
np.savetxt(file_prefix + '_vibfreq.dat',tabfreq,'%d','\t\t\t',header="Freq (cm-1)\t\tDegen",footer="\n# Molecule's energy ="+str(scf_energy))
