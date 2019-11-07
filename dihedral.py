#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

# USAGE:
# ./rgyr.py system_descripter pdb_file traj_file 

# calculate radius of gyration for a solute;  

# PREAMBLE:

import numpy as np
import MDAnalysis
import sys
import math

def computePbcDist2(r1,r2,box):
        dist2 = 0

        for j in range(0,3):
                temp = r1[j]-r2[j]
                if temp < -box[j]/2.0:
                        temp += box[j]
                elif temp > box[j]/2.0:
                        temp -= box[j]
                dist2 += temp*temp

        return dist2;

def computeDihedral(r1, r2, r3, r4, box):
	vec1=np.zeros(3, dtype=np.float)
	vec2=np.zeros(3, dtype=np.float)
	vec3=np.zeros(3, dtype=np.float)
 
	for j in range(0,3):
		# Define the 3 vectors within the dihedral
		vec1[j] = r2[j]-r1[j]
		vec2[j] = r3[j]-r2[j]
		vec3[j] = r4[j]-r3[j]

		# Pbc
		if vec1[j] < -box[j]/2.0:
                        vec1[j] += box[j]
		elif vec1[j] > box[j]/2.0:
                        vec1[j] -= box[j]
		if vec2[j] < -box[j]/2.0:
                        vec2[j] += box[j]
		elif vec2[j] > box[j]/2.0:
                        vec2[j] -= box[j]
		if vec3[j] < -box[j]/2.0:
                        vec3[j] += box[j]
		elif vec3[j] > box[j]/2.0:
                        vec3[j] -= box[j]

	# Calculate angle between 2 planes
	n1_v = np.cross(vec1, vec2)
	n2_v = np.cross(vec2, vec3)
	m_v = np.cross(n1_v, vec2)
	
	n1_norm = n1_v/np.linalg.norm(n1_v)
	n2_norm = n2_v/np.linalg.norm(n2_v)
	m_norm = m_v/np.linalg.norm(m_v)
	
	x = np.dot(n1_norm, n2_norm)
	y = np.dot(m_norm, n2_norm)
	dih = -np.arctan2(y,x) 

	return dih;

system = sys.argv[1]
pdb = sys.argv[2]
traj = sys.argv[3]

res_list = [2,3,4,5,6,12,11,10,9,8,15,16,17,18,19,25,24,23,22,21]
omega_sel=[]
phi_sel=[]
psi_sel=[]


u = MDAnalysis.Universe('%s' %(pdb), '%s' %(traj))

for k in range(len(res_list)):
	if 1 <= res_list[k] <= 5:
		i = res_list[k]
		j = i - 1
		l = i + 1
	if 9 <= res_list[k] <= 13:
		i = res_list[k]
		j = i + 1
		l = i - 1		
	if 14 <= res_list[k] <= 18:
		i = res_list[k]
		j = i - 1
		l = i + 1
	if 22 <= res_list[k] <= 26:
		i = res_list[k]
		j = i + 1
		l = i - 1
	
	if i == 6 or i == 19:
		omega_sel1 = []
		omega_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
		omega_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
		omega_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
		omega_sel1.append(u.select_atoms("resid %s and name CA" %(j)))  
		omega_sel.append(omega_sel1)
		
		phi_sel1 = []
		phi_sel1.append(u.select_atoms("resid %s and name CA" %(l)))  
		phi_sel1.append(u.select_atoms("resid %s and name N" %(i)))  
		phi_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
		phi_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
		phi_sel.append(phi_sel1)
		
		psi_sel1 = []
		psi_sel1.append(u.select_atoms("resid %s and name N" %(i)))  
		psi_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
		psi_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
		psi_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
		psi_sel.append(psi_sel1)

	if i == 8 or i == 21:
		omega_sel1 = []
		omega_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
		omega_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
		omega_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
		omega_sel1.append(u.select_atoms("resid %s and name CA" %(j)))  
		omega_sel.append(omega_sel1)
		
		phi_sel1 = []
		phi_sel1.append(u.select_atoms("resid %s and name CB" %(l)))  
		phi_sel1.append(u.select_atoms("resid %s and name N" %(i)))  
		phi_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
		phi_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
		phi_sel.append(phi_sel1)
		
		psi_sel1 = []
		psi_sel1.append(u.select_atoms("resid %s and name N" %(i)))  
		psi_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
		psi_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
		psi_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
		psi_sel.append(psi_sel1)
		
	else:
		omega_sel1 = []
		omega_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
		omega_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
		omega_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
		omega_sel1.append(u.select_atoms("resid %s and name CA" %(j)))  
		omega_sel.append(omega_sel1)
		
		phi_sel1 = []
		phi_sel1.append(u.select_atoms("resid %s and name C" %(l)))  
		phi_sel1.append(u.select_atoms("resid %s and name N" %(i)))  
		phi_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
		phi_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
		phi_sel.append(phi_sel1)
		
		psi_sel1 = []
		psi_sel1.append(u.select_atoms("resid %s and name N" %(i)))  
		psi_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
		psi_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
		psi_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
		psi_sel.append(psi_sel1)


	print k,phi_sel[k][0],phi_sel[k][1],phi_sel[k][2],phi_sel[k][3]
	print k,psi_sel[k][0],psi_sel[k][1],psi_sel[k][2],psi_sel[k][3]
	print k,omega_sel[k][0],omega_sel[k][1],omega_sel[k][2],omega_sel[k][3]
#
#	if res_list[k] =

#	phi_sel.append(u.select_atoms("resid %s and name C" %(i)), u.select_atoms("resid %s and name CA" %(i)), u.select_atoms("resid %s and name N" %(i)), u.select_atoms("resid %s and name C" %(j)))
#	psi_sel.append(u.select_atoms("resid %s and name N" %(i)), u.select_atoms("resid %s and name C" %(j)), u.select_atoms("resid %s and name CA" %(j)), u.select_atoms("resid %s and name N" %(j)))
	
sel1=u.select_atoms('resid 1:13 and resname PDI')
sel2=u.select_atoms('resid 14:26 and resname PDI')
sel3=u.select_atoms('resid 11').residues[0].phi_selection()
sel4=u.select_atoms('resid 11').residues[0].psi_selection()
sel5=u.select_atoms('resid 11').residues[0].omega_selection()
#print sel3, sel4, sel5
# SUBROUTINES:

# MAIN PROGRAM:

out = open('%s.dihedral.dat' %(system), 'w')

phi_data = np.zeros(len(res_list), dtype=np.float)
psi_data = np.zeros(len(res_list), dtype=np.float)
omega_data = np.zeros(len(res_list), dtype=np.float)

for ts in u.trajectory:
	dist = math.sqrt(computePbcDist2(sel1.center_of_mass(),sel2.center_of_mass(),u.dimensions[:3]))
	out.write('%10.6f' %(dist))
	for i in range(len(res_list)):
		#print i,phi_sel[i][0],phi_sel[i][1],phi_sel[i][2],phi_sel[i][3]
		phi_data[i]=computeDihedral(phi_sel[i][0].center_of_mass(), phi_sel[i][1].center_of_mass(), phi_sel[i][2].center_of_mass(), phi_sel[i][3].center_of_mass(), u.dimensions[:3])
		omega_data[i]=computeDihedral(omega_sel[i][0].center_of_mass(), omega_sel[i][1].center_of_mass(), omega_sel[i][2].center_of_mass(), omega_sel[i][3].center_of_mass(), u.dimensions[:3])
		psi_data[i]=computeDihedral(psi_sel[i][0].center_of_mass(), psi_sel[i][1].center_of_mass(), psi_sel[i][2].center_of_mass(), psi_sel[i][3].center_of_mass(), u.dimensions[:3])

		out.write(' %10.6f %10.6f %10.6f' %(phi_data[i]*180/np.pi, omega_data[i]*180/np.pi, psi_data[i]*180/np.pi))
		
#		if i == 5:
#			dih_phi = sel3.dihedral.value()
#			dih_psi = sel4.dihedral.value()
#			dih_omega = sel5.dihedral.value()
#			print dih_phi, phi_data[i]*180/np.pi
#			print dih_psi, psi_data[i]*180/np.pi
#			print dih_omega, omega_data[i]*180/np.pi
			
		
	out.write('\n')
out.close()

