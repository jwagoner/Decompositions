#!/usr/bin/python

import sys, string, math, os
import numpy as np
from numericalCumulants import *


def getCBRenergies(nC,nStates,bGoodBinder):
	"""Populate arrays of energies for a CBR-like process

	We will assume a ligand, solvent, and 2 residues (0 and 1).  
 	The extension to a large protein is easy because residue 0 can
	actually be assumed to contain however many components.  
	Assuming that nC = 5, the components are as follows:
	0:  Res0-solvent (+ whatever else we want as reference)
	1:  Res0-ligand
	2:  ligand-solvent
	3:  Res1-solvent
	4:  Res1-ligand
	We will break nStates in half.  The first half corresponds to
	solvent in the binding site.  The second half is solvent evacuated
	from the binding site.  The ligand is always in the binding site,
	ala the CBR formalism
	"""

	nS0 = nStates/2
	nS1 = nStates-nS0

	# Energy distibutions for these 2 state decompositions
	# vv good binder
	if (bGoodBinder): mu0 = np.array( [-3.0 , -12.0 , 60.0 , -4.0 , -14.0] )
	# vv Res1 is a bad binder compared to the Res1-solvent
        else: mu0 = np.array( [-3.0 , -12.0 , 60.0 , -17.0 , -14.0] )
	sigma0 = np.array( [1.5 , 1.5 , 4.0 , 1.5 , 1.5] )
	""" Note the differences.  The below state lowers the residue-solvent
		interactions.  It also removes the huge ligand-solvent
		penalty since the solvent is evacuated
	"""
        mu1 = np.array( [-1.0 , -9.0 , -2.5 , -1.0 , -14.0] )
        sigma1 = np.array( [1.5 , 1.5 , 1.0 , 1.5 , 1.5] )


	energies = np.zeros( (nC,nStates) )
	for i in range(nC):
		
		e0 = np.random.normal(mu0[i],sigma0[i],nS0)
                e1 = np.random.normal(mu1[i],sigma1[i],nS1)
		energies[i] = np.hstack( (e0,e1) )

	return energies


def main():
	

	nC = 5 # 5 components
	nStates = 500
	bGoodBinder = False

	energies = getCBRenergies(nC,nStates,bGoodBinder)

	# Let's get G_tot-G_03 (free energy of binding)
	Ltot = np.array( [np.ones(nC)] )
	Ztot = getPartitionFunction(energies,Ltot)
	L03 = np.array( [ [1.0 , 0.0 , 0.0 , 1.0 , 0.0 ] ] )
	Z03 = getPartitionFunction(energies,L03)
	Gtot_G03 = -math.log(Ztot/Z03)

        # Let's get G_012-G_0 (free energy of binding
	#	without residue 1)
        L012 = np.array( [[1.0 , 1.0 , 1.0 , 0.0 , 0.0 ]] )
        Z012 = getPartitionFunction(energies,L012)
        L0 = np.array( [[1.0 , 0.0 , 0.0 , 0.0 , 0.0 ]] )
        Z0 = getPartitionFunction(energies,L0)
        G012_G0 = -math.log(Z012/Z0)

	# Component 1
        L234 = np.array( [[0.0 , 0.0 , 1.0 , 1.0 , 1.0 ]] )
        Z234 = getPartitionFunction(energies,L234)
        L3 = np.array( [[0.0 , 0.0 , 0.0 , 1.0 , 0.0 ]] )
        Z3 = getPartitionFunction(energies,L3)
        G234_G3 = -math.log(Z234/Z3)

	# Other endpoints
        L0123 = np.array( [[1.0 , 1.0 , 1.0 , 1.0 , 0.0 ]] )
        Z0123 = getPartitionFunction(energies,L0123)
        L0134 = np.array( [[1.0 , 1.0 , 0.0 , 1.0 , 1.0 ]] )
        Z0134 = getPartitionFunction(energies,L0134)
        L0234 = np.array( [[1.0 , 0.0 , 1.0 , 1.0 , 1.0 ]] )
        Z0234 = getPartitionFunction(energies,L0234)
        
	Gtot = -math.log(Ztot)
	G0123 = -math.log(Z0123)
	G0134 = -math.log(Z0134)
	G0234 = -math.log(Z0234)

	# ddG
	ddG = Gtot_G03-G012_G0

	# new the brady-like component for residue 1
	spectators = np.array( [[1.0 , 0.0 , 0.0 , 1.0 , 0.0 ]] )
        scaled = np.array( [[0.0 , 1.0 , 1.0 , 0.0 , 1.0 ]] )
        bradycomp = np.array( [[0.0 , 0.0 , 0.0 , 0.0 , 1.0 ]] )
	dG4 = bradyTI(energies,spectators,scaled,bradycomp,1000)

        bradycomp = np.array( [[0.0 , 0.0 , 1.0 , 0.0 , 0.0 ]] )
        dG2 = bradyTI(energies,spectators,scaled,bradycomp,1000)

        bradycomp = np.array( [[0.0 , 1.0 , 0.0 , 0.0 , 0.0 ]] )
        dG1 = bradyTI(energies,spectators,scaled,bradycomp,1000)

	print "COMPARISONS for two components (ddG ; Brady ; dG)"
	print "C4:",ddG,";",dG4,";",Gtot-G0123
	print "C2:",Gtot_G03-G234_G3,";",dG1,";",Gtot-G0234
	print "Total binding free energy",Gtot_G03
	print "TI components",dG1,dG2,dG4

if __name__ == "__main__":
        main()



