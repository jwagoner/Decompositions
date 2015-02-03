"""Manipulate numerical examples of free energy processes using cumulants"""

import sys, string, random, math, os
import numpy as np


def getPartitionFunction(energies,lambdas):
	"""Return Z given arrays of energys and a set of lambda values"""

	# lambdas shape should be 1 x ncomponents
	# energies shape should be ncomponents x nstates
	expEnergies = np.exp(-1*np.dot(lambdas,energies))
	Z = expEnergies.sum()
	return Z

def bradyTI(energies,spectators,scaled,bradycomp,npts):
	"""Perform Brady-like TI on the component indexed by bradies

	Args:
	energies:  ncomponents x nstates
	spectators:  1 x ncomponents.  set to 0 for non participators or
		components that are scaled.  set to 1 for spectators
	scaled: 1 x ncomponents.  set to 0 for non participators and
		spectators.  set to 1 for all components that are scaled
	bradycomp:  1 x nocomponents.  set to 1 for the component we
		are performing TI over.  0 otherwise.
	npts:  number of points for the TI numerical integration  
	"""


	gs = 1./npts
  	dHdl_sum = 0.0
	
	for i in range(npts):
		lam = i*gs

		pStates = np.exp(-1*np.dot(spectators,energies)-1*lam*np.dot(scaled,energies))
		Z = pStates.sum()
	
		pStates = (1./Z)*pStates

		# we assume a linear scaling of energies
		dHdl = np.dot(bradycomp,energies)

		pStates_t = pStates.reshape(-1,1)

		dHdl_tot = np.dot(dHdl,pStates_t)

		dHdl_sum = dHdl_sum+dHdl_tot*gs
	
	return dHdl_sum





