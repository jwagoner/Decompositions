#!/usr/bin/python

import sys, string, math

import cumulant_base as cb


def cancelBradyDerivatives(myC):
	"""Demonstrate that all derivative components cancel in the ways we expect"""

        # define spectator and scaled components
        myC.compScale[0]=False
        myC.compScale[1]=True
        myC.compScale[2]=True

	bSpec = [] # labels spectators
        bInt = [] # labels the free energy component of interest
        for i in range(3):
		bSpec.append(False)
		bInt.append(False)
	bSpec[0] = True
	bInt[2] = True

	print "Calculating cumulants up to order",myC.nOrders
	myC.cumulants = cb.getCumulants(myC.nC,myC.nOrders)
	print "Done."

	# perform the Brady splittings
        for i in range(myC.nOrders):
                myC.cumulants[i] =  cb.splitCumulantM1(myC.cumulants[i],bSpec,bInt,False)

        myC.getDerivatives()
	'''assume that f(\lambda) is identical and separable 
		for all scaled components'''
	cb.separateDerivatives(myC.derivatives,myC.nOrders)

	# combine, rearrange, simplify all derivatives
	myC.fullExpansion = myC.simplifyDerivatives()

	# now let's play around and cancel all terms
	# first, separate the full expansion into separate strings organized
	# by powers of beta
	myC.betaPowers(myC.fullExpansion)

        for i in range(myC.nOrders+1):
                cbi = cb.cancelTermsString(myC.betaList[i],True)
                myC.betaListCombined[i] = cbi
                if (i < myC.nOrders): print "\n\n\n\nCOMBINED",i,myC.betaListCombined[i]


def main():

	nC = 3 # number of components
	nOrders = 4 # number of cumulant orders to test
	myC = cb.Cumulants(nC,nOrders)

	cancelBradyDerivatives(myC)
	

if __name__ == "__main__":
	main()





