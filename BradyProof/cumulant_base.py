"""Symbolically manipulate cumulants and their derivatives

Words are stored as follows:
M2_h2d0h1d1
First two characters denote the moment order (2nd moment here)
This is followed by underscore and four character segments
Each four character segment denotes the hamiltonian component (2nd 
(character) and the derivative order (4th character).  This doesn't 
allow us to go over order 9 for moments, components, or derivatives.
Every independant string is an average quantity, taken over the 
ensemble containing all components turned on (the lambda = 1 state). 
Defining the ensemble is only necessary for derivatives. 
"""

import sys, string, math, os

# defines our cutoff for 'zero'
FSMALL = 10E-10

def getIthMoment(nC,iOrder):
	"""Return moment of arbitrary order

	Args:
	nC: number of Hamiltonian components
	iOrder: order of the moment
	"""

	nwords = int(math.pow(nC,iOrder))
	tot = ""
	for i in range(nwords):
		mystr = "M"+str(iOrder)+"_"

		for j in range(iOrder):
			jmod = int(math.pow(nC,j+1))
			jderv = int(math.pow(nC,j))
			jndx = (i%jmod)/jderv
			jword = "h"+str(jndx)+"d0"
			mystr = mystr+jword
	
		tot = tot+mystr+"+"

	tot = tot[:-1]
        return tot






# multiply two strings, throw the prefix on front
def multiplyStrings(str1,str2,mypre):
	"""Return the product of two strings and a prefactor"""

	
	words1 = []
	words2 = []
	tot = ""

	# First, break up each string into its components
	word = ""
	for i in range(len(str1)):
		ch1 = str1[i]
		if (ch1 == "+"):
			words1.append(word)
			word = ""
		else:
			word = word+ch1
	words1.append(word)
	word = ""
        for i in range(len(str2)):
                ch1 = str2[i]
                if (ch1 == "+"):
                        words2.append(word)
                        word = ""                               
                else:
                        word = word+ch1
        words2.append(word)				

	nw1 = len(words1)
	nw2 = len(words2)
	print nw1,"x",nw2,"terms"
	# Multiply components together
	for i in range(nw1):
		word1_tot = words1[i]
		# look for a prefactor, like beta or 3, etc.
		bPos = False
		ii=0
		while (not bPos):
			ch1 = word1_tot[ii]
			if (ch1 == "M"): bPos = True
			ii = ii+1
		p1 = word1_tot[:ii-1]
		word1 = word1_tot[ii-1:]

		for j in range(nw2):

			word2_tot = words2[j]
                       	bPos = False
                       	jj=0
                	while (not bPos):
                	        ch1 = word2_tot[jj]
                	        if (ch1 == "M"): bPos = True
                	        jj = jj+1
                	p2 = word2_tot[:jj-1]
                	word2 = word2_tot[jj-1:]

			if (p1 == "" and p2 == ""):
				pre = ""
			elif (p1 != "" and p2 != ""):
				pre = p1+p2
			elif (p1 == ""):
				pre = p2
			else: pre = p1

			mterm = mypre+pre+word1+"*"+word2
			tot = tot+mterm+"+"
	# leave off the last "+" sign
	return tot[:-1]

def getCumulants(nC,nOrders):
        """Return list of cumulants indexed by order"""

        moments = []
        cumulants = []
	fac = 1
        for i in range(nOrders):
		# get the ith-moment
                moments.append(getIthMoment(nC,i+1))
		# put on the appropriate prefactor
                fac = fac*(i+1)
		prefactor = "[1/"+str(fac)+"]*beta"+str(i+1)
		moments[i] = addPrefactor(prefactor,moments[i])		

	# combine all moments into one big string
	momentExpansion = combineStrings(moments,nOrders)
	# now, let taylor expand the logarithm
	# put each term into a gigant string
	# the first term is just the moment expansion
	cumTermI = momentExpansion
	cumTermI = simplifyPrefactors(cumTermI)
	cumTermI = arrangeWords(cumTermI)
	cumTermI = cancelTermsString(cumTermI,False)
	momentExpansion = cumTermI
	cumulantExpansion = cumTermI
	for i in range(nOrders-1):
		print "Expanding the logarithm, term",i+1
		# i = 0 corresponds to n = 2 in the logarithm expansion
		cumTermI = multiplyStrings(cumTermI,momentExpansion,"")
		cumTermI = cutTerms(cumTermI,nOrders)
		cumTermI = simplifyPrefactors(cumTermI)
 		cumTermI = arrangeWords(cumTermI)
		cumTermI = cancelTermsString(cumTermI,False)
		if (i%2==0): prefactor = "[-1/"+str(i+2)+"]"
		else: prefactor = "[1/"+str(i+2)+"]"
		cumTermTemp = addPrefactor(prefactor,cumTermI)
		cumulantExpansion = cumulantExpansion+"+"+cumTermTemp 

	# divide through by beta and clean up the prefactors
	cumulantExpansion = addPrefactor("beta-1",cumulantExpansion)
 	cumulantExpansion = simplifyPrefactors(cumulantExpansion)

	for i in range(nOrders):
		# look through the string for the ith cumulant
		cumulants.append(separateCumulantFromString(cumulantExpansion,i+1))

        return cumulants


def separateCumulantFromString(str1,iOrder):
        words1 = []
        tot = ""

        # First, break up each string into its components
        word = ""
        for i in range(len(str1)):
                ch1 = str1[i]
                if (ch1 == "+"):
                        words1.append(word)
                        word = ""
                else:
                        word = word+ch1
        words1.append(word)

        nw1 = len(words1)
        for i in range(nw1):
                w1 = words1[i]
		tn = 0
		for j in range(len(w1)):
			ch1 = w1[j]
			if (ch1 == "M"):
				tn = tn+int(w1[j+1])

		if (tn == iOrder):
                	tot = tot+w1+"+"
        tot = tot[:-1]
        return tot

def cutTerms(str1,nOrder):
	"""Cut all terms more than the given order

	Args:
	str1: string we're cutting
	nOrder: highest allowed order. 
	"""

        words1 = []
        tot = ""

        # First, break up each string into its components
        word = ""
        for i in range(len(str1)):
                ch1 = str1[i]
                if (ch1 == "+"):
                        words1.append(word)
                        word = ""
                else:
                        word = word+ch1
        words1.append(word)

        nw1 = len(words1)
        for i in range(nw1):
                w1 = words1[i]
                tn = 0
                for j in range(len(w1)):
                        ch1 = w1[j]
                        if (ch1 == "M"):
                                tn = tn+int(w1[j+1])

                if (tn <= nOrder):
                        tot = tot+w1+"+"
        tot = tot[:-1]
        return tot




def addPrefactor(pre,str1):
        words1 = []
        tot = ""

        # First, break up each string into its components
        word = ""
        for i in range(len(str1)):
                ch1 = str1[i]
                if (ch1 == "+"):
                        words1.append(word)
                        word = ""
                else:
                        word = word+ch1
        words1.append(word)

        nw1 = len(words1)
        for i in range(nw1):
                word1_temp = words1[i]
		newword = pre+"*"+word1_temp
		tot = tot+newword+"+"
	tot = tot[:-1]
	return tot

		


def getFreeEnergy(nOrders,cumulants):
	# all we do here is add the appropriate prefactor
	# c1 can stay as-is
	c2pre = "[1/2]*beta"
	c3pre = "[1/6]*beta*beta"
	c4pre = "[1/24]*beta*beta*beta"
        c5pre = "[1/120]*beta*beta*beta*beta"
	for i in range(1,nOrders):
		if (i == 1):
			c2 = addPrefactor(c2pre,cumulants[1])
			cumulants[1] = c2
		elif (i == 2):
                	c3 = addPrefactor(c3pre,cumulants[2])
                	cumulants[2] = c3
		elif (i == 3):
			c4 = addPrefactor(c4pre,cumulants[3])
                        cumulants[3] = c4
                elif (i == 4):
                        c5 = addPrefactor(c5pre,cumulants[4])
                        cumulants[4] = c5
		else: 
			print "getFreeEnergy() not ready for order",i
			sys.exit()
	return cumulants

def simplifyFraction(frac1):

        num = ""
        den = ""
        bDen = False
	if (frac1[0] != "[" or frac1[-1] != "]"):
		emsg = "simplfyFraction: why are you feeding me %s?" %(frac1)
		sys.exit(emsg)

        for kk in range(1,len(frac1)-1):
                ch1 = frac1[kk]
                if (ch1 == "/"):
                        bDen = True
                        continue
                if (bDen):  den = den+ch1
                else: num = num+ch1

        if (num == 0):
                emsg = "simplifyFraction: Don't feed me zeros, bro"
                sys.exit(emsg)

	num = int(num)
	den = int(den)

        if (num < 0): bNeg = True
        else: bNeg = False

        anum = abs(num)
        aden = abs(den)

        # get the greatest common divisor, euclidean algorithm
        while (anum > 0):
                tmp = anum;
                anum = aden % anum;
                aden = tmp;
        # aden is now the gcd

        # check
        if (num%aden != 0 or den%aden != 0):
                emsg = "Broken Euclidean algorithm!"
                sys.exit(emsg)

        new_num = num/aden
        new_den = den/aden

        #if (bNeg): p1 = "[-"
        #else: 
        p1 = "["
        p1 = p1+str(new_num)+"/"+str(new_den)+"]"

        return p1


# This function combines numbers, powers of beta, and separable derivatives
# It does not touch derivatives that have not yet been separated from 
# the Hamiltonian

def simplifySinglePrefactor(p1):
	Nneg = 0
	Nbeta = 0
	numbers = []
	Nnumbers = 0
	Drvs = []
	NDrvs = 0
	drvOrder = 5 # below, we increase this if we find higher orders
	oldOrder = drvOrder
	DCounts = []
	for i in range(drvOrder):  DCounts.append(0)	

        # First, break up each prefactor into its terms
        term = ""
	terms = []
        for i in range(len(p1)):
               ch1 = p1[i]
               if (ch1 == "*"):
                        terms.append(term)
                        term = ""
               else:
                        term = term+ch1
        if (term != ""): terms.append(term)

	nterms = len(terms)
	# look through each term
	for i in range(nterms):
		t1 = terms[i]
		if (t1 == "[-1]"): Nneg = Nneg+1
		elif (t1 == "beta"): Nbeta = Nbeta+1
                elif (len(t1) > 4 and t1[:4] == "beta"): Nbeta = Nbeta+int(t1[4:])
		elif (t1[0] == "D"):
			Drvs.append(t1)
			NDrvs = NDrvs+1
			if (int(t1[1]) >= drvOrder):
				drvOrder = int(t1[1])+1
				for i in range(oldOrder,drvOrder):
					DCounts.append(0)
				oldOrder = drvOrder
		else:
			numbers.append(t1)
			Nnumbers = Nnumbers+1

	# how many betas?
	if (Nbeta == 0): bterm = ""
	elif (Nbeta < 0): sys.exit("Error: negative beta term")
	else: bterm = "beta"+str(Nbeta)+"*"
		
	# and the multiplicative factor?
	numerator = 1
	denominator = 1

	for i in range(Nnumbers):	
		num = numbers[i]
		bDen = False
		s_num = ""
		s_den = ""
		for ii in range(1,len(num)-1):
			ch1 = num[ii]
			if (ch1 == "/"): 
				bDen = True
				continue
			if (bDen):  s_den = s_den+ch1
			else: s_num = s_num+ch1
		numerator = numerator*int(s_num)
		denominator = denominator*int(s_den)

        # are we negative?
	if (numerator < 0): Nneg = Nneg+1
	numerator = abs(numerator)
        if ((Nneg%2)==0): bNegative = False
        else: bNegative = True


	if (Nnumbers == 0 and bNegative):
		newMult = "[-1/1]*"
	elif (Nnumbers == 0): newMult = ""
	elif (bNegative):
		newMult = "[-"+str(numerator)+"/"+str(denominator)+"]"
	else:
		newMult = "["+str(numerator)+"/"+str(denominator)+"]"
	# simplify with the Euclidean algorithm
	if (newMult != ""): 
		newMult = simplifyFraction(newMult)+"*"

	dterm = ""
	for i in range(NDrvs):
		mydrv = Drvs[i]
		order = int(mydrv[1])
		oval = DCounts[order]
		nval = oval+1
		DCounts[order] = nval
	for i in range(1,drvOrder):
		mcount = DCounts[i]
		if (mcount > 0):
			mt = "D"+str(i)+"_"+str(mcount)+"*"
			dterm = dterm+mt

	pre = newMult+bterm+dterm
	#print "OLD",p1,"NEW",pre
	return pre

def simplifyPrefactors(str1):
        words1 = []
        tot = ""

        # First, break up each string into its components
        word = ""
        for i in range(len(str1)):
                ch1 = str1[i]
                if (ch1 == "+"):

			# let's process our word
			# look for the prefactor
                        bPos = False
                       	ii=0
                       	while (not bPos):
                               	ch1 = word[ii]
                        	if (ch1 == "M"): bPos = True
                        	ii = ii+1
                	p1 = word[:ii-1]
			if (p1 != ""): p1 = p1[:-1]
                	suffix = word[ii-1:]
			newp1 = simplifySinglePrefactor(p1)
			newword = newp1+suffix
			tot = tot+newword+"+"

                        word = ""
                else:
                        word = word+ch1


        bPos = False
        ii=0
        while (not bPos):
                ch1 = word[ii]
                if (ch1 == "M"): bPos = True
                ii = ii+1
        p1 = word[:ii-1]
        if (p1 != ""): p1 = p1[:-1]
        suffix = word[ii-1:]
        newp1 = simplifySinglePrefactor(p1)
        newword = newp1+suffix
        tot = tot+newword

	return tot



def combineStrings(strings,num):

	tot = ""
	for i in range(num):
		tot = tot+strings[i]+"+"
	tot = tot[:-1]	
	return tot


def arrangeSingleWord(word):
        tot = []
        nC = int(word[1]) # the total number of components being averaged
        components = []
        for i in range(nC):
                comp = word[3+4*i:7+4*i]
                components.append(comp)		
		
	# order the components using an N! sort. 
	# need to rewrite this if we go to high orders
		
	crossed = []
	ncrossed = 0

	newWord = word[:3]
	for i in range(nC):
		# we're looking for rank i			
		iword = ""
		ndx = ""
		for j in range(nC):
				
			# check to see if we've cross j off the list
			bCrossed = False
			for k in range(ncrossed):			
				if (j == crossed[k]): 
					bCrossed = True
					break
			if (bCrossed): continue

			jword = components[j]			

			if (iword == ""): 
				iword = jword	
				i_d = iword[3]
				i_c = iword[1]
				ndx = j
				continue
			# first, derivative number
			j_d = jword[3]
			j_c = jword[1]

			if (j_d > i_d): continue
			if (j_d < i_d): 
				ndx = j
				iword = jword
				i_d = j_d
				i_c = j_c
				continue
			if (j_d == i_d):
				# now, sort by component number
				if (j_c >= i_c): continue
				if (j_c < i_c):
					i_d = j_d
					i_c = j_c
					iword = jword
					ndx = j
		newWord = newWord+iword
		ncrossed = ncrossed+1
		crossed.append(int(ndx))
	return newWord	

# pretty straightforward.  Prue the digits from 'word', make it into a number
def getRank(word):
	iC = int(word[1])
	
	snum = str(iC)
	for i in range(iC*2):
		digit = word[4+i*2]
		snum = snum+digit
	return int(snum)	

def arrangeWordProduct(word2_list):

	nw = len(word2_list)
	
	crossed = []
	ncrossed = 0
	tot = ""
	
	for i in range(nw):
		# i is the rank

		iword = ""
		ndx = ""
		inum = ""
		# look through the list of words
		for j in range(nw):

                        # check to see if we've cross j off the list
                        bCrossed = False
                        for k in range(ncrossed):
                                if (j == crossed[k]):
                                        bCrossed = True
                                        break
                        if (bCrossed): continue

			jword = word2_list[j]
			jnum = getRank(jword)				

			if (iword == ""):
				iword = jword
				inum = jnum
				ndx = j

			if (jnum < inum):
				inum = jnum
				iword = jword
				ndx = j

                ncrossed = ncrossed+1
                crossed.append(int(ndx))

		tot = tot+iword+"*"
	tot = tot[:-1]
	return tot

def arrangeWords(str1):
        words1 = []
        tot = ""

        # First, break up each string into its components
        word = ""
        for i in range(len(str1)):
                ch1 = str1[i]
                if (ch1 == "+"):
                        words1.append(word)
			word = ""
                else:
                        word = word+ch1
        words1.append(word)

        nw1 = len(words1)
        for i in range(nw1):
                word1_temp = words1[i]
                # second, look for a prefactor, like beta or 3, etc.
                bPos = False
                ii=0
                while (not bPos):
                        ch1 = word1_temp[ii]
                        if (ch1 == "M"): bPos = True
                        ii = ii+1
                p1 = word1_temp[:ii-1]
                word1_tot = word1_temp[ii-1:]
                # third step is to break up word1_tot into multiplicative 
                                # components
                w1 = ""
                word1_list = []
                for ii in range(len(word1_tot)):
                        ch1 = word1_tot[ii]
                        if (ch1 == "*"):
                                word1_list.append(w1)
                                w1 = ""
                        else: w1 = w1+ch1

                word1_list.append(w1)

                nlist = len(word1_list)

		word2_list = []
                for j in range(nlist):
                        jword = word1_list[j]
			# Fix this single word
			new_jword = arrangeSingleWord(jword)
			word2_list.append(new_jword)

		# now, we need to order word2_list.  this is harder
		newBigWord = arrangeWordProduct(word2_list)
		tot = tot+p1+newBigWord+"+"
	tot = tot[:-1]
	return tot

def splitCumulantM1(str1,bSpectator,bInterest,bSpectatorSplitting):
        words1 = []
        tot = ""

        # First, break up each string into its components
        word = ""
        for i in range(len(str1)):
                ch1 = str1[i]
                if (ch1 == "+"):
                        words1.append(word)
                        word = ""
                else:
                        word = word+ch1
        words1.append(word)

        nw1 = len(words1)
        for i in range(nw1):
                word1_temp = words1[i]
                # second, look for a prefactor, like beta or 3, etc.
                bPos = False
                ii=0
                while (not bPos):
                        ch1 = word1_temp[ii]
                        if (ch1 == "M"): bPos = True
                        ii = ii+1
                p1 = word1_temp[:ii-1]
                word1_tot = word1_temp[ii-1:]

                # third step is to break up word1_tot into multiplicative 
                        # components

                w1 = ""
                word1_list = []
                for ii in range(len(word1_tot)):
                        ch1 = word1_tot[ii]
                        if (ch1 == "*"):
                                word1_list.append(w1)
                                w1 = ""
                        else: w1 = w1+ch1

                word1_list.append(w1)
                nlist = len(word1_list)

		den = 0
		num = 0
		bAuto = False
                for j in range(nlist):
                        jword = word1_list[j]
			jn = int(jword[1])
			# don't include this word at all if has any spectators
			# Otherwise, take a simple fraction where the numerator 
			# is the # component of interest and the denominator
			# is the total non-spectator components in this word
			for jj in range(jn):
				jnum = int(jword[4+4*jj])
				if (bSpectator[jnum]): 
					bAuto = True
					break	
				else:
					den = den+1
				if (bInterest[jnum]):
					num = num+1
			if (bAuto): break
		if (bAuto and bSpectatorSplitting):
			newP = "[1/1]*"
                        newWord = newP+p1+word1_tot
                        tot = tot+newWord+"+"
		if (den > 0 and num > 0 and not bAuto):
			newP = "["+str(num)+"/"+str(den)+"]*"
			newWord = newP+p1+word1_tot
			tot = tot+newWord+"+"
	tot = tot[:-1]
	return tot 

def splitCumulantWeighted(str1,bSpectator,bInterest,bSpectatorSplitting,weights):
	"""Takes a list of weights that define cumulant coefficients"""

        words1 = []
        tot = ""

        # First, break up each string into its components
        word = ""
        for i in range(len(str1)):
                ch1 = str1[i]
                if (ch1 == "+"):
                        words1.append(word)
                        word = ""
                else:
                        word = word+ch1
        words1.append(word)

        nw1 = len(words1)
        for i in range(nw1):
                word1_temp = words1[i]
                # second, look for a prefactor, like beta or 3, etc.
                bPos = False
                ii=0
                while (not bPos):
                        ch1 = word1_temp[ii]
                        if (ch1 == "M"): bPos = True
                        ii = ii+1
                p1 = word1_temp[:ii-1]
                word1_tot = word1_temp[ii-1:]

                # third step is to break up word1_tot into multiplicative 
                        # components

                w1 = ""
                word1_list = []
                for ii in range(len(word1_tot)):
                        ch1 = word1_tot[ii]
                        if (ch1 == "*"):
                                word1_list.append(w1)
                                w1 = ""
                        else: w1 = w1+ch1

                word1_list.append(w1)
                nlist = len(word1_list)

                den = 0
                num = 0
                bAuto = False
                for j in range(nlist):
                        jword = word1_list[j]
                        jn = int(jword[1])
                        # don't include this word at all if has any spectators
                        # Otherwise, take a simple fraction where the numerator 
                        # is the # component of interest and the denominator
                        # is the total non-spectator components in this word
                        for jj in range(jn):
                                jnum = int(jword[4+4*jj])
                                if (bSpectator[jnum]):
                                        bAuto = True
                                        break
                                else:
                                        den = den+weights[jnum]
                                if (bInterest[jnum]):
                                        num = num+weights[jnum]
                        if (bAuto): break
                if (bAuto and bSpectatorSplitting):
                        newP = "[1/1]*"
                        newWord = newP+p1+word1_tot
                        tot = tot+newWord+"+"
                if (den > 0 and num > 0 and not bAuto):
                        newP = "["+str(num)+"/"+str(den)+"]*"
                        newWord = newP+p1+word1_tot
                        tot = tot+newWord+"+"
        tot = tot[:-1]
        return tot



def splitBrady(cumulants,nOrders):
		
	for i in range(nOrders):
		ci = splitCumulantM1(cumulants[i])
		cumulants[i] = ci
	return cumulants


def separateDerivativeString(str1):
        words1 = []
        tot = ""


        # First, break up each string into its components
        word = ""
        for i in range(len(str1)):
                ch1 = str1[i]
                if (ch1 == "+"):
                        words1.append(word)
                        word = ""
                else:
                        word = word+ch1
        words1.append(word)

        nw1 = len(words1)
        for i in range(nw1):
                word1_temp = words1[i]
                # second, look for a prefactor, like beta or 3, etc.
                bPos = False
                ii=0
                while (not bPos):
                        ch1 = word1_temp[ii]
                        if (ch1 == "M"): bPos = True
                        ii = ii+1
                p1 = word1_temp[:ii-1]
                word1_tot = word1_temp[ii-1:]

                # third step is to break up word1_tot into multiplicative 
                        # components

                word1_list = []
		newWord = ""
		bSkip = False
		myP = ""
                for ii in range(len(word1_tot)):

			if (bSkip):
				bSkip = False
				continue

                        ch1 = word1_tot[ii]
			newWord = newWord+ch1
			if (ch1 == "d"):
				bSkip = True
				newWord = newWord+"0"

				myP = myP+"D"+word1_tot[ii+1]+"*"

		newWordF = myP+p1+newWord
                tot = tot+newWordF+"+"
        tot = tot[:-1]
        return tot



def separateDerivatives(derivatives,nOrders):
        for i in range(nOrders):
                di = separateDerivativeString(derivatives[i])
                derivatives[i] = di
	return derivatives


def addNumbers(fac1,fac2):
        num1 = ""
        den1 = ""
        bDen = False
        for kk in range(1,len(fac1)-1):
                ch1 = fac1[kk]
                if (ch1 == "/"):
                        bDen = True
                        continue
                if (bDen):  den1 = den1+ch1
                else: num1 = num1+ch1

        num2 = ""
        den2 = ""
        bDen = False
        for kk in range(1,len(fac2)-1):
               ch1 = fac2[kk]
               if (ch1 == "/"):
                        bDen = True
                        continue
               if (bDen): den2 = den2+ch1
               else: num2 = num2+ch1
        num1 = int(num1)
        num2 = int(num2)
        den1 = int(den1)
        den2 = int(den2)

        if (den1 != den2):
                den = den1*den2
                num = num2*den1+num1*den2
        else:
                num = num2+num1
                den = den2
	temp = "["+str(num)+"/"+str(den)+"]"
        number = simplifyFraction(temp)
        return number




def cancelTermsString(str1,bSimplify):
        words1 = []
        tot = ""
        # First, break up each string into its components
        word = ""
        for i in range(len(str1)):
                ch1 = str1[i]
                if (ch1 == "+"):
                        words1.append(word)
                        word = ""
                else:
                        word = word+ch1
        words1.append(word)

        nw1 = len(words1)
	searchList = [] # this is our new list of words that we'll use to cancel terms
        for i in range(nw1):
                word1_temp = words1[i]
                # second, pull out the prefactor
                bPos = False
                ii=0
                while (not bPos):
                        ch1 = word1_temp[ii]
                        if (ch1 == "M"): bPos = True
                        ii = ii+1
                p1 = word1_temp[:ii-1]
                word1_tot = word1_temp[ii-1:]
		# let's simplify word1_tot for viewing purposes
		word1_tot = word1_tot+"*"
		wn = len(word1_tot)
		word1_nums = ""
		w1 = ""
		for ii in range(wn):
                        ch1 = word1_tot[ii]
                        if (ch1 == "*"):
				myn = int(w1[1])
				numword = "M"
				for kk in range(myn):
					numword = numword+w1[4+4*kk]
				word1_nums = word1_nums+numword
				w1 = ""
			else:
				w1 = w1+ch1


                w1 = ""
                p2 = ""
                bpow = 0
		number = "[1/1]"
                for ii in range(len(p1)):
                        ch1 = p1[ii]
                        if (ch1 == "*"):
                                if (w1[0] == "["):
		                     	s_num = ""
                			s_den = ""
					bDen = False	
                			for kk in range(1,len(w1)-1):
                        			ch1 = w1[kk]
                        			if (ch1 == "/"):
                                			bDen = True
                                			continue
                        			if (bDen):  s_den = s_den+ch1
                        			else: s_num = s_num+ch1
					number = float(s_num)/float(s_den)
					number = w1
                                else:
                                        p2 = p2+w1+"*"
                                w1 = ""
                        else: w1 = w1+ch1
		if (bSimplify): newWord = p2+word1_nums
		# vv we artificially added a "*" to the end of word1_tot.  Remove it. 
		else: newWord = p2+word1_tot[:-1]
		wordInfo = number,newWord
		searchList.append(wordInfo)

	for i in range(nw1):

		Iinfo = searchList[i]
		inum = Iinfo[0]
		iword = Iinfo[1]
		for j in range(i+1,nw1):

			Jinfo = searchList[j]
			jnum = Jinfo[0]
			jword = Jinfo[1]
			if (jnum == "[0/1]"): continue

			if (iword == jword):
				inum = addNumbers(inum,jnum)
				newIinfo = inum,iword
				searchList[i] = newIinfo

				newJnum = "[0/1]"
				newJinfo = newJnum,jword
				searchList[j] = newJinfo


	# now put all the words back together
	for i in range(nw1):
		Iinfo = searchList[i]
                inum = Iinfo[0]
                iword = Iinfo[1]
		#if (inum == 0 or math.fabs(inum) < FSMALL): continue
		if (inum[1] != "-" and int(inum[1]) == 0): continue

		newWord = inum+"*"+iword
		tot = tot+newWord+"+"
	if (tot != ""): tot = tot[:-1]
		
	return tot


class Cumulants:
	def __init__(self,nC,nOrders):
		self.nC = nC
	
		self.nOrders = nOrders
		if (self.nOrders > 9):
			emsg = "We can't handle orders > 9, need to recode all word analyses to do so"
			sys.exit(emsg)

		self.moments = []
		self.cumulants = []
		self.derivatives = []
		self.betaList = []
		self.fullExpansion = ""
		self.betaListCombined = []
		self.bPrint = False


		for i in range(self.nOrders+1):
			self.betaList.append("")
			self.betaListCombined.append("")

		self.compScale = []	# True if this component is scaled.  False if 
			# it is a spectator.  Defaults to True
		for i in range(self.nC):
			self.compScale.append(True)

	def getSingleWordDerivative(self,word,prefix):

		# "word" has no prefix, so we don't need to worry about that yet.  
		# we just put the prefix back on at the end

		tot = []
		nC = int(word[1]) # the total number of components being averaged
		components = []
		for i in range(nC):
			comp = word[3+4*i:7+4*i]
			components.append(comp)

		# First, take the derivative of the argument directly
		drv1 = ""
		for i in range(nC):

			# first, create a suffix of all components not i
			isuffix = ""
			for ii in range(nC):
				if (i != ii):
					isuffix = isuffix+components[ii]

			comp = components[i]
			cnum = int(comp[1])
			dnum = int(comp[3])
			if (self.compScale[cnum]):
				newcomp = comp[:3]+str(dnum+1)

				newword = prefix+"M"+str(nC)+"_"+newcomp+isuffix
				tot.append(newword)


		# next, take the derivative of the exponential in the numerator
		for i in range(self.nC):
			if (self.compScale[i]):
				# derivative of the exponential in the numerator
				newword1 = prefix+"[-1]*beta*M"+str(nC+1)+"_"+word[3:]+"h"+str(i)+"d1"
                                # derivative of the exponential in the denominator
				newword2 = prefix+"beta*"+word+"*M1_h"+str(i)+"d1"
				tot.append(newword1)
				tot.append(newword2)


		return tot


	def getStringDerivative(self,str1):
                words1 = []
                tot = ""

                # First, break up each string into its components
                word = ""
                for i in range(len(str1)):
                        ch1 = str1[i]
                        if (ch1 == "+"):
                                words1.append(word)
                                word = ""
                        else:
                                word = word+ch1
                words1.append(word)

                nw1 = len(words1)
                for i in range(nw1):
                        word1_temp = words1[i]
                        # second, look for a prefactor, like beta or 3, etc.
                        bPos = False
                        ii=0
                        while (not bPos):
                                ch1 = word1_temp[ii]
                                if (ch1 == "M"): bPos = True
                                ii = ii+1
                        p1 = word1_temp[:ii-1]
                        word1_tot = word1_temp[ii-1:]
			
			# third step is to break up word1_tot into multiplicative 
				# components

			w1 = ""
			word1_list = []
			for ii in range(len(word1_tot)):
				ch1 = word1_tot[ii]
				if (ch1 == "*"):
					word1_list.append(w1)
					w1 = ""
				else: w1 = w1+ch1

			word1_list.append(w1)
			nlist = len(word1_list)

			for j in range(nlist):
				jword = word1_list[j]

				# jsuffix handles the product rule by throwing everything that
				# we aren't differentiating onto the end of the term
				jsuffix = ""
				for jj in range(nlist):
					jjword = word1_list[jj]
					if (j != jj):
						jsuffix = jsuffix+jjword+"*"	
				# remove the last * symbol
				jsuffix = "*"+jsuffix
				jsuffix = jsuffix[:-1]


				jderv = self.getSingleWordDerivative(jword,p1)
				# jderv is a list of terms that are to be added, then multiplied by
					# jsuffix.  p1 is placed on the front of each term in jderv
					# in the above function
				nj = len(jderv)
				for k in range(nj):
					kword = jderv[k]
					
					ktot = kword+jsuffix
					
					tot = tot+ktot+"+"
		
		tot = tot[:-1]
		return tot



	def getDerivatives(self):
		self.derivatives = []
		for i in range(self.nOrders):
			#self.bPrint = True
			di = self.getStringDerivative(self.cumulants[i])
			self.derivatives.append(di)

		#print "DERIVATIVES"
		#print self.derivatives[0]
		#print ""





	def simplifyDerivatives(self):
		tot = combineStrings(self.derivatives,self.nOrders)

		# simplify prefactors
		tot = simplifyPrefactors(tot)

		# arrange terms.  This orders the multiplication of words.  It also orders
			# the configuration of the words themselves
		tot = arrangeWords(tot)

		return tot
		#print ""
		#print "Simplified"
		#print tot

	def betaPowers(self,str1):
                words1 = []
                tot = ""

                # First, break up each string into its components
                word = ""
                for i in range(len(str1)):
                        ch1 = str1[i]
                        if (ch1 == "+"):
                                words1.append(word)
                                word = ""
                        else:
                                word = word+ch1
                words1.append(word)

                nw1 = len(words1)
                for i in range(nw1):
                        word1_temp = words1[i]
                        # second, look for a prefactor, like beta or 3, etc.
                        bPos = False
                        ii=0
                        while (not bPos):
                                ch1 = word1_temp[ii]
                                if (ch1 == "M"): bPos = True
                                ii = ii+1
                        p1 = word1_temp[:ii-1]
                        word1_tot = word1_temp[ii-1:]


			# calculate the powers of beta
			w1 = ""
			p2 = ""
			bpow = 0
                        for ii in range(len(p1)):
                                ch1 = p1[ii]
                                if (ch1 == "*"):
					if (len(w1) == 5 and w1[:4] == "beta"):
						bpow = int(w1[4])
                                        else: 
						p2 = p2+w1+"*"
                                        w1 = ""
                                else: w1 = w1+ch1
			## p1 ends on "*", so the above covers all prefactor elements
			#if (len(w1) == 5 and w1[:4] == "beta"):
           	        #     	bpow = int(w1[4])
                        #else: 
                        #        p2 = p2+w1+"*"

			newWord = p2+word1_tot
			ostr = self.betaList[bpow]
			if (ostr == ""):  nstr = newWord
			else: nstr = ostr+"+"+newWord
			self.betaList[bpow] = nstr




