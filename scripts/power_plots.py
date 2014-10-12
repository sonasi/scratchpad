#!/usr/bin/env python
"""
This script runs analysis on a GWAS/PheWAS dataset to determine the power and other 
properties.  Much of the analysis is based the paper:

	Experimental Designs for Robust Detection of Effects in Genome-Wide Case-Control Studies
	(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3241427)

"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from prettytable import PrettyTable
import matplotlib.mlab as mlab


def power(MAF, OR, PREV):

	# Define the Bayes Factor, sample size,  prior sample size, and case/control ratio
	B = 1e6 
	n = 300000
# 	MAF = 0.05

	a = 1
	c = 0.1
	c = PREV
# 	OR = 2.

	N = np.array([[100, 50],[1000,10]], np.float)
	N = np.array([[(1-c)*n, (1-c)*n*MAF],[c*n,OR*c*n*MAF]], np.float)

	n = np.sum(N)
	loud = False


	# Define the test statistic for a given Bayes Factor
	Z_B2 = (n+a)/n**2 * math.log( (n+a)/a * B**2 )
	
	if loud:
		print Z_B2
	Z_B = math.sqrt(Z_B2)
	sigma_n2 = 1/N[0][0] + 1/N[0][1] + 1/N[1][0] + 1/N[1][1] 
	sigma_12 = n * sigma_n2
	sigma_1 = math.sqrt(sigma_12)

	etahat = math.log(  (N[0][0] * N[1][1]) /  (N[1][0] * N[0][1])  )  

	P = 1 - stats.norm.cdf(math.sqrt(float(n))*(Z_B-etahat/sigma_1))


	x = PrettyTable(["", "X=0", "X=1"])
	x.add_row(["Control", N[0][0], N[0][1]] )
	x.add_row(["Cases", N[1][0], N[1][1]] )

	if loud:
		print "n=", n
		print "log(OR)=", etahat
		print "Contingency Table:"
		print x
		print "OR:", 1.2
		print "Z_B:", Z_B
		print "Power", P

	return P

power2 = np.vectorize(power)

# print "Autism sample:", power()


delta1 = 0.0005
delta2 = 0.1
x = np.arange(0.0001, 0.01, delta1)
y = np.arange(2.0, 5.0, delta2)
X, Y = np.meshgrid(x, y)
Z1 = power2(X, Y, 0.05)
# difference of Gaussians
Z = Z1


# Create a simple contour plot with labels using default colors.  The
# inline argument to clabel will control whether the labels are draw
# over the line segments of the contour, removing the lines beneath
# the label
plt.figure()
CS = plt.contour(X, Y, Z)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('Statistical Power')
plt.xlabel('MAF')
plt.ylabel('Odds Ratio')
plt.show()
