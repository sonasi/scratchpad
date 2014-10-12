import scipy.integrate as integrate
import scipy.stats as stats
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.special as special
import bernoulli as ber
from mpmath import *
import exact_stats as exact

#N = [[19780.0, 627.0], [51.0, 20.0]]
#N = [[8169.0, 913.0], [278.0, 22.0]]
# N = [[18963.0,8169.0], [751.0, 278.0]]

gamma = 0.57721566490153286060651209008240243104215933593992
PI = 3.14159265


tot = 30000.
prev = 0.01
maf = 0.01
mod = 2.

N = [[(1-prev)*(1-maf)*tot, (1-prev)*(maf)*tot], [(prev)*(1-maf)*tot, (mod)*(prev)*(maf)*tot]]
#N = [[int((1-prev)*(1-maf)*tot), int((1-prev)*(maf)*tot)], [int((prev)*(1-maf)*tot), int((mod)*(prev)*(maf)*tot)]]
#N_null = [ [ N[0][0], N[0][1] ], [ N[1][0], N[1][0] * N[0][1] / N[0][0] ] ]

N = [[25633, 78],[99, 5]]
N = [[24370, 72], [1052, 2]]
N = [[7,4], [1, 6]]

N = [[24381, 70], [762, 1]]
N = [[25603, 125], [24, 4]]
N = [[23102.0, 29.0], [1201.0, 2.4544824899355007]]
N = [[23101, 28], [1200, 1]]

N = [[18336, 71], [551, 1]]
N = [[21987, 17], [1912, 1]]
N = [[23305, 42], [4173, 1]]
N = [[21764, 2443], [65, 1]]
N = [[53796, 10], [305, 43]]
N = [[36052, 22876], [309, 27]]
N = [[34145, 138], [18677, 80]]
N = [[34145, 18677], [138, 80]]

N_null = [[N[0][0], N[0][1]], [N[1][0], float(N[1][0])*float(N[0][1])/float(N[0][0])]]




print N_null

N_flat = [item for sublist in N for item in sublist]
N_flat_null = [item for sublist in N_null for item in sublist]

C = [1.,-1.,-1.,1.]

print N, N_flat
print N_null

print "Contingency Table:", N
print "metric:", C

OR =  float(N[0][0]) * float(N[1][1]) /  float(N[1][0] * N[0][1])
lOR = math.log(OR)

print "OR:", OR
print "log(OR):", lOR

C = [1.,-1.,-1.,1.]

C_prior = [0.5,0.5,0.5,0.5]


a_prior = 1.

a_ = [n+a_prior for n in N_flat]
a_null = [n+c for n, c in zip(N_flat_null, C_prior)]


mu    = np.sum(  [  c*(math.log(j-0.5+1./(24.*j))) for c, j in zip(C, a_) ]  )
sigma = math.sqrt(	np.sum(  [  1./(j-0.5+1./(24.*j)) for c, j in zip(C, a_) ]  )	)
print "Bloch/Watson  {:.4e}, sigma {:.4e}\n".format(mu, sigma)


mu    = np.sum(  [  c*(math.log(j)) for c, j in zip(C, a_) ]  )
sigma = math.sqrt(	np.sum(  [  1./(j) for c, j in zip(C, a_) ]  )	)
print "Naive mu  {:.4e}, sigma {:.4e}\n".format(mu, sigma)


mu    = np.sum(  [  c*(math.log(j) - 1./(2.*j)) for c, j in zip(C, a_) ]  )
sigma = math.sqrt(	np.sum(  [  1./(j) for c, j in zip(C, a_) ]  )	)
print "Normal mu  {:.4e}, sigma {:.4e}\n".format(mu, sigma)

mu = exact.b_l(1,a_)
sigma = math.sqrt(exact.b_l(2,a_))
print "Exact mu (b_l)  {:.4e}, sigma {:.4e}\n".format(mu, sigma )


print "O(min_a) = ", 1./np.min(N)**3


size = 100000

fig = plt.figure()
ax = fig.add_subplot(111)



# A = np.random.poisson(N[0][0], size=size)
# B = np.random.poisson(N[1][1], size=size)
# C = np.random.poisson(N[1][0], size=size)
# D = np.random.poisson(N[0][1], size=size)

# S = np.random.multinomial(tot, [N[0][0]/tot, N[0][1]/tot, N[1][0]/tot, N[1][1]/tot], size=size)
# A = [x[0] for x in S]
# B = [x[3] for x in S]
# C = [x[2] for x in S]
# D = [x[1] for x in S]

S = np.random.dirichlet(  a_, size=size)
A = [tot*x[0] for x in S]
B = [tot*x[3] for x in S]
C = [tot*x[2] for x in S]
D = [tot*x[1] for x in S]

# A = np.random.gamma(N[0][0]+0.5,1., size=size)
# B = np.random.gamma(N[1][1]+0.5,1., size=size)
# C = np.random.gamma(N[1][0]+0.5,1., size=size)
# D = np.random.gamma(N[0][1]+0.5,1., size=size)

import collections
dict_ = collections.Counter()

x = []
for i in range(len(A)):
		val = 0
		a = A[i]
		b = B[i]
		c = C[i] 
		d = D[i]
		if  B[i] == 0:
			val = 0
		if  D[i] == 0:
			val = 0
		if (A[i] * b * C[i] * d > 0):	
			val = math.log(  float(A[i])*float(B[i]) / (float(C[i])*float(D[i]))  )   #/    (sigma*math.sqrt(tot))
# 		if (val < 0.005):
# 			print a, b, c, d, val
		dict_[b]+=1
		x.append(val)
	
	

n, bins, patches = ax.hist(x, 50, normed=1, facecolor='green', alpha=0.75)

print 'bins:', bins

bincenters = 0.5 * (bins[1:]+bins[:-1])

m = 6
ks = 4


u = [exact.pdf( phi, a_null) for phi in [bin for bin in bincenters]]
l = ax.plot(bincenters, u, 'b', linewidth=1)

v = [exact.cdf( phi, a_null) for phi in [bin for bin in bincenters]]
l = ax.plot(bincenters, v, 'c', linewidth=1)

q = [exact.pdf( phi, a_) for phi in [bin for bin in bincenters]]
l = ax.plot(bincenters, q, 'r--', linewidth=1)

w = [exact.cdf( phi, a_) for phi in [bin for bin in bincenters]]
l = ax.plot(bincenters, w, 'k.', linewidth=1)
# 
# r = [exact.pdf4( phi, a_, m=ks) for phi in [bin for bin in bincenters]]
# l = ax.plot(bincenters, r, 'y--', linewidth=1)

# z = mlab.normpdf( bincenters, mu, sigma) 
# l = ax.plot(bincenters, z, 'b--', linewidth=1)

print phi
print a_
print a_null
print "Kurtosis:", stats.kurtosis(x)

testval = 1.
print "\nTest Values (phi=",testval,"), a=", a_,":\n"

iters = 1000
print 'Latorre (numpy)\t', exact.pdf(testval, a_)

print iters, testval, a_
print 'Latorre ({})\t{}'.format(iters, exact.pdf2(testval, a_))
# print 'Edgewor ({})\t{}'.format(ks, exact.pdf3(testval, a_,ks))
# print 'Gram-C. ({})\t{}'.format(ks, exact.pdf4(testval, a_,ks))

print ''
# print 'pval (Edgeworth)\t', 1.-exact.pval2(testval, a_)
# print 'pval (Latorre)\t', exact.pval(testval, a_)

print ''
print 'cdf (Latorre)\t {}'.format(exact.cdf(testval, a_))
print 'cdf (Edgewor)\t {}'.format(exact.cdf2(testval, a_,m=5))
print 'cdf (Normals)\t {}'.format(exact.cdf_norm(testval, a_))
# print 'cdf (Numeric)\t {}'.format(exact.cdf_num(testval, a_))

print ''




ax.set_xlabel('')
ax.set_ylabel('Probability')
ax.grid(True)


def hist_pdf(x, n, bins):
	tot = 0
	for k in range(len(n)):
		if (x > bins[k]) and (x < bins[k+1]):
			return n[k]
	print "bin not found\n"

# mod = 1.
P_exp = a_[0]*a_[3]/((a_[2])*(a_[1]))
Y = math.log(P_exp)


X = mu

ax.axvline(x=mu, ymin=0, ymax=100)

ax.axvline(x=mu, ymin=0, ymax=100)
ax.axvline(x=exact.b_l(1,a_), ymin=0, ymax=100, color='r')


fisher_res = stats.fisher_exact(N)

print "Fisher p-value:{}".format(fisher_res)

plt.savefig("logOR_null_test2.png")
plt.savefig("logOR_null_test2.pdf")