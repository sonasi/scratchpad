import scipy.stats as stats
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

N = [[19780.0, 627.0], [51.0, 20.0]]
#N = [[8169.0, 913.0], [278.0, 22.0]]
#N = [[18963.0,8169.0], [751.0, 278.0]]

tot = 30000
prev = 0.0031
maf = 0.025

N = [[(1-prev)*(1-maf)*tot, (1-prev)*(maf)*tot],[(prev)*(1-maf)*tot, (prev)*(maf)*tot]]
OR =  N[0][0] * N[1][1] /  (N[1][0] * N[0][1])
mu = math.log(OR)

sigma  = math.sqrt(1./N[0][0] + 1./N[1][1] + 1./N[1][0] + 1./N[0][1])
sigma2 = math.sqrt(1./N[0][0] + 1./N[1][0] + 1./N[0][1])




mod = 2
sigma_test  = math.sqrt(1./N[0][0] + 1./(mod*N[1][1]) + 1./N[1][0] + 1./N[0][1])
test_OR = math.log(  N[0][0] * mod*N[1][1] /  (N[1][0] * N[0][1])  ) / (  sigma_test * math.sqrt( N[0][0] + mod*N[1][1] + N[1][0] + N[0][1] )  )  




print "Contingency Table:", N
print "Odds Ratio",  N[0][0] * mod*N[1][1] /  (N[1][0] * N[0][1]) 
print "Log OR", mu

print "sigma", sigma
print "sigma2", sigma2

size = 500000

fig = plt.figure()
ax = fig.add_subplot(111)



# A = np.random.poisson(N[0][0], size=size)
# B = np.random.poisson(N[1][1], size=size)
# C = np.random.poisson(N[1][0], size=size)
# D = np.random.poisson(N[0][1], size=size)

A = np.random.negative_binomial(N[0][0]+0.5,0.5, size=size)
B = np.random.negative_binomial(N[1][1]+0.5,0.5, size=size)
C = np.random.negative_binomial(N[1][0]+0.5,0.5, size=size)
D = np.random.negative_binomial(N[0][1]+0.5,0.5, size=size)



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
			val = math.log(  float(A[i])*float(b) / (float(C[i])*float(d))  )   /    (sigma*math.sqrt(tot))
# 		if (val < 0.005):
# 			print a, b, c, d, val
		dict_[b]+=1
		x.append(val)
	
print dict_

n, bins, patches = ax.hist(x, 100, normed=1, facecolor='green', alpha=0.75)

bincenters = 0.5 * (bins[1:]+bins[:-1])

z = mlab.normpdf( bincenters, mu*math.sqrt(tot)/(sigma), 1./math.sqrt(tot))







test_val = 2.
y = [0]*len(bincenters)
k = 2.


for k in range(1,15):
# 	norm_const = stats.poisson.pmf(k, N[1][1])
	norm_const = stats.nbinom.pmf(k, N[1][1]+0.5, 0.5)
	if k == 0:
		y1 = [0]*len(bincenters)
		y1[0] = size * norm_const
		y = [a+b for a,b in zip(y,y1)]
	else:
		mu2 = math.log(N[0][0] * k /  (N[1][0] * N[0][1]))/(sigma*math.sqrt(tot))
		y1 = mlab.normpdf( bincenters, mu2, 1.5*(sigma2/sigma)*1./math.sqrt(tot)) 
		y = [(a+b*norm_const) for a,b in zip(y,y1)]






def exact_pdf(X, N):
	print "\n Exact PDF", X, N

	sigma = math.sqrt(1./N[0][0] + 1./N[1][1] + 1./N[1][0] + 1./N[0][1])
	sigma2 = math.sqrt(1./N[0][0] + 1./N[1][0] + 1./N[0][1])
	max = int(math.ceil(X)+2)
	
	for k in range(1,10000):
		mu2 = math.log(N[0][0] * k /  (N[1][0] * N[0][1]))/(sigma*math.sqrt(tot))
		if X < mu2:
			max = k
			break
	
	print int(max)
	
	tot_val = 0
	for k in range(1, max):
		mu2 = math.log(N[0][0] * k /  (N[1][0] * N[0][1]))/(sigma*math.sqrt(tot))
		norm_const = stats.poisson.pmf(k, N[1][1])
# 		norm_const = stats.nbinom.pmf(k, N[1][1]+0.5, 0.5)

		val = stats.norm.pdf(X, loc = mu2, scale=1.5*(sigma2/sigma)*1./math.sqrt(tot) )
		tot_val += val * norm_const

	print "f(",X,"):", tot_val
	print "-----------------\n"

	return tot_val	



def exact_cdf(X, N):
	print "\n Exact CDF", X, N
	import decimal
	
	sigma = math.sqrt(1./N[0][0] + 1./N[1][1] + 1./N[1][0] + 1./N[0][1])
	sigma2 = math.sqrt(1./N[0][0] + 1./N[1][0] + 1./N[0][1])
	max = int(math.ceil(X)+2)
	tot = np.sum(N)

	for k in range(1,10000):
		mu2 = math.log(N[0][0] * k /  (N[1][0] * N[0][1]))/(sigma*math.sqrt(tot))
		if X < mu2:
			max = k
			break
	print int(max)

	
	tot_val = stats.poisson.pmf(0, N[1][1])
	for k in range(1, max):
		mu2 = math.log(N[0][0] * k /  (N[1][0] * N[0][1]))/(sigma*math.sqrt(tot))
		norm_const = stats.poisson.pmf(k, N[1][1])
# 		norm_const = stats.nbinom.pmf(k, N[1][1]+0.5, 0.5)

		val = stats.norm.cdf(X, loc = mu2, scale=1.5*(sigma2/sigma)*1./math.sqrt(tot) )
		tot_val += val * norm_const
	print "P-val(",X,"):", 1-tot_val
	print "-----------------\n"

	
	return tot_val	



# y = mlab.normpdf( bincenters, mu*math.sqrt(tot)/(sigma), 1./math.sqrt(tot))

print "Poisson Norm. Const:  Poisson({:.2f}, {:.2f}) = {}".format(test_val, N[1][1], norm_const)


l = ax.plot(bincenters, y, 'r--', linewidth=1)
l = ax.plot(bincenters, z, 'b--', linewidth=1)

ax.set_xlabel('Smarts')
ax.set_ylabel('Probability')
ax.grid(True)

X = 0.034
X = test_OR

ax.axvline(x=X, ymin=0, ymax=100)

def hist_pdf(x, y, z, n, bins):
	print "Histogram (empiracal) PDF", x
	tot = 0
	for k in range(len(n)):
		tot = tot + y[k]
# 		print tot
		if (x > bins[k]) and (x < bins[k+1]):
			print n[k]
			return n[k], y[k], z[k]
	print "bin not found\n"
	

# pdf_vec = np.vectorize(exact_pdf)
# Z1 = pdf_vec(bincenters, N)
# print Z1

plt.savefig("logOR_null_test.png")


print "Test Val:", test_OR

print hist_pdf(X, y, z, n, bins)
print exact_pdf(X, N)
print exact_cdf(X, N)


def prior(l):
	X = 8.
	return stats.poisson.pmf(X, l)/math.sqrt(X)



print "Fischer Exact:"
print [[N[0][0], N[0][1]], [N[1][0], mod*N[1][1]]]
oddsratio, pvalue = stats.fisher_exact([[N[0][0], N[0][1]], [N[1][0], mod*N[1][1]]])
print "OR: {}, p-val: {}".format(oddsratio, pvalue)
#n, bins, patches = ax.hist(x, 100, normed=1, facecolor='green', alpha=0.75, cumulative=True)


# plt.show()