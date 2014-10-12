
import scipy.stats as stats
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab



fig = plt.figure()
ax = fig.add_subplot(111)

r = 400.
sigma = math.sqrt(r)

bincenters = np.array([x for x in range(int(r-3*math.sqrt(r)), int(r+3*math.sqrt(r)) )])

z = mlab.normpdf( bincenters, r, sigma) #*math.sqrt(tot)/(sigma), 1.4/math.sqrt(tot))
y = stats.nbinom.pmf( bincenters, r+0.5, 0.5 )
x = stats.poisson.pmf( bincenters, r )

l = ax.plot(bincenters, z, 'b--', linewidth=1)
l = ax.plot(bincenters, y, 'r--', linewidth=1)
l = ax.plot(bincenters, x, 'g--', linewidth=1)

ax.set_xlabel('')
ax.set_ylabel('Probability')
ax.grid(True)


plt.savefig("logOR_null_test2.png")
