
import numpy as np
import sys
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects.packages import importr

class Skat_Wrapper:

	def __init__(self):
		
		self.skat = importr("SKAT", robject_translations = {"Beta.Weights": "Beta_Weights2", "to.period": "to_period2"})

	def p_value_cont(self, X, y, Z):

		form = rpy2.robjects.Formula('y~X')

		env = form.environment
		env['y'] = y
		env['X'] = X
		env['Z'] = Z

		obj = self.skat.SKAT_Null_Model(form, out_type="C")
		res = self.skat.SKAT_CommonRare(Z, obj)

		return res[0][0]

	def p_value_bin(self, X, y, Z):

		form = rpy2.robjects.Formula('y~X')

		env = form.environment
		env['y'] = y
		env['X'] = X
		env['Z'] = Z

		obj = self.skat.SKAT_Null_Model(form, out_type="D")
		res = self.skat.SKAT_CommonRare(Z, obj)

		return res[0][0]

	

skat = importr("SKAT", robject_translations = {"Beta.Weights": "Beta_Weights2", "to.period": "to_period2"})

ro.r('library(SKAT)')
ro.r('data(SKAT.example)')
ro.r('names(SKAT.example)')

data = ro.r('attach(SKAT.example)')
skat_inst = ro.r('obj<-SKAT_Null_Model(y.c ~ X, out_type="C")')
res = ro.r('SKAT(Z, obj)$p.value')



print 'Result =',  res[0]


X = np.array(data['X'])
y_b = np.array(data['y.b'])
y_c = np.array(data['y.c'])
Z = np.array(data['Z'])


form = rpy2.robjects.Formula('y_c~X')

env = form.environment
env['y_c'] = y_c
env['X'] = X
env['Z'] = Z

obj = skat.SKAT_Null_Model(form, out_type="C")
res = skat.SKAT_CommonRare(Z, obj)
print 'p-value1', res[0][0]

sw = Skat_Wrapper()
print 'p-value2', sw.p_value_cont(X, y_c, Z)
print 'p-value3', sw.p_value_cont(X, y_b, Z)



a1 = np.array([[1., 4.], [0., 0. ], [0., 1.], [0., 0.]])
a2 = np.array([[1., 0., 0., 0.], [0., 0., 0., 1.], [0., 0., 1., 0.], [0., 2., 0., 0.]])
v = np.array([1, 0, 0, 0])

print 'p-value4', sw.p_value_cont(a1, v, a2)



m = rpy2.robjects.r.matrix(rpy2.robjects.IntVector(range(10)), nrow=5)

v_ = rpy2.robjects.vectors.IntVector(v)
m1 = rpy2.robjects.r.matrix(a1, nrow=4)
m2 = rpy2.robjects.r.matrix(a2, nrow=4)

form = rpy2.robjects.Formula('v_~m1')
env = form.environment
env['v_'] = v_
env['m1']=m1
env['m2']=m2

obj = skat.SKAT_Null_Model(form, out_type="D")
res3 = skat.SKAT_CommonRare(m2, obj)

print 'p-value 3', res3[0][0]


# df_pheno = pandas.io.parsers.read_csv("/home/vagrant/data/labwas/data/pheno_Ws_no_lab_count_req/pheno_{}_Exome_demos_Ws_Age18.txt".format('%SAT'), dtype=str, sep='\t', names=['FID', 'IID', 'Mean', 'Median', 'Min', 'Max', 'Ever_high', 'Ever_low', 'First'])

