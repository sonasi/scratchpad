import numpy as np
import pandas as pd
from prettytable import PrettyTable
import math
import scipy.stats as stats
import statsmodels.api as sm
import plinkio_example as pl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import exact_posterior as exact
from sklearn.linear_model import LogisticRegression

df_demos = pd.read_csv("../data/demos_Ws_view.tsv", '\t')

def hist_pdf(x, y, n, bins):
	tot = 0
	for k in range(len(n)):
		tot = tot + y[k]
		if (x > bins[k]) and (x < bins[k+1]):
			print "PDF Value (asymtotic):", n[k], "(", y[k]/np.sum(y),")"
			return n[k]

	print "bin not found"
	return -100

def exact_pdf(X, N):
	print "\n Exact CDF", X, N

	sigma = math.sqrt(1./N[0][0] + 1./N[1][1] + 1./N[1][0] + 1./N[0][1])
	sigma2 = math.sqrt(1./N[0][0] + 1./N[1][0] + 1./N[0][1])
	max = int(math.ceil(X)+2)
	tot = np.sum(N)
	
	for k in range(1,10000):
		mu2 = math.log(N[0][0] * k /  (N[1][0] * N[0][1]))/(sigma*math.sqrt(tot))
		if X < mu2:
			max = k + 10
			break
	print max
	tot_val = 0
	for k in range(1, max):
		mu2 = math.log(N[0][0] * k /  (N[1][0] * N[0][1]))/(sigma*math.sqrt(tot))
		norm_const = stats.poisson.pmf(k, N[1][1])
		val = stats.norm.pdf(X, loc = mu2, scale=(sigma2/sigma)*1./math.sqrt(tot) )
		tot_val += val * norm_const
	
	return tot_val	

def exact_cdf(X, N):

	print "\n Exact CDF", X, N
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
		val = stats.norm.cdf(X, loc = mu2, scale=(sigma2/sigma)*1./math.sqrt(tot) )
		tot_val += val * norm_const


	print "pval:", 1-tot_val
	print "-----------------\n"
	
	return tot_val	

def logodds_exact(X, N):

	size = 1000000

	OR =  N[0][0] * N[1][1] /  (N[1][0] * N[0][1])
	mu = math.log(OR)

	N112 = N[1][0]*N[0][1]/N[0][0]
	sigma1 = math.sqrt(1./N[0][0] + 1./N[1][1] + 1./N[1][0] + 1./N112)
	sigma2 = math.sqrt(1./N[0][0] + 1./N[1][1] + 1./N[1][0] + 1./N[0][1])
	sigma3 = math.sqrt(1./N[0][0] + 1./N[1][0] + 1./N[0][1])

	N_null = [ [N[0][0],N[0][1]], [N[1][0], N112]  ]

# 	A = np.random.poisson(N[0][0], size=size)
# 	B = np.random.poisson(N112, size=size)
# 	C = np.random.poisson(N[1][0], size=size)
# 	D = np.random.poisson(N[0][1], size=size)
# 	E = np.random.poisson(N[1][1], size=size)


	A = np.random.negative_binomial(N[0][0]+0.5,0.5, size=size)
	B = np.random.negative_binomial(N112 + 0.5 ,0.5, size=size)
	C = np.random.negative_binomial(N[1][0]+0.5,0.5, size=size)
	D = np.random.negative_binomial(N[0][1]+0.5,0.5, size=size)
	E = np.random.negative_binomial(N[1][1]+0.5,0.5, size=size)



	tot = np.sum(N)


	z = []
	for i in range(len(A)):
		val = 0
		a = A[i]
		b = E[i]
		c = C[i] 
		d = D[i]
		if  E[i] == 0:
			val = 0
		if  D[i] == 0:
			val = 0
		if (A[i] * b * C[i] * d > 0):	
			val = math.log(float(A[i])*float(b)/(float(C[i])*float(d)))/(sigma2*math.sqrt(tot))
		z.append(val)


	print "Contingency Table:", N
	print "Null Contingency Table:", N_null

	print "Odds Ratio", OR
	print "Log OR", mu
	print "sigma1", sigma1
	print "sigma2", sigma2

	fig = plt.figure()
	ax = fig.add_subplot(111)

	print N



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
			val = math.log(float(A[i])*float(b)/(float(C[i])*float(d)))/(sigma1*math.sqrt(tot))
		x.append(val)

 	n, bins, patches = ax.hist(z, 100, normed=1, facecolor='red', alpha=0.75)
	n, bins, patches = ax.hist(x, 100, normed=1, facecolor='green', alpha=0.75)

	bincenters = 0.5*(bins[1:]+bins[:-1])

	y = mlab.normpdf( bincenters, 0, 1./math.sqrt(tot))
	l = ax.plot(bincenters, y, 'r--', linewidth=1)

	


	ax.set_xlabel('log(OR)')
	ax.set_ylabel('Probability')
	ax.grid(True)

  	ax.axvline(x=X, ymin=0, ymax=100)

	print "test Value:", X
	plt.savefig("logOR_null.png")

	
	res = hist_pdf(X, y, n, bins)


	print exact_pdf(X, N_null)

	print 1-exact_cdf(X, N_null)

	if res < -1:
		return 1-exact_cdf(X, N_null)
	else:
		return res		




def get_data(SNP, pheno):
	return pl.get_SNP_samples(SNP), pl.pheno_from_file(pheno)
	
	
	
	


def logistic_regression(df_pheno, df_snp, df_control_snp=False, loud=False, drop_covariates=False):
	'''
	Perform logistic regression
	
	input: pd dataframes of individuals: df_pheno, df_snp
	keyword terms: loud - be verbose
	returns: results of logit
	'''
	
	df_logit = pd.merge(df_snp, df_pheno, on="FID")
	
	cols_to_keep = ['pheno', 'min_age', 'Male', 'SNP']

	if drop_covariates:
		cols_to_keep = ['pheno', 'SNP']	
	
	if isinstance(df_control_snp, type(df_snp)):
		df_control_snp = df_control_snp.rename(columns={'SNP': 'SNP_control'})
		df_logit = pd.merge(df_logit, df_control_snp, on="FID")
		cols_to_keep.append('SNP_control')

	if loud:
		print "\nTesting SNP, {}, and phenotype, {}".format(df_snp.snp_name, df_pheno.pheno_name)

	# Shift the phenotype variable
	df_logit['pheno'] = df_logit[df_logit.columns[2]].map({2: 1, 1: 0, -9:-9})
# 	df_logit['SNP'] = df_logit['SNP'].map({0: 0, 1: 1, 3:2})


	df_logit = df_logit[df_logit['SNP'] < 3]
	if loud:
		print "\nContingency Table:"
	 
	cross = pd.crosstab(df_logit['pheno'], df_logit['SNP'], rownames=['pheno'])

	if len(df_logit[df_logit['SNP']==2]) == 0:
		cross[2] = 0


	x = PrettyTable(["Allele", "X=0", "X=1", "X=2"])
	x.add_row(["Control", cross[0][0], cross[1][0], cross[2][0]] )
	x.add_row(["Cases",  cross[0][1],  cross[1][1], cross[2][1]] )
	
	if loud:
		print x

	# Merge with Demographic info
	df_logit_clean = df_logit[df_logit['pheno'] > -1]

	# Do like plink: if either of the homozygote numbers is zero drop that column:
# 	if cross[2][0] == 0 or cross[2][1] == 0:
	df_logit_clean = pd.merge(df_demos[['FID', 'min_age', 'gender']], df_logit_clean, on="FID")

	print cols_to_keep

	# Make the gender column a numeral
	df_logit_clean['Male'] = df_logit_clean['gender'].map({'M': 1,'F': 0})
	
	# Dummy Variables (is this necessary or correct?)
	dummy_ranks = pd.get_dummies(df_logit_clean['SNP'], prefix='allele')
	
	# More rearranging
	data = df_logit_clean[cols_to_keep]   #.join(dummy_ranks.ix[:, 'allele_1':])
	data['intercept'] = 1.0
# 	data['exp'] = math.e


	train_cols = data.columns[1:]	
	print train_cols
#   	train_cols = ['SNP', 'intercept']

	if loud:
		print "\nTraining on columns:", train_cols, "\n"
	
	logit = sm.Logit(data['pheno'], data[train_cols], missing='drop')
	result = logit.fit()

	if loud:
		print "\nLLR:\t\t{:.2f}".format(result.llr)
		print "\nDOF:\t\t{:.2f}".format(result.df_model)
		print "LLR p-value:\t{:.2e}".format(result.llr_pvalue)
	
	
	params = result.params
	pvalues = result.pvalues
	
	conf = result.conf_int()
	conf['OR'] = np.exp(params)
	conf['P'] = pvalues
	conf.columns = ['2.5%', '97.5%', 'OR', 'P']	
	
	x = PrettyTable()
	x.add_column("Var", train_cols)
	x.add_column("OR", [ "{:.2f}".format(res) for res in np.exp(params)])	
	x.add_column("P",  [ "{:.1e}".format(res) for res in result.pvalues])	

# 	B_obs = math.exp(result.llf-result.llnull)

	if loud:
		print x
# 		print "Bayes Factor:", B_obs

	print "{:.2f}\t{:.2e}".format(  np.exp(result.params['SNP']), result.pvalues['SNP'])
	return (  np.exp(result.params['SNP']), result.pvalues['SNP'] )
	
#  	clf_l1_LR = LogisticRegression(penalty='l1', tol=0.01)
#  	clf_l1_LR.fit(X, y)
# 	
# 	coef_l1_LR = clf_l1_LR.coef_.ravel()
# 	print clf_l1_LR.coef_
# 	print "params: ",clf_l1_LR.get_params(deep=True)



def cdf(x, scale=1.):
        val = 1-0.5*(1+math.erf(x/math.sqrt(2*scale**2)))
        return val



# Read in only the phenotypes of interest
def compute_stats(df_pheno, df_snp, loud=False):

	S_case = df_pheno[df_pheno[df_pheno.columns[1]]==2].FID.tolist()
	S_cont = df_pheno[df_pheno[df_pheno.columns[1]]==1].FID.tolist()

	G_0 = df_snp[df_snp['SNP']==0].FID.tolist() 
	G_1 = df_snp[df_snp['SNP']==1].FID.tolist() 
	G_2 = df_snp[df_snp['SNP']==2].FID.tolist() 

	set(S_cont).intersection(G_0)

	N = np.array([
				 [len(set(S_cont).intersection(G_0)), len(set(S_cont).intersection(G_1)), len(set(S_cont).intersection(G_2))],
				 [len(set(S_case).intersection(G_0)), len(set(S_case).intersection(G_1)), len(set(S_case).intersection(G_2))]
				 ]
				 , np.float)

	x = PrettyTable(["Allele:", "X=0", "X=1", "X=2"])
	x.add_row(["Control", N[0][0], N[0][1], N[0][2]] )
	x.add_row(["Cases", N[1][0], N[1][1], N[1][2]] )
	
	if loud:
		print "\nTesting SNP, {}, and phenotype, {}\n".format(df_snp.snp_name, df_pheno.pheno_name)
		print "Contingency Table:"
	 	print x
	OR1 =  (N[0][0] * N[1][1]) /  (N[1][0] * N[0][1])  
	OR2 =  (N[0][0] * N[1][2]) /  (N[1][0] * N[0][2]) 

	n = np.sum(N)
	c = (N[0][0] + N[0][1] + N[0][2])/n
	Nmin = 0

	if N[0][2]*N[1][2] == 0:
		OR =  (N[0][0] * N[1][1]) /  (N[1][0] * N[0][1])  
		sigma_n2 = (1/N[0][0] + 1/N[0][1] + 1/N[1][0] + 1/N[1][1])
		etahat = math.log( OR )
		sigma_null =  (  1/N[0][0] + 1/N[0][1]  +  (1-c)/N[0][0] + (1-c)/N[0][1]  )
		N = [[N[0][0], N[0][1]], [N[1][0], N[1][1]]]
		Nmin = np.min(N)
	else:
		OR = math.sqrt( N[0][0] * N[1][2] / ( N[0][2] * N[1][0]) )
		etahat = 0.5 * math.log( N[0][0] * N[1][2] / ( N[0][2] * N[1][0]) )
		sigma_n2 = 0.5 * (1/N[0][0] + 1/N[0][2] + 1/N[1][0] + 1/N[1][2])
		sigma_null =  0.5 * (  1/N[0][0] + 1/N[0][2]  +  (1-c)/N[0][0] + (1-c)/N[0][2]  )
		N = [[N[0][0], N[0][2]], [N[1][0], N[1][2]]]
		Nmin = np.min(N)
#  	etahat = math.fabs(etahat)
	



	
	
	
	

# 	OR = (N[0][0] * N[1][1]) /  (N[1][0] * N[0][1])
# 	OR = (N[0][1]/N[0][0] +  N[0][2]/N[0][1] ) / (N[1][1]/N[1][0] +  N[1][2]/N[1][1] )
	
	
# 	      etahat = math.log((N[0][0] * N[1][1]) /  (N[1][0] * N[0][1])  )  


# 	print "n=\t\t{:.3f}".format(n)
# 	print "OR=\t\t{:.3f}".format(OR)
# 	print "Log(OR)=\t{:.3f}".format(etahat)
	
	
	a = 1
	B = 1e6 

	# Define the test statistic for a given Bayes Factor
	Z_B2 = (n+a)/n**2 * math.log( (n+a)/a * B**2 )
	Z_B = math.sqrt(Z_B2)
	
	# Eqtn.'s A48, Ball [1]
	sigma_n = math.sqrt(sigma_n2)
	sigma_12 = n * sigma_n2
	sigma_1 = math.sqrt(sigma_12)

# 	print ""
# 

	Z_n = etahat / sigma_1
	
# 	print Z_n*math.sqrt(n)
	

# 	
# 	print "P-Value: cdf({:.2},scale={:.2}) \t{:.2e}".format( Z_n*math.sqrt(n), math.sqrt(1./n), P_value)
# 	print "P-Value2:\t{:.2e}".format(P_value2)
# 	print "P-Value3:\t{:.2e}".format(P_value3)
# 
# 	print ""
# 		
# 	print "sigma_1:\t{:.3f}".format(sigma_1)

	N112 = N[1][0]*N[0][1]/N[0][0]
	N_null = [ [N[0][0],N[0][1]], [N[1][0], N112]  ]

	print "L(H0)", exact_pdf(Z_n, N_null)
	print "L(H1)", exact_pdf(Z_n, N)

	B_obs = math.sqrt(a) / math.sqrt(n+a) \
  		  * math.exp(  n**2 * Z_n**2 / (2.* (n + a) )  )
  		   		  
  	pi = stats.norm.pdf(0, loc = 0, scale=1./a)
	g  = stats.norm.pdf(0, loc = Z_n* math.sqrt(n/(n+a)), scale=math.sqrt(1./(n+a)))
	B_acc = pi/g	


	Power = 1 - stats.norm.cdf( math.sqrt(float(n))*(Z_B-etahat/sigma_1) )
	P_val = 2 * ( 1 - stats.norm.cdf(  math.sqrt(n) * math.fabs(Z_n) ) )
	

	# Compute the "chi-squared llr" p-value
	llr = 2*math.log(B_obs)
	P_val2 = stats.chisqprob(  math.fabs(llr), 1 )
	P_val3 = logodds_exact( Z_n, N)
	
	
# 	print "Bayes Factor\t{:.2f}".format(B_obs)
# 	print "Power:\t\t{:.3f}".format(Power)
#
	if loud:
	 	print "etahat:\t{:.3f}".format(etahat)

	 	print "Z_n\t\t{:.3f} \nsigma_norm\t{:.5f}\n".format(Z_n, math.sqrt(1./n) )

		print "pi, {}".format(pi)
		print "g, {}".format(g)
	
		print "B_acc, {}\n".format(B_acc)

	 	print "P-Value: \t{:.2e}".format(P_val)
	 	print "P-Value2:\t{:.2e}".format(P_val2)
	 	print "P-Value3:\t{:.2e}".format(P_val3)

	 	print "llr:\t{:.2f}".format(llr)

		print "\nResults:"
		print "SNP\t\tpheno\t\tOR\tsigma_1\tZ_n\tP_val\t\tP_val2\t\tB_obs\t\tPower\tminN"
		print "{}\t{}\t{:.3f}\t{:.2f}\t{:.2f}\t{:.2e}\t{:.2e}\t{:.2e}\t{:.3f}\t{}".format( df_snp.snp_name, df_pheno.pheno_name, OR, sigma_1, Z_n, P_val, P_val2, B_obs, Power, Nmin )

	return [df_snp.snp_name, df_pheno.pheno_name, OR, sigma_1, Z_n, P_val, P_val2, B_obs, Power, Nmin]



def test_functions():
	df_crosscheck = pd.read_csv("../data/Exome_36k_Ws_GAE3_view_view.csv", dtype=str)
	df_groups = df_crosscheck.groupby('SNP')

	for name, group in df_groups:
		df_snp = pl.get_SNP_samples(name)
	
		for i, row in group.iterrows():
			df_pheno = pl.pheno_from_file(str(row['jd_code']))
# 			logistic_regression(df_pheno, df_snp)
			res = compute_stats(df_pheno, df_snp)	
			print res			
