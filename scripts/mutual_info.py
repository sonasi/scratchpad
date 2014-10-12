#!/usr/bin/env python

'''
File: mutual_info.py
Created: 09-Aug-2014 Joseph P. Bochenek
1.0
Description: Mutual Information 
'''


print(__doc__)

import time as time

from plinkio import plinkfile
from plinkio import cplinkio

import collections
import numpy as np
import pandas
from prettytable import PrettyTable
import math
import scipy.stats as stats

import plinkio_example as pl

def multiple_snp_sets(*args):
	
	
	
	return p1


# def mutual_information(*args, **kwargs ):
def mutual_information(df_snp_1, df_snp_2 ):
	'''
	Compute the mutual information of two SNPs, or two snips conditional on a phenotype:
	
	http://en.m.wikipedia.org/wiki/Mutual_information
	

	Input:
	as many dataframes as you want representing genotypes
	
	Return:
	Tuple with the mutual information and the entropy of df_snp_1
	'''
	
	allele_types = [0, 1, 2]
	total = len( set(df_snp_1.FID.tolist())  )
	
	MI_sum = 0
	Entropy_sum = 0
	
	for allele1 in allele_types:	

		p_1 = float( len( set(df_snp_1[df_snp_1['SNP']==allele1].FID.tolist()) ) ) \
			  / total
		
		Entropy_sum += -1 * p_1 * math.log( p_1, 2)

		for allele2 in allele_types:	

			p_j = float( len( set(df_snp_1[df_snp_1['SNP']==allele1].FID.tolist() ) \
					  &  set( df_snp_2[df_snp_2['SNP']==allele2].FID.tolist())))  \
				  / total


			p_2 = float( len( set(df_snp_2[df_snp_2['SNP']==allele2].FID.tolist() ) ) ) \
				  / total	

			if not(p_j == 0):			
				MI_sum += p_j * math.log( p_j / (p_1*p_2), 2 )
	
	return (MI_sum, Entropy_sum)



def combine_snps_name(args):

	snps = []
	for arg in args:	
		print "Processing SNP,", arg
		snps.append(pl.get_SNP_samples(arg))
	
	return combine_snps(snps)

def combine_snps(args):

	if len(args) == 1:
		return args[0]

	df_ref = args[0]
	df_cat = pandas.concat( args, axis=1 )
	combined = []

	for row in df_cat['SNP'].iterrows():
		combined.append( np.max(row[1]) )

	df_comb = pandas.DataFrame( zip(df_ref.FID, combined), index=df_ref.index, columns = ['FID', 'SNP'] ) 
	return df_comb
	
	

def mi_venn_snps(df_snp_1, df_snp_2, df_snp_3, file='venn_snp.png'):
	import matplotlib_venn
	from matplotlib import pyplot as plt

	n_A = len( set(df_snp_1[df_snp_1['SNP'].isin([1,2])].FID.tolist() ) ) 
	n_B = len( set(df_snp_2[df_snp_2['SNP'].isin([1,2])].FID.tolist() ) )
	n_C = len( set(df_snp_3[df_snp_3['SNP'].isin([1,2])].FID.tolist() ) ) 
	n_AB =  len( set(df_snp_1[df_snp_1['SNP'].isin([1,2])].FID.tolist() ) \
			   & set(df_snp_2[df_snp_2['SNP'].isin([1,2])].FID.tolist() ) )
	n_AC =  len( set(df_snp_1[df_snp_1['SNP'].isin([1,2])].FID.tolist() ) \
			   & set(df_snp_3[df_snp_3['SNP'].isin([1,2])].FID.tolist() ) )
	n_BC =  len( set(df_snp_2[df_snp_2['SNP'].isin([1,2])].FID.tolist() ) \
			   & set(df_snp_3[df_snp_3['SNP'].isin([1,2])].FID.tolist() ) )
	n_ABC = len( set(df_snp_1[df_snp_1['SNP'].isin([1,2])].FID.tolist() ) \
			   & set(df_snp_2[df_snp_2['SNP'].isin([1,2])].FID.tolist() ) \
			   & set(df_snp_3[df_snp_3['SNP'].isin([1,2])].FID.tolist() ) )
	n_Abc = n_A - n_AC - n_AB + n_ABC
	n_aBc = n_B - n_BC - n_AB + n_ABC
	n_ABc = n_AB - n_ABC
	n_abC = n_C - n_BC - n_AC + n_ABC
	n_AbC = n_AC - n_ABC
	n_aBC = n_BC - n_ABC
	n_ABC = n_ABC
	
	print n_Abc, n_aBc, n_ABc, n_abC, n_AbC, n_aBC, n_ABC
	
	#  (Abc, aBc, ABc, abC, AbC, aBC, ABC)

	matplotlib_venn.venn3(subsets = (n_Abc, n_aBc, n_ABc, n_abC, n_AbC, n_aBC, n_ABC), set_labels = ('SNP_1', 'SNP_2', 'SNP_3'))

	plt.savefig(file)
	plt.close()



def mi_venn(df_snp_1, df_snp_2, df_pheno, file='venn.png'):
	import matplotlib_venn
	from matplotlib import pyplot as plt
	
	pheno = df_pheno.pheno_name

	pheno_total = float(len(set(df_pheno[df_pheno['pheno_'+pheno]==2].FID.tolist()) )  )
	
	snp_codes = [1,2]
	
	print "Venn Diagram"
	
	
	p_1 = float( len( set(df_snp_1[df_snp_1['SNP'].isin([1,2])].FID.tolist())   \
				& set(df_pheno[df_pheno['pheno_'+pheno]==2].FID.tolist()) )  )

	p_2 = float( len( set(df_snp_2[df_snp_2['SNP'].isin([1,2])].FID.tolist())   \
				& set(df_pheno[df_pheno['pheno_'+pheno]==2].FID.tolist()) )  )
				

	n_ABC = float( len( set(df_snp_1[df_snp_1['SNP'].isin([1,2])].FID.tolist() ) \
						  &  	  set(df_snp_2[df_snp_2['SNP'].isin([1,2])].FID.tolist() ) \
						  &  	  set(df_pheno[df_pheno['pheno_'+pheno]==2].FID.tolist())
						  )  )
	
	n_AB =  len( set(df_snp_1[df_snp_1['SNP'].isin([1,2])].FID.tolist() ) \
				&  	   set(df_snp_2[df_snp_2['SNP'].isin([1,2])].FID.tolist() ) )
	
	n_AC =  len( set(df_snp_1[df_snp_1['SNP'].isin([1,2])].FID.tolist() ) \
						  &  	  set(df_pheno[df_pheno['pheno_'+pheno]==2].FID.tolist())
						  )  
	
	n_BC =  len( set(df_snp_2[df_snp_2['SNP'].isin([1,2])].FID.tolist() ) \
						  &  	  set(df_pheno[df_pheno['pheno_'+pheno]==2].FID.tolist())
						  )
	n_A = len( set(df_snp_1[df_snp_1['SNP'].isin([1,2])].FID.tolist() ) ) 
	n_B = len( set(df_snp_2[df_snp_2['SNP'].isin([1,2])].FID.tolist() ) )
	
	n_Abc = n_A - n_AC - n_AB + n_ABC
	n_aBc = n_B - n_BC - n_AB + n_ABC
	n_ABc = n_AB - n_ABC
	n_abC = pheno_total - n_BC - n_AC + n_ABC
	n_AbC = n_AC - n_ABC
	n_aBC = n_BC - n_ABC
	n_ABC = n_ABC

	print "n_pheno: {}\t n_SNP_common: {}\t n_SNP_rare: {}\n".format(pheno_total, n_A, n_B)


	print n_Abc, n_aBc, n_ABc, n_abC, n_AbC, n_aBC, n_ABC
	
	#  (Abc, aBc, ABc, abC, AbC, aBC, ABC)

	matplotlib_venn.venn3(subsets = (n_Abc, n_aBc, n_ABc, n_abC, n_AbC, n_aBC, n_ABC), set_labels = ('SNP_common', 'SNP_rare', 'Pheno'))

	plt.savefig(file)
	plt.close()





def conditional_mutual_information(df_snp_1, df_snp_2, df_pheno ):
	'''
	Compute the mutual information of two SNPs, or two snips conditional on a phenotype:
	
	http://en.m.wikipedia.org/wiki/Mutual_information
	

	Input:
	as many dataframes as you want representing genotypes
	
	Return:
	Tuple with the mutual information and the entropy of df_snp_1
	'''
	
	allele_types = [0, 1, 2]
	pheno_vals = [2]
	pheno = df_pheno.pheno_name
	
	total = float(len(set(df_pheno[df_pheno['pheno_'+pheno]>-1].FID.tolist()) )  )
	
	print "Disease Prevalence: {}".format(compute_prevalance(df_pheno))
	
	

	MI_sum = 0
	Entropy_sum = 0
	p_sum = 0

	for pheno_val in pheno_vals:

		p_12 = 0	

		pheno_total = float(len(set(df_pheno[df_pheno['pheno_'+pheno]==pheno_val].FID.tolist()) )  )
		
		p_z =  pheno_total \
			/ total
		

		for allele1 in allele_types:	

			p_1 = float( len( set(df_snp_1[df_snp_1['SNP']==allele1].FID.tolist())   \
				& set(df_pheno[df_pheno['pheno_'+pheno]==pheno_val].FID.tolist()) )  ) \
				  / pheno_total
		
			if not(p_1 == 0):	
				Entropy_sum += -1 * p_z * p_1 * math.log( p_1, 2)

			for allele2 in allele_types:	
				
 				print "Pheno:", pheno_val, "Allele1:", allele1, "Allele2:", allele2
				
				p_j = float( len( set(df_snp_1[df_snp_1['SNP']==allele1].FID.tolist() ) \
						  &  	  set(df_snp_2[df_snp_2['SNP']==allele2].FID.tolist() ) \
						  &  	  set(df_pheno[df_pheno['pheno_'+pheno]==pheno_val].FID.tolist())
						  ))  \
					  / pheno_total


				p_2 = float( len( set(df_snp_2[df_snp_2['SNP']==allele2].FID.tolist() ) 
								& set(df_pheno[df_pheno['pheno_'+pheno]==pheno_val].FID.tolist()) ) ) \
					  / pheno_total	
			
 				print "\t", 'p_z', p_z, 'n_12', p_j*pheno_total, "p_1", p_1, "p_2", p_2, "p_12", p_j, "p_1*p_2", p_1*p_2 
				
				if not(p_j == 0):			
					MI_sum += p_z * p_j * math.log( p_j / (p_1*p_2), 2 )
					print "\t MI_i", p_z * p_j * math.log( p_j / (p_1*p_2), 2 )

					print "\t", p_j/p_1
				
				if allele1 == allele2:
					if not((p_1*pheno_total) == 0):
						p_12 += p_1 * p_j*pheno_total / (p_1*pheno_total)
				
		p_sum += p_z

		print 'p_z', p_z, 'p_12', p_12, '\n'
	
	print 'pheno_total', pheno_total, 'p_sum', p_sum
		
	return (MI_sum, Entropy_sum, MI_sum/Entropy_sum)


def combine_snps(args):
	if len(args) == 1:
		return args[0]

	df_ref = args[0]
	df_cat = pandas.concat( args, axis=1 )
	combined = []

	for row in df_cat['SNP'].iterrows():
		combined.append( np.max(row[1]) )

	df_comb = pandas.DataFrame( zip(df_ref.FID, combined), index=df_ref.index, columns = ['FID', 'SNP'] ) 
	return df_comb
	
	

def compute_prevalance(df_pheno):
	pheno = df_pheno.columns.tolist()[1]
	has_pheno = float(len(set(df_pheno[df_pheno[pheno]==2].FID.tolist()) )  )		
	no_pheno = float(len(set(df_pheno[df_pheno[pheno]==1].FID.tolist()) )  )		

	return has_pheno/(has_pheno+no_pheno)

	

def p12(df_common, snp_set, df_pheno):

	pheno = df_pheno.columns.tolist()[1]
	
	df_snps_cmb = combine_snps(snp_set)
	
	common_exp = len(  set(df_common[df_common['SNP']>0].FID.tolist())  & set(df_pheno[df_pheno[pheno]>1].FID.tolist()) )
	
	common_and_rare_exp = len(   
						set(df_snps_cmb[df_snps_cmb['SNP']>0].FID.tolist())  &
						set(df_common[df_common['SNP']>0].FID.tolist())  &
						set(df_pheno[df_pheno[pheno]>1].FID.tolist())
						)
						
	return float(common_and_rare_exp) / common_exp

