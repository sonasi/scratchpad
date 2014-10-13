#!/usr/bin/env python

'''
File: plinkio_example.py
Created: 22-May-2014 Joseph P. Bochenek
1.0
Description: Demonstrates how to read a plink binary (.bed) file using 
			 libplinkio, and has a few functions for using SNPs and phenotypes
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
from operator import itemgetter


# Initialize plink library and read phenotype file
infile = "/home/vagrant/work/DBMI/data/2_Exome_36K_Ws" 
plink_file = plinkfile.open( infile )

if not plink_file.one_locus_per_row( ):
	 print( "This script requires that snps are rows and samples columns." )
	 exit( 1 )

sample_list = plink_file.get_samples( )
locus_list = plink_file.get_loci( )

# Phenotype  file
phenotypefile = "../data/pheno_Exome_36k_MCC_Ws_MIN2.txt"
# phenotypefile = "../data/pheno_autism_exclutions2.txt"
df_full_pheno = pandas.read_csv(phenotypefile, '\t')



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

		for allele1 in allele_types:	

			p_j = float( len( set(df_snp_1[df_snp_1['SNP']==allele1].FID.tolist() ) \
					  &  set( df_snp_2[df_snp_2['SNP']==allele2].FID.tolist())))  \
				  / total


			p_2 = float( len( set(df_snp_2[df_snp_2['SNP']==allele2].FID.tolist() ) ) ) \
				  / total	
			
			MI_sum += p_j * math.log( p_j / (p_1*p_2), 2 )
	
	return (MI_sum, Entropy_sum)





def combine_snps(args):

	if len(args) < 2:
		return args[0]

	df_ref = args[0]
	df_cat = pandas.concat( args, axis=1 )
	combined = []

	for row in df_cat['SNP'].iterrows():
		combined.append( np.max(row[1]) )

	df_comb = pandas.DataFrame( zip(df_ref.FID, combined), index=df_ref.index, columns = ['FID', 'SNP'] ) 
	return df_comb


def cont_table(snp_name, pheno_name, verbose=False):
	'''
	Make contingency table from SNP name and JD Code 
	
	input:
		snp_name
		pheno_name
	'''
	 
	df_snp = get_SNP_samples(snp_name)
	df_pheno = pheno_from_file(pheno_name)
	return make_cont_table(df_snp, df_pheno, verbose=verbose)
	

def make_cont_table(df_snp, df_pheno, verbose=False):

	pheno = df_pheno.columns.tolist()[1]

	if verbose:
		print "Phenotype: {}\tprev.: {:.3e}".format(pheno, compute_prevalance(df_pheno))

	try:
		maf = df_snp.maf
	except:
		maf_ = [ len(df_snp[df_snp['SNP']==x].SNP.tolist()) for x in range(4)  ]
		maf = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))
		
	try:
		snp_name = df_snp.snp_name
	except:
		snp_name = 'unnamed'

	if verbose:		
		print "SNP: {}\t MAF:{:.3e}".format(snp_name, maf)


	df_logit = pandas.merge(df_snp, df_pheno, on="FID")			
	cross = pandas.crosstab(df_logit[pheno], df_logit['SNP'], rownames=[pheno])
	
	if verbose:
		print cross

	try:
		cross = cross.drop(3, 1)
	except:
		if verbose:
			print 'no error column'
		
	# Discard exclusions
	z = [cross[x].tolist()[1:] for x in cross]
	a = [ [x[0] for x in z], [x[1] for x in z] ]

	if len(a[0])==3:
		x = PrettyTable(["Allele:", "X=0", "X=1", "X=2"])
		x.add_row(["Control", a[0][0], a[0][1], a[0][2]] )
		x.add_row(["Cases", a[1][0], a[1][1], a[1][2]] )
	if len(a[0])==2:
		x = PrettyTable(["Allele:", "X=0", "X=1"])
		x.add_row(["Control", a[0][0], a[0][1] ] )
		x.add_row(["Cases", a[1][0], a[1][1] ] )
	
	return a




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
	pheno_vals = [1, 2]
	pheno = df_pheno.pheno_name
	
	total = float(len(set(df_pheno[df_pheno['pheno_'+pheno]>-1].FID.tolist()) )  )
	
	print "Disease Prevalence: {}".format(compute_prevalance(df_pheno))
	
	MI_sum = 0
	Entropy_sum = 0
	
	for pheno_val in pheno_vals:
	
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
						  &  	  set( df_snp_2[df_snp_2['SNP']==allele2].FID.tolist() ) \
						  &  	  set(df_pheno[df_pheno['pheno_'+pheno]==pheno_val].FID.tolist())
						  ))  \
					  / pheno_total


				p_2 = float( len( set(df_snp_2[df_snp_2['SNP']==allele2].FID.tolist() ) 
								& set(df_pheno[df_pheno['pheno_'+pheno]==pheno_val].FID.tolist()) ) ) \
					  / pheno_total	
			
				print "\t", p_1, p_2, p_j
				
				if not(p_j == 0):			
					MI_sum += p_z * p_j * math.log( p_j / (p_1*p_2), 2 )

				print MI_sum
	
	return (MI_sum, Entropy_sum, 2**MI_sum/2**Entropy_sum)




prev_map = dict()

def compute_prevalance(df_pheno):
	pheno = df_pheno.columns.tolist()[1]
	has_pheno = float(len(set(df_pheno[df_pheno[pheno]==2].FID.tolist()) )  )		
	no_pheno = float(len(set(df_pheno[df_pheno[pheno]==1].FID.tolist()) )  )		

	return has_pheno/(has_pheno+no_pheno)





def pheno_from_file(code):
	'''
	Retreives a phenotype from the phenotype file and returns a dataframe
	
	input: a phenotype code, e.g. '297.1'
	returns: a dataframe with the individuals in the phenotype and their phenotype codes
	'''

	st = time.time()

	cols_to_keep = ['FID', 'pheno_' + code]
# 	df_pheno = pandas.concat([x.ix[:, cols_to_keep] for x in pandas.read_csv(phenotypefile, '\t', chunksize=200, na_values=[-9], na_filter=False)])
	df_pheno = df_full_pheno[cols_to_keep]
	
	pheno_stats = collections.Counter()
	for index, thing in df_pheno.iterrows():
		pheno_stats[thing['pheno_' + code]]+=1

# 	print "Has disease:", pheno_stats[2]
# 	print "Doesn't have disease:", pheno_stats[1]
# 	print "Excluded:", pheno_stats[-9]
# 
# 	print "\nTotal:", len(df_pheno)

	df_pheno.pheno_name = code

# 	print "Execution time: {:.2f}s".format(time.time() - st)

	return df_pheno
	

	 
def pheno_from_db(code):
	'''
	Load Phenotype using an SD query (not finished)
	'''
	query6 = "select * from ICD_CODES where IND_SEQ in (select IND_SEQ from ICD_CODES where CODE like \'714.0%\' group by IND_SEQ having count(*) > 1)"
	df_RA_codes = pandas.io.sql.frame_query(query6, conn)




def get_SNP_pos(SNP):
	'''
	Utility function to get the position of a SNP in the binary file
	'''
	pos = 0
	for locus in locus_list:
		if SNP == locus.name:
			SNP_pos = pos
			return SNP_pos
		pos+=1
	raise ValueError("Error, snp, {} not found in {}.".format(SNP , infile))

	

snp_array_list = dict()
	
def store_snp_sample(snp_list):
	'''
	Plink can't handle exm-rs names in the SNP filter, so let's just filter the SNPs ourselves, stores a 
	list of snparray objects for use with get_SNP_samples_subset
	
	be careful about memory
	'''
	
	positions = []
	snp_list_found = []
	
	for snp_name in snp_list:
		snp_name = snp_name.replace('rs', 'exm-rs')
		try:
			positions.append( get_SNP_pos( snp_name ) )
			snp_list_found.append( snp_name )
		except:
			print "SNP {} not found".format(snp_name)

	snps = zip(positions, snp_list_found)
	snps_sorted = sorted(snps, key=itemgetter(0))
	
	print snps_sorted
		
	cplinkio.reset_row(plink_file.handle)
	
	list_counter = 0

	i = 0
	while plink_file:
		try:
			iter = plink_file.next()
			if i == snps_sorted[list_counter][0]:	
	# 				print "SNP",locus_list[i],"found at", pos
				snp_array_list[ snps_sorted[list_counter][1].replace('exm-rs', 'rs') ] = iter
				list_counter+=1
			i+=1
		except:
			print "done."
			break

	return snp_array_list
	
	



def get_SNP_samples_fast(SNP):

	df_snp = pandas.DataFrame(zip([int(x.iid) for x in sample_list], [x for x in snp_array_list[SNP] ]), columns=['FID', 'SNP'])
	df_snp.snp_name = SNP
# 	print "SNP", SNP, "found.  Allele counts:", snp.allele_counts()

# 	print "Elapsed time: {:.2f}s".format(time.time() - st)

	return df_snp
	
	
	
def get_SNP(SNP):
	'''
	Get the Genotype information for a SNP from the BED file
	
	Description: The BED file is not indexed so we have to find the position of the 
	SNP, then iterate to find the position of the genotype

	Returns: A function containing the Genotype information for a SNP
	'''

	pos = get_SNP_pos( SNP )
# 	plink_file = plinkfile.open( infile )
	cplinkio.reset_row(plink_file.handle)
	
	if pos > -1:
		SNP_info = 0
	
		i = 0
		while plink_file:
			try:
				iter = plink_file.next()
				if i == pos:	
	# 				print "SNP",locus_list[i],"found at", pos
					return iter
				i+=1
			except:
				print "done."
				break
	
	raise ValueError("Error, snp, {} not found in {}.".format(SNP , infile))
	


def get_SNP_samples(SNP):
	'''
	Return the individuals with a SNP and their allele count
	
	input: SNP name, e.g. 'exm12345'
	output: Dataframe with SNP list
	'''

	st = time.time()

	snp = get_SNP(SNP)
	if snp < 0:
		raise ValueError("Error, snp, {} not found in {}.".format(SNP , infile))

	df_snp = pandas.DataFrame(zip([int(x.iid) for x in sample_list], [x for x in snp]), columns=['FID', 'SNP'])
	df_snp.snp_name = SNP
	
	maf_ = snp.allele_counts()
	maf = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))
	df_snp.maf = maf
	
# 	print "SNP", SNP, "found.  Allele counts:", snp.allele_counts()
# 	print "Elapsed time: {:.2f}s".format(time.time() - st)

	return df_snp
	

def compute_MAF(SNP):
	'''
	Compute the MAF frequency for a snp

	Description: this info is precomputed in the binary (e.g. get_SNP('exm1234').count), this function is 
		really just for testing purposes

	input: SNP name, e.g. 'exm12345'
	output: contingency table
	'''
	st = time.time()
	
	maf_ = collections.Counter()
	
	snp = get_SNP(SNP)
	
	print "From plink file:", snp.allele_counts()
	
	if snp < 0:
		raise ValueError("Error, snp, {} not found in {}.".format(SNP , infile))

	for sample, genotype in zip( sample_list, snp ):
		maf_[genotype]+=1
#  		print( "Individual {0} has genotype {1} for snp {2}.".format( sample.iid, genotype, SNP ) )	

	print "computed:", maf_
	print "native:", snp.allele_counts()

	print "Minor Allele Frequency:", (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))

	print "Elapsed time: {:.2f}s".format(time.time() - st)

	print "Total", np.sum(maf_.values())
	print maf_
	return dict(zip(maf_.keys(), maf_.values()))
	

