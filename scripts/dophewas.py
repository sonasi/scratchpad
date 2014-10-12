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
import scipy as sp

sp.disp=False





import plinkio_example as pl

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
	from operator import itemgetter
	snps_sorted = sorted(snps, key=itemgetter(0))
	
	print snps_sorted
		
	cplinkio.reset_row(pl.plink_file.handle)
	
	list_counter = 0
		
	i = 0
	while pl.plink_file:
		try:
			iter = pl.plink_file.next()
			if i == snps_sorted[list_counter][0]:	
	# 				print "SNP",locus_list[i],"found at", pos
				snp_array_list[ snps_sorted[list_counter][1].replace('exm-rs', 'rs') ] = iter
				list_counter+=1
			i+=1
		except:
			print "done."
			break

	return snp_array_list
	
	
import exact_stats as exact
reload(exact)
import genotype_math as gm


def compute_maf(df_snp):
	maf_ = []
	maf_.append(len(df_snp[df_snp['SNP']==0]))
	maf_.append(len(df_snp[df_snp['SNP']==1]))
	maf_.append(len(df_snp[df_snp['SNP']==2]))
	maf_.append(len(df_snp[df_snp['SNP']==3]))
	maf = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))
	return maf
	
def phewas_top_results():

	import MySQLdb
	
# 	conn = MySQLdb.connect (host = '127.0.0.1',
# 							   user = "root",
# 							   passwd = "faf&7_Ew",
# 							   port = 9871)
#                            
# 	cur = conn.cursor()
#     cur.execute("use GWAxPheWAS")             
# 	query = 'select SNP, jd_code, jd_string, odds_ratio, P, MAF_W from Exome_36k_Ws_GAE3_view order by P asc limit 2000'
# 	df_top_pvals = pandas.io.sql.frame_query(query, conn)

	df_top_pvals = pandas.read_pickle('df_top_pvals.pkl')
	snp_list = df_top_pvals.drop_duplicates('SNP').SNP.tolist()
	
	if not(pl.snp_array_list):
		snps_list = list(set(snp_list))
		pl.store_snp_sample(snps_list)
	
	ofile = open('phewas_top_results_full3.txt', 'w')
	ofile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('SNP', 'pheno', 'maf', 'PheWAS_MAF_W', 'prev', 'table', 'minN', 'jd_string', 'PheWAS_OR', 'PheWAS_P', 'logit_OR', 'logit_p-val', 'fisher_OR', 'fisher_pval', 'exact_OR', 'exact_p-val', 'pval_norm', 'bayes factor'))

	grps = df_top_pvals.groupby(['SNP'])

	for snp_name, grp in grps:

		try:
			df_snp = pl.get_SNP_samples_fast( snp_name )
		except:
			print 'SNP', snp_name, 'not found'
			continue
			
		maf = compute_maf(df_snp)

		for iter, item in grp.iterrows():
			pheno_name = 'pheno_' + item['jd_code']

			df_pheno = pandas.DataFrame( zip(pl.df_full_pheno['FID'],   pl.df_full_pheno[pheno_name]), columns=['FID', pheno_name] )

			df_logit = pandas.merge(df_snp, df_pheno, on="FID")			
			cross = pandas.crosstab(df_logit[pheno_name], df_logit['SNP'], rownames=[pheno_name])


			prev = pl.compute_prevalance(df_pheno)
			filename = 'plots/{}_{}'.format(snp_name, pheno_name)

			try:
				cross = cross.drop(3, 1)
			except:
				pass
				
			# Discard exclusions
			z = [cross[x].tolist()[1:] for x in cross]
			a = [ [x[0] for x in z], [x[1] for x in z] ]
	
			print a
			
			result = [0.,0.]		
			try:
				result = gm.logistic_regression(df_pheno, df_snp, loud=False)
			except:
				pass
			result2 = exact.pval_null( a )
			result3 = exact.Bayes_factor( a )
			fisher_res = exact.fisher_exact( a )


			pval_norm = exact.pval_norm( a )
			minN = np.min(exact.normalize_table(a))						



			print snp_name, pheno_name, maf, item['MAF_W'], prev, a, minN, item['jd_string'], item['odds_ratio'], item['P'], result[0], result[1], fisher_res[0], fisher_res[1], result2[0], result2[1], pval_norm, result3[0]
			res = '{}\t{}\t{:.3e}\t{}\t{:.3e}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.3e}\t{:.2f}\t{:.3e}\t{:.2f}\t{:.3e}\t{:.3e}\t{:.2f}\n'.format( snp_name, pheno_name, maf, item['MAF_W'], prev, a, minN, item['jd_string'], item['odds_ratio'], item['P'], result[0], result[1], fisher_res[0], fisher_res[1], result2[0], result2[1], pval_norm[1], result3[0])
			ofile.write( res )

# 			exact.test_plot_with_toys(a, filename=filename)
			



def run_phewas(mpi=1):
	
	i = 0
	j=0
	
	nloci = len(pl.locus_list)
	nphen = len(pl.df_full_pheno.columns)
	nsamp = len(pl.sample_list)
	
	print 'PheWAS on {} loci and {} phenotypes'.format( nloci, nphen )
	print 'Total Tests: {:.2e}'.format(  nloci * nphen )
	
	ofile = open('phewas_results_{}.txt'.format(mpi), 'w')
	ofile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('ind', 'SNP', 'pheno', 'maf', 'prev', 'table', 'minN', 'logit_OR', 'logit_p-val', 'exact_OR', 'exact_p-val', 'bayes factor'))

	pl.cplinkio.reset_row(pl.plink_file.handle)

	while pl.plink_file:
		i+=1
		if i%mpi:
			continue


		iter = pl.plink_file.next()
		maf_ = iter.allele_counts()
		maf = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))
		
		
		print maf_


		
		name = pl.locus_list[i].name
		

# 		if maf_[1] + maf_[2] < 5:
# 			print "Skipping ", name, " only ", maf_[1] + maf_[2], "cases."
# 			continue

		k = 0
				
		df_snp = pandas.DataFrame(zip([int(x.iid) for x in pl.sample_list], [x for x in iter ]), columns=['FID', 'SNP'])
		npheno = len(pl.df_full_pheno.columns[2:])

		for column in pl.df_full_pheno.columns[2:]:

			k+=1
			print "Processing {}/{} phenotypes for {}/{} SNPs.  Current SNP, {} (maf={}).".format(k, npheno, i, len(pl.locus_list), name, maf)
			
			df_pheno = pandas.DataFrame( zip(pl.df_full_pheno['FID'],   pl.df_full_pheno[column]), columns=['FID', column] )			

			prev = pl.compute_prevalance(df_pheno)
			
# 			if maf*prev*nsamp < 1:
# 				print "skipping combination, combined freq = ", maf*prev*nsamp
# 				continue
				
			df_logit = pandas.merge(df_snp, df_pheno, on="FID")			
			cross = pandas.crosstab(df_logit[column], df_logit['SNP'], rownames=[column])
						
			try:
				cross = cross.drop(3, 1)
			except:
				pass
				
			# Discard exclusions
			z = [cross[x].tolist()[1:] for x in cross]
			a = [ [x[0] for x in z], [x[1] for x in z] ]

			result = [0,0]
			result2 = [0,0]
			result3 = [0,0]
			
			j+=1
			
			if not(np.sum(a[1][1:]) == 0) or (a[1][0]==0):
				try:
					result2 = exact.pval_null( a )	
					print result2	
				except:
					pass
				try:
					result3 = exact.Bayes_factor( a )	
					print result3
				except:
					pass

				try:
					result = gm.logistic_regression(df_pheno, df_snp, loud=False)
				except:
					pass
				
				minN = np.min(a)						
				res = '{}\t{}\t{}\t{:.3e}\t{:.3e}\t{}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.2f}\n'.format(j, name, column, maf, prev, a, minN, result[0], result[1], result2[0], result2[1], result3[0])
				ofile.write( res )




def check_gene(snp, pheno):

	df_pheno = pl.pheno_from_file(pheno)
	pl.cplinkio.reset_row(pl.plink_file.handle)
	
	
	snp_file_pos = pl.get_SNP_pos(snp)
	snp_bp_pos = pl.locus_list[snp_file_pos].bp_position
	
	prev_or = 1
	snp_list_full = []
	snp_list_names = []
	p_value_prev = 1
	
	lfile = open('sd_test_log.txt', 'w')
		
	for locus in pl.locus_list:
		iter = pl.plink_file.next()
# 		print locus.name, iter.allele_counts()
		
		if locus.chromosome == 16:
			if abs(locus.bp_position - snp_bp_pos) < 20000:
			

				print '\n', locus.name, locus.bp_position - snp_bp_pos				
				df_snp = pandas.DataFrame(zip([int(x.iid) for x in pl.sample_list], [x for x in iter ]), columns=['FID', 'SNP'])
			
				maf_ = [ len(df_snp[df_snp['SNP']==x].SNP.tolist()) for x in range(4)  ]
				maf = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))

			
				snp_list_test = []
				snp_list_test += snp_list_full
				snp_list_test.append(df_snp)
				df_snp_cmb = pl.combine_snps(snp_list_test)

				a = pl.make_cont_table(df_snp, df_pheno)
				print exact.pval_null(a)
				odds_ratio, p_value = exact.pval_null(a)

				print '\t', a
				print '\t this SNP:', odds_ratio, p_value

				a = pl.make_cont_table(df_snp_cmb, df_pheno)
				odds_ratio2, p_value2 = exact.pval_null(a)
				print '\t', a
				print '\t combination ({}): {}\t{}'.format(snp_list_names, odds_ratio2, p_value2)
				
				if (p_value < p_value_prev): # and (odds_ratio > prev_or)
					prev_or = odds_ratio
					p_value_prev = p_value
					snp_list_full.append(df_snp)
					snp_list_names.append(locus.name)
					print "\tADDING SNP:", snp_list_names



				pheno = df_pheno.columns.tolist()[1]
				lfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(locus.name, pheno, (locus.bp_position - snp_bp_pos), odds_ratio, p_value, odds_ratio2, p_value2   )    )
				
				df_logit = pandas.merge(df_snp, df_pheno, on="FID")			
				cross = pandas.crosstab(df_logit[pheno], df_logit['SNP'], rownames=[pheno])
	
				verbose=False		
				try:
					cross = cross.drop(3, 1)
				except:
					if verbose:
						print 'no error column'
				
				# Discard exclusions
				z = [cross[x].tolist()[1:] for x in cross]
				a = [ [x[0] for x in z], [x[1] for x in z] ]
				print a




def sa_by_region(snp_name, pheno, d_bp = 20000):

	df_pheno = pl.pheno_from_file(pheno)
	pheno = df_pheno.columns.tolist()[1]

	pl.cplinkio.reset_row(pl.plink_file.handle)
	
	snp_file_pos = pl.get_SNP_pos(snp_name)
	snp_bp_pos = pl.locus_list[snp_file_pos].bp_position
	chromosome = pl.locus_list[snp_file_pos].chromosome

	snp_list = []
	
	df_snp_common = pl.get_SNP_samples(snp_name)
	maf_ = [ len(df_snp_common[df_snp_common['SNP']==x].SNP.tolist()) for x in range(4)  ]
	maf_common = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))
	
	print 'Searching for SA with', snp_name, '\tmaf:', maf_common


	a1 = pl.make_cont_table(df_snp_common, df_pheno)
	print a1
	print exact.pval_null(a1)
	odds_ratio1, p_value1 = exact.pval_null(a1)
	print 'odds ratio: {}\t p-value: {}\n'.format( odds_ratio1, p_value1)


	pl.cplinkio.reset_row(pl.plink_file.handle)


	print 'Starting scan'
	lfile = open('sd_test_log.txt', 'w')
	
	snp_list_desc = []



	for locus in pl.locus_list:

		iter = pl.plink_file.next()
		maf_ = iter.allele_counts()
		maf = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))



		if locus.chromosome == chromosome:
			if abs(locus.bp_position - snp_bp_pos) < d_bp:

				df_snp = pandas.DataFrame(zip([int(x.iid) for x in pl.sample_list], [x for x in iter ]), columns=['FID', 'SNP'])
			
# 				maf_ = [ len(df_snp[df_snp['SNP']==x].SNP.tolist()) for x in range(4)  ]
# 				maf_ = iter.allele_count()
				
				a = pl.make_cont_table(df_snp, df_pheno)
				if (maf > 0.05) or (maf > 0.4*maf_common):
					continue
				if a==False:
					odds_ratio, p_value = -1, -1
				else:
					odds_ratio, p_value = exact.pval_null(a)

				print 'Testing SNP: {}, maf: {},  d-bp: {}, cont. table: {}, OR {}, p-val {}'.format( locus.name, maf, locus.bp_position - snp_bp_pos, a, odds_ratio, p_value )
	


				if p_value < 0.05 and p_value > 0:
					print 'adding SNP'
					snp_list.append(df_snp)
					snp_list_desc.append([locus.name, maf, odds_ratio, p_value, a])
			
	
	
	
	df_snp_comb = pl.combine_snps(snp_list)
	maf_ = [ len(df_snp_comb[df_snp_comb['SNP']==x].SNP.tolist()) for x in range(4)  ]
	maf_common = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))

	overlap_matrix = np.array([[0,0,0],[0,0,0]])

	for i in [1, 2]:
		for j in [0, 1, 2]:
			common_and_rare_exp = len(   set(df_snp_comb[df_snp_comb['SNP']==j].FID.tolist())   & \
					 set(df_snp_common[df_snp_common['SNP']==j].FID.tolist())  & \
					 set(df_pheno[df_pheno[pheno]==i].FID.tolist()) \
				 )
			overlap_matrix[i-1][j] = common_and_rare_exp


	a2 = pl.make_cont_table(df_snp_comb, df_pheno)


	odds_ratio, p_value = exact.pval_null(a2)	


	print "Rare SNP List:"
	for desc in snp_list_desc:
		print '\t'.join([str(x) for x in desc])
	
	print exact.normalize_table(a1) 
	print exact.normalize_table(a2)

	print 'Common Matrix'
	print a1
	print 'Combined Matrix'
	print a2
	print 'Overlap Matrix'
	print overlap_matrix




	om = overlap_matrix
	remainder_matrix = [ [ (a1[0][0] + om[0][1] + om[0][2]) , (a1[0][1] - om[0][1]) , (a1[0][2] - om[0][2]) ], [ (a1[1][0] + om[1][1] + om[1][2]) , (a1[1][1] - om[1][1]) , (a1[1][2] - om[1][2]) ] ]

	print 'Common SNPs odds ratio: {}\t p-value: {}\t maf_cmb: {}\n'.format( odds_ratio1, p_value1, maf_common)
	print 'Combinde SNPs odds ratio: {}\t p-value: {}\t maf_cmb: {}\n'.format( odds_ratio, p_value, maf_)
	print 'Remainder OR, p-val', exact.pval_null(remainder_matrix)
	
	return 
	

	