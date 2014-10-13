#!/usr/bin/env python

'''
File: synth_assoc_analysis.py
Created: 24-Aug-2014 Joseph P. Bochenek
1.0
Description: Scan Exome dataset for synthetic association by collapsing all causal rare variants in the same 
region as a common variant in to a 'superSNP'.  That is, all low MAF variants in the same gene(s) as the common variant which have 
some association with the same phenotype are combined inclusively.  The method is similar to Combined Collapsing (
(http://varianttools.sourceforge.net/Association/CMC)
'''

from plinkio import plinkfile
import pandas as pd
import plinkio_example as pl
import exact_stats as exact
import numpy as np
import math 
import re
from itertools import izip
import bisect
from scipy.stats import linregress
import genotype_math as gm
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
import rpy2.robjects.numpy2ri
from string import find

rpy2.robjects.numpy2ri.activate()

non_decimal = re.compile(r'[^\d.]+')

query = "select * from (select catalog_UP_brief.Chr, catalog_UP_brief.SNP, catalog_traits.jd_code, Trait, Reported_Gene, Mapped_gene, p_value, OR_or_beta_text, CI_95, Intergenic from catalog_UP_brief join catalog_traits on catalog_traits.trait = catalog_UP_brief.Disease_trait) as t1 join GWAxPheWAS.Exome_36k_Ws_GAE3_view on GWAxPheWAS.Exome_36k_Ws_GAE3_view.SNP = t1.SNP and GWAxPheWAS.Exome_36k_Ws_GAE3_view.jd_code = t1.jd_code"

CHR = 16
d_bp = 25000

import utils
data_dir = utils.data_dir

# df_genemap = pd.read_csv('Homo_sapiens.GRCh37.73.gtf.gz', compression='gzip', sep='\t', names=['id', 'processed_transcript', 'exon', 'start', 'end', 'other', 'strand', 'other2', 'gene_id', 'transcript_id', 'exon_number', 'gene_name', 'gene_biotype', 'pseudogene'], index_col=False)
# df_genes = df_genemap[df_genemap['gene_biotype']=='protein_coding'].groupby('gene_name')
# df_gene_reg = df_genes['start'].agg([np.min, np.max])
# df_cat_assoc_map = pd.read_pickle('~/work/Rare_gene_NHGS_project/scripts/df_cat_assoc_map.pkl')
# select  catalog_UP_brief.Chr, count( catalog_UP_brief.Chr) from catalog_UP_brief join catalog_traits on catalog_traits.trait = catalog_UP_brief.Disease_trait where catalog_traits.match_available = 1 group by catalog_UP_brief.Chr
# Read catalog studies


# df_genemap = pd.read_csv('/home/vagrant/data/Homo_sapiens.GRCh37.73.gtf.gz', compression='gzip', sep='\t', names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'exon_number', 'gene_name', 'gene_biotype', 'pseudogene'], index_col=False)
# df_genemap = df_genemap.drop_duplicates(['seqname', 'start', 'end', 'gene_name'])

def load():
	global df_genes
	global df_demos
	df_genemap = pd.read_pickle('../data/df_genemap.pkl')
	df_genes = df_genemap[df_genemap['gene_biotype']=='protein_coding'].groupby(['gene_name', 'seqname'])
	df_demos = pd.read_csv("../data/demos_Ws_view.tsv", '\t')

# df_gene_reg = df_genes['start'].agg([np.min, np.max])
# df_genes.to_pickle('df_genes.pkl')
# df_gene_reg = pd.read_pickle('df_gene_map_cleaned.pkl')
# df_genes = pd.read_pickle('df_genes.pkl')

class Skat_Wrapper:

	def __init__(self):
		
		self.skat = importr("SKAT", robject_translations = {"Beta.Weights": "Beta_Weights2", "to.period": "to_period2"})

	def p_value(self, X, y, Z, mode="D"):

		form = Formula('y~X')

		env = form.environment
		env['y'] = y
		env['X'] = X
		env['Z'] = Z

		obj = self.skat.SKAT_Null_Model(form, out_type=mode)
		res = self.skat.SKAT(Z, obj, r_corr=0.)

		return res[0][0]

def inrange( x, start, end ):
	if x > end or x < start:
		return False
	else:
		return True

def run_continuous(df, snps, snps_cmn, chr, method='gene'):
	'''
	Scan exome using lab phenotypes using linear regression.
	
	arguments:
		df - pd dataframe of 
		snps - table of snps in chromosome
		snps_cmn - table of common snps in df
		chr - do one chromosome at a time (memory managment)
	
	'''

	lfile = open('sd_labs_genebased_chr_{}.txt'.format(chr), 'w')

	data =	[ 'SNP', 'Reported_Gene', 'Chr', 'Lab_name', 'Study', 'Trait', 'p-value_reported (Catalog)', 'beta_reported (Catalog)', 'MAF_reported (catalog)' , 'beta_common (Exome)', 'p-value_common (Exome)', 'MAF_common (exome)', 'beta_combined', 'p-value_combined', 'MAF_combined', 'SNPs_name', 'SNPs_MAF', 'SNPS_beta' , 'SNPS_p-val' ]

 	lfile.write(
 	'\t'.join([str(x) for x in data]) + '\n'
 	)

	print "Continuous Phenotypes, Chromosome", chr, len(df)
	
	positions = []
	for key, val in snps.iteritems():
		bisect.insort(positions, int(val[0]))
	position_to_snp = dict( izip([int(val[0]) for val in snps.values()], snps.keys()))
	k=0
	for iter, row in df.iterrows():	
		k+=1
		lab_name, snp_name, chr = row['Lab_name'], row['SNP.1'], row['Chr']
		p_val_reported, OR_reported, MAF_reported = row['p_value'], row['OR_or_beta'], row['Risk_Allele_Frequency']
		study, trait = row['Study'], row['Disease_Trait']

		print '\n Starting scan pair', k, '/', len(df), '--', lab_name, snp_name

		try:		
			df_labs = pd.io.parsers.read_csv("../data/pheno_Ws_no_lab_count_req/pheno_{}_Exome_demos_Ws_Age18.txt".format(lab_name), sep='\t', names=['FID', 'IID', 'Mean', 'Median', 'Min', 'Max', 'Ever_high', 'Ever_low', 'First'])	
		except:
			print 'The data for the lab,', lab_name, ', doesn\'t seem to exist.'
			continue
			
		labs_vals = df_labs.Mean[1:].values
		
		inExome = True
		if (str(row['SNP'])=='nan'):
			inExome = False
		
		snp_list_desc = []
		snps_selected = []

		if (method=='gene'):
			
			reported_genes = [gene_str.strip() for gene_str in row['Reported_Gene'].split(',')]
 			print 'Reported Genes', reported_genes
			
			for gene in reported_genes:	
				try:
					gene_definitions = df_genes.get_group((gene,chr))
					for ind, row in gene_definitions.iterrows():
						min = row['start'] - 50000
						max = row['end'] + 50000
						
						if not(inrange( snp_bp_pos, min, max )):
							continue
					
# 					min = df_gene_reg['amin'][gene]
# 					max = df_gene_reg['amax'][gene]

# 					print '\t (min, max)', min, max
						left = bisect.bisect_left(positions, min)
						right = bisect.bisect_right(positions, max)

						snps_selected += [position_to_snp[pos] for pos in positions[left:right]]
					
				except:
					'Gene', gene, 'not found.'
				
			snps_selected = list(set(snps_selected))

		else:
			for name, locus in snps.iteritems():

				bp_position = locus[0]

				if abs(bp_position - snp_bp_pos) < d_bp:
					snps_selected.append(name)

		
		snp_bp_pos = 0

		if inExome:
			snp_name_common, maf_common = row['SNP'], float(row['MAF_W'])
			reported_beta  = row['mean_Beta']
			reported_p_val  = row['mean_P']

			if find(snp_name_common, 'rs') > -1:
				snp_name_common = 'exm-'+snp_name_common
			locus = snps_cmn[snp_name_common]
			print locus 
			df_snp_common = pd.DataFrame(zip([x.iid for x in pl.sample_list], [x for x in locus[2] ]), columns=['FID', 'SNP'])
			df_pheno_geno = pd.merge(df_snp_common, df_labs, on='FID')
			
			print 'common SNP', reported_beta, reported_p_val
			print 'reported SNP', OR_reported, p_val_reported,  MAF_reported
 			print 'Linear Regression (numpy)', np.polyfit( [int(x) for x in df_pheno_geno.SNP], [float(x) for x in df_pheno_geno.Mean], 1 )
			beta_c, intercept_c, r_value_c, p_value_c, stderr_c = linregress([int(x) for x in df_pheno_geno.SNP], [float(x) for x in df_pheno_geno.Mean])
			print 'Linear Regression (scipy)', beta_c, p_value_c
			snp_bp_pos = locus[0]

# 		print snps_selected
# 		print 'inExome', inExome
		else:
			snp_name_common, maf_common = '-', '-'
			beta_c, intercept_c, r_value_c, p_value_c, stderr_c = '-', '-', '-', '-', '-'
			
		snp_list = []

		print snps_selected

		for name in snps_selected:				

			locus = snps[name]
			
			
		
			bp_position, maf, iter = int(locus[0]), float(locus[1]), locus[2]
			
			if inExome:
				maf_compare = maf_common
				beta_compare = beta_c
				if (maf > 0.05) or (maf > 0.4*maf_common):
					continue			
			else:
				maf_compare = MAF_reported
				beta_compare = OR_reported
				if (maf > 0.05) or (maf > 0.4*float(MAF_reported)):
					continue			


			df_snp = pd.DataFrame(zip([x.iid for x in pl.sample_list], [x for x in locus[2] ]), columns=['FID', 'SNP'])
			df_pheno_geno = pd.merge(df_snp, df_labs, on='FID')
			beta, intercept, r_value, p_value, stderr = linregress([int(x) for x in df_pheno_geno.SNP], [float(x) for x in df_pheno_geno.Mean])

			# only include if they are in the right direction.

			print 'Testing SNP: {}\t pos: {}\t maf: {}/{}, \t beta: {}/{} \t p-val:{}'.format( name, bp_position - snp_bp_pos, maf, maf_compare, beta, beta_compare, p_value )

			if math.log(maf_compare / 0.5) > 0:
				print 'Flipping the sign for major/minor Allele freq.   '
				
				
			if beta * beta_compare * math.log(maf_compare / 0.5) > 0:
				continue

			print beta, beta_compare, p_value


			if p_value < 0.05 and p_value > 0:
				print 'adding SNP'
				snp_list.append(df_snp)
				snp_list_desc.append([name, maf, beta, p_value])

		if len(snp_list) == 0:
			continue
	
		df_snp_comb = pl.combine_snps(snp_list)
		maf_ = [ len(df_snp_comb[df_snp_comb['SNP']==x].SNP.tolist()) for x in range(4)  ]
		maf_comb = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))
		
		df_pheno_geno = pd.merge(df_snp, df_labs, on='FID')
		beta_n, intercept_n, r_value_n, p_value_n, stderr_n = linregress([int(x) for x in df_pheno_geno.SNP], [float(x) for x in df_pheno_geno.Mean])


		print "\n------------------\nRare SNP List:"
		print 'SNP\tMAF\tOR\tp-val\tmatrix'
		for desc in snp_list_desc:
			print '\t'.join([str(x) for x in desc])
		print

		print '\nCommon SNPs odds ratio: {}\t p-value: {}\t maf_cmb: {}'.format( beta_c, p_value_c, maf_common)
		print 'Combinde SNPs odds ratio: {}\t p-value: {}\t maf_cmb: {}'.format( beta_n, p_value_n, maf_comb)
 		data = [snp_name, row['Reported_Gene'], chr, lab_name, study, trait, p_val_reported, OR_reported, MAF_reported, beta_c, p_value_c, maf_common, beta_n, p_value_n, maf_comb, ','.join([str(z[0]) for z in snp_list_desc]), ','.join([str(z[1]) for z in snp_list_desc]), ','.join([str(z[2]) for z in snp_list_desc]), ','.join([str(z[3]) for z in snp_list_desc]) ]
 		lfile.write( '\t'.join( [str(x) for x in data] ) + '\n') 
 	lfile.close()

def getdata2(chr, df_catalog):

	snps_cmn = dict()
	snp_list_cmn = ['exm-'+str(x) for x in df_catalog.SNP.tolist()]

	infile ="{}/2_Exome_36K_Ws_chr{}".format(data_dir, chr) 
	plink_file = plinkfile.open( infile )
	locus_list = plink_file.get_loci( )
	total = 0
	
	for locus in locus_list:
		total+=1

		iter = plink_file.next()

		if (locus.name in snp_list_cmn) and (locus.chromosome == chr):
			maf_ = iter.allele_counts()
			maf = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))
			if not(total%10000): print total, maf

			snps_cmn[locus.name] = [int(locus.bp_position), maf, iter]

	return snps_cmn

def getdata(chr):

	total = 0
	sel = 1
	infile = "{}/2_Exome_36K_Ws_chr{}".format(data_dir, chr) 

	print infile
	plink_file = plinkfile.open( infile )
# 	sample_list = plink_file.get_samples( )
	locus_list = plink_file.get_loci( )

# 	pl.cplinkio.reset_row(pl.plink_file.handle)
	snps = dict()
	for locus in locus_list:

		total+=1
		if not(total%10000): print total, sel

		if (locus.chromosome == chr):
	
			iter = plink_file.next()
			maf_ = iter.allele_counts()
			maf = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))
			if not(total%10000): print total, sel, maf
			if not(sel%1000): print total, sel, maf
			sel+=1
			snps[locus.name] = [int(locus.bp_position), maf, iter]
			
	return snps


def loop_pairs(df_catalog):

	for chr in range(1, 24):
		print chr
		

def dopair_lab( gene, chr, lab_name, snps, snp_cmn=False ):

	try:	
		df_labs = pd.io.parsers.read_csv("../data/labwas/data/pheno_Ws_no_lab_count_req/pheno_{}_Exome_demos_Ws_Age18.txt".format(lab_name), sep='\t', names=['FID', 'IID', 'Mean', 'Median', 'Min', 'Max', 'Ever_high', 'Ever_low', 'First'], skiprows=1)
		df_labs['FID'] = df_labs['FID'].astype(int)
	except:
		print 'The data for the lab,', lab_name, ', doesn\'t seem to exist.'
		return 
	sw = Skat_Wrapper()
	snp_list_all = []

	positions = []
	for key, val in snps.iteritems():
		bisect.insort(positions, int(val[0]))
	position_to_snp = dict( izip([int(val[0]) for val in snps.values()], snps.keys()))

	snps_list = []
	snps_selected = []

	try:
		print 'Trying', gene, chr
		gene_definitions = df_genes.get_group((gene,chr))

		print "Found", len(gene_definitions), "snps"
		
		for ind1, row1 in gene_definitions.iterrows():
			min = row1['start'] - 40000
			max = row1['end'] + 40000

			print "row", ind1, min, max
		
# 			if not(inrange( snp_bp_pos, min, max )):
# 				print "\tSNP not in gene", snp_bp_pos, min, max
# 				continue
	
			left = bisect.bisect_left(positions, min)
			right = bisect.bisect_right(positions, max)

			print "SNPS:", [position_to_snp[pos] for pos in positions[left:right]]

			snps_selected += [position_to_snp[pos] for pos in positions[left:right]]

		snps_selected = list(set(snps_selected))
		print 'selected', len(snps_selected), 'snps'

	except:
		'Gene', gene, 'not found.'
		

	

	for name in snps_selected:				

		locus = snps[name]
	
		bp_position, maf, iter = locus[0], locus[1], locus[2]

		print bp_position, maf

		if (maf > 0.05):
			continue			
				
		df_snp = pd.DataFrame(zip([int(x.iid) for x in pl.sample_list], [x for x in iter ]), columns=['FID', 'SNP'])

		df_pheno_geno = pd.merge(df_snp, df_labs, on='FID')
		beta, intercept, r_value, p_value, stderr = linregress([int(x) for x in df_pheno_geno.SNP], [float(x) for x in df_pheno_geno.Mean])

		print name, beta, intercept, r_value, p_value, stderr
		
		snp_list_all.append(df_snp)
		
		if p_value < 0.005:
			snps_list.append(df_snp)

	
	skat_pval = -1

	# run skat
	current = snp_list_all[0].rename(columns={'SNP':'SNP_1'})
	if len(snp_list_all) > 1:
		for i, frame in enumerate(snp_list_all[1:], 2):
			current = current.merge(frame, on='FID').rename(
									 columns={'SNP':'SNP_%d' % i})
	current = pd.merge(current, df_labs[['FID','Median']], on='FID')
	
	current = pd.merge(df_demos[['FID', 'min_age', 'gender']], current, on="FID")
	current['gender'] = current['gender'].map({'M': 1,'F': 0})
	
	y = np.array([x-1 for x in current[current.columns[-1]].tolist()])
	X = np.array([tuple[2:4] for tuple in current.itertuples() ])
	Z = np.array([tuple[4:-1] for tuple in current.itertuples() ])
	skat_pval = sw.p_value(X, y, Z, mode="C") 

	print "SKAT:", skat_pval 

	print len(snps_list), 'rare SNPs with MAF < 0.05'

	df_snp_comb = pl.combine_snps(snps_list)
	df_pheno_geno = pd.merge(df_snp_comb, df_labs, on='FID')
	res_t = linregress([int(x) for x in df_pheno_geno.SNP], [float(x) for x in df_pheno_geno.Mean])

	print "Regression Combined:", res_t[0], res_t[3],

	# Get the common one now
	if snp_cmn:
		locus = snps['exm-'+snp_cmn]
		bp_position, maf, iter = locus[0], locus[1], locus[2]
		df_snp_cmn = pd.DataFrame(zip([int(x.iid) for x in pl.sample_list], [x for x in iter ]), columns=['FID', 'SNP'])
		df_pheno_geno = pd.merge(df_snp_cmn, df_labs, on='FID')
		res_c = linregress([int(x) for x in df_pheno_geno.SNP], [float(x) for x in df_pheno_geno.Mean])
		print "\nRegression Common:", res_c[0], res_c[3]

	if snp_cmn:
		snp_rem_temp = pd.merge(df_snp_cmn.rename(columns={'SNP':'SNP_cmn'}), df_snp_comb.rename(columns={'SNP':'SNP_cmb'}), on='FID')
		Z = []
		for x, y in zip(snp_rem_temp.SNP_cmn.tolist(), snp_rem_temp.SNP_cmb.tolist()):
			z = x
			if y <= x:
				z = x-y
			Z.append(z)
		snp_rem = pd.DataFrame(zip(snp_rem_temp.FID.tolist(), Z) , columns=['FID', 'SNP'])
		print "\nCombined minor alleles:", len(df_snp_comb[df_snp_comb['SNP']==1]) + 2*len(df_snp_comb[df_snp_comb['SNP']==2])
		print "Remainder minor alleles:", len(snp_rem[snp_rem['SNP']==1]) + 2*len(snp_rem[snp_rem['SNP']==2])
		df_pheno_geno = pd.merge(snp_rem, df_labs, on='FID')
		res_r = linregress([int(x) for x in df_pheno_geno.SNP], [float(x) for x in df_pheno_geno.Mean])
		print "Regression Remainder:", res_r[0], res_r[3]
		
	return current




def runloop(df_catalog, snps, snps_cmn, chr, method='gene'):
	'''
	Scan exome using categorical phenotypes using fisher or exact test.
	
	arguments:
		df - pd dataframe of 
		snps - table of snps in chromosome
		snps_cmn - table of common snps in df
		chr - do one chromosome at a time (memory managment)
	
	'''
	lfile = open('sd_cond_genebased_chr_{}_v3.txt'.format(chr), 'w')
	data =['SNP_name', 'Phenotype', 'Chr', 'Disease_Trait', 'Sample Size (Reported)', 'p-val_rep', 'OR_rep', 'MAF_rep', 'Gene_rep', 'OR_common', 'p-val_common', 'MAF_common', 'OR_combined', 'p-val_combined', 'MAF_combinded', 'OR_remainder', 'p-val_remainder', 'P_12', 'SKAT_p-val',  'Logit_no-covariates', 'Logit_conditional', 'Logit_nonconditional', 'SNPS_rare', 'MAFs_rare_snps', 'p-vals_rare_snps', 'ORs_rare_snps',  'A_common', 'A_cmb', 'A_rem' ]
	lfile.write( '\t'.join( [str(x) for x in data] ) + '\n' )


	positions = []
	for key, val in snps.iteritems():
		bisect.insort(positions, int(val[0]))
	position_to_snp = dict( izip([int(val[0]) for val in snps.values()], snps.keys()))

	k = 0
	print 'Searching', len(df_catalog), 'pairs in chromosome', chr
	sw = Skat_Wrapper()

	for iter, row in df_catalog.iterrows():	
		k+=1
		pheno = row['jd_code']
		snp_name = row['SNP']
		snp_cmn_pos = row['Chr_pos']
	
		print '\n------\nStarting scan pair', k, '/', len(df_catalog), '--', pheno, snp_name

		try:
			df_pheno = pl.pheno_from_file(str(pheno))
		except:
			continue
			
		pheno = df_pheno.columns.tolist()[1]
	
		snp_bp_pos = snps_cmn['exm-'+snp_name][0]
		chromosome = chr

		snp_list = []

		iter = snps_cmn['exm-'+snp_name][2]	
		df_snp_common = pd.DataFrame(zip([int(x.iid) for x in pl.sample_list], [x for x in iter ]), columns=['FID', 'SNP'])
		maf_common = snps_cmn['exm-'+snp_name][1]
	
		print 'Searching for SA with', snp_name, '\tmaf:', maf_common

		a1 = pl.make_cont_table(df_snp_common, df_pheno)
		print a1
		print exact.pval_null(a1)
		odds_ratio1, p_value1 = exact.pval_null(a1)
		print 'odds ratio: {}\t p-value: {}\n'.format( odds_ratio1, p_value1)


		snp_list_desc = []

		snp_list_desc = []
		snps_selected = []

		if (method=='gene'):

			reported_genes = [gene_str.strip() for gene_str in row['Reported_Gene'].split(',')]
			mapped_genes = [gene_str.strip() for gene_str in row['Mapped_gene'].split(' - ')]
			
 			print 'Reported Genes, Mapped_genes', reported_genes, mapped_genes
 			
 			gene_list = list(set(reported_genes + mapped_genes))
			
# 			for gene in reported_genes:	
# 				try:
# 					min = df_gene_reg['amin'][gene]
# 					max = df_gene_reg['amax'][gene]
# 
#  					print '\t (min, max)', min, max
# 					left = bisect.bisect_left(positions, min)
# 					right = bisect.bisect_right(positions, max)
# 
# 					snps_selected += [position_to_snp[pos] for pos in positions[left:right]]
# 					
# 				except:
# 					'Gene', gene, 'not found.'




			for gene in gene_list:	
				try:
					print 'Trying', gene, chr
					gene_definitions = df_genes.get_group((gene,chr))
					
					for ind1, row1 in gene_definitions.iterrows():
						min = row1['start'] - 40000
						max = row1['end'] + 40000
					
						if not(inrange( snp_bp_pos, min, max )):
							print "\tSNP not in gene", snp_bp_pos, min, max
							continue
				
						left = bisect.bisect_left(positions, min)
						right = bisect.bisect_right(positions, max)

						print "SNPS:", [position_to_snp[pos] for pos in positions[left:right]]

						snps_selected += [position_to_snp[pos] for pos in positions[left:right]]
				except:
					'Gene', gene, 'not found.'
					

				
			snps_selected = list(set(snps_selected))





		else:
			for name, locus in snps.iteritems():

				bp_position = locus[0]

				if abs(bp_position - snp_bp_pos) < d_bp:
					snps_selected.append(name)

# 		print 'Selected:', snps_selected

		for name in snps_selected:				

			locus = snps[name]
		
			bp_position, maf, iter = locus[0], locus[1], locus[2]
			
			if (maf > 0.05) or (maf > 0.4*maf_common):
				continue			
					
			df_snp = pd.DataFrame(zip([int(x.iid) for x in pl.sample_list], [x for x in iter ]), columns=['FID', 'SNP'])
	
# 				maf_ = [ len(df_snp[df_snp['SNP']==x].SNP.tolist()) for x in range(4)  ]
# 				maf_ = iter.allele_count()
		
			a = pl.make_cont_table(df_snp, df_pheno)

			if a==False:
				odds_ratio, p_value = -1, -1
			else:
# 					odds_ratio, p_value = exact.pval_null(a)
				odds_ratio, p_value = exact.fisher_exact(a)	

			print 'Testing SNP: {}, maf: {},  d-bp: {}, cont. table: {}, OR {}, p-val {}'.format( name, maf, bp_position - snp_bp_pos, a, odds_ratio, p_value )

			# only include if they are in the right direction.
			if (odds_ratio1 <= 0) or (odds_ratio <= 0):
				continue
			if math.log(odds_ratio1) * math.log(odds_ratio) * math.log(maf_common / 0.5) > 0:
				continue

			if p_value < 0.05 and p_value > 0:
				print 'adding SNP'
				snp_list.append(df_snp)
				snp_list_desc.append([name, maf, odds_ratio, p_value, a])
	

# 		pd.merge(snp_list, on='FID')

# 		
		

		if len(snp_list) == 0:
			continue

		skat_pval = -1

		# run skat
		if len(snp_list) > 1:
			current = snp_list[0].rename(columns={'SNP':'SNP_1'})
			for i, frame in enumerate(snp_list[1:], 2):
				current = current.merge(frame, on='FID').rename(
										 columns={'SNP':'SNP_%d' % i})
			current = pd.merge(current, df_pheno.drop(df_pheno[df_pheno[df_pheno.columns[-1]]==-9].index.tolist()), on='FID')
# 			current = current.drop(current[current['pheno_446.4']==-9].index.tolist())
			
			current = pd.merge(df_demos[['FID', 'min_age', 'gender']], current, on="FID")
			current['gender'] = current['gender'].map({'M': 1,'F': 0})
			
			y = np.array([x-1 for x in current[current.columns[-1]].tolist()])
			X = np.array([tuple[2:4] for tuple in current.itertuples() ])
			Z = np.array([tuple[4:-1] for tuple in current.itertuples() ])
			try:
				skat_pval = sw.p_value_bin(X, y, Z) 
			except:		
				skat_pval = -1

				
	
		df_snp_comb = pl.combine_snps(snp_list)
		maf_ = [ len(df_snp_comb[df_snp_comb['SNP']==x].SNP.tolist()) for x in range(4)  ]
		maf_comb = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))

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


		print "\n------------------\nRare SNP List:"
		print 'SNP\tMAF\tOR\tp-val\tmatrix'
		for desc in snp_list_desc:
			print '\t'.join([str(x) for x in desc])
		print
	
		print exact.normalize_table(a1) 
		print exact.normalize_table(a2)

		print 'Common Matrix'
		print a1
		print 'Combined Matrix'
		print a2
		print 'Overlap Matrix'
		print [list(x) for x in overlap_matrix]
		
		P_12 = float(np.sum(overlap_matrix[1][1:]))/float(np.sum(a1[1][1:]))

		om = overlap_matrix
		remainder_matrix = [ [ (a1[0][0] + om[0][1] + om[0][2]) , (a1[0][1] - om[0][1]) , (a1[0][2] - om[0][2]) ], [ (a1[1][0] + om[1][1] + om[1][2]) , (a1[1][1] - om[1][1]) , (a1[1][2] - om[1][2]) ] ]
		print 'Remainder Matrix'
		print remainder_matrix

		print '\nCommon SNPs odds ratio: {}\t p-value: {}\t maf_cmb: {}'.format( odds_ratio1, p_value1, maf_common)
		print 'Combinde SNPs odds ratio: {}\t p-value: {}\t maf_cmb: {}'.format( odds_ratio, p_value, maf_comb)
		or2, pval2 = exact.pval_null(remainder_matrix)
		print 'Remainder  odds ratio: {}\t p-value: {}'.format(or2, pval2)

		try:
			lr_nocovar = gm.logistic_regression(df_pheno, df_snp_common, drop_covariates=True)
			lr_cond = gm.logistic_regression(df_pheno, df_snp_common, df_snp_comb)
			lr_noncond = gm.logistic_regression(df_pheno, df_snp_common)
		except:
			lr_nocovar = -1
			lr_cond = -1
			lr_noncond = -1

		data = [snp_name, pheno, chr, row['Disease_Trait'], row['Initial_Sample_Size'], row['p_value'], row['OR_or_beta_text'], row['Risk_Allele_Frequency'], ','.join(gene_list), odds_ratio1, p_value1, maf_common, odds_ratio, p_value, maf_comb, or2, pval2, P_12, skat_pval, lr_nocovar, lr_cond, lr_noncond, ','.join([str(z[0]) for z in snp_list_desc]), ','.join([str(z[1]) for z in snp_list_desc]), ','.join([str(z[2]) for z in snp_list_desc]), ','.join([str(z[3]) for z in snp_list_desc]), a1, a2, remainder_matrix ]
		lfile.write( '\t'.join( [str(x) for x in data] ) + '\n') 
	lfile.close()


def runloop_nocond(df_catalog_full, snps, chr, method='gene'):
	
	lfile = open('sd_nocond_genebased_chr_{}_v3.txt'.format(chr), 'w')

	data =['SNP_name', 'Phenotype', 'Chr', 'Disease_Trait', 'Sample Size', 'p-val_rep', 'OR_rep', 'MAF_rep', 'Gene_rep', 'OR_combined', 'p-val_combined', 'MAF_combined', 'P_12_max', 'SKAT_p-val',  'Logit_no-covariates', 'Logit', 'SNPS_rare', 'MAFs_rare_snps', 'p-vals_rare_snps', 'ORs_rare_snps',  'A_cmb' ]

	lfile.write( '\t'.join( [str(x) for x in data] ) + '\n' )


	positions = []
	for key, val in snps.iteritems():
		bisect.insort(positions, int(val[0]))
	position_to_snp = dict( izip([int(val[0]) for val in snps.values()], snps.keys()))
	
	sw = Skat_Wrapper()	

	k = 0
	print 'Searching', len(df_catalog_full), 'pairs in chromosome', chr

	for iter, row in df_catalog_full.iterrows():	

# 		############
# 		if iter > 20:
# 			break
# 		#############

		k+=1
		pheno = row['jd_code']
		snp_name = row['SNP']
	
		print 'Starting scan pair', k, '/', len(df_catalog_full), '--', pheno, snp_name
		try:
			df_pheno = pl.pheno_from_file(str(pheno))
		except:
			continue
			
		pheno = df_pheno.columns.tolist()[1]
	
		snp_bp_pos = float(row['Chr_pos'])
		chromosome = chr

		snp_list = []
		maf_common = float(non_decimal.sub('', row['Risk_Allele_Frequency']))
		OR_common = float(row['OR_or_beta'])

		print 'Searching for SA with', snp_name, '\tmaf:', maf_common

		odds_ratio1, p_value1 = OR_common, float(row['p_value'])
		print 'odds ratio: {}\t p-value: {}\n'.format( odds_ratio1, p_value1)


	
		snp_list_desc = []
		snps_selected = []

		if (method=='gene'):
			
			reported_genes = [gene_str.strip() for gene_str in row['Reported_Gene'].split(',')]
			mapped_genes = [gene_str.strip() for gene_str in row['Mapped_gene'].split(' - ')]
			
 			print 'Reported Genes, Mapped_genes', reported_genes, mapped_genes
 			
 			gene_list = list(set(reported_genes + mapped_genes))


			for gene in gene_list:	
				try:
					print 'Trying', gene, chr
					gene_definitions = df_genes.get_group((gene,chr))
					
					for ind1, row1 in gene_definitions.iterrows():
						min = row1['start'] - 40000
						max = row1['end'] + 40000
						if not(inrange( snp_bp_pos, min, max )):
							continue

						left = bisect.bisect_left(positions, min)
						right = bisect.bisect_right(positions, max)

						snps_selected += [position_to_snp[pos] for pos in positions[left:right]]
				except:
					'Gene', gene, 'not found.'
					

			snps_selected = list(set(snps_selected))

		else:
			for name, locus in snps.iteritems():

				bp_position = locus[0]

				if abs(bp_position - snp_bp_pos) < d_bp:
					snps_selected.append(name)

# 		print 'Selected:', snps_selected

		for name in snps_selected:				

				locus = snps[name]
				bp_position, maf, iter = locus[0], locus[1], locus[2]
			
				if (maf > 0.05) or (maf > 0.4*maf_common):
					continue

				print name, maf

				df_snp = pd.DataFrame(zip([int(x.iid) for x in pl.sample_list], [x for x in iter ]), columns=['FID', 'SNP'])

				a = pl.make_cont_table(df_snp, df_pheno)
				print a
				if a==False:
					odds_ratio, p_value = -1, -1
				else:
# 					odds_ratio, p_value = exact.pval_null(a)
					odds_ratio, p_value = exact.fisher_exact(a)	

				print 'Testing SNP: {}, maf: {},  d-bp: {}, cont. table: {}, OR {}, p-val {}'.format( name, maf, bp_position - snp_bp_pos, a, odds_ratio, p_value )

				# only include if they are in the right direction.
				if (odds_ratio1 <= 0) or (odds_ratio <= 0):
					continue
				if math.log(odds_ratio1) * math.log(odds_ratio) * math.log(maf_common / 0.5) > 0:
					continue

				if p_value < 0.05 and p_value > 0:
					print 'adding SNP'
					snp_list.append(df_snp)
					snp_list_desc.append([name, maf, odds_ratio, p_value, a])
		


		if len(snp_list) == 0:
			continue
	
		skat_pval = -1

		# run skat
		if len(snp_list) > 1:
			current = snp_list[0].rename(columns={'SNP':'SNP_1'})
			for i, frame in enumerate(snp_list[1:], 2):
				current = current.merge(frame, on='FID').rename(
										 columns={'SNP':'SNP_%d' % i})
			current = pd.merge(current, df_pheno.drop(df_pheno[df_pheno[df_pheno.columns[-1]]==-9].index.tolist()), on='FID')
# 			current = current.drop(current[current['pheno_446.4']==-9].index.tolist())
			
			current = pd.merge(df_demos[['FID', 'min_age', 'gender']], current, on="FID")
			current['gender'] = current['gender'].map({'M': 1,'F': 0})
			
			y = np.array([x-1 for x in current[current.columns[-1]].tolist()])
			X = np.array([tuple[2:4] for tuple in current.itertuples() ])
			Z = np.array([tuple[4:-1] for tuple in current.itertuples() ])
			try:
				skat_pval = sw.p_value_bin(X, y, Z) 
			except:		
				skat_pval = -1

	
		df_snp_comb = pl.combine_snps(snp_list)
		maf_ = [ len(df_snp_comb[df_snp_comb['SNP']==x].SNP.tolist()) for x in range(4)  ]
		maf_cmb = (2*maf_[2]+maf_[1])/(2.*float(maf_[0]+maf_[1]+maf_[2]))



		a2 = pl.make_cont_table(df_snp_comb, df_pheno)
		odds_ratio, p_value = exact.pval_null(a2)	
		
		if maf_common > 0.5:
			P_12 = ( maf_cmb * odds_ratio ) / ( (1-OR_common) / maf_common )
		else:
			P_12 = ( maf_cmb * odds_ratio ) / ( OR_common * maf_common )

		try:
			lr_nocovar = gm.logistic_regression(df_pheno, df_snp_common, drop_covariates=True)
			lr = gm.logistic_regression(df_pheno, df_snp_common)
		except:
			lr_nocovar = -1
			lr = -1

		print "\n------------------\nRare SNP List:"
		print 'SNP\tMAF\tOR\tp-val\tmatrix'
		for desc in snp_list_desc:
			print '\t'.join([str(x) for x in desc])
		print

		print 'Combined Matrix'
		print a2

		print '\nCommon SNPs odds ratio: {}\t p-value: {}\t maf_common: {}'.format( odds_ratio1, p_value1, maf_common)
		print 'Comb SNPs odds ratio: {}\t p-value: {}\t maf_cmb: {}'.format( odds_ratio, p_value, maf_cmb)

		print "\n------------------\n\n"

		data = [snp_name, pheno, chr, row['Disease_Trait'], row['Initial_Sample_Size'], row['p_value'], OR_common, maf_common, ','.join(gene_list),  odds_ratio, p_value, maf_cmb, P_12, skat_pval, lr, lr_nocovar,   ','.join([str(z[0]) for z in snp_list_desc]), ','.join([str(z[1]) for z in snp_list_desc]), ','.join([str(z[2]) for z in snp_list_desc]), ','.join([str(z[3]) for z in snp_list_desc]), a2 ]
		lfile.write( '\t'.join( [str(x) for x in data] ) + '\n')

	lfile.close()



def main():
# 	for chr in range(16, 23):
#   	df_catalog = pd.io.parsers.read_csv("PheWAS_translate_echo.csv", dtype=str)
 	df_catalog_full = pd.io.parsers.read_csv("../data/PheWAS_translate_full.csv", dtype=str)

# 	df_catalog_lab = pd.io.parsers.read_csv('PheWAS_translate_noncond_continuous.csv', dtype=str)
# 	df_catalog_lab = df_catalog_lab.drop_duplicates(['SNP', 'Disease_Trait', 'OR_or_beta'])

	# select * from catalog_UP_brief join Lab_view on lower(catalog_UP_brief.p_value_text) like CONCAT('%', lower(Lab_view.Lab_name), '%')  where continuous_trait = 1
# select * from catalog_UP_brief join Lab_view on lower(catalog_UP_brief.p_value_text) like CONCAT('%', lower(Lab_view.Lab_name), '%')  where continuous_trait = 1
# 
# select unique(SNP) from catalog_UP_brief join Lab_view on lower(catalog_UP_brief.p_value_text) like CONCAT('%', lower(Lab_view.Lab_name), '%')  where continuous_trait = 1
# 
# select * from LabWAS.Exome_36k_adult_Wnh_GA_view where SNP in (select SNP from catalog_UP_brief join Lab_view on lower(catalog_UP_brief.p_value_text) like CONCAT('%', lower(Lab_view.Lab_name), '%')  where continuous_trait = 1) limit 100

#  	df_catalog = df_catalog.drop_duplicates(['SNP', 'jd_code'])
  	df_catalog_full = df_catalog_full.drop_duplicates(['SNP', 'jd_code'])


	for chr in range(1,19):
		print 'Processing Chromosome', chr
#   		snps_cmn = getdata2(chr, df_catalog)
 		snps = getdata(chr)
# 		df_catalog_lab_chr = df_catalog_lab[df_catalog_lab['Chr']==str(chr)]
		
# 		print len(df_catalog_lab_chr)
# 		run_continuous(df_catalog_lab_chr, snps, snps_cmn, chr)

# 		df_catalog_ = df_catalog[df_catalog['Chr']==str(chr)]
#    		df_catalog_ = df_catalog[df_catalog['Chr']==str(chr)]
		df_catalog_full_ = df_catalog_full[df_catalog_full['Chr']==str(chr)]
		runloop_nocond(df_catalog_full_, snps, chr)
#   		runloop(df_catalog_, snps, snps_cmn, chr)


	
if __name__ == "__main__":
    main()
