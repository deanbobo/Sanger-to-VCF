#!/usr/bin/env python

import os
import argparse
from argparse import RawTextHelpFormatter
import re
import pandas as pd

parser = argparse.ArgumentParser(description='This script will take two fasta sequences and align them to the reference using BLAST. Then will use BLAST\'s BTOP string to make a VCF file. You must specify two fasta files using -fa1 and -fa2. The -sample parameter will put the sample ID in the VCF header. You can change the blast database using -db parameter.\n', formatter_class=RawTextHelpFormatter)
parser.add_argument('-sample', metavar='<sample_name>', help='sample name for header in VCF', required=True)
parser.add_argument('-fa1', metavar='<fa1>',help='First fasta file (1 seqeunce per file)', required=True)
parser.add_argument('-fa2', metavar='<fa1>',help='Second fasta file (1 seqeunce per file)', required=True)
parser.add_argument('-verbose', action='store_true', help='prints more information [verbose mode]')
parser.add_argument('-out',metavar='<out>', help='Output file. Note: ".vcf" will be appended automatically', required=True)
parser.add_argument('-db', metavar='<db>', required=True)

args = parser.parse_args()

sample=args.sample
fa1=args.fa1
fa2=args.fa2
out=args.out
db=args.db

print('using reference: '+str(db))
#blast the two sequences

cmd1='blastn -query '+fa1+' -db '+db+' -outfmt "6 qseqid sseqid sstart btop send" -parse_deflines | head -n 1 > fa1.aln.tmp'
print('running command: '+str(cmd1))
os.system(cmd1)
fa1_blast = open('fa1.aln.tmp').read().rstrip().split('\t')
os.remove('fa1.aln.tmp')

cmd2='blastn -query '+fa2+' -db '+db+' -outfmt "6 qseqid sseqid sstart btop send" -parse_deflines | head -n 1 > fa2.aln.tmp'
print('running command: '+str(cmd2))
os.system(cmd2)
fa2_blast = open('fa2.aln.tmp').read().rstrip().split('\t')
os.remove('fa2.aln.tmp')

# print fa1_blast
# print fa2_blast

fa1_chr = fa1_blast[1]
fa1_start = int(fa1_blast[2])
fa1_btop_string = fa1_blast[3]
fa1_stop = int(fa1_blast[4])

fa2_chr = fa2_blast[1]
fa2_start = int(fa2_blast[2])
fa2_btop_string = fa2_blast[3]
fa2_stop = int(fa2_blast[4])

if fa1_chr != fa2_chr:
	print('Cannot continue. The sequences did not align to the same chromosome.\n')
	exit(0)

aln_start=min(fa1_start, fa2_start)
aln_stop=max(fa1_stop,fa2_stop)
aln_len= aln_stop - aln_start

fa1_btop = filter(None, re.split('(\d+)', fa1_btop_string))
fa2_btop = filter(None, re.split('(\d+)', fa2_btop_string))


#create pandas array
df = pd.DataFrame(columns=['cov','ref','alt1','alt2'])#range(aln_start,aln_stop+1))

for i in range(aln_start,aln_stop+1):
	df.loc[i]=''

df.loc[:,'cov']=0 #df.loc[rows,'columns'] = value

#BLAST BTOP string goes QUERY|REFERENCE

fa1_pos=fa1_start
for element in fa1_btop:

	if element.isdigit():  #for exact matches
		#print element
		element=int(element)
		for i in range(fa1_pos,fa1_pos+element): #they are matches
			df['cov'][i]+=1
		fa1_pos += element #increment fa1_pos

	else: #for mismatches and gaps
		#var_list = [element[i:i+2] for i in range(0, len(element), 2)]

		for j in [element[j:j+2] for j in range(0, len(element), 2)]:
			#print i
			
			if re.match('[ACTGN]-',j):#it is an insertion  #insertions look like [ACTG]-
				pass #go to next iteration
				#do not increment position
			
			elif re.match('-[ACTGN]',j): #it is a deletion #deletions look like    -[ACTG]
				fa1_pos+=1 #increment the position
				#print "it is a deletion, increment"
			
			elif re.match('N',j): #if it is an N
				fa1_pos+=1 #increment position
				#print "there is an N, increment position"
			
			else: #it is a mismatch 
				#print "it is a mismatch, add to coverage, increment position"
				df['cov'][fa1_pos]=1 #add 1 to coverage
				df['ref'][fa1_pos]=j[1:2] #add ref allele
				df['alt1'][fa1_pos]=j[0:1] #add alt allele

				fa1_pos+=1
fa2_pos=fa2_start
for element in fa2_btop:

	if element.isdigit():  #for exact matches
		element=int(element)
		for i in range(fa2_pos,fa2_pos+element): #they are matches
			df['cov'][i]+=1
		fa2_pos += element #increment fa2_pos


	else: #for mismatches and gaps

		for j in [element[j:j+2] for j in range(0, len(element), 2)]:
			
			if re.match('[ACTGN]-',j):#it is an insertion  #insertions look like [ACTG]-
				pass #go to next iteration
				#do not increment position
			
			elif re.match('-[ACTGN]',j): #it is a deletion #deletions look like    -[ACTG]
				fa2_pos+=1 #increment the position

			
			elif re.match('N',j): #if it is an N
				fa2_pos+=1 #increment position
			
			else: #it it is a mismatch, add to coverage, increment position
				if df['ref'][fa2_pos] and df['ref'][fa2_pos]!=j[1:2]:
					print('Error: the reference allele identified at position ',fa2_pos,'  differ between the two alignments. Cannot continue.')
					exit(1)
				else:
					df['cov'][fa2_pos]+=1 #add 1 to coverage
					df['ref'][fa2_pos]=j[1:2] #add ref allele
					df['alt2'][fa2_pos]=j[0:1] #add alt allele

				fa2_pos+=1



#make the VCF
vcfout = open(out+'.vcf','w')
vcfout.write('##fileformat=VCFv4.2\n')
vcfout.write('##reference='+db+'\n')
vcfout.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">\n')
vcfout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
vcfout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')

bedout = open(out+'.bed','w')
bedout.write('#chr\tstart\tstop\n')


bedpos=''

for index, row in df.iterrows():  #iterate through dataframe (df)
	if row['cov']==2:  #if it is covered by both sequences
		
		#ADD TO BED  #THE LAST RANGE IS WRITTEN AFTER FOR LOOP
		if bedpos=='': #starting
			bedpos=index
			bedout.write(str(fa1_chr)+'\t'+str(index-1)+'\t')
		elif index == bedpos+1: #we're contiguous 
			bedpos+=1
		else: #we're not contiguous
			bedout.write(str(bedpos)+'\n') #close the previous range
			bedout.write(str(fa1_chr)+'\t'+str(index-1)+'\t')  #start new range
			bedpos=index
		
			
		
		#MAKE VCF
		if (row['alt1'] or row['alt2']):
			GT='0/0'
			AC=0
			altstr=''
			if row['alt1'] == row['alt2']:
				altstr=row['alt1']
				GT='1/1'
				AC=2
			elif row['alt1'] and row['alt2']:
				altstr=row['alt1']+','+row['alt2']
				GT='1/2'
				AC=2
			else:
				altstr=''.join(filter(None,[row['alt1'],row['alt2']]))
				GT='0/1'
				AC=1			
			outlist = [fa1_chr,index,'.',row['ref'],altstr,'50','PASS','AC='+str(AC),'GT',GT]
			vcfout.write('\t'.join([str(i) for i in outlist])+'\n')

#last range of bed
bedout.write(str(bedpos)+'\n')

covered_regions=df[df['cov']==2].count()['cov']
uncovered_regions=df[df['cov']<2].count()['cov']
print('Result:\tChr:'+str(fa1_chr)+'\tThe two sequnces covered '+str(covered_regions)+' bases. '+str(uncovered_regions)+' bases were not covered by both Sanger sequences.\n')

bedout.close()
vcfout.close()
