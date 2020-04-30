
######################################################################
# A. Pendleton
# 2019-11-13
# Formatting input files for the Jupyter notebook that will allow interaction with 
#    expression and coexpression data 
######################################################################


from optparse import  OptionParser
import sys
import os
import numpy as np

##############################################################################

USAGE = """

python format_website_inputs.py	--outdir <Output directory> --samples <Sample, condition, experiment file> 


		
#Input descriptions
Required inputs: 

	outdir ==  Output directory - where all files will be written
	
	samples == Sample file - this file will need to be generated/populated manually *BEFORE* running this script. 
			Each row details the sample (e.g. Leaf_01), the condition (e.g. High Light), and experimental group (e.g. Light). 
			Must be tab delimited. 


Optional inputs:


	#ANNOTATION DATA
	
	annotations == Functional annotation file. Each gene has a row. 
	
	annotationColumns == Which columns in the annotations file have the annotations you wish to link to your gene of interest? 
				(example = "2,5,6" will extract the 2nd, 5th, and 6th columns 
				-OR- 
				"all" will extract all columns in file) 
	
	
	#HOMOLOG DATA
	
	
	homologs == File with homologs linked to genes in your species of interest 
			(e.g. column numbers with Arabidopsis, maize, yeast, or human homologs) 
	
	homologColumns == Which columns in the homologs file have the homologs you wish to link to your gene of interest.
				(example = "3,4" will extract the 3rd and 4th columns 
				-OR- 
				"all" will extract all columns in file) 
	
	
	#DIFFERENTIAL EXPRESSION DATA
	Note: All files require the same format/structure. *MUST* have header lines
	See description for de_fileList below:
	
	
	de_file_list == Differential expression file - tab-delimited file where:
							column 1 = condition
							column 2 = experiment
							column 3 = path to relevant differential expression output file 
							         with FDR, pvalues, log fold change (e.g. PATH/drought_de_output.txt)
	de_gene_column == [integer] column in the input differential expression files where
			the gene ID can be found (e.g. 4 will equal 4th column)
			
	de_pvalue_column == [integer] column in the input differential expression files where
			the test *p-value* can be found (e.g. 4 will equal 4th column)

	de_FDR_column == [integer] column in the input differential expression files where
			the test *FDR* can be found (e.g. 4 will equal 4th column)
			
	de_logFC_column == [integer] column in the input differential expression files where
			the test *log fold change (FC)* can be found (e.g. 4 will equal 4th column)

	#GLOBAL EXPRESSION DATA 
	Note: *MUST* have header line with sample information. We recommend the VST normalized
	   expression matrix that is produced by the coexpression network pipeline from the 
	   Wisecaver lab. But any gene expression matrix can work as input as long as the 
	   file structure described below is maintained:	
	
	exp_matrix = Multi column matrix file with gene expression values. *TAB-DELIMITED*
				row 1 = header line with sample names.
				column 1 = gene ID
				column 2-N (where N equals number of samples + 1) = expression value that 
					corresponds to the correct sample identifier
	
"""

parser = OptionParser(USAGE)
parser.add_option('--outdir',dest='outdir', help = 'Output directory - where all files will be written')
parser.add_option('--samples',dest='samples', help = 'Tab delmited sample file with a minimum of three columns that include: Sample ID, Condition, Experimental group. ')
parser.add_option('--annotations',dest='annotations', help = 'Annotation file not provided', default= '')
parser.add_option('--annotationColumns',dest='annotationColumns', help = ' Which columns in the annotations file have the annotation(s) you wish to link to your gene of interest?', default= '')
parser.add_option('--homologs',dest='homologs', help = 'Homolog file not provided', default= '')
parser.add_option('--homologColumns',dest='homologColumns', help = ' Which columns in the homolog file have the homolog(s) you wish to link to your gene of interest?', default= '')
parser.add_option('--de_file_list',dest='de_file_list', help = 'Differential expression file list', default= '')
parser.add_option('--de_pvalue_column',dest='de_pvalue_column', help = 'Which column in the DE file has the p-value?', default= '')
parser.add_option('--de_gene_column',dest='de_gene_column', help = 'Which column in the DE file has the gene ID?', default= '')
parser.add_option('--de_FDR_column',dest='de_FDR_column', help = 'Which column in the DE file has the FDR value?', default= '')
parser.add_option('--de_logFC_column',dest='de_logFC_column', help = 'Which column in the DE file has the log fold change value?', default= '')
parser.add_option('--exp_matrix',dest='exp_matrix', help = 'File that has normalized expression values for all samples being assessed', default= '')

(options, args) = parser.parse_args()

if options.outdir is None:
    parser.error('Output directory path not given.')
if options.samples is None:
    parser.error('Sample file not given.')

###############################################################################


#Define output directory variable
outdir = options.outdir 
if '/' != outdir[-1]:
	outdir = outdir + '/'
print(outdir)


################
### Process sample file
################

print('\n##SAMPLES ')
sampleFile = options.samples
print('\n' + 'Processing the sample file: \n %s ' % sampleFile +  '\n')


samples, conditions, experiments = [], [], []
sampleDict = {}

lineCount = 0
for line in open(sampleFile, 'r'):
	lineCount += 1 #track line number, so we can ignore header
	
	#skip header
	if '#' in line[0] or lineCount == 1:
		continue 
	
	line=line.rstrip().split('\t') #split line by tab
	sample, condition, experiment = line[0], line[1], line[2]
	
	#append each data point to the lists to track number of unique sample IDs, 
	samples.append(samples)
	conditions.append(condition)
	experiments.append(experiment)
	 
	#Store in dictionary
	if sample not in sampleDict.keys():
		sampleDict[sample] = {}
		sampleDict[sample]['Condition'] = condition
		sampleDict[sample]['Experiment'] = experiment


#Check that there are no redundant sample identifiers
if len(sampleDict.keys()) != len(samples): #len(np.unique(sampleDict.keys())):
	print('ERROR: Are there redundant sample identifiers in your samples?')
	print('%i == number of samples' % len(sampleDict.keys()))
	print('%i == numnber of unique sample identifiers' % len(np.unique(samples))) #len(np.unique(sampleDict.keys())))
	print('exiting...')
	sys.exit()
else:
	print('%i samples in input' % len(sampleDict.keys()))
print('%i conditions in input' % len(np.unique(conditions)))
print('%i experiments in input' % len(np.unique(experiments)))


print(sampleDict.keys())


####################################################################
### Process annotation file
####################################################################
"""
print('\n##ANNOTATIONS ')

#Inputs:
annotFile = options.annotations
annotFileColumns = options.annotationColumns

#Outputs
out_annotfile = outdir + 'Genes_annotations.txt'
out_annotFile = open(out_annotfile, 'w')
print('Will write the extracted annotation information to:\n %s' % out_annotfile, '\n')

if annotFile != '':
	print('\n' + 'Processing the annotation file: \n %s ' % annotFile +  '\n')
	print('Will extract columns: %s ' % annotFileColumns)


	######## Determine columns to extract
	#If no columns were provided, report error and exit
	if annotFileColumns == '':
		print('ERROR: You must provide columns to extract from annotation file!')
		print('Exiting...')
		sys.exit()
	#In case the user put quotation marks, get rid of them
	annotFileColumns = str(annotFileColumns.replace('"', '').replace("'", ""))

	columnsToExtract = [] # set equal to zero
	
	#If user wants all columns, then define columnsExtract as equal to all columns
	#    (excluding first column == geneID) to end of headerline
	if annotFileColumns == 'all':
		for line in open(annotFile, 'r'):
			for i in range(1, len(line.split('\t'))):
				columnsToExtract.append(i)
			break
		print(columnsToExtract)
		
	else:
		annotFileColumns = annotFileColumns.replace('"','') #replace space if given by user 
		annotFileColumns = annotFileColumns.replace(' ','') #replace space if given by user 
		#If only one column given
		if ',' not in annotFileColumns:
			columnsToExtract.append(int(annotFileColumns))
		#If more than one column given, then split by comma and append each entry
		else:
			for i in annotFileColumns.split(','):
				columnsToExtract.append(int(i))
			print(columnsToExtract)	
		
	#### Now that we've determined the number of columns to extract, now lets extract 
	####    the data from the annotation file
	lineCount = 0 
	for line in open(annotFile, 'r'):
		lineCount += 1
		
		#skip header line
		if '#' == line[0] or lineCount ==1:
			line=line.rstrip().split('\t')
			#Write out "GeneID" as first column's header
			out_annotFile.write('#GeneID\t')
			for c in columnsToExtract[1:]:
				out_annotFile.write(line[c-1] + '\t') #You want column -1, as this "index-speak" in python
			out_annotFile.write('\n')
			continue
		
		line=line.rstrip().split('\t') #Strip and split line
		
		for c in columnsToExtract:
			#if for some reason this data point is missing in the in file, simply 
			#   write out a tab
			if len(line) < c:
				out_annotFile.write('\t')
			#if the data is there, then add it 
			else:
				out_annotFile.write(line[c-1] + '\t') #You want column -1, as this "index-speak" in python
		
		#Now finish with writing out a new line
		out_annotFile.write('\n')
	out_annotFile.close()
	
else:
	print('\nNo annotations file to parse')




####################################################################
### Process homolog file
####################################################################

print('\n##HOMOLOGS ')

#Inputs:
homologFile = options.homologs
homologFileColumns = options.homologColumns

#Outputs
out_homologfile = outdir + 'Genes_homologs.txt'
out_homologFile = open(out_homologfile, 'w')
print('Will write the extracted annotation information to:\n %s' % out_homologfile, '\n')

if homologFile != '':
	print('\n' + 'Processing the annotation file: \n %s ' % homologFile +  '\n')
	print('Will extract columns: %s ' % homologFileColumns)


	######## Determine columns to extract
	#If no columns were provided, report error and exit
	if homologFileColumns == '':
		print('ERROR: You must provide columns to extract from annotation file!')
		print('Exiting...')
		sys.exit()
	#In case the user put quotation marks, get rid of them
	homologFileColumns = str(homologFileColumns.replace('"', '').replace("'", ""))

	columnsToExtract = [] # set equal to zero
	
	#If user wants all columns, then define columnsExtract as equal to all columns
	#    (excluding first column == geneID) to end of headerline
	if homologFileColumns == 'all':
		for line in open(homologFile, 'r'):
			for i in range(1, len(line.split('\t'))):
				columnsToExtract.append(i)
			break
		print(columnsToExtract)
		
	else:
		homologFileColumns = homologFileColumns.replace('"','') #replace space if given by user 
		homologFileColumns = homologFileColumns.replace(' ','') #replace space if given by user 
		#If only one column given
		if ',' not in homologFileColumns:
			columnsToExtract.append(int(homologFileColumns))
		#If more than one column given, then split by comma and append each entry
		else:
			for i in homologFileColumns.split(','):
				columnsToExtract.append(int(i))
			print(columnsToExtract)	
		
	#### Now that we've determined the number of columns to extract, now lets extract 
	####    the data from the annotation file
	lineCount = 0 
	for line in open(homologFile, 'r'):
		lineCount += 1
		
		#skip header line
		if '#' == line[0] or lineCount ==1:
			line=line.rstrip().split('\t')
			#Write out "GeneID" as first column's header
			out_homologFile.write('#GeneID\t')
			for c in columnsToExtract[1:]:
				out_homologFile.write(line[c-1] + '\t') #You want column -1, as this "index-speak" in python
			out_homologFile.write('\n')
			continue
		
		line=line.rstrip().split('\t') #Strip and split line
		
		for c in columnsToExtract:
			#if for some reason this data point is missing in the in file, simply 
			#   write out a tab
			if len(line) < c:
				out_homologFile.write('\t')
			#if the data is there, then add it 
			else:
				out_homologFile.write(line[c-1] + '\t') #You want column -1, as this "index-speak" in python
		
		#Now finish with writing out a new line
		out_homologFile.write('\n')
	out_homologFile.close()
	
else:
	print('\nNo homologs file to parse')
"""



####################################################################
#### Differential Expression Files (e.g. EdgeR, DeSeq2, etc. outputs)
####################################################################

print('\n\n' + '#### DIFFERENTIAL EXPRESSION\n')

#File with condition being tested for differential expression, the experiment this 
#    test would fall under (e.g. Drought), and the path to the consistently 
# 	formatted DE output files (e.g. outputs from DeSeq or EdgeR) 

de_fileList = options.de_file_list

#Track the conditions tested for differential expression
de_conditions = []

#Columns with stats that we need
geneColumn = options.de_gene_column
pvalueColumn = options.de_pvalue_column
FDRColumn = options.de_FDR_column
logFCColumn = options.de_logFC_column


#Store all de values in a dictionary where the key is the gene
deDict = {}

#If not an empty variable, then start parsing
if de_fileList != '':

	#Let's first check to make sure that the user has provided a value for 
	#  pvalue, fdr, and log fc columns to extract
	if geneColumn == '':
		print('ERROR: Differential expression file provided but no indication of which column has the gene ID')
		print('Exiting...')
		sys.exit()
	if pvalueColumn == '':
		print('ERROR: Differential expression file provided but no indication of which column has the pvalue')
		print('Exiting...')
		sys.exit()
	if FDRColumn == '':
		print('ERROR: Differential expression file provided but no indication of which column has the FDR')
		print('Exiting...')
		sys.exit()
	if logFCColumn == '':
		print('ERROR: Differential expression file provided but no indication of which column has the log fold change value')
		print('Exiting...')
		sys.exit()
		
	#Now parse the input file to determine the paths to those we need for writing out
	lineCount = 0
	for line in open(de_fileList, 'r'):
		lineCount += 1
		
		if '#' in line[0] or lineCount == 1: #skip header
			continue
			
		line = line.rstrip().split('\t') #strip and split tab by line
		
		#Determine condition
		condition = line[0]
		de_conditions.append(condition) #for tracking 
		
		if condition not in conditions: #make sure the condition is in our condition list
			print('ERROR: the condition below is not in our global conditions list:')
			print(condition + '\n', line)
			print("Exiting...")
			sys.exit()
		#Determine experiment
		experiment = line[1]
		if experiment not in experiments: #make sure the condition is in our condition list
			print('ERROR: the experiment below is not in our global experiment list:')
			print(experiment + '\n', line)
			print("Exiting...")
			sys.exit()	
			
		#Determine DE file with results
		defile = line[2]
		
		#Read through the de file to extract out the data we want	
		de_lineCount = 0	

		for line in open(defile, 'r'):
			de_lineCount += 1
			#Check to see if we have headers in the file. Can be checked by seeing if the
			#   column specified by user to be the pvalue is a float or not 
			if de_lineCount == 1:
				continue
				
			line = line.rstrip().split('\t') # split the file by tab
			geneID = line[int(geneColumn) - 1] #get gene ID
			pvalue = line[int(pvalueColumn) - 1] #get pvalue
			FDR = line[int(FDRColumn) - 1] #get FDR
			logFC = line[int(logFCColumn) - 1] #get log fold change
			
			#add geneID to dictionary for storing all data if not already in there
			if geneID not in deDict.keys():
				deDict[geneID] = {}
			
			#now create another key linked to that gene for the condition linked to this file
			# as determined by the de_file_list FileCache
			deDict[geneID][condition] = ['', pvalue, FDR, logFC] #keep index =0 empty for now	
	
#else, skip as there's no file to parse
else:
	print('\nNo differential expression description/path file to parse')


if len(deDict.keys()) > 0:
	#Report how many genes' differential expression data was stored:	
	print('%i genes stored in differential expression dictionary deDict' % len(deDict.keys()))

	print('The following differential expression conditions were processed:')
	print('\t'.join(map(str, de_conditions)))
elif de_fileList != '':
	print('Why were no genes stored as keys in the deDict?')
	print('Exiting...')
	sys.exit()
	

################################################
#### Expression matrix from coexpression network pipeline
################################################

print('\n\n\n', '####### MEAN EXPRESSION CALCULATIONS', '\n')

#Get file and column extraction data from user provided information
expressionFile =  options.exp_matrix #OLD -- '/depot/jwisecav/data/pendlea/setaria/APL_Data/ME034_Assembly/v1.0/9_CoexpressionNetworks/52Conditions_ME034v1.0/gene_counts_normalized_vst_transformed.matrix'#'/depot/jwisecav/data/pendlea/coexpression_assessments/development/52_conditions_newpipeline/Setaria_A10_normalized_vst_transformed.matrix'
print('Parsing gene expression data from expression matrix file:\n '  , expressionFile)

#To track total expression data, values and conditions:
sample2Column = {}
expDict = {}
condition_expDict = {}
sample_expDict = {}
exp_conditions = [] #Keep track of which conditions we even quantified expression on 


if expressionFile != '':
	lineCount = 0 #For tracking line numbers (distinguishing headers from data lines)
	
	#Parse the user-provided file that links each sample to its expression quantification fiel
	for line in open(expressionFile, 'r'):
		lineCount += 1
		line=line.rstrip().split('\t') #strip and split tab by line

		##HEADER LINE
		#parse header to get sample identifiers
		if lineCount == 1:
			for i in range(0,len(line)):
				sample = line[i]
				#Check that the sample is in the dictionary
				#If not, it may have been ignored for QA/QC reasons, but still alert the user
				if sample not in sampleDict.keys():
					print('ERROR: Sample not in sample list provided in samples.txt:\n %s \nSkipping...' % sample)
					continue
				sample2Column[i+1] = sample
			continue

		#DATA LINES		
		#Now parse the data lines which are lines 2 and on		
		geneID = line[0]
		
		#if gene ID not already in the expDict dictionary, then add
		if geneID not in sample_expDict.keys():
			#Create empty keys for each gene ID 
			sample_expDict[geneID] = {}
			condition_expDict[geneID] = {}
			
		#Create empty values for each sample we want to store the exp. data for
		for s in sampleDict.keys():
			sample_expDict[geneID][s] = ''
					
			#Now create empty arrays linked to each possible condition 
			for c in conditions:
				condition_expDict[geneID][c] = []

		#Get the expression values
		for i in range(1,len(line)):
			#Get the expression Value
			expValue = float(line[i])
		
			#Determine which sample it belongs to based on which column you're reading in
			sample = sample2Column[i]
			
			#Keep track of the expression across each condition as we will want to 
			#   calculate the average expression value for plotting in the jupyter notebook
			#Determine the condition that this sample is linked to:
			Condition = sampleDict[sample]['Condition']
			exp_conditions.append(Condition)
			condition_expDict[geneID][Condition].append(expValue)
		
			#Now also track per sample
			sample_expDict[geneID][sample] = expValue
		#print(sample_expDict[geneID]) ##REMOVE

#else, skip as there's no file to parse
else:
	print('\nNo quantification of expression description/path file to parse')


#Now report how many genes expression data was stored:
if len(expDict.keys()) > 0: #If there's even a populated dictionary, report gene count
	print('%i genes stored in expression dictionary expDict' % len(expDict.keys()))


#If the dictionaries of either of 
if len(condition_expDict.keys()) > 0 and len(deDict.keys()) > 0:
	for geneID in condition_expDict.keys():
		for condition in de_conditions: #expDict[geneID].keys():
			#If the gene isn't in the differential expression dictionary (deDict) at all,
			#   then skip it
			if geneID not in deDict.keys():
				print('gene ID missing from deDict:' % geneID)
				continue
			
			else: #meaning the gene IS in the dictionary:
				if len(condition_expDict[geneID][condition]) == 0 or condition_expDict[geneID][condition] == '':
					averageExpression = 0.0
				if len(condition_expDict[geneID][condition]) == 1:
					averageExpression = condition_expDict[geneID][condition][0]
				if len(condition_expDict[geneID][condition]) > 1:
					averageExpression = np.mean(condition_expDict[geneID][condition])
			#Hashed out below for troubleshooting
			#print(condition, 'averageExpression = ', averageExpression)
			#print(condition, 'before', deDict[geneID][condition])
			deDict[geneID][condition][0] = averageExpression
			#print(condition, 'after', deDict[geneID][condition])

"""################################################
####WRITE OUT THE OVERALL EXPRESSION FILE
################################################

#Output file name to be written to
total_expfile = outdir + 'Genes_expression.txt'
total_expFile = open(total_expfile, 'w')
print('Writing total expression data to:\n%s' % total_expfile + '\n\n')

#Write header line first
total_expFile.write('#GeneID')
for sample in sampleDict.keys():
	total_expFile.write('\t%s' % sample) #tab delimited, write out conditions to headerline
total_expFile.write('\n') #Write new line once header line is done being written out


#Now write out the expression data
for geneID in sample_expDict.keys():
	outLine = [] #Line to ultimately write out
	outLine.append(geneID)  #append the gene ID
	for sample in sampleDict.keys(): #For each condition we'll append the gene expression value in that condition
		outLine.append(sample_expDict[geneID][sample])
	#Now write the outline to the outfile
	total_expFile.write('\t'.join(map(str,outLine)) + '\n')
#print('\t\tDone writing the expression file...')
total_expFile.close()
"""


################################################
####WRITE OUT THE GLOBAL DIFFERENTIAL EXPRESSION FILE
################################################

#Output file name to be written to
total_DEfile = outdir + 'DE_experiments.txt'
total_DEFile = open(total_DEfile, 'w')
print('Writing total DE data to:\n%s' % total_DEfile, '\n\n')

print('Just below new lines... deDict')
print(geneID,deDict[geneID])

#Write headerline
total_DEFile.write('#GeneID')
for condition in de_conditions:
	total_DEFile.write('\t%s' % condition) #tab delimited, write out conditions to headerline
total_DEFile.write('\n') #Write new line once header line is done being written out

print('de_conditions = ', de_conditions)

#Now add the gene expression data
for geneID in sample_expDict.keys():
	lineOut = []
	#To first column write the gene ID
	lineOut.append(geneID)
	#print('here geneid = ', geneID)##REMOVE
	#for each condition, get the expression stats

	for condition in de_conditions:
		#averageExpression, pvalue, FDR, logFC = float(deDict[geneID][condition][0]), float(deDict[geneID][condition][1]), float(deDict[geneID][condition][2]), float(deDict[geneID][condition][3])
		averageExpression = float(deDict[geneID][condition][0])
		pvalue, FDR, logFC = float(deDict[geneID][condition][1]), float(deDict[geneID][condition][2]), float(deDict[geneID][condition][3])

		lineOut.append('%f,%f,%f,%f' % (averageExpression, pvalue, FDR, logFC))
	total_DEFile.write('\t'.join(map(str,lineOut)) + '\n')
#print('\t\tDone writing to differential expression outfile...')
total_DEFile.close()	





################################################
#### WRITING NETWORK MODULE FILES
################################################

# .... There should be no need to do this as they come out of the coexp pipe
#    in the proper format

