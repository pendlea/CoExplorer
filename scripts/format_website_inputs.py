
######################################################################
# A. Pendleton
# 2019-11-13
# Formatting input files for the Jupyter notebook that will allow interaction with 
#    expression and coexpression data 
######################################################################


from optparse import  OptionParser
import sys
import os
import pandas as pd
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
condition_list = []

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

	#Track in condition_list
	condition_list.append([condition, sample])
	
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
print('conditions: ')
print(np.unique(conditions))
print('\n\n' + '%i experiments in input' % len(np.unique(experiments)))
print('experiments: ')
print(np.unique(experiments))


#Make condition_list a dataframe
condition_df = pd.DataFrame(condition_list, columns=['Condition', 'Sample'])
print(condition_df)


####################################################################
### Process annotation file
####################################################################

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


#If homolog file is provided, let's parse it
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

	#Use provided columns to produce string of which columns to extract
	extractString = []
	for c in homologFileColumns.split(','):
		#Subtract one from the column number to become the index
		extractString.append(int(c)-1)
	print('Extracting columns from homolog file: ',  extractString)


	#Read in the homolog File with pandas
	homolog_data_table = pd.read_csv(homologFile, delimiter='\t')


	#Rewrite to the subsetted table with columns to extract only
	homolog_data_table = homolog_data_table.iloc[:, extractString]

	#Ensure each comma has a space after it for viewing in website tables
	homolog_data_table = homolog_data_table.replace({',': ', '}, regex=True)

	#Rename first column
	homolog_data_table = homolog_data_table.rename(columns={ homolog_data_table.columns[0]: "#GeneID" })

	#Write to outfile
	homolog_data_table.to_csv(out_homologfile, sep='\t', index=False)


else:
	print('\nNo homologs file to parse')




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

#Create empty dataframe to track global DE results
global_de_df = pd.DataFrame()

#Track how many samples/condition DE that has been processed
sampleCount = 0 

#Clear
genes = []

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
	
		if '#' in line[0] or lineCount == 1: #skip header or marked out lines
			continue
		
		line = line.rstrip().split('\t') #strip and split tab by line
	
		#Add one to sample count
		sampleCount += 1
	
		#Determine condition
		condition = line[0]
		de_conditions.append(condition) #for tracking 
	
		if condition not in conditions: #make sure the condition is in our condition list
			print('ERROR: the condition below is not in our global conditions list:')
			print(condition + '\n', line)
			print('Please check that conditions in your DE file list match that of your sample/condition file ')
			print("Exiting...")
			sys.exit()
		#Determine experiment
		experiment = line[1]
		if experiment not in experiments: #make sure the condition is in our condition list
			print('ERROR: the experiment below is not in our global experiment list:')
			print(experiment + '\n', line)
			print('Please check that experiments in your DE file list match that of your sample/condition file ')
			print("Exiting...")
			sys.exit()	
		
		#Determine DE file with results
		defile = line[2]
	
		#Reset for tracking with each new file
		genes, de_list = [], [] #clear both
	
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
			pvalue = float(line[int(pvalueColumn) - 1]) #get pvalue
			FDR = float(line[int(FDRColumn) - 1]) #get FDR
			logFC = float(line[int(logFCColumn) - 1]) #get log fold change
		
			#Now add information to the gene
			genes.append(geneID)
			de_list.append('%f,%f,%f' % (pvalue, FDR, logFC))
		
		#Store as a dataframe
		de_df = pd.DataFrame(de_list, columns=[condition], index=[genes])
				
		#Merge with global df after each new condition is processed
		global_de_df = pd.concat([global_de_df, de_df], axis=1)

	#Name the index of the final dataframe to be the genes
	global_de_df.index = genes

	#Sort the dataframe by gene so we can merge faster with the expression data down the road
	global_de_df = global_de_df.sort_index()

	#Get sneak peak of the data frame
	global_de_df.head()

#else, skip as there's no file to parse
else:
	print('\nNo differential expression description/path file to parse')


if len(genes) > 0:
	#Report how many genes' differential expression data was stored:	
	print('%i genes stored in differential expression dataframe' % len(genes))

	print('The following differential expression conditions were processed:')
	print('\t'.join(map(str, de_conditions)))


################################################
#### Total expression profiles (e.g. Kallisto, Cufflinks, Salmon, etc. per-gene expression
####     quantification files. Must be tab delimited)
################################################

print('\n\n\n', '####### MEAN EXPRESSION CALCULATIONS', '\n')

#Get file and column extraction data from user provided information
#expressionFile = '/depot/jwisecav/data/pendlea/coexpression_assessments/development/52_conditions_newpipeline/Setaria_A10_normalized_vst_transformed.matrix'
expressionFile = options.exp_matrix


#If the user did provide an expression file, then lets process it!
if expressionFile != '':

	#Read in gene expression matrix using pandas
	exp_data_table = pd.read_csv(expressionFile, delimiter='\t')


	#Create an empty dataframe to store the global (i.e. all conditions) gene average expression levels
	global_mean_df = pd.DataFrame()

	for condition in np.unique(de_conditions):     
		#For each condition in the differential expression analysis, 
		#.  let's pull out the samples that represent that condition from the expression matrix
		samples = condition_df['Sample'].loc[condition_df['Condition'] == condition]
		
		#Get average for the samples within that condition
		mean = exp_data_table[exp_data_table.columns.intersection(samples)].mean(axis=1)
		#Store as dataframe, column name = condition
		mean_df = pd.DataFrame(mean, columns=[condition])
   
		#For each condition processed, add to the global_mean dataframe
		global_mean_df = pd.concat([global_mean_df, mean_df], axis=1)
	
	#Check out final matrix
	global_mean_df.head()

	#Sort the dataframe by gene name (i.e. index) for merging with de dataframe
	global_mean_df = global_mean_df.sort_index()

	#Now report how many genes expression data was stored:
	print('%i genes stored in expression dataframe' % len(global_mean_df.index))

	
#else, skip as there's no file to parse
else:
	
	print('\nNo quantification of expression description/path file to parse')





################################################
##### MERGE THE DE AND EXPRESSION DATAFRAMES
################################################

#If the user did provide an expression file, then lets process it!
if expressionFile != '' and de_fileList != '':
	global_mean_df_index = global_mean_df.index

	print('%i genes in global gene expression dataframe' % len(global_mean_df.index))

	global_de_df_index = global_de_df.index
	print('%i genes in global differential expression dataframe' % len(global_de_df.index))

	#Find intersection between two so 
	len(global_mean_df_index.intersection(global_de_df_index))

	#Only get rows from diff exp dataframe that are in expression dataframe
	subsetted_global_de_df = global_de_df.loc[global_mean_df_index]

	### Merge the de and expression average matrices together so that the cells are merged results, comma delmited
	merged_final_df = global_mean_df.astype(str).add(',').add(subsetted_global_de_df.astype(str))

	#Move the geneID from index to its own column
	merged_final_df = merged_final_df.reset_index()

	#Rename from 'index' to '#geneID'
	merged_final_df.rename(columns={'index':'#GeneID'}, inplace=True)

	#Print length of the final dataframe
	print('%i genes in merged final dataframe' % len(merged_final_df.index))




################################################
####WRITE OUT THE OVERALL EXPRESSION FILE
################################################

print('##### WRITE OUT GENE EXPRESSION FILE')

total_expfile = outdir + 'Genes_expression.txt'

if expressionFile != '':
	print('Writing total DE data to:\n%s' % total_expfile, '\n\n')

	#Move the geneID from index to its own column
	exp_data_table = exp_data_table.reset_index()

	#Rename from 'index' to '#geneID'
	exp_data_table.rename(columns={'index':'#GeneID'}, inplace=True)

	#Write out using pandas df.to_csv function
	exp_data_table.to_csv(total_expfile, sep='\t', index=False)

else:
	print('No expression file provided, only empty Genes_expression.txt being created')

################################################
####WRITE OUT THE GLOBAL DIFFERENTIAL EXPRESSION FILE
################################################

if expressionFile != '' and de_fileList != '':

	print('##### WRITE OUT DIFFERENTIAL EXPRESSION FILE')

	#Output file name to be written to
	total_DEfile = outdir + 'DE_experiments.txt'
	print('Writing total DE data to:\n%s' % total_DEfile, '\n\n')

	#Write out using pandas df.to_csv function
	merged_final_df.to_csv(total_DEfile, sep='\t', index=False)





################################################
#### WRITING NETWORK MODULE FILES
################################################

# .... There should be no need to do this as they come out of the coexp pipe
#    in the proper format

