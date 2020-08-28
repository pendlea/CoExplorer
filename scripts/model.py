# model.py - Storage access for coexp notebook
# rcampbel@purdue.edu - November 2019

import os
import datetime
import csv
import tempfile
import subprocess
from numpy import mean
import base64
import sys

class Model:

    DATA_DIR       = 'data'
    TEMP_ROOT      = 'temp'
    SAMPLES        = 'samples.txt'
    HOMOLOGS       = 'Genes_homologs.txt'
    ANNOTATIONS    = 'Genes_annotations.txt'
    EXPERIMENTS    = 'DE_experiments.txt'
    EXPRESSIONS    = 'Genes_expression.txt'
    NET_DIR        = 'Networks'
    MODULE         = 'Module_Summary.txt'
    ABC_EXT        = '.abc'

    HEADER_FLAG    = '#'
    GENE_ID        = 'GeneID'
    EMPTY_GENE_IDS = []

    SAMP_COL       = 0
    COND_COL       = 1
    EXP_COL        = 2

    COEXP          = '_Coexp_'
    IMAGE_ERROR    = 'ERROR at open image file: '

    OUTPUT_FORMAT_HEADER = 'Fields per gene ID (row) and condition (column): <avg_exp,p_value,FDR,fold_chg>'

    GREP_CMD = 'grep'

    def __init__(self):
        self.view      = None
        self.ctrl      = None
        self.root      = self.DATA_DIR
        self.exps      = [] # experiment list from samples file - filled when view calls get_samples()
        self.samp_dict = {} # sample dictionary - filled when view calls get_samples()
        self.cond_dict = {} # condition dictionary for samples - filled when view calls get_samples()

        # Dict where key is gene ID and value is experiment data for gene/condition
        self.filter_results = {}

        # Dict where key is gene ID in filter results and value is dictionary of annotations for that gene
        self.filter_results_annos = {}

        # String containing module data for download
        self.module_download_data = ''

        # Read headers out of certain files to create field name lists

        # Get conditions from header in experiments file
        with open(os.path.join(self.root,self.EXPERIMENTS),'rt') as handle:
            self.cond = next(csv.reader(handle,delimiter='\t'),None)[1:] # Skip "#GeneID"

        # Get samples from header in expression file (used by plotting)
        with open(os.path.join(self.root,self.EXPRESSIONS),'rt') as handle:
            self.samples = next(csv.reader(handle,delimiter='\t'),None)[1:] # Skip "#GeneID"

        # Get annotation fields from header in annotations file (used by plotting)
        with open(os.path.join(self.root,self.ANNOTATIONS),'rt') as handle:
            self.anno = next(csv.reader(handle,delimiter='\t'),None) # Do NOT skip "#GeneID"

    def intro(self,view,ctrl):
        self.view = view
        self.ctrl = ctrl

    def get_annos(self,genes):
        '''Get functional annotations of genes'''
        ret = {}

        for line in self.grep_lookup(genes,os.path.join(self.root,self.ANNOTATIONS)):
            annos = {}

            for i in range(1,len(self.anno)):  # Assum gene ID is in 1st col, so skip it
                annos[self.anno[i]] = line[i]

            ret[line[0]] = annos

        return ret

    def get_samples(self):
        '''Provide data for samples tab & save data for use when plotting'''
        with open(os.path.join(self.root,self.SAMPLES),'rt') as handle:

            # Get all sample data - including header

            reader    = csv.reader(handle,delimiter='\t')
            header    = next(reader,None)
            header[0] = header[0][1:] # Remove '#'
            data      = list(reader)  # Read all of remaining file

            exps  = set()

            for row in data:
                sample     = row[self.SAMP_COL]
                condition  = row[self.COND_COL]
                experiment = row[self.EXP_COL]

                # Retain dictionaries to be used by plotting NOTE: Side effect
                self.samp_dict[sample   ] = [condition,experiment] # Tracks how samples are linked to conditions and experiments
                self.cond_dict[condition] = experiment

                exps.add(experiment)

            self.exps = sorted(list(exps)) # Retain list of experiments NOTE: Side effect

            return (header,data)

    def get_expression(self,geneList):
        '''Provide expression data for given gene list'''

        # Build dictionary to store expression data for genes in gene list

        exp_dict      = {}
        samp_exp_dict = {} # Store expression for ALL samples in one dictionary

        for line in self.grep_lookup(geneList,os.path.join(self.root,self.EXPRESSIONS)):

            # Build dictionary where dict[gene][sample]=[[sample,expression],[sample,expression],[sample,expression]...]

            gene_id           = line[0]  # Assume gene ID always in 1st column
            exp_dict[gene_id] = {}

            for condition in self.cond_dict.keys(): # Set up cleared lists to track expression levels per condition
                experiment = self.cond_dict[condition]

                if experiment not in exp_dict[gene_id].keys():
                    exp_dict[gene_id][experiment] = {}

                exp_dict[gene_id][experiment][condition] = []

            samp_exp_dict[gene_id] = {}

            # Create new key per sample ID
            for sample_id in self.samples:
                sample_col_num   = self.samples.index(sample_id)   # Find sample column number
                expression_value = float(line[sample_col_num+1])   # Extract expression level from same column
                condition        = self.samp_dict[sample_id][0]    # Find condition & experiment to which sample belongs
                experiment       = self.samp_dict[sample_id][1]

                # Store in dictionary
                exp_dict[gene_id][experiment][condition].append(expression_value)

                # Populate dictionary to store samples' expression
                samp_exp_dict[gene_id][sample_id] = expression_value

        self.ctrl.debug('Experimental expression values for %i genes stored' % len(exp_dict.keys()))

        # Calc average expression value, per gene, per condition

        mean_expression_dict = {}

        for gene_id in exp_dict.keys():
            mean_expression_dict[gene_id] = {} # Create an average dictionary for this gene

            for condition in self.cond_dict.keys():
                experiment = self.cond_dict[condition] # Get condition's experimental group

                # Create experiment key
                if experiment not in mean_expression_dict[gene_id].keys():
                    mean_expression_dict[gene_id][experiment] = {}

                # Get array of total expression values for this gene under given condition
                all_expression_values = exp_dict[gene_id][experiment][condition]

                # # Store expression value as calculated mean, minimum, and maximum expression
                mean_expression_dict[gene_id][experiment][condition] = [
                    mean(all_expression_values)
                    ,min(all_expression_values)
                    ,max(all_expression_values)
                ]

        self.ctrl.debug('Average expression levels stored for %i genes' % len(mean_expression_dict.keys()))
        return (mean_expression_dict,samp_exp_dict,self.samp_dict)

    def get_networks(self):
        return os.listdir(os.path.join(self.root,self.NET_DIR))

    def get_module_data(self,network,output_widget):
        '''Get module data based on selected network and add calculated data'''
        qry_gen     = list(self.filter_results) # Create list of valid gene IDs (dictionay's keys)
        qry_gen_set = set(qry_gen)
        disp_data   = [] # Create list of module records (lists with all fields as strings) for display
        value_data  = [] # Create list of tupeles (module,recovered) corresponding to disp_data

        # Read thru module summary file, selecting data, and write to module output file

        # NOTE Assumes dir name is exact same as select option.
        #      If they differ, then use NET_PREFIX to match and fix up.
        with open(os.path.join(self.root,self.NET_DIR,network,self.MODULE),'r') as summ_in:

            # Add header line to download string
            self.module_download_data = self.HEADER_FLAG + '\t'.join(self.view.get_module_export_header()) + '\n'

            # Start reading module file
            reader = csv.reader(summ_in,delimiter='\t')
            next(reader,None) # Skip header

            for line in reader:
                mod_gen_set = set(line[3].split())
                recovered   = qry_gen_set.intersection(mod_gen_set)

                if not recovered:
                    continue # Only show lines with at least one recovered gene

                num_qry_mod        = len(qry_gen_set.intersection(mod_gen_set)) # Count query genes in module
                num_both           = len(qry_gen_set.union(       mod_gen_set)) # Count superset of query and mod genes

                output_fields = [
                    #___formating____ _______________value________________    ___________description__________
                                      line[0]                                 # Module ID
                    ,                 line[1]                                 # Quality
                    ,                 line[2]                                 # p-value
                    ,'{:d}'.format(   len(list(mod_gen_set))                ) # Number of genes in module
                    ,'{:d}'.format(   num_qry_mod                           ) # Number of query genes in module
                    ,'{:.4f}'.format( len(qry_gen) / num_both               ) # Jaccard index
                    ,'{:.4f}'.format( num_qry_mod / len(qry_gen) * 100.0    ) # "% Query Genes Recovered"
                    ,','.join(        recovered                             ) # "Recovered Query Genes"
                    ,','.join(        qry_gen_set.difference(mod_gen_set)   ) # Missing query genes
                    ,','.join(        mod_gen_set                           ) # All genes in module
                ]

                # Add to download string
                self.module_download_data += '\t'.join(output_fields) + '\n'

                # Add to return values (for select widget)
                disp_data.append(output_fields[:-2]) # Data user will see. Hide last col, tho is used in export/download
                value_data.append((line[0],list(recovered))) # Data used if line is selected for network plotting

        return (disp_data,value_data)

    def clear_filter_results(self):
        self.filter_results = {}

    def clear_module_download(self):
        self.module_download_data = ''

    def search(self,target_ids,tpm_thresh,pval_thresh,fdr_thresh):
        '''Find data for given gene(s) and condition tests, save to class properites'''

        # Get data for given gene IDs - FILTER 1: Valid gene?
        # Note special case (not target_ids): Here, there's no limit on genes to consider
        for line in self.grep_lookup(target_ids,os.path.join(self.DATA_DIR,self.EXPERIMENTS),(not target_ids)):
            try:
                valid_row  = False # Current line passed all tests?
                parsed_row = []   # All parsed data from current line

                # Run thru records for ea. condition
                # NOTE: Could optimize by building list of non-clear cond tests and only checking those conditions
                for i,csv_data in enumerate(line[1:]): # skip 1st col continaing gene ID
                    record   = csv_data.split(',') # Parsed values for this gene and condition
                    tpm      = float(record[0])
                    min_test = self.view.filter_conditons[i][0].icon # icon = "value" of 3-state button

                    # FILTER 2: Correct TPM value?
                    if ( min_test == self.view.FILT_CLR                       or
                            min_test == self.view.FILT_POS and tpm >= tpm_thresh or
                            min_test == self.view.FILT_NEG and tpm <  tpm_thresh
                        ):
                        pval       = float(record[1])
                        fdr        = float(record[2])
                        fchg       = float(record[3])
                        over_test  = self.view.filter_conditons[i][1].icon  # icon = "value" of 3-state button
                        under_test = self.view.filter_conditons[i][2].icon  #   "       "

                        # FILTER 3: Diff expression: correct p-value, false discovery rate, and fold change values?
                        if  (
                            ( (over_test  == self.view.FILT_CLR                                                                 ) or
                              (over_test  == self.view.FILT_POS and     (pval <= pval_thresh and fdr <= fdr_thresh and fchg > 0)) or
                              (over_test  == self.view.FILT_NEG and not (pval <= pval_thresh and fdr <= fdr_thresh and fchg > 0))    )
                                 and
                            ( (under_test == self.view.FILT_CLR                                                                 ) or
                              (under_test == self.view.FILT_POS and     (pval <= pval_thresh and fdr <= fdr_thresh and fchg < 0)) or
                              (under_test == self.view.FILT_NEG and not (pval <= pval_thresh and fdr <= fdr_thresh and fchg < 0))    )
                            ):
                            parsed_row.append(record) # Accumulate parsed data for gene
                            valid_row = True
                        else:
                            # Failed a diff exp test, don't save this gene, skip remaining conditions
                            valid_row = False
                            break

                    else:
                        # Failed min exp test, don't save this gene
                        valid_row = False
                        break

                # Gene's data passed all tests? Save the parsed data
                if valid_row:
                    self.filter_results[line[0].strip()] = parsed_row # Results indexed by gened ID (line[0])

            except:
                self.ctrl.debug('search(): Error at "'+line[0]+'", exp #'+str(i))
                raise

    def add_annos(self):
        '''Append annotation data to each gene record in search results'''
        self.filter_results_annos = {}

        try:
            # Build text file with search terms
            with tempfile.NamedTemporaryFile(mode='w+t',delete=False) as tmpfile:

                latch = False

                for gene_id in self.filter_results.keys():

                    if latch:
                        tmpfile.write('\n')
                    else:
                        latch = True

                    tmpfile.write(gene_id)

            # Run grep with search terms file
            raw = subprocess.check_output([self.GREP_CMD,'-F','-f',tmpfile.name,os.path.join(self.root,self.ANNOTATIONS)])

            # Get results
            for line in raw.decode('utf-8').split('\n'): # NOTE relies on utf8 encoded

                # Skip headers
                if not line.startswith(self.HEADER_FLAG) and not line.strip() == '':
                    annos = {}
                    line  = line.split('\t')

                    # Creat dictionary
                    for i in range(1,len(self.anno)):  # Assum gene ID is in 1st col, so skip it
                        annos[self.anno[i]] = line[i]

                    self.filter_results_annos[line[0]] = annos

            os.remove(tmpfile.name)

        except Exception as e:
            self.ctrl.debug('add_annos() exception: "'+str(e)+'"')

    def write_filtered_data(self):
        '''Dump filtered gene data to output widget for export'''

        data =  self.HEADER_FLAG + self.HEADER_FLAG + self.OUTPUT_FORMAT_HEADER + '\n'
        data += self.HEADER_FLAG + self.GENE_ID + '\t' + '\t'.join(self.cond) + '\t' + '\t'.join(self.anno[1:]) + '\n'

        try:
            for gene_id,parsed_row in self.filter_results.items():
                line = gene_id

                for record in parsed_row:
                    line += '\t' + ','.join(record)

                for key,value in self.filter_results_annos[gene_id].items():
                    line += '\t' + value

                data += line + '\n'
        except:
            self.ctrl.debug('write_filtered_data() - EXCEPTION:'+sys.exc_info()[0])
            raise

        return data

    def translate_genes(self,genes,funcs):
        '''Translate given gene IDs into native IDs'''

        terms        = [] # List of grep-specific search terms
        homo_results = [] # Results from homologs search
        anno_results = [] # Results from annotations search
        ret          = [] # Return value

        # Use homologs file to translate genes as needed
        homo_results = self.grep_reverse(genes,os.path.join(self.root,self.HOMOLOGS))

        # Use annotations file to translate functions as needed
        anno_results = self.grep_reverse(funcs,os.path.join(self.root,self.ANNOTATIONS))

        self.ctrl.debug('Translation results: '+str(homo_results)+' & '+str(anno_results))

        # Provide the union of results of both grep runs - avoid union with empty set
        if homo_results and anno_results:
            ret = list(set(homo_results) & set(anno_results))
        elif not homo_results and not anno_results:
            ret = []
        elif homo_results:
            ret = homo_results
        else:
            ret = anno_results

        return ret

    def grep_args(self,terms):
        '''Convert list of search terms to list of args for grep command'''
        args = []

        for term in terms:
            args.append('-e')
            args.append(term)

        return args

    def grep_reverse(self,search_terms,search_file):
        '''Grep through given file using given search term list, return list of index (1st col) values'''
        found = []
        args  = ['-i'] + self.grep_args(search_terms)

        if args:

            try:
                raw = subprocess.check_output([self.GREP_CMD] + args + [search_file])

                # Pull index (1st col) out of results
                for line in raw.decode('utf-8').split('\n'): # NOTE relies on encoded
                    first = line.split('\t')[0]

                    # NOTE Should we verify that search term IS in index (1st col)?

                    if not first.startswith(self.HEADER_FLAG) and not first == '':
                        found.append(first)
            except Exception as e:
                self.ctrl.debug('grep_reverse() exception: "'+str(e)+'"')
                found = []

        return found

    def grep_lookup(self,search_terms,search_file,get_all=False):
        '''Grep through given file using given search term, return list of lines w/parsed fields'''
        found = []

        if get_all:
            self.ctrl.debug('grep_lookup(): get_all = True')
            args = ['-F', '-e', '''''']
        else:
            args = self.grep_args(search_terms)

        if args:

            try:
                raw = subprocess.check_output([self.GREP_CMD] + args + [search_file])
                #self.ctrl.debug('grep_lookup(): raw = "'+str(raw)+'"')

                # Pull index (1st col) out of results
                for line in raw.decode('utf-8').split('\n'): # NOTE relies on encoded

                    # Shoudl we verify that search term IS NOT in index (1st col)?

                    if not line.startswith(self.HEADER_FLAG) and not line.strip() == '':
                        found.append(line.split('\t'))

            except Exception as e:
                self.ctrl.debug('grep_lookup() exception: "'+str(e)+'"')
                found = []

        return found


