# view.py - User interface for coexp notebook
# rcampbel@purdue.edu - November 2019

import sys
import ipywidgets as ui
import urllib

class View:

    TABS           = ['Welcome','Samples','Filter','Plot by Experiment','Plot by Co-expression','Plot Differential Expression']
    WELCOME1_TITLE = 'Using This Tool'
    WELCOME1_TEXT  = '''
    <p>In the <b>Samples</b> tab above, you can review the available datasets and conditions</p>
    <p>In the <b>Filter</b> tab, you can search for genes and functional terms of interest. You can also filter by
    minimum and differential expression. Once you applied your filter(s), the data can be viewed and dowloaded.</p>
    <p>Plots of the data can be seen on the <b>Plot by Experiment</b> and <b>Plot by Co-expression</b> tabs.</p>
    '''
    WELCOME2_TITLE = 'Methods'
    WELCOME2_TEXT  = '''
    <b>Gene Differential Expression Analysis</b>
    <p>
    Three prime end sequencing data was obtained for 397 Setaria viridis cultivar ME034V samples from leaf, sheath, 
    and root tissue. Each library was subsequently trimmed and quality filtered using Trimmomatic v0.36 (Bolger, Lohse, 
    and Usadel 2014) to remove the TruSeq three prime adapter, trim the first three bp and last six bp (due to drops in 
    quality scores), perform sliding window quality filtration (4 bp step size, average quality score of 20), filter 
    reads with final whole-read quality scores less than twenty, and filter reads with final minimum read lengths 
    less than 50 bp. Libraries of poor quality, as determined by FASTQC (Andrews 2010), or with less than 1.5 million 
    reads were eliminated from the downstream analysis. Following these filtrations, 354 samples remained for processing 
    through differential expression and co-expression network analyses (see table under the ‘Samples’ tab for the full 
    list of samples and their descriptions)
    </p>
    <p>Quantification of expression was performed using Kallisto (--single --single-overhang -l
    200 -s 30 -t 5; Bray et al. 2016) using version 2.0 of the ME034V gene annotation set
    (unpublished), an improvement on the ME034V version 1.0 (Thielen and Pendleton et al. 2020).
    Differential expression analysis was performed using the EdgeR R package (Robinson, McCarthy,
    and Smyth 2010) in paired tests of test versus control conditions.</p>
    <b>Co-Expression Network Construction and Analysis</b>
    <p>The expression results (CPMs; counts per million) from Kallisto were used as input
    for the Mutual Ranks to Modules co-expression network pipeline
    (<a href="https://github.rcac.purdue.edu/jwisecav/coexp-pipe" target="_blank">coexp-pipe</a>).
    Genes without expression support (CPM > 0) in at least three libraries or with less than
    ten reads mapped across all libraries were excluded from the gene expression matrix.
    The expression matrix was then transformed using  the Variance Stabilized Transformation
    (VST) as implemented in the R package DeSeq2 (Love, Huber, and Anders 2014). Pairwise
    Pearson’s correlation coefficients were calculated between all possible gene pairs,
    followed by calculation of mutual rank (MR) scores. MR scores were transformed to network
    edge weights using the exponential decay function e^(-(MR-1/x)); three different networks
    were constructed with x set to 5, 10, and 25, respectively. Finally, co-expressed modules
    are delineated using ClusterONE (Nepusz, Yu, and Paccanaro 2012). Average expression values
    (post VST) are plotted in the ‘Plot by Experiment’ tab by experiment.</p>
    <b>Works Cited</b>
    <p>Andrews, S. (2010) "FastQC: a quality control tool for high throughput sequence data".</p>
    <p>Bray, N.L., Pimentel, H., Melsted, P., and Pachter, L. (2016) “Near-optimal probabilistic
    RNA-seq quantification.” Nature Biotechnology 34: 525–527.</p>
    <p>Love, M.I., Huber, W., & Anders, S. (2014). “Moderated estimation of fold change and
    dispersion for RNA-seq data with DESeq2.” Genome biology, 15(12): 550.</p>
    <p>Nepusz, T., Yu, H., & Paccanaro, A. (2012) “Detecting overlapping protein complexes in
    protein-protein interaction networks.” Nature Methods, 9: 471.</p>
    <p>Robinson, M.D., McCarthy, D.J., & Smyth, G.K. (2010) “edgeR: a Bioconductor package
    for differential expression
    analysis of digital gene expression data." Bioinformatics 26(1): 139-140.</p>
    <p>Thielen, P.M., Pendleton, A.L., Player, R.A., Bowden, K.V., Lawton, T.J., & Wisecaver,
    J.H. (2020). “Reference genome for the highly transformable Setaria viridis cultivar
    ME034V.” bioRxiv.</p>
    '''
    SAMPLES1_TITLE = 'Samples'
    FILTER1_TITLE  = 'Filter by Gene ID'
    #FILTER1_TEXT   = '<b>Provide either Setaria or Arabidopsis gene IDs</b>'
    FILTER1_TEXT   = '<b>Provide Arabidopsis thaliana gene IDs (eg: AT4G13615). </b>'
    FILTER1_TEXT_2 = 'You can also search for OrthoFinder-identified orthologs to the following species/genomes: Setaria viridis ME034V (e.g. Svm5G0050410), Setaria viridis A10 (eg: Sevir.5G400800), S. italica (eg: Seita.1G019400), sorghum (eg: Sobic.006G004700), rice (eg: Os03g27570), corn (eg: Zm00001d024902), Arabidopsis (eg: AT2G46040). Search with single gene, or comma-separated list of genes.'
    FILTER1_HINT   = 'e.g. Sevir.5G400800, AT1G44575'
    FILTER2_TITLE  = 'Filter by Function'
    #FILTER2_TEXT   = '<b>Provide either GO ID or KEGG PATHWAY ID</b>'
    FILTER2_TEXT   = '<b>Provide keyword (eg: cytochrome p450), GO ID (eg: GO:0009522), Interpro ID (eg: IPR013087), KEGG ID (eg: K01835), or E.C. code (eg: EC:2.1.1.112).<b>'
    FILTER2_HINT   = 'e.g. GO:0009522, KOG0190'
    FILTER3_TITLE  = 'Filter by Expression'
    FILTER4_TITLE  = 'Filtered Gene List'
    FILTER5_TEXT   = '<b>Set minimally expressed threshold</b>'
    FILTER6_TEXT   = 'TPM'
    FILTER7_TEXT   = '<b>Set differentially expressed thresholds</b>'
    FILTER8_TEXT   = 'P-value'
    FILTER9_TEXT   = 'FDR '
    DEFAULT_TPM    = '1'
    DEFAULT_PVAL   = '0.05'
    DEFAULT_FDR    = '0.05'
    FILTER10_TEXT  = '<b>Select filters</b>'
    FILTER11_TEXT  = 'Minimally expressed'
    FILTER12_TEXT  = 'Over expressed'
    FILTER13_TEXT  = 'Under expressed'
    FILTER14_TEXT  = 'Filtered Gene List'
    FILTER15_TEXT  = 'Show'
    FILTER16_TEXT  = 'entries'
    FILTER17_TEXT  = 'Gene ID'
    FILTER22_TEXT  = 'GO terms'
    FILTER23_TEXT  = 'KEGG pathways'
    FILTER24_TITLE = 'Export'
    FILTER_APPLY   = 'Apply Filter'
    FILTER_DOWNLD  = 'Download&nbspSelected&nbspData'
    FILTER_PROG    = '<br>Filtering data...'
    FILTER_PROG_ALL= '<br>Assembling ALL data! WARNING: This may take some time...'
    FILTER_REFEXP  = 'Create Download Link'

    PLOTEX1_TITLE  = 'Plot Filtered Data'
    PLOTEX1_TEXT   = 'Plotting is restricted to 100 genes.'
    PLOTEX2_TITLE  = 'Plot Line Graphs'
    PLOTEX2_TEXT   = '<b>Select experiment</b>'
    PLOTEX3_TITLE  = 'Plot Heatmaps'
    PLOTEX3_TEXT   = '<b>Select experiment</b>'
    PLOTCO1_TITLE  = 'Select Network Graph'

    PLOTCO1_TEXT_A = '''
    <p>This tab displays the co-expressed gene sets (i.e., modules) from global gene co-expression networks. Modules listed below contain one or more query genes (specified in the Filter tab).
    </p>
    <p>Global co-expression networks were created using the Mutual Rank (MR) score first described by <a href="https://www.ncbi.nlm.nih.gov/pubmed/19767600/" target="_blank">Obayashi & Kinoshita</a>. Exponential decay functions were used to create MR-transformed edge weights as described in here: https://github.rcac.purdue.edu/jwisecav/coexp-pipe and here: [Wisecaver et al. 2017 Plant Cell](http://www.plantcell.org/content/29/5/944).
    </p>
    <p>We impose an edge weight cutoff of ≥ 0.01 and vary the exponential decay denominator, thereby varying the number of edges (ie gene associations) retained in the networks.
    </p>
    '''
    PLOTCO1_TEXT_B = '''
    </p>
    <p>Using 5 as the decay constant creates the fastest rate of decay (Fig 1: light blue curve); MR scores ≤ approx 27 are included as edge weights. Using 25 as the decay constant creates the slowest rate of decay (Fig 1: light green curve); MR scores ≤ approx 130 are included. Using 10 as the decay constant puts you somewhere in between (dark blue curve).
    </p>
    <p>Different rates of decay result in different average module sizes. Fast rate of decay = small modules; Slow rate of decay = large modules. <b>This makes different networks perform better for different metabolic pathways.</b> If your pathway of interest consists of fewer genes, a faster rate of decay may work better for you, because the modules will be smaller and may contain less false positives. If your pathway of interest consists of more genes, or if you want to recover a large list of gene candidates, a slower rate of decay may work better for you, because the modules will be larger.
    </p>
    '''
    PLOTCO1_LABEL  = '''<b>Select a decay denominator</b>'''

    PLOTCO2_TITLE  = 'Network Graph of Selected Module'
    PLOTCO3_TEXT   = '<br><b>Select a module</b>'

    PLOTDE1_TEXT   = 'Select gene(s) to plot. Use Ctrl (Cmd) key while selecting additional genes. Use Shift key to select a range of genes.'
    PLOTDE2_TEXT   = 'Create Plot(s)'
    PLOTDE3_TEXT   = 'Select Gene(s) and Plot Differential Expression'
    PLOTDE4_TEXT   = 'Genes and Annotations'

    PLOT_DOWNLOAD  = 'Download Selected Plot'
    MODULE_HEADER  = [
        ["          ","       ","           ","Genes in","Query Genes","Jaccard","% Query Genes","Recovered                      ","Missing Query                                       ","All Genes                                       "],
        ["Module    ","Quality","p-value    ","Module  ","in Module  ","Index  ","Recovered    ","Query Genes                    ","Genes                                               ","in Module                                       "]
    ]
    MODULE_PROG    = '(Working...)'
    NO_MODULE_DATA = '(Filtered gene list is empty.)'
    EMPTY_LIST_MSG = '''<br>(There's no data to display. Press "Apply Filter" to search using filter criteria specified above.)'''
    ALL            = 'All'
    EMPTY          = ''      # For dropdown item

    MODULE_DOWNLD  = 'Download&nbspSelected&nbspData'

    LO10 = ui.Layout(width='10%')
    LO15 = ui.Layout(width='15%')
    LO20 = ui.Layout(width='20%')
    LO25 = ui.Layout(width='25%')

    # Support matrix of custom, 3-state buttons....
    FILT_CLR               = 'XXX'     # Button value
    FILT_POS               = 'check'   #  "
    FILT_NEG               = 'times'   #  "
    FILTER_CONDS_TITLE     = u'<br><b>Set expression tests per condition</b> (\u2714: "must be", \u2717: "must NOT be", clear: any value accepted)'
    FILTER_CONDS_CON_HDR   = '<b>Conditon</b>'
    FILTER_CONDS_MIN_HDR1  = 'Minimally'
    FILTER_CONDS_MIN_HDR2  = 'Expressed'
    FILTER_CONDS_DIF_HDR1  = 'Differentially Expressed'
    FILTER_CONDS_DIF_HDR2O = 'Over'
    FILTER_CONDS_DIF_HDR2U = 'Under'
    FILTER_CONDS_COL_WIDTH = '100px'
    FILTER_CONDS_COL_DBLWI = '200px'

    def __init__(self):
        self.tabs    = None # Main UI container

    def intro(self,model,ctrl):
        self.model = model
        self.ctrl  = ctrl

    def display(self,debug=False):
        '''Build and show notebook user interface'''
        self.build()

        if debug:
            self.debug_output = ui.Output(layout={'border': '1px solid black'})
            display(ui.VBox([self.tabs,self.section('Debug',[self.debug_output])]))
        else:
            display(self.tabs)

    def debug(self,text):
        with self.debug_output:
            print(text)

    def section(self,title,contents):
        '''Create a collapsible widget container'''

        if type(contents) == str:
            contents = [ui.HTML(value=contents)]

        ret = ui.Accordion(children=tuple([ui.VBox(contents)]))
        ret.set_title(0,title)
        return ret

    def columnize(self,guide,data,sep='  '):
        '''Pad or truncate columns in text lines based on guide'''
        rows = []

        #Sort by column with query count, in reverse order
        #data.sort(key=lambda x: x[4], reverse=True)
        
        
        for r in range(len(data)):
            row = ''

            for c in range(len(data[r])): # guide)):
                text      = data[r][c]
                col_width = len(guide[c])
                diff      = col_width - len(data[r][c])

                if diff < 0:
                    row += text[:diff]         + sep
                else:
                    row += text + (' ' * diff) + sep

            rows.append(row)

        return rows

    def build(self):
        '''Create user interface widgets'''

        tabs = [];

        # Tab 1: Welcome ======================================================
        content = []
        content.append(self.section(self.WELCOME1_TITLE,self.WELCOME1_TEXT))
        content.append(self.section(self.WELCOME2_TITLE,self.WELCOME2_TEXT))
        tabs.append(ui.VBox(content))

        # Tab 2: Samples ======================================================

        content = []

        # CSS (stle) and start the table
        table = '''<style>
                        .sample_cell {padding-right: 32px     }
                        .sample_even {background   : White    }
                        .sample_odd  {background   : Gainsboro}
                    </style><table>'''

        samp_hdrs,samp_data = self.model.get_samples()

        for item in samp_hdrs:
            table += '<th class="sample_cell">' + item + '</th>'

        for i,line in enumerate(samp_data):

            if i % 2 == 0:
                table += '<tr class="sample_even">'
            else:
                table += '<tr class="sample_odd">'

            for item in line:
                table += '<td class="sample_cell">' + item + '</td>'

            table += '</tr>'

        table += '</table>'

        widgets = []
        widgets.append(ui.HTML(value=table))
        content.append(self.section(self.SAMPLES1_TITLE,widgets))

        tabs.append(ui.VBox(content))

        # Tab 3: Filter ======================================================

        self.filter_txt_gene    = ui.Text(description='',value='',placeholder=self.FILTER1_HINT  ) # Test: Sevir.9G584100/5G108900
        self.filter_txt_func    = ui.Text(description='',value='',placeholder=self.FILTER2_HINT  )
        self.filter_txt_tpm     = ui.Text(description=self.FILTER6_TEXT,value=self.DEFAULT_TPM ,placeholder='NA')
        self.filter_txt_pval    = ui.Text(description=self.FILTER8_TEXT,value=self.DEFAULT_PVAL,placeholder='NA')
        self.filter_txt_fdr     = ui.Text(description=self.FILTER9_TEXT,value=self.DEFAULT_FDR ,placeholder='NA')
        self.filter_conditons   = [] # Created in loop below
        self.filter_btn_apply   = ui.Button(description=self.FILTER_APPLY,icon='filter',layout=self.LO20)
        self.filter_ddn_ndisp   = ui.Dropdown(options=['25','50','100',self.ALL],layout=self.LO10)
        self.filter_html_output = ui.HTML(self.EMPTY_LIST_MSG)
        self.filter_btn_refexp  = ui.Button(description=self.FILTER_REFEXP,icon='download',layout=self.LO20)

        content = []

        widgets = []
        widgets.append(ui.HTML(value=self.FILTER1_TEXT))
        widgets.append(ui.HTML(value=self.FILTER1_TEXT_2))
        widgets.append(self.filter_txt_gene)
        content.append(self.section(self.FILTER1_TITLE,widgets))

        widgets = []
        widgets.append(ui.HTML(value=self.FILTER2_TEXT))
        widgets.append(self.filter_txt_func)
        content.append(self.section(self.FILTER2_TITLE,widgets))

        widgets = []
        widgets.append(ui.HTML(value=self.FILTER5_TEXT))
        widgets.append(self.filter_txt_tpm)
        widgets.append(ui.HTML(value=self.FILTER7_TEXT))
        widgets.append(self.filter_txt_pval)
        widgets.append(self.filter_txt_fdr )

        # Conditions ("Set expressions tests..."): matrix of 3-state buttons

        # Title
        widgets.append(ui.HTML(value=self.FILTER_CONDS_TITLE))

        # Column headers (2 rows of headers)

        widgets.append(ui.HTML('<style> .ctrhdr {text-align: center; font-weight: bold;} </style>'))
        div = '<div class="ctrhdr">'
        vid = '</div>'

        row = []
        row.append(ui.Label(value=''                         ,layout=self.LO15                                   ))
        row.append(ui.HTML(div+self.FILTER_CONDS_MIN_HDR1+vid,layout=ui.Layout(width=self.FILTER_CONDS_COL_WIDTH)))
        row.append(ui.HTML(div+self.FILTER_CONDS_DIF_HDR1+vid,layout=ui.Layout(width=self.FILTER_CONDS_COL_DBLWI)))
        widgets.append(ui.HBox(row))

        row = []
        row.append(ui.HTML(    self.FILTER_CONDS_CON_HDR      ,layout=self.LO15                                   ))
        row.append(ui.HTML(div+self.FILTER_CONDS_MIN_HDR2 +vid,layout=ui.Layout(width=self.FILTER_CONDS_COL_WIDTH)))
        row.append(ui.HTML(div+self.FILTER_CONDS_DIF_HDR2O+vid,layout=ui.Layout(width=self.FILTER_CONDS_COL_WIDTH)))
        row.append(ui.HTML(div+self.FILTER_CONDS_DIF_HDR2U+vid,layout=ui.Layout(width=self.FILTER_CONDS_COL_WIDTH)))
        widgets.append(ui.HBox(row))

        # Buttons
        for item in self.model.cond:
            row = []
            row.append(ui.Label(value=item         ,layout=self.LO15                                   ))
            row.append(ui.Button(icon=self.FILT_CLR,layout=ui.Layout(width=self.FILTER_CONDS_COL_WIDTH)))
            row.append(ui.Button(icon=self.FILT_CLR,layout=ui.Layout(width=self.FILTER_CONDS_COL_WIDTH)))
            row.append(ui.Button(icon=self.FILT_CLR,layout=ui.Layout(width=self.FILTER_CONDS_COL_WIDTH)))
            self.filter_conditons.append((row[1],row[2],row[3]))
            widgets.append(ui.HBox(row))

        content.append(self.section(self.FILTER3_TITLE,widgets))

        # "Filtered gene list"

        widgets = []

        row = []
        row.append(self.filter_btn_apply)
        row.append(ui.HTML('<div style="text-align: right;">'+self.FILTER15_TEXT+'</div>',layout=self.LO15))
        row.append(self.filter_ddn_ndisp)
        row.append(ui.HTML('<div style="text-align: left;">' +self.FILTER16_TEXT+'</div>',layout=self.LO10))
        widgets.append(ui.HBox(row))

        widgets.append(ui.HBox([self.filter_html_output],layout={'width':'90vw'})) # Cause horiz. scrollbar (was: 1200)

        content.append(self.section(self.FILTER4_TITLE,widgets))

        # "Export"

        self.filter_out_export = ui.Output(layout={'border': '1px solid black'})
        row = self.section(self.FILTER24_TITLE,[ui.VBox([self.filter_btn_refexp,self.filter_out_export])])
        #row.selected_index = None
        content.append(row)

        tabs.append(ui.VBox(content))

        # Tab 4: Plot-experiment ======================================================

        content = []
        content.append(self.section(self.PLOTEX1_TITLE,self.PLOTEX1_TEXT))

        # Plot line graphs ---

        self.plotex_ddn_selex_lg = ui.Dropdown(options=[self.EMPTY]+self.model.exps,value=None,disabled=True)

        widgets = []

        row = []
        row.append(ui.HTML(value=self.PLOTEX2_TEXT))
        row.append(ui.Label(value='',layout=ui.Layout(width='60%'))) # Cheat: spacer
        widgets.append(ui.HBox(row))

        widgets.append(self.plotex_ddn_selex_lg)
        widgets.append(self.ctrl.plotter.line_plot) # Use widget from plotter
        content.append(self.section(self.PLOTEX2_TITLE,widgets))

        # Plot heatmaps ---

        self.plotex_ddn_selex_hm = ui.Dropdown(options=[self.EMPTY,self.ALL]+self.model.exps,value=None,disabled=True)
        self.plotex_img_dispp_hm = ui.Output()

        self.ctrl.plotter.out_plot_msg(self.plotex_img_dispp_hm,self.ctrl.plotter.HEAT_INIT_TITLE)

        widgets = []

        row = []
        row.append(ui.HTML(value=self.PLOTEX3_TEXT))
        row.append(ui.Label(value='',layout=ui.Layout(width='60%'))) # Cheat: spacer
        widgets.append(ui.HBox(row))

        widgets.append(self.plotex_ddn_selex_hm)
        widgets.append(self.plotex_img_dispp_hm)
        content.append(self.section(self.PLOTEX3_TITLE,widgets))

        tabs.append(ui.VBox(content))

        # Tab 5: Plot-coexpression ======================================================

        # Get list of networks and sort them using inline function that parses out numeric portion of network name
        networks = sorted(self.model.get_networks(),key=lambda x: int(x.split('_')[1]))

        # Build display/value pairs as list of tuples for use in dropdown menu

        disp_vals = []

        for i,net in enumerate(networks):
            #disp_vals.append((net.replace('_',' '),net) )
            disp_vals.append((net.split('_')[1],net) )

        self.plotco_ddn_netw  = ui.Dropdown(options=[(self.EMPTY,self.EMPTY)]+disp_vals,value=None,disabled=True)

        # Download button
        self.plotco_btn_modu = ui.Button(description=self.FILTER_REFEXP,icon='download',layout=self.LO20)

        # Create a fixed width data selection widet
        self.plotco_sel_modu = ui.Select(rows=10,options=[],value=None,layout={'width':'99%'})
        self.plotco_sel_modu.add_class('selmono') # Use JavaScript to spec fixed-width font (see custom CSS)

        # Module selection # NOTE If using latest ver of Jupyter, consider using grid widet instead of fixed font select widget

        # Create a fixed width header widget that will match data widget
        header = [self.MODULE_HEADER[0][:-2],self.MODULE_HEADER[1][:-2]] # Hide last 2 cols, tho they're used in export/download
        header = self.columnize(self.MODULE_HEADER[0][:-2],header)
        header = ui.Select(options=header,disabled=True,value=None,layout=ui.Layout(height='3em',width='99%'))
        header.add_class('selmono') # Use JavaScript to spec fixed-width font (see custom CSS)

        content = []

        # Image
        with open('images/DecayRates.png','rb') as handle:
            image = handle.read()

        widgets = []

        widgets.append(ui.HTML(self.PLOTCO1_TEXT_A))
        widgets.append(ui.Image(value=image,format='png',width=500,height=259))
        widgets.append(ui.HTML(self.PLOTCO1_TEXT_B))

        content.append(self.section(self.PLOTEX1_TITLE,widgets))

        # Select Network Graph

        widgets = []
        widgets.append(ui.HTML(value=self.PLOTCO1_LABEL))
        widgets.append(self.plotco_ddn_netw)

        row = []
        row.append(ui.HTML(value=self.PLOTCO3_TEXT))
        row.append(ui.Label(value='',layout=ui.Layout(width='60%'))) # Cheat: spacer
        widgets.append(ui.HBox(row))

        widgets.append(ui.HBox([ui.VBox([header,self.plotco_sel_modu],layout={'width':'100%'})])) # was: ,layout={'width':'90vw'}

        content.append(self.section(self.PLOTCO1_TITLE,widgets))

        # "Export"
        self.plotco_out_export = ui.Output(layout={'border': '1px solid black'})
        row = self.section(self.FILTER24_TITLE,[ui.VBox([self.plotco_btn_modu,self.plotco_out_export])])
        row.selected_index = None
        content.append(row)

        # Network Graph of Selected Module NOTE Uses widget created by plotter
        content.append(self.section(self.PLOTCO2_TITLE,[self.ctrl.plotter.net_plot]))

        tabs.append(ui.VBox(content))

        # Tab 6: Plot-differentialy  ======================================================

        content = []
        content.append(self.section(self.PLOTEX1_TITLE,self.PLOTDE1_TEXT))

        self.plotdif_sel_title = ui.HTML(value=self.PLOTDE4_TEXT,layout={'width':'99%'})
        self.plotdif_sel_genes = ui.SelectMultiple(rows=10,options=[],value=[],layout={'width':'99%'})
        self.plotdif_img_disp  = ui.Output()
        self.plotdif_btn_plot  = ui.Button(description=self.PLOTDE2_TEXT,icon='bar-chart',layout=self.LO20)

        self.ctrl.plotter.out_plot_msg(self.plotdif_img_disp,self.ctrl.plotter.DIF_INIT_TITLE)

        widgets = []
        widgets.append(self.plotdif_sel_title)
        widgets.append(self.plotdif_sel_genes)
        widgets.append(self.plotdif_btn_plot)
        widgets.append(ui.Label(value='',layout=ui.Layout(width='60%'))) # Cheat: spacer
        widgets.append(self.plotdif_img_disp)
        widgets.append(ui.Label(value='',layout=ui.Layout(width='60%'))) # Cheat: spacer

        content.append(self.section(self.PLOTDE3_TEXT,widgets))

        tabs.append(ui.VBox(content))

        # Tab bar ======================================================

        self.tabs = ui.Tab()

        for i in range(len(self.TABS)):
            self.tabs.set_title(i,self.TABS[i])

        self.tabs.children = tuple(tabs)

    def update_filtered_gene_list(self):
        '''Update filtered genes list (and diff exp selection) with new data '''

        # Calc output line limit
        if self.filter_ddn_ndisp.value == self.ALL:
            limit = sys.maxsize
        else:
            limit = int(self.filter_ddn_ndisp.value)

        # CSS (style)

        output = '''<style>
                    .op th {
                        padding    : 3px;
                        border     : 1px solid black;
                        font-size  : 6px !important;
                        text-align : center;
                        line-height: 14px;
                        background-color: lightgray;
                    }
                    .op td {
                        padding    : 3px;
                        border     : 1px solid black;
                        font-size  : 6px !important;
                        text-align : left;
                        line-height: 12px;
                    }
                    </style>'''

        # Table start and header start
        output += '<br><table class="op"><tr><th>'+self.FILTER17_TEXT+'</th>'

        # Column headers
        for anno in self.model.anno[1:]:  # Skip first header since its for gene ID
            output         += '<th class="op">'+anno+'</th>'

        output += '</tr>'

        try:
            # Also build new options for gene selection in diff. exp. plotting
            plotdif_options = []

            # Build table rows
            for count,(gene_id,annos) in enumerate(self.model.filter_results_annos.items()):
                output       += '<tr><td class="op">'+gene_id+'</td>'
                plotdif_line  = gene_id + ': '

                for key,value in annos.items():
                    output       += '<td class="op">'+value+'</td>'
                    plotdif_line += value + ' '

                output += '</tr>'

                plotdif_options.append((plotdif_line,gene_id)) # [('One',1),('Two',2)]

                if count+1 >= limit:
                    break

            output += '</table>' # End table

            self.filter_html_output.value  = output             # Update filtered gene list (search results)
            self.plotdif_sel_genes.options = plotdif_options    # Update diff exp gene menu
        except:
            self.ctrl.debug('update_filtered_gene_list(): EXCEPTION')
            raise

    def set_plot_status(self,enable):
        '''Change status of plot-related widgets based on availability of filter results'''

        pltr = self.ctrl.plotter
        pltr.clear_plots(self.plotex_img_dispp_hm,self.plotdif_img_disp)

        self.plotdif_sel_genes.value = []

        if enable:
            self.plotex_ddn_selex_lg.disabled = False
            self.plotex_ddn_selex_hm.disabled = False
            self.plotco_ddn_netw.disabled     = False
            self.plotdif_sel_genes.disabled    = False

            pltr.line_plot.layout.title = pltr.LINE_PROMPT_TITLE
            pltr.net_plot.layout.title  = pltr.NET_PROMPT_TITLE

            pltr.out_plot_msg(self.plotex_img_dispp_hm,pltr.HEAT_PROMPT_TITLE)
            pltr.out_plot_msg(self.plotdif_img_disp   ,pltr.DIF_PROMPT_TITLE )
        else:
            self.plotex_ddn_selex_lg.disabled = True
            self.plotex_ddn_selex_hm.disabled = True
            self.plotco_ddn_netw.disabled     = True
            self.plotdif_sel_genes.disabled    = True

            pltr.line_plot.layout.title = pltr.LINE_INIT_TITLE
            pltr.net_plot.layout.title  = pltr.NET_INIT_TITLE

            pltr.out_plot_msg(self.plotex_img_dispp_hm,pltr.HEAT_INIT_TITLE)
            pltr.out_plot_msg(self.plotdif_img_disp   ,pltr.DIF_INIT_TITLE )
            self.ctrl.set_module_data([self.NO_MODULE_DATA],[(None,None)])

    def get_module_export_header(self):
        '''Generate module output header for export'''
        ret = []

        for i in range(len(self.MODULE_HEADER[1])):
            pre   = self.MODULE_HEADER[0][i].strip()
            title = self.MODULE_HEADER[1][i].strip()

            if not pre == '':
                title = pre + ' ' + title

            ret.append(title)

        return ret

    def output_data_link(self,output_widget,data_str):
        '''Create data URI link to download data'''

        pre  = '<a download="coexplorer.csv" target="_blank" href="data:text/csv;charset=utf-8,'
        post = '">Download</a>'

        with output_widget:
            display(ui.HTML(pre+urllib.parse.quote(data_str)+post))



