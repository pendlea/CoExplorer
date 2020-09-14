# plotter.py - Co-expression plots
# A. Pendleton, November 2019

import networkx as nx  # Plotting networks
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

# For heatmap plots
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns
from IPython.display import clear_output

class Plotter:

    PLOT_LIMIT          =  100

    COLORS              =  [
        "#63b598", "#ce7d78", "#ea9e70", "#a48a9e", "#c6e1e8", "#648177", "#0d5ac1",
        "#f205e6", "#1c0365", "#14a9ad", "#4ca2f9", "#a4e43f", "#d298e2", "#6119d0",
        "#d2737d", "#c0a43c", "#f2510e", "#651be6", "#79806e", "#61da5e", "#cd2f00",
        "#9348af", "#01ac53", "#c5a4fb", "#996635", "#b11573", "#4bb473", "#75d89e",
        "#2f3f94", "#2f7b99", "#da967d", "#34891f", "#b0d87b", "#ca4751", "#7e50a8",
        "#c4d647", "#e0eeb8", "#11dec1", "#289812", "#566ca0", "#ffdbe1", "#2f1179",
        "#935b6d", "#916988", "#513d98", "#aead3a", "#9e6d71", "#4b5bdc", "#0cd36d",
        "#250662", "#cb5bea", "#228916", "#ac3e1b", "#df514a", "#539397", "#880977",
        "#f697c1", "#ba96ce", "#679c9d", "#c6c42c", "#5d2c52", "#48b41b", "#e1cf3b",
        "#5be4f0", "#57c4d8", "#a4d17a", "#225b8e", "#be608b", "#96b00c", "#088baf",
        "#f158bf", "#e145ba", "#ee91e3", "#05d371", "#5426e0", "#4834d0", "#802234",
        "#6749e8", "#0971f0", "#8fb413", "#b2b4f0", "#c3c89d", "#c9a941", "#41d158",
        "#fb21a3", "#51aed9", "#5bb32d", "#807fbc", "#21538e", "#89d534", "#d36647",
        "#7fb411", "#0023b8", "#3b8c2a", "#986b53", "#f50422", "#983f7a", "#ea24a3",
        "#79352c", "#521250"
    ]

    GENE_ID_LBL         = 'Gene ID'

    LINE_INIT_TITLE     = '(No filter results.)'
    LINE_PROMPT_TITLE   = '(No experiment selected.)'
    LINE_UPDATE_TITLE   = 'Updating plot...'
    LINE_PREFIX_TITLE   = 'Average Gene Expression - '
    LINE_THRESH_TOP     = 50  # Above which need to use complex lines (dot, dash, etc)
    LINE_THRESH_MID     = 20  # Above which need to use semi-complex lines (eg, dot)
    LINE_AXIS_TITLE_X   = 'Experimental Conditions'
    LINE_AXIS_TITLE_Y   = 'Expression'
    LINE_GRID_COLOR     = 'whitesmoke'
    LINE_BACKGROUND     = 'white'

    HEAT_INIT_TITLE     = '(No filter results.)'
    HEAT_PROMPT_TITLE   = '(No experiment selected.)'
    HEAT_UPDATE_TITLE   = 'Updating plot...'
    HEAT_PREFIX_TITLE   = 'Heatmap - '
    HEAT_FONT_SIZE      = 10
    HEAT_FIG_SIZE       = (22,24)
    HEAT_LESS_THAN_TWO  = '(Less than two results.)'

    NET_INIT_TITLE      = '(No filter results.)'
    NET_PROMPT_TITLE    = '(No module selected.)'
    NET_UPDATE_TITLE    = 'Updating plot...'
    NET_PREFIX_TITLE    = 'Network '
    NET_TITLE_FONT_SIZE = 16
    NET_MARKER_SIZE     = 22
    NET_EDGE_WITH       = 0.5
    NET_EDGE_COLOR      = '#888'
    NET_BACKGROUND      = 'rgba(0,0,0,0)'
    NET_HEIGHT          = 600
    NET_WIDTH           = 1200
    NET_CREDIT          = "Python code: <a href='https://plot.ly/ipython-notebooks/network-graphs/'> https://plot.ly/ipython-notebooks/network-graphs/</a>"
    NET_MARGIN          = dict(b=20,l=5,r=5,t=40)
    NET_ANNO_X          = 0.005
    NET_ANNO_Y          = -0.002
    NET_GENE_COLOR_EMPH = 'dodgerblue'
    NET_GENE_COLOR_NORM = 'lightgray'

    DIF_INIT_TITLE      = '(No filter results.)'
    DIF_PROMPT_TITLE    = '(No gene selected.)'
    DIF_UPDATE_TITLE    = 'Updating plot...'
    DIF_PREFIX_TITLE    = 'Differential Expression:'
    DIF_FONT_SIZE       = 10
    DIF_FIG_SIZE        = (22,24)


    def __init__(self,ctrl):
        self.ctrl = ctrl
        init_notebook_mode(connected=False)

        # Line plot widget =========================================================

        # Create data as list of empty Scatter objects

        data        = []
        line_styles = ['solid','dot','dash']

        for i in range(self.PLOT_LIMIT):

            if   i > self.LINE_THRESH_TOP: line_style = line_styles[i % len(line_styles)  ]
            elif i > self.LINE_THRESH_MID: line_style = line_styles[i % len(line_styles)-1]
            else                         : line_style = line_styles[0]

            scatter = go.Scatter(
                x     = []
                ,y    = []
                ,name = ''
                ,mode = 'lines+markers'
                ,line = dict(
                    width  = 4
                    ,dash  = line_style
                    ,color = self.COLORS[i]
                )
            )

            data.append(scatter)

        self.line_plot = go.FigureWidget(
            data    = data
            ,layout = go.Layout(
                title         = self.LINE_INIT_TITLE
                ,yaxis        = dict(
                    title      = self.LINE_AXIS_TITLE_Y
                    ,showgrid  = True
                    ,gridcolor = self.LINE_GRID_COLOR
                )
                ,xaxis        = dict(
                    title      = self.LINE_AXIS_TITLE_X
                    ,showgrid  = True
                    ,gridcolor = self.LINE_GRID_COLOR
                )
                ,plot_bgcolor = self.LINE_BACKGROUND
                ,dragmode     = 'select'
            )
        )

        # Network plot widget =========================================================

        # Prep edge traces for later use
        edge_trace = go.Scatter3d(
            x         = [],
            y         = [],
            z         = [],
            line      = dict(width=self.NET_EDGE_WITH,color=self.NET_EDGE_COLOR),
            hoverinfo = 'none',
            mode      = 'lines')

        # Prep node traces for later use
        node_trace = go.Scatter3d(
            x          = []
            ,y         = []
            ,z         = []
            ,text      = []
            ,mode      = 'markers'
            ,hoverinfo = 'text'
            ,marker    = dict(
                color  = []
                ,size  = self.NET_MARKER_SIZE
            )
        )

        axis = dict(
            showbackground  = False
            ,showline       = False
            ,zeroline       = False
            ,showgrid       = False
            ,showticklabels = False
            ,title          = ''
        )

        # Create network widget
        self.net_plot = go.FigureWidget(
            data        = [edge_trace,node_trace]
            ,layout     = go.Layout(
                title        = self.NET_INIT_TITLE
                ,titlefont   = dict(size=self.NET_TITLE_FONT_SIZE)
                ,showlegend  = False
                ,hovermode   = 'closest'
                ,margin      = self.NET_MARGIN
                ,annotations = [ dict(
                        text       = self.NET_CREDIT
                        ,showarrow = False
                        ,xref      = 'paper', yref = 'paper'
                        ,x         = self.NET_ANNO_X
                        ,y         = self.NET_ANNO_X
                ) ]
                ,plot_bgcolor  = self.NET_BACKGROUND
                ,scene         = dict(
                                    xaxis  = dict(axis)
                                    ,yaxis = dict(axis)
                                    ,zaxis = dict(axis)
                )
                ,height = self.NET_HEIGHT
                ,width  = self.NET_WIDTH
            )
        )

    def out_plot_msg(self,output_widget,text):
        '''Replace current plot output with message'''
        with output_widget:
            clear_output(wait=True)
            print(text)

    def limit_num_genes(self,gene_list):
        '''Truncate list of genes if needed'''
        self.ctrl.debug('Plot limit: %i max, %i genes in list' % (self.PLOT_LIMIT,len(gene_list)))
        return gene_list[:min(self.PLOT_LIMIT,len(gene_list))]

    def build_annotation_text(self,annos,gene_id):
        '''Return string with formatted annotation data for given gene'''
        text          = self.GENE_ID_LBL + ': ' + gene_id + '<br>'

        if any(key == gene_id for key in annos):

            for key,value in annos[gene_id].items():

                if not value == '':
                    if len(value) > 50:
                        value = value[0:47] + '...'
                    text += key              + ': ' + value   + '<br>'
        else:
            text += '(no data)<br>'

        return text

    def draw_line_plot(self,experiement,gene_list):
        '''Expression line plots for gene list generated from sample filtration steps'''
        self.line_plot.layout.title = self.LINE_UPDATE_TITLE
        gene_list                   = self.limit_num_genes(gene_list)
        avg_exp, _, _               = self.ctrl.model.get_expression(gene_list)

        for i,gene_id in enumerate(gene_list):
            x,y = [],[]

            # Plot average expression level for each condition
            for condition in avg_exp[gene_id][experiement].keys():
                x.append(condition)
                y.append(float(avg_exp[gene_id][experiement][condition][0])) # average expression

            self.line_plot.data[i].x     = x
            self.line_plot.data[i].y     = y
            self.line_plot.data[i].name  = gene_id

        self.line_plot.layout.title = self.LINE_PREFIX_TITLE + experiement

    def clear_plots(self,heatmap_out_widget,differential_out_widget):
        '''Erase all data from all plot displays'''

        # Clear line plot
        for i in range(self.PLOT_LIMIT):
            self.line_plot.data[i].x    = []
            self.line_plot.data[i].y    = []
            self.line_plot.data[i].name = ''

        # Clear heatmap
        with heatmap_out_widget:
            clear_output(wait=True)

        # Clear differential
        with differential_out_widget:
            clear_output(wait=True)

        # clear network plot
        edge_trace = self.net_plot.data[0]
        node_trace = self.net_plot.data[1]
        edge_trace.x    = []
        edge_trace.y    = []
        edge_trace.z    = []
        node_trace.x    = []
        node_trace.y    = []
        node_trace.z    = []
        node_trace.text = []

    def draw_heatmap_plot(self,experiment,gene_list,out_widget):
        '''Create new plot widget, fill with new line plot, append to output'''
        # Dendrograms are determined by distance matrix calculations, color intensity of cells are z-score normalized gene expression matrix data
        self.out_plot_msg(out_widget,self.HEAT_UPDATE_TITLE) # NOTE Also clears plot

        gene_list                    = self.limit_num_genes(gene_list)
        _ ,sample_expression,samples = self.ctrl.model.get_expression(gene_list)

        with out_widget:
            experiment_expression = {}

            # Plot average expression level for each condition
            for gene in gene_list:
                experiment_expression[gene] = {}      # Create empty key in dictionary for this gene
                x, y                        = [], []
                index                       = 0

                # If you're wanting to plot ALL conditions/experiments in one heatmap:
                if experiment == self.ctrl.view.ALL:

                    for sample in samples.keys():
                        experiment_expression[gene][sample] = sample_expression[gene][sample]
                else:
                    for sample in sample_expression[gene].keys():

                        if samples[sample][1] == experiment:
                            experiment_expression[gene][sample] = sample_expression[gene][sample]

            df             = pd.DataFrame(experiment_expression)                         # Store expression data in pandas dataframe for plotting
            df[df.columns] = df[df.columns].astype(float)                                # Make all columns have float data so they can be plotted
            df             = df.transpose()                                              # Transpose for viewing
            df_norm_col    = (df-df.mean())/df.std()                                     # Normalize by column by zscore in dataframe
            df             = df_norm_col.dropna(axis='columns')                          # Remove any 'na' that may arise due to z-score standardizations
            sns_plot       = sns.clustermap(df,figsize=self.HEAT_FIG_SIZE,xticklabels=1) # Draw clustered heat map - Zscore (z_score=0/1) normalize cols=1, or rows=0

            # Plotting layout
            plt.setp(sns_plot.ax_heatmap.yaxis.get_majorticklabels()                      ,rotation=0 ,fontsize=self.HEAT_FONT_SIZE)
            sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(),rotation=90,fontsize=self.HEAT_FONT_SIZE)
            sns_plot.ax_heatmap.set_yticklabels(sns_plot.ax_heatmap.get_ymajorticklabels(),rotation=0 ,fontsize=self.HEAT_FONT_SIZE)

            clear_output(wait=True)
            print(self.HEAT_PREFIX_TITLE + experiment)
            plt.show()

    def draw_differential_plot(self,genes,results,conds,out_widget):
        '''Create new plot widget, fill with new diff. exp. plot, append to output'''

        try:
            self.out_plot_msg(out_widget,self.DIF_UPDATE_TITLE) # NOTE Also clears plot

            with out_widget:
                clear_output(wait=True)

                for gene in genes:
                    print(self.DIF_PREFIX_TITLE+' '+gene)

                    fchanges    = []
                    palette     = {}

                    # Prep: For this gene 1) build list of fold change values 2) Build color palette
                    for i in range(len(conds)):
                        cond        = conds[i]
                        pvalue      = float(results[gene][i][1])
                        fchange     = float(results[gene][i][3])

                        fchanges.append(fchange)

                        # Color palette based on fold changes (up vs down) and significance (p-values)
                        if   pvalue < 0.01 and fchange < 0: palette[cond] = '#0072b2' # Sig down
                        elif pvalue < 0.05 and fchange < 0: palette[cond] = '#56b4e9' # "
                        elif pvalue < 0.01 and fchange > 0: palette[cond] = '#f77f00' # Sig up
                        elif pvalue < 0.05 and fchange > 0: palette[cond] = '#fcbf49' # "
                        else:                               palette[cond] = 'grey'    # Non sig

                    # Create plot for this gene

                    fig,ax = plt.subplots(figsize=(12,8))
                    g      = sns.barplot(x=fchanges,y=conds,palette=palette)  # horizontal bar chart

                    # Set xaxis limits based on abs value of max fold change for gene
                    max_fchange = max(fchanges)
                    #g.set(xlim=(-1*max_fchange-(0.35*max_fchange),max_fchange+(0.35*max_fchange)))
                    g.set(xlim=(-1*max_fchange-1,max_fchange+1))

                    # Axis labels
                    g.set_xlabel('Fold Change',fontsize=14,fontweight='bold')
                    g.set_ylabel('Condition'  ,fontsize=14,fontweight='bold')

                    # For legend, shrink current axis's height by 10% on the bottom
                    box = ax.get_position()
                    ax.set_position([box.x0,box.y0+box.height*0.1,box.width,box.height*0.9])
                    
                    # Add vertical line along x=0
                    ax.axvline(x=0, color='black')
                    
                    # Create color samples for legend
                    color_lines = [ Line2D([0],[0],color='#0072b2' ,lw=4),
                                    Line2D([0],[0],color='#56b4e9',lw=4),
                                    Line2D([0],[0],color='#f77f00',lw=4),
                                    Line2D([0],[0],color='#fcbf49',lw=4),
                                    Line2D([0],[0],color='gray'   ,lw=4)]
            
                    # Create legend, place below plot
                    ax.legend(color_lines,
                            [    'Highly underexpressed (p < 0.01)'
                                ,'Underexpressed (p < 0.05)'
                                ,'Highly overexpressed (p < 0.01)'
                                ,'Overexpressed (p < 0.05)'
                                ,'No DE (p >= 0.05)']
                            ,loc            = 'upper center'
                            ,bbox_to_anchor = (0.5, -0.1)
                            ,fancybox       = True
                            ,shadow         = True
                            ,ncol           = 3
                    )

                    # Plot
                    plt.tight_layout()
                    plt.show()
        except:
            self.ctrl.debug('draw_differential_plot(): EXCEPTION')
            raise


    def draw_network_plot(self,module_id,abc_path,gene_list):
        '''Draw coexpression networks for modules that contain genes of interest & genes with desired expression patterns'''
        # module_id: "N010M00910"
        # abc_path : "./Networks/Network_10/N010M00910.abc"
        # gene_list: "['Sevir.5G108900',....]"    NOTE These are only the _recovered_ genes

        self.ctrl.debug('Net plot: ' + module_id+', '+abc_path+', and '+str(gene_list))

        self.net_plot.layout.title = self.NET_UPDATE_TITLE

        # Get access to traces within widget data NOTE Changing node/edge_trace sends updates to plot displayed in browser
        edge_trace = self.net_plot.data[0]
        node_trace = self.net_plot.data[1]

        # Clear plot of any previous data
        edge_trace.x    = []
        edge_trace.y    = []
        edge_trace.z    = []
        node_trace.x    = []
        node_trace.y    = []
        node_trace.z    = []
        node_trace.text = []

        # Create networkX graph
        G   = nx.read_weighted_edgelist(abc_path)
        pos = nx.spring_layout(G,dim=3) # This is a randomized layout of the network

        # Edge traces
        for edge in G.edges():
            x0,y0,z0         = pos[edge[0]]
            x1,y1,z1         = pos[edge[1]]
            edge_trace['x'] += tuple([x0,x1,None])
            edge_trace['y'] += tuple([y0,y1,None])
            edge_trace['z'] += tuple([z0,z1,None])

        # Add xy positions of each node/gene
        for node in G.nodes():
            x,y,z            = pos[node]
            node_trace['x'] += tuple([x])
            node_trace['y'] += tuple([y])
            node_trace['z'] += tuple([z])

        # Highlight selected genes

        colors = []
        annos  = self.ctrl.model.get_annos(G.nodes) # Data for all genes

        for gene_id in G.nodes():

            if gene_id in gene_list:
                colors.append(self.NET_GENE_COLOR_EMPH)
            else:
                colors.append(self.NET_GENE_COLOR_NORM)

            node_trace['text'] += tuple([self.build_annotation_text(annos,gene_id)])

        node_trace['marker']['color'] = colors
        self.net_plot.layout.title = self.NET_PREFIX_TITLE + module_id

