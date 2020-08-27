# controller.py - Central logic for coexp notebook
# rcampbel@purdue.edu - November 2019

import os

# Avoids warning: "numpy.dtype size changed, may indicate binary incompatibility"
import warnings
warnings.filterwarnings('ignore')

from scripts import plotter

class Controller:

    VALUE = 'value' # for observe calls

    def __init__(self):
        self.debugging     = False # NOTE Change to False to hide debug output
        self.debug_buffer  = []
        self.display_ready = False

    def intro(self,model,view):
        self.model = model
        self.view  = view

    def start(self):
        self.plotter = plotter.Plotter(self) # Send self so plotter can call debug()
        self.view.display(debug=self.debugging)
        self.display_ready = True
        self.observe()

    def debug(self,text):

        if self.debugging:
            self.debug_buffer.append(text)

            if self.display_ready:

                for line in self.debug_buffer:
                    self.view.debug('DEBUG: '+line)

                self.debug_buffer = []

    def observe(self):
        '''Connect widgets to controller callbacks'''

        self.view.filter_btn_apply.on_click(self.apply_filter)
        self.view.filter_ddn_ndisp.observe(self.ndisp_changed,self.VALUE)

        self.view.plotex_ddn_selex_lg.observe(self.linegraph_experiment_selected,self.VALUE)
        self.view.plotex_ddn_selex_hm.observe(self.heatmap_experiment_selected  ,self.VALUE)

        self.view.plotco_ddn_netw.observe(self.network_selected,self.VALUE)
        self.view.plotco_sel_modu.observe(self.module_selected ,self.VALUE)

        # Only monitor the min expressed filter widgets - to change diff expressed widgets as needed
        for item in self.view.filter_conditons:

            for button in item:
                button.on_click(self.three_state_pressed)

        self.view.filter_btn_refexp.on_click(self.fill_results_export)
        self.view.plotco_btn_modu.on_click(self.fill_module_export)

        self.view.plotdif_btn_plot.on_click(self.diff_plot)

    def fill_results_export(self,change):
        # Generate output file
        self.view.filter_out_export.clear_output()

        if self.model.filter_results:
            self.view.output_data_link(self.view.filter_out_export,self.model.write_filtered_data())

    def fill_module_export(self,change):
        # Generate output file
        self.view.plotco_out_export.clear_output()

        if self.model.filter_results:
            self.view.output_data_link(self.view.plotco_out_export,self.model.module_download_data)

    def three_state_pressed(self,button):
        '''React to user pressing a 3-state button in fitler'''

        # Change state of button that was pressed

        if button.icon   == self.view.FILT_CLR:
            button.icon   = self.view.FILT_POS

        elif button.icon == self.view.FILT_POS:
            button.icon   = self.view.FILT_NEG

        elif button.icon == self.view.FILT_NEG:
            button.icon   = self.view.FILT_CLR

        # Adust diff exp buttons
        for item in self.view.filter_conditons:

            # Which min exp button was pressed (i.e. on which row)?
            if button == item[0]:

                # Set corresponding buttons accordingly
                if button.icon  == self.view.FILT_NEG:
                    item[1].icon = self.view.FILT_CLR
                    item[2].icon = self.view.FILT_CLR
                    item[1].disabled = True
                    item[2].disabled = True
                else:
                    item[1].disabled = False
                    item[2].disabled = False

                break

    def set_image(self,widget,img_data):
        '''Update image widget with new image'''

        # Expecting image_data = (image_bytes,(width,height))
        widget.width  = img_data[1][0]
        widget.height = img_data[1][1]
        widget.value  = img_data[0]     # Triggers update to disp img

    def parse(self,text):
        '''Isolate tokens from within CSV text field'''
        text   = text.replace(';',',')        # Support semicolon as delimiter by converting to commas
        tokens = text.split(',')              # Separate out individual tokens
        tokens = list(map(str.strip,tokens))  # Create list of stripped tokens
        return list(filter(None,tokens))      # Remove empty strings

    def apply_filter(self,change):
        '''React to apply filter button press'''

        self.view.filter_html_output.value = self.view.FILTER_PROG

        # Clear export widgets
        self.view.filter_out_export.clear_output()
        self.view.plotco_out_export.clear_output()

        # Reset some plot-control widgets (others done in set_plot_status())
        self.view.plotex_ddn_selex_lg.value = self.view.EMPTY
        self.view.plotex_ddn_selex_hm.value = self.view.EMPTY
        self.view.plotco_ddn_netw.value     = self.view.EMPTY
        self.set_module_data([''],[(None,None)])
        self.model.clear_module_download()

        # Get IDs from UI
        gene_ids   = self.parse(self.view.filter_txt_gene.value)
        func_ids   = self.parse(self.view.filter_txt_func.value)
        self.debug('Preped gene IDs:'+str(gene_ids))
        self.debug('Preped func IDs:'+str(func_ids))

        # Translate to native.
        target_ids = self.model.translate_genes(gene_ids,func_ids)
        self.debug('Searching for '+str(len(target_ids))+' IDs; first few: '+str(target_ids[:10]))

        perform_search = True

        if not target_ids:

            if len(gene_ids) > 0 or len(func_ids) > 0: # Did user specify gene(s) or funtion(s)?
                perform_search = False
                self.debug('Translation failed.')
            else:
                self.view.filter_html_output.value = self.view.FILTER_PROG_ALL
                self.debug('WARNING: Considering ALL genes.')

        self.model.clear_filter_results() # New search attempt so reset

        if perform_search: # Either terms were left empty (user wants full results) or at least one Sevir ID is available

            # Get thresholds from from UI
            tpm_thresh  = float(self.view.filter_txt_tpm.value )
            pval_thresh = float(self.view.filter_txt_pval.value)
            fdr_thresh  = float(self.view.filter_txt_fdr.value )

            # Search for valid data, (results stored in model)
            self.model.search(target_ids,tpm_thresh,pval_thresh,fdr_thresh)

            # Get annotation data (stored in model) for search results
            self.model.add_annos()

        # Refresh output widgets
        self.refresh_filter_output()

    def ndisp_changed(self,change):
        self.refresh_filter_output()

    def refresh_filter_output(self):
        # Enable or disable controls based on filter results
        if self.model.filter_results:
            self.view.update_filtered_gene_list()    # Update output table in filter tab
            self.view.set_plot_status(enable=True)
        else:
            self.view.filter_html_output.value = self.view.EMPTY_LIST_MSG
            self.view.filter_btn_downd.update()  # Disable download
            self.view.set_plot_status(enable=False)

    def linegraph_experiment_selected(self,change):
        '''React to experiment selection for line graph'''
        experiment = change['owner'].value

        if experiment != self.view.EMPTY:
            gene_list   = list(self.model.filter_results.keys())
            self.plotter.draw_line_plot(experiment,gene_list)

    def heatmap_experiment_selected(self,change):
        '''React to experiment selection for heatmap'''
        experiment = change['owner'].value

        if experiment != self.view.EMPTY:
            gene_list = list(self.model.filter_results.keys())

            if len(gene_list) > 1:
                self.plotter.draw_heatmap_plot(experiment,gene_list,self.view.plotex_img_dispp_hm)
            else:
                self.plotter.out_plot_msg(self.view.plotex_img_dispp_hm,self.plotter.HEAT_LESS_THAN_TWO)

    def network_selected(self,change):
        '''React to network selection'''
        network = change['owner'].value
        self.view.plotco_out_export.clear_output()

        if network != self.view.EMPTY:
            self.set_module_data([self.view.MODULE_PROG],[(None,None)]) # Tell user working on it (may take a while)

            # Any filter results available?
            if self.model.filter_results:
                # Get module data and format it
                display,values = self.model.get_module_data(network,self.view.plotco_out_export)
                display        = self.view.columnize(self.view.MODULE_HEADER[0],display)
            else:
                display = [self.view.NO_MODULE_DATA]
                values  = [(None,None)]
                self.ctrl.debug('No module data: "'+str(module_data)+'"')

            self.set_module_data(display,values)

    def set_module_data(self,display,values):
        '''Suspend callback, update module widget w/new data, reinstate callback'''
        self.view.plotco_sel_modu.unobserve(self.module_selected,self.VALUE)
        self.view.plotco_sel_modu.options = zip(display,values) # Single list of tuples, each with corresponding items from both lists
        self.view.plotco_sel_modu.value   = None
        self.view.plotco_sel_modu.observe(self.module_selected,self.VALUE)

    def module_selected(self,change):
        '''React to module selection'''
        line,recovered = change['new']
        module         = line.split()[0]
        self.debug('module_selected(): module='+module+', recovered='+str(recovered)+', line='+line)

        # Build path to ABC file: "./Networks/Network_10/N010M00910.abc"
        abc_path = os.path.join(
            self.model.DATA_DIR               # '.'
            ,self.model.NET_DIR               # 'Networks'
            ,self.view.plotco_ddn_netw.value  # 'Network_10'
            ,module + self.model.ABC_EXT      # 'N010M00910' + '.abc'
        )

        self.plotter.draw_network_plot(module,abc_path,recovered)

    def diff_plot(self,change):
        '''React to gene selection for diff exp plotting'''
        try:
            self.plotter.draw_differential_plot(
                self.view.plotdif_sel_genes.value
                ,self.model.filter_results
                ,self.model.cond
                ,self.view.plotdif_img_disp
            )
        except:
            self.ctrl.debug('diff_plot(): EXCEPTION')
            raise
