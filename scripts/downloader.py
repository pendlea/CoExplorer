# downloaded.py - Custom download button for coexp notebook
# rcampbel@purdue.edu - January 2020

import os
import ipywidgets

class Downloader:

    def __init__(self,text,width=None,path=None):
        '''Create link widget that looks like ipywidgets button'''

        if width:
            self.width = width
        else:
            self.width = '25%'

        self.text   = text
        self.widget = ipywidgets.HTML() # <- NOTE: Outside code must access ".widget" (e.g. for display)
        self.update(path) # Sets .path

    def update(self,new_path=None):
        '''Fill content of widget'''

        if not new_path:
            new_path = ''

        self.path = new_path

        if self.path == '':
            status =  'pointer-events: none;' # Disabled
            color  =  'gray'
        else:
            status =  ''                      # Endabled
            color  =  'black'

        # Replace HTML contents  # TODO Find way to make Firefox download dialog appear (or open in new tab)
        self.widget.value =  '''
            <a  href     = "'''+self.path+'''"
                download = "'''+os.path.basename(self.path)+'''"
                target   = "_blank"
                style    =
                    " width          : '''+self.width+''';
                    color            : '''+color+''';
                    background-color : #EEEEEE;
                    text-decoration  : none;
                    '''+status+'''
                    padding-left     : 10px;
                    padding-right    : 10px;
                    padding-top      : 5px;
                    padding-bottom   : 5px;">
                    <i class="fa fa-download"></i>&nbsp'''+self.text+'</a>'

