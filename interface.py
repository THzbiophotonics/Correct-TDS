#!/usr/bin/env python
# coding: utf-8



#!/usr/bin/python
# -*- coding: latin-1 -*-

import os, sys
from Correct_TDS import fitf as TDS
from Correct_TDS.fitc import Controler
import Correct_TDS.constants as csts
import subprocess
import numpy as np
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtGui
from PyQt5 import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import sip
import shutil
import h5py
import time
from scipy import signal
from pathlib import Path as path_
import multiprocessing
from datetime import datetime

plt.rcParams.update({'font.size': 13})

# ========================== progress bar stylesheet ========================= #
StyleSheet = '''
#RedProgressBar {
    text-align: center;
}
#RedProgressBar::chunk {
    background-color: #F44336;
}
#GreenProgressBar {
    min-height: 12px;
    max-height: 12px;
    border-radius: 6px;
}
#GreenProgressBar::chunk {
    border-radius: 6px;
    background-color: #009688;
}
#BlueProgressBar {
    border: 2px solid #2196F3;
    border-radius: 5px;
    background-color: #E0E0E0;
}
#BlueProgressBar::chunk {
    background-color: #2196F3;
    width: 10px; 
    margin: 0.5px;
}
'''


try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    size = comm.Get_size()
except:
    print('mpi4py is required for parallelization')
    myrank=0

ROOT_DIR = path_(__file__).parent

def deleteLayout(layout):
    if layout is not None:
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
            else:
                deleteLayout(item.layout())
        sip.delete(layout)

def to_sup(s):
    """Convert a string of digit to supescript"""
    sups = {u'0': u'\u2070',
            u'1': u'\u00b9',
            u'2': u'\u00b2',
            u'3': u'\u00b3',
            u'4': u'\u2074',
            u'5': u'\u2075',
            u'6': u'\u2076',
            u'7': u'\u2077',
            u'8': u'\u2078',
            u'9': u'\u2079'}

    return ''.join(sups.get(char, char) for char in s)


class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'
   
   

graph_option_2=None
preview = 1
apply_window = 0


class MyTableWidget(QWidget):

    def __init__(self, parent,controler):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)

        # Initialize tab screen
        self.tab = Initialisation_tab(self,controler)

        # Add tabs to widget
        self.layout.addWidget(self.tab)
        self.setLayout(self.layout)


###############################################################################
###############################################################################
#######################   Initialisation tab   ################################
###############################################################################
###############################################################################

class Initialisation_tab(QWidget):
    def __init__(self, parent, controler):
        super().__init__(parent)
        self.controler = controler
        self.setMinimumSize(640, 480)
        hlayout  = QHBoxLayout()
        vlayout = QVBoxLayout()
        
        self.init_param_widget = InitParam_handler(self, controler)
        self.init_param_widget.refresh()

        vlayout.addWidget(self.init_param_widget)
        hlayout.addLayout(vlayout)
        self.setLayout(hlayout)
        

        
class InitParam_handler(QWidget):
    def __init__(self, parent, controler):
        super().__init__(parent)
        self.parent = parent
        self.controler = controler

    def refresh(self):
        init_instance = InitParamWidget(self.parent,self.controler)
        try:
            deleteLayout(self.layout())
            self.layout().deleteLater()
        except:
            main_layout = QVBoxLayout()
            main_layout.addWidget(init_instance)
            self.setLayout(main_layout)


class InitParamWidget(QWidget):
    def __init__(self, parent, controler):
        super().__init__(parent)
        self.parent = parent
        self.controler = controler
        self.controler.addClient(self)
        self.path_data = None
        self.path_data_ref = None
        
        self.dialog = QDialog()
        self.dialog.ui = Ui_Dialog()
        self.dialog.ui.setupUi(self.dialog, controler)
        
        self.dialog_match = QDialog()
        self.dialog_match.ui = Ui_Dialog_match()
        
        # ================= change: add a ref loading bool to control ================ #
        self.ref_loaded = False
        # ---------------------------------------------------------------------------- #
        label_width=1500
        text_box_width=150
        text_box_height=22
        
        # We create the text associated to the text box

        self.label_path_data = QLabel('Select data (hdf5) \u00b9')
        self.label_path_data.setAlignment(Qt.AlignVCenter)
        self.label_path_data.resize(200, 100)
        self.label_path_data.resize(self.label_path_data.sizeHint())
        self.label_path_data.setMaximumHeight(text_box_height)
        
        self.button_ask_path_data = QPushButton('browse')
        self.button_ask_path_data.resize(200, 100)
        self.button_ask_path_data.resize(self.button_ask_path_data.sizeHint())
        self.button_ask_path_data.setMaximumHeight(text_box_height)
        self.button_ask_path_data.clicked.connect(self.get_path_data)
        
        
        self.label_path_without_sample = QLabel('Select data without sample (optional) \u00b2')
        self.label_path_without_sample.setAlignment(Qt.AlignVCenter)
        self.label_path_without_sample.resize(200, 100)
        self.label_path_without_sample.resize(self.label_path_without_sample.sizeHint())
        self.label_path_without_sample.setMaximumHeight(text_box_height)
        
        
        self.button_ask_path_without_sample = QPushButton('browse')
        self.button_ask_path_without_sample.resize(200, 100)
        self.button_ask_path_without_sample.resize(self.button_ask_path_without_sample.sizeHint())
        self.button_ask_path_without_sample.clicked.connect(self.get_path_data_ref)
        
        self.label_match_sample_ref = QLabel('Match sample and reference files')
        self.label_match_sample_ref.setAlignment(Qt.AlignVCenter)
        self.label_match_sample_ref.resize(200, 100)
        self.label_match_sample_ref.resize(self.label_match_sample_ref.sizeHint())
        self.label_match_sample_ref.setMaximumHeight(text_box_height)
        
        self.button_match_sample_ref = QPushButton('match')
        self.button_match_sample_ref.resize(200, 100)
        self.button_match_sample_ref.resize(self.button_match_sample_ref.sizeHint())
        self.button_match_sample_ref.setMaximumHeight(text_box_height)
        self.button_match_sample_ref.clicked.connect(self.open_dialog_match)
        
        
        self.label_data_length = QLabel('Set part of data to analyze? (optional)')
        self.label_data_length.setAlignment(Qt.AlignVCenter)
        self.label_data_length.resize(200, 100)
        self.label_data_length.resize(self.label_path_data.sizeHint())
        self.label_data_length.setMaximumHeight(text_box_height)
        
        self.button_ask_data_length = QPushButton('set') #unset after
        self.button_ask_data_length.resize(200, 100)
        self.button_ask_data_length.resize(self.button_ask_path_data.sizeHint())
        self.button_ask_data_length.setMaximumHeight(text_box_height)
        self.button_ask_data_length.clicked.connect(self.open_dialog)
        

        self.button_data = QPushButton('Submit data')

        self.button_parameters = QPushButton('Submit parameters')
        self.button_parameters.setMaximumHeight(22)
        #TODO
        # self.button_parameters.clicked.connect(self.on_click_param)
        self.button_parameters.pressed.connect(self.pressed_loading)
        self.button_parameters.setEnabled(False)

        # We create a button to extract the information from the text boxes
        self.button = QPushButton('Submit / Preview')
        self.button.pressed.connect(self.pressed_loading)
        self.button.setEnabled(False)
        #TODO
        #TOVERIFY
        # self.button.clicked.connect(self.on_click)
        # self.button.clicked.connect(self.on_click_print)
        self.button.setMaximumHeight(text_box_height)
        
        # Filter or not filter
        self.LFfilter_label = QLabel('Filter low frequencies?\t      ')
        self.LFfilter_choice = QComboBox()
        self.LFfilter_choice.addItems(['No','Yes'])
        self.LFfilter_choice.setMaximumWidth(text_box_width)
        self.LFfilter_choice.setMaximumHeight(text_box_height)
        self.LFfilter_choice.setCurrentIndex(1)
        
        
        self.HFfilter_label = QLabel('Filter high frequencies?\t      ')
        self.HFfilter_choice = QComboBox()
        self.HFfilter_choice.addItems(['No','Yes'])
        self.HFfilter_choice.setMaximumWidth(text_box_width)
        self.HFfilter_choice.setMaximumHeight(text_box_height)
        
        self.label_start = QLabel('\tStart (Hz)')
        self.label_end   = QLabel('\tEnd (Hz)')
        self.label_sharp = QLabel('Sharpness of frequency filter \u00b3')
        self.start_box = QLineEdit()
        self.end_box   = QLineEdit()
        self.sharp_box = QLineEdit()
        self.start_box.setMaximumWidth(text_box_width)
        self.start_box.setMaximumHeight(text_box_height)
        self.start_box.setText("0.18e12")
        self.end_box.setMaximumWidth(text_box_width)
        self.end_box.setMaximumHeight(text_box_height)
        self.end_box.setText("6e12")
        self.sharp_box.setMaximumWidth(text_box_width-85)
        self.sharp_box.setMaximumHeight(text_box_height)
        self.sharp_box.setText("2")
        
        # Super resolution
        self.label_super = QLabel("            Super resolution ")
        self.options_super = QComboBox()
        self.options_super.addItems(['No','Yes'])
        self.options_super.setMinimumWidth(text_box_width-75)
        self.options_super.setMaximumWidth(text_box_width)
        self.options_super.setMaximumHeight(text_box_height)
        self.options_super.setCurrentIndex(0)
    
        
        # Delay
        self.label_delay = QLabel("Correct delay?")
        self.label_delay.setMaximumWidth(label_width)
        self.label_delay.setMaximumHeight(text_box_height)
        self.options_delay = QComboBox()
        self.options_delay.addItems(['No','Yes'])
        self.options_delay.setMaximumWidth(text_box_width-75)
        self.options_delay.setMaximumHeight(text_box_height)
        self.options_delay.setCurrentIndex(1)
        self.delayvalue_label = QLabel("Delay absolute value")
        self.delay_limit_box = QLineEdit()
        self.delay_limit_box.setMaximumWidth(text_box_width-24)
        self.delay_limit_box.setMaximumHeight(text_box_height)
        self.delay_limit_box.setText("10e-12")
        
        # Leftover noise
        self.label_leftover = QLabel("Correct amplitude noise?")
        self.label_leftover.setMaximumWidth(label_width)
        self.label_leftover.setMaximumHeight(text_box_height)
        self.options_leftover = QComboBox()
        self.options_leftover.addItems(['No','Yes'])
        self.options_leftover.setMaximumWidth(text_box_width-75)
        self.options_leftover.setMaximumHeight(text_box_height)
        self.options_leftover.setCurrentIndex(1)
        self.leftovervaluea_label = QLabel("Absolute value of         a")
        self.leftovera_limit_box = QLineEdit()
        self.leftovera_limit_box.setMaximumWidth(text_box_width-24)
        self.leftovera_limit_box.setMaximumHeight(text_box_height)
        self.leftovera_limit_box.setText("10e-2")
        #self.leftovervaluec_label = QLabel(" c")
        #self.leftoverc_limit_box = QLineEdit()
        #self.leftoverc_limit_box.setMaximumWidth(text_box_width-24)  
        # Dilatation
        self.label_dilatation = QLabel("Correct dilatation?")
        self.label_dilatation.setMaximumWidth(label_width)
        self.label_dilatation.setMaximumHeight(text_box_height)
        self.options_dilatation = QComboBox()
        self.options_dilatation.addItems(['No','Yes'])
        self.options_dilatation.setMaximumWidth(text_box_width-75)
        self.options_dilatation.setMaximumHeight(text_box_height)
        self.dilatationvaluea_label = QLabel("Absolute value of         \u03B1    ")
        self.dilatationa_limit_box = QLineEdit()
        self.dilatationa_limit_box.setMaximumWidth(text_box_width-24)
        self.dilatationa_limit_box.setMaximumHeight(text_box_height)
        self.dilatationa_limit_box.setText("10e-3")
        self.dilatationvalueb_label = QLabel("b")
        self.dilatationb_limit_box = QLineEdit()
        self.dilatationb_limit_box.setMaximumWidth(text_box_width-24) 
        self.dilatationb_limit_box.setMaximumHeight(text_box_height)
        self.dilatationb_limit_box.setText("10e-12")
        
        #periodic sampling
        self.label_periodic_sampling = QLabel("Correct periodic sampling?")
        self.label_periodic_sampling.setMaximumWidth(label_width)
        self.label_periodic_sampling.setMaximumHeight(text_box_height)
        self.options_periodic_sampling = QComboBox()
        self.options_periodic_sampling.addItems(['No','Yes'])
        self.options_periodic_sampling.setMaximumWidth(text_box_width-75)
        self.options_periodic_sampling.setMaximumHeight(text_box_height)
        self.options_periodic_sampling.setCurrentIndex(1)
        self.periodic_sampling_freq_label = QLabel("Frequency [THz]")
        self.periodic_sampling_freq_limit_box = QLineEdit()
        self.periodic_sampling_freq_limit_box.setText("7.5")
        self.periodic_sampling_freq_limit_box.setMaximumWidth(text_box_width-24)
        self.periodic_sampling_freq_limit_box.setMaximumHeight(text_box_height)
        


        # Organisation layout
        self.hlayout6=QHBoxLayout()
        self.hlayout7=QHBoxLayout()
        #TOVERIFY:
        self.hlayout99=QHBoxLayout()
        
        self.hlayout8=QHBoxLayout()
        self.hlayout9=QHBoxLayout()
        self.hlayout10=QHBoxLayout()
        self.hlayout11=QHBoxLayout()
        self.hlayout12=QHBoxLayout()
        self.hlayout13=QHBoxLayout()
        self.hlayout14=QHBoxLayout()
        self.hlayout17=QHBoxLayout()
        self.hlayout18=QHBoxLayout()
        self.hlayout19=QHBoxLayout()
        self.hlayout20=QHBoxLayout()
        self.hlayout21=QHBoxLayout()
        self.hlayout22=QHBoxLayout()
        self.hlayout23=QHBoxLayout()
        self.vlayoutmain=QVBoxLayout()
        
        
        self.hlayout6.addWidget(self.label_path_data,20)
        self.hlayout6.addWidget(self.button_ask_path_data,17)


        self.hlayout7.addWidget(self.label_path_without_sample,20)
        self.hlayout7.addWidget(self.button_ask_path_without_sample,17)
        
        self.hlayout99.addWidget(self.label_match_sample_ref,20)
        self.hlayout99.addWidget(self.button_match_sample_ref,17)
        
        self.hlayout11.addWidget(self.label_data_length,20)
        self.hlayout11.addWidget(self.button_ask_data_length,17)
        
        self.hlayout8.addWidget(self.button_data)
        
        self.hlayout9.addWidget(self.label_delay,0)
        self.hlayout9.addWidget(self.options_delay,1)

        self.hlayout9.addWidget(self.delayvalue_label,0)
        self.hlayout9.addWidget(self.delay_limit_box,1)
        
        self.hlayout12.addWidget(self.label_leftover,0)
        self.hlayout12.addWidget(self.options_leftover,0)
        self.hlayout12.addWidget(self.leftovervaluea_label,0)
        self.hlayout12.addWidget(self.leftovera_limit_box,1)
        #self.hlayout14.addWidget(self.leftovervaluec_label,0)
        #self.hlayout14.addWidget(self.leftoverc_limit_box,1)
        
        self.hlayout13.addWidget(self.label_dilatation,1)
        self.hlayout13.addWidget(self.options_dilatation,1)
        self.hlayout13.addWidget(self.dilatationvaluea_label,0)
        self.hlayout13.addWidget(self.dilatationa_limit_box,0)
        #self.hlayout13.addWidget(self.dilatationvalueb_label,0)
        #self.hlayout13.addWidget(self.dilatationb_limit_box,0)
        
        self.hlayout14.addWidget(self.label_periodic_sampling,0)
        self.hlayout14.addWidget(self.options_periodic_sampling,0)
        self.hlayout14.addWidget(self.periodic_sampling_freq_label,0)
        self.hlayout14.addWidget(self.periodic_sampling_freq_limit_box,1)
        
        #TOADD:
        # self.hlayout10.addWidget(self.button_parameters)
        
        self.hlayout17.addWidget(self.LFfilter_label,1)
        self.hlayout17.addWidget(self.LFfilter_choice,0)
        self.hlayout17.addWidget(self.label_start,1)
        self.hlayout17.addWidget(self.start_box,0)
        
        self.hlayout18.addWidget(self.HFfilter_label,1)
        self.hlayout18.addWidget(self.HFfilter_choice,0)
        self.hlayout18.addWidget(self.label_end,1)
        self.hlayout18.addWidget(self.end_box,0)

        self.hlayout19.addWidget(self.label_sharp,1)
        self.hlayout19.addWidget(self.sharp_box,0)
        self.hlayout19.addWidget(self.label_super,1)
        self.hlayout19.addWidget(self.options_super,0)

        
        sub_layoutv2 = QVBoxLayout()
        sub_layoutv3 = QVBoxLayout()
        
        sub_layoutv2.addLayout(self.hlayout6)
        sub_layoutv2.addLayout(self.hlayout7)
        sub_layoutv2.addLayout(self.hlayout99)
        sub_layoutv2.addLayout(self.hlayout11)
        
        sub_layoutv3.addLayout(self.hlayout9)
        sub_layoutv3.addLayout(self.hlayout13)
        sub_layoutv3.addLayout(self.hlayout12)
        sub_layoutv3.addLayout(self.hlayout14)
        sub_layoutv3.addLayout(self.hlayout10)

        # sub_layoutv2.addLayout(self.hlayout23)
        init_group = QGroupBox()
        init_group.setTitle('Data Initialization')
        init_group.setLayout(sub_layoutv2)
        

        self.vlayoutmain.addWidget(init_group)
        
        init_group2 = QGroupBox()
        init_group2.setTitle('Correction parameters')
        init_group2.setLayout(sub_layoutv3)
        
        sub_layoutv = QVBoxLayout()
        #sub_layoutv.addLayout(self.hlayout21)
        sub_layoutv.addLayout(self.hlayout22)

        sub_layoutv.addLayout(self.hlayout17)
        sub_layoutv.addLayout(self.hlayout18)
        sub_layoutv.addLayout(self.hlayout19)
        #sub_layoutv.addLayout(self.hlayout20)
        filter_group = QGroupBox()
        filter_group.setTitle('Filters')
        filter_group.setLayout(sub_layoutv)
        self.vlayoutmain.addWidget(filter_group)
        #TOADD
        # self.vlayoutmain.addWidget(self.button)
        self.vlayoutmain.addWidget(init_group2)

        
        self.action_widget = Action_handler(self,controler)
        self.action_widget.refresh()
        self.save_param = Saving_parameters(self,controler)
        #self.save_param.refresh()
        self.graph_widget = Graphs_optimisation(self,controler)
        
        
        self.vlayoutmain.addWidget(self.action_widget)
        self.vlayoutmain.addWidget(self.save_param)
        
        self.text_box = TextBoxWidget(self, controler)
        self.vlayoutmain.addWidget(self.text_box,0)
        
        main_layout = QHBoxLayout()
        main_layout.addLayout(self.vlayoutmain,0)
        main_layout.addWidget(self.graph_widget,1)
        self.setLayout(main_layout)

    def pressed_loading1(self):
        self.controler.loading_text()
        
    def open_dialog(self):
        self.dialog.exec_()
    
    def open_dialog_match(self):
        if csts.refs:
            self.dialog_match.ui.setupUi(self.dialog_match, self.controler)
        # self.button_match_sample_ref.setEnabled(False)
            self.dialog_match.exec_()
        else:
            self.dialog_match.ui.show_info_noref_messagebox()
    
    # def on_click_print(self):
    #     print(self.LFfilter_choice.currentIndex())
    #     print(self.HFfilter_choice.currentIndex())
    #     print(self.start_box.text())
    #     print(self.end_box.text())
    #     print(self.sharp_box.text())
    #     print(self.options_super.currentIndex())
        

    def on_click(self,data_,ref_):
        # ============================================================================ #
        #                   change funtion to remove preview of graph                  #
        # ============================================================================ #
        global graph_option_2, preview
        try:
            Lfiltering_index = self.LFfilter_choice.currentIndex()
            Hfiltering_index = self.HFfilter_choice.currentIndex()
            cutstart = float(self.start_box.text())
            cutend   = float(self.end_box.text())
            cutsharp = float(self.sharp_box.text())
            modesuper = self.options_super.currentIndex()
            # print(f" super = {modesuper}")
            
            trace_start = 0
            trace_end = -1
            time_start = 0
            time_end = -1
            self.button_ask_data_length.setText("set")
            
            if self.dialog.ui.length_initialized:
                trace_start = self.dialog.ui.trace_start
                trace_end = self.dialog.ui.trace_end
                time_start = self.dialog.ui.time_start
                time_end = self.dialog.ui.time_end
                self.button_ask_data_length.setText(str(trace_start)+"-"+str(trace_end))

            try:
                # self.controler.choices_ini(self.path_data, self.path_data_ref, trace_start, trace_end, time_start, time_end,
                # ============================================================================ #
                #            change: modify function to select sample and ref files            #
                # ============================================================================ #
                # self.controler.choices_ini(data_, self.path_data_ref, trace_start, trace_end, time_start, time_end,Lfiltering_index, Hfiltering_index, cutstart,cutend, cutsharp, modesuper, apply_window)
                csts.modesuper = modesuper
                self.controler.choices_ini(data_, ref_, trace_start, trace_end, time_start, time_end,Lfiltering_index, Hfiltering_index, cutstart,cutend, cutsharp, modesuper, apply_window)
                # ============================= change : comment ============================= #
                # graph_option_2='Pulse (E_field)'
                preview = 1
                # ============================= change : comment ============================= #
                # self.graph_widget.refresh()
                # ================ change: add text to mention number of files =============== #
                self.controler.refreshAll3(" Data initialization done | "+str(self.controler.data.numberOfTrace)+ " time traces loaded between ["+ str(int(self.controler.data.time[0])) + " , " + str(int(self.controler.data.time[-1]))+ "] ps" + "; for " + str(len(csts.files)) + " files")
            except Exception as e:
                # ---------------------------------------------------------------------------- #
                print(e)
                self.controler.error_message_path3()
                return(0)
        except Exception as e:
            print(e)
            self.controler.refreshAll3("Invalid parameters, please enter real values only")
        self.controler.initialised=1
        self.controler.optim_succeed = 0
        print(f"submit files finished")

    def on_click_param(self):
        global preview, graph_option_2
        if not self.controler.initialised:
            self.controler.refreshAll3("Please submit initialization data first")
            return(0)
        try:
            fit_delay = 0
            delay_guess = 0
            delay_limit = 0
            fit_leftover_noise = 0
            leftover_guess = np.zeros(2)
            leftover_limit = np.zeros(2)
            fit_periodic_sampling = 0
            periodic_sampling_freq_limit = 0
            fit_dilatation = 0
            dilatation_limit = np.zeros(2)
            dilatationmax_guess = np.zeros(2)

            if self.options_delay.currentIndex():
                fit_delay = 1
                delay_limit = float(self.delay_limit_box.text())
                if delay_limit <=0:
                    raise ValueError
            if self.options_leftover.currentIndex():
                fit_leftover_noise = 1
                leftover_limit[0] = float(self.leftovera_limit_box.text())
                #leftover_limit[1] = float(self.leftoverc_limit_box.text())
                leftover_limit[1] = 10e-100
                
                if leftover_limit[0]<=0 or leftover_limit[1]<=0:
                    raise ValueError

            if self.options_dilatation.currentIndex():
                fit_dilatation = 1
                dilatation_limit[0] = float(self.dilatationa_limit_box.text())
                dilatation_limit[1] = float(self.dilatationb_limit_box.text())
                
                if dilatation_limit[0]<=0 or dilatation_limit[1]<=0:
                    raise ValueError
            if self.options_periodic_sampling.currentIndex():
                fit_periodic_sampling = 1
                periodic_sampling_freq_limit = float(self.periodic_sampling_freq_limit_box.text())
                if periodic_sampling_freq_limit < 0 or periodic_sampling_freq_limit > self.controler.myglobalparameters.freq[-1]*1e-12:
                    raise ValueError
                    
            try:
                self.controler.choices_ini_param( fit_delay, delay_guess, delay_limit, fit_dilatation, dilatation_limit, dilatationmax_guess,
                                               fit_periodic_sampling, periodic_sampling_freq_limit,
                                              fit_leftover_noise ,leftover_guess, leftover_limit)
                self.controler.refreshAll3(" Parameters initialization done")
                print(f"submit parameters done")
                self.controler.initialised_param = 1
                self.controler.optim_succeed = 0
            except Exception as e:
                print(e)
                return(0)
        except:
            self.controler.refreshAll3("Invalid parameters, please enter real positive values only and valid frequency range")
        # ======================= change: remove graph preview ======================= #
        # try:
            # preview = 1
            # self.graph_widget.refresh()
        # except:
            # print("unknown error")
        # ---------------------------------------------------------------------------- #
            
    def get_path_data(self):
        # # options = QFileDialog.Options()
        # # options |= QFileDialog.DontUseNativeDialog
        # # fileName, _ = QFileDialog.getOpenFileName(self,"Data (hdf5 file)", options=options, filter="Hdf5 File(*.h5)")
        # =========================== initialize csts.files ========================== #
        csts.files = []
        # ============================================================================ #
        #                     change function to get multiple files                    #
        #    store the files in as a list in another python file called constants.py   #
        # ============================================================================ #
        # ===================== change : change to get filenames ===================== #
        # fileName, _ = QFileDialog.getOpenFileName(self,"Data (hdf5 file)", filter="Hdf5 File(*.h5)")
        fileNames, _ = QFileDialog.getOpenFileNames(parent=self,caption=f"Select multiple h5 files",directory=f"{csts.init_directory}", filter="Hdf5 File(*.h5)")
        # ================== change: fill files list with filenames ================== #
        csts.files = fileNames
        # print(fileNames)
        try:
            name = path_(fileNames[0]) #make sure files are loaded 
            self.path_data=path_(fileNames[0]).parent # a verifier 
            file_name = os.path.basename(fileNames[0])
            csts.init_directory = self.path_data # to make default directory
            # name=os.path.basename(fileName)
            if name:
                # ======================== change: change butoon text ======================== #
                self.button_ask_path_data.setText(f"loadded {file_name}") # replace "len(csts.files)" by "file_name"
                # self.button_ask_path_data.setText(name)
                self.controler.refreshAll3(f"loaded files:\n{[path_(csts.files[i]).name for i in range(len(csts.files))]}")
        # ---------------------------------------------------------------------------- #
            else:
                self.button_ask_path_data.setText("browse")
            self.controler.initialised = 0
        except:
            self.button_ask_path_data.setText("No files were loaded")
            self.controler.error_message_path3()
            

    def get_path_data_ref(self):
        # options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        # fileName, _ = QFileDialog.getOpenFileName(self,"Reference Data (hdf5 file)", options=options,  filter="Hdf5 File(*.h5)")
        # =========================== initialize csts.refs =========================== #
        csts.refs = []
        # fileName, _ = QFileDialog.getOpenFileName(self,"Reference Data (hdf5 file)",  filter="Hdf5 File(*.h5)")
        # ==================================== new =================================== #
        fileNames, _ = QFileDialog.getOpenFileNames(parent=self,caption=f"Select multiple reference h5 files",directory=f"{csts.init_directory}", filter="Hdf5 File(*.h5)")
        # ================ cahnge : fill refs list with ref filenames ================ #
        csts.refs = fileNames
        try:
            names = [path_(fileNames[i]) for i in range(len(csts.refs))] # make sure files are loaded
            self.path_data_ref = path_(fileNames[0]).parent
            ref_name = os.path.basename(fileNames[0])
            # self.path_data_ref=fileName
            # name=os.path.basename(fileName)
            # if name:
            # ================================== change ================================== #
            if names:
                # self.button_ask_path_without_sample.setText(name)
                self.button_ask_path_without_sample.setText(f"loaded {ref_name} refs") # replace len(csts.refs) by ref_name
                self.controler.refreshAll3(f"loaded refs:\n{[path_(csts.refs[i]).name for i in range(len(csts.refs))]}")
            # ---------------------------------------------------------------------------- #
            else:
                self.button_ask_path_without_sample.setText("browse")
            self.controler.initialised = 0
            self.ref_loaded = True
            
        except:
            self.button_ask_path_without_sample.setText("No files were loaded")
            self.controler.error_message_path3()


    def pressed_loading(self):
        self.controler.loading_text3()
        
    def refresh(self):# /!\ Optimization_tab is not a client of the controler, 
    #this is not called by the controler, only when action_choice is changed.
        deleteLayout(self.action_widget.layout())
        self.action_widget.refresh()
        self.controler.refreshAll3('')
        

class Ui_Dialog(object):
    def setupUi(self, Dialog, controler):
        self.controler = controler
        
        self.length_initialized = 0
        self.trace_start = 0
        self.trace_end = -1
        self.time_start = 0
        self.time_end = -1
        
        self.dialog = Dialog
        self.dialog.resize(400, 126)
        self.dialog.setWindowTitle("Set length of data to analyze (optional)")
        
        self.label_data_length = QLabel('Number of time traces (0,...,N-1)   ')
        
        self.label_data_length_start = QLabel('  start    ')
        self.length_start_limit_box = QLineEdit()
        
        self.label_data_length_end = QLabel('    end      ')
        self.length_end_limit_box = QLineEdit()
        
        self.label_time_length = QLabel('Length of time trace (optional)   ')
        
        self.label_time_length_start = QLabel('start (ps)')
        self.time_start_limit_box = QLineEdit()
        
        self.label_time_length_end = QLabel('end (ps)')
        self.time_end_limit_box = QLineEdit()
        
        self.button_submit_length = QPushButton('Submit')
        self.button_submit_length.clicked.connect(self.action)
        
        self.hlayout0=QHBoxLayout()
        self.hlayout1=QHBoxLayout()
        self.hlayout2=QHBoxLayout()
        
        self.hlayout0.addWidget(self.label_data_length,1)
        self.hlayout0.addWidget(self.label_data_length_start,1)
        self.hlayout0.addWidget(self.length_start_limit_box,1)
        self.hlayout0.addWidget(self.label_data_length_end,1)
        self.hlayout0.addWidget(self.length_end_limit_box,0)
        self.hlayout1.addWidget(self.label_time_length,1)
        self.hlayout1.addWidget(self.label_time_length_start,1)
        self.hlayout1.addWidget(self.time_start_limit_box,0)
        self.hlayout1.addWidget(self.label_time_length_end,1)
        self.hlayout1.addWidget(self.time_end_limit_box,1)
        self.hlayout2.addWidget(self.button_submit_length)
        
        self.sub_layoutv10 = QVBoxLayout(self.dialog)
        self.sub_layoutv10.addLayout(self.hlayout0)
        self.sub_layoutv10.addLayout(self.hlayout1)
        self.sub_layoutv10.addLayout(self.hlayout2)

    def action(self):

        self.length_initialized = 0
        try:
            if self.length_start_limit_box.text() or self.length_end_limit_box.text():
                trace_start = int(self.length_start_limit_box.text())
                trace_end = int(self.length_end_limit_box.text())

                if trace_start  < 0 or trace_start > trace_end:
                    raise ValueError
                else:
                    self.trace_start = trace_start
                
                if trace_end  < 0:
                        raise ValueError
                else:
                    self.trace_end = trace_end
                    
                self.length_initialized = 1
                
                    
            if self.time_start_limit_box.text() or self.time_end_limit_box.text():
                time_start = float(self.time_start_limit_box.text())
                time_end = float(self.time_end_limit_box.text())
                #TODO
                #conditions for time

                self.length_initialized = 1
            
            self.dialog.close()
        except Exception as e:
            print(e)
            self.controler.refreshAll3("Invalid values, please enter positive values and trace start number must be less than or equal to trace end")
            return(0)

class Ui_Dialog_match(object):
    def setupUi(self, Dialog, controler):
        self.controler = controler
        
        
        self.lists_ordered = 0
        # ======================= keep the same for the moment ======================= #
        self.dialog = Dialog
        self.dialog.resize(400, 126)
        # ---------------------------------------------------------------------------- #
        self.dialog.setWindowTitle("match samples to references")
        
        i = 0
        i2 = 0
        
        self.label_files = QLabel(f"Files")
        
        self.labels_files_text = [f"file {i}" for i in range(len(csts.files))]
        self.labels_files = [QLabel(f"file {label}") for label in self.labels_files_text ]
        # print(self.labels_files)
        
        
        self.comboboxes_files = [QComboBox() for i in range(len(csts.files))]
        for combobox in self.comboboxes_files:
            combobox.addItems([path_(csts.files[i]).name for i in range(len(csts.files))])
            combobox.setCurrentIndex(i)
            i+=1
        
        # print(self.comboboxes_files)
        
        self.label_refs = QLabel(f"Reference Files")
        
        self.labels_refs_text = [f"ref {i}" for i in range(len(csts.files))]
        self.labels_refs = [QLabel(f"file {label}") for label in self.labels_refs_text ]
        
        self.comboboxes_refs = [QComboBox() for i in range(len(csts.files))]
        for combobox in self.comboboxes_refs:
            combobox.addItems([path_(csts.refs[i]).name for i in range(len(csts.refs))])
            combobox.setCurrentIndex(i2)
            i2+=1
        
        
        self.button_match = QPushButton(f"Match")
        self.button_match.resize(100, 50)
        self.button_match.clicked.connect(self.action_match)
        
        
        if not csts.files:
            self.show_warning_messagebox()
        
        
        try:
            if csts.refs:
                if self.comboboxes_refs[-1].currentText() == "" :
                    self.show_critical_messagebox()
            else:
                self.show_info_noref_messagebox()
        except:
            pass

        self.button_match.setEnabled(False)
        self.verify_last_combo()
        
        for combobox in self.comboboxes_refs:
            combobox.activated.connect(self.verify_last_combo)
        
        
        self.vlayout = QVBoxLayout()
        self.hlayout = QHBoxLayout()
        self.hlayout_button  = QHBoxLayout()
        
        
        for i in range(len(csts.files)):
            hlayout = QHBoxLayout()
            hlayout.addWidget(self.labels_files[i])
            hlayout.addWidget(self.comboboxes_files[i])
            hlayout.addWidget(self.labels_refs[i])
            hlayout.addWidget(self.comboboxes_refs[i])
            self.vlayout.addLayout(hlayout)
            
        # print(self.vlayout)
            
            # self.vlayout_files.addWidget(self.labels_files[i])
            # self.vlayout_files.addWidget(self.comboboxes_files[i])
            # self.vlayout_files.addWidget(self.labels_refs[i])
            # self.vlayout_files.addWidget(self.comboboxes_refs[i])
        
        # for i in range(len(csts.files)):
        
        self.hlayout_button.addWidget(self.button_match)
        
        
        
        self.general_layout = QVBoxLayout(self.dialog)
        self.general_layout.addLayout(self.vlayout)
        self.general_layout.addLayout(self.hlayout_button)
        
        
        self.general_layout.deleteLater()

    def action_match(self):
        self.lists_ordered = 0
        # self.length_initialized = 0
        # print(csts.files)
        try:
            if self.comboboxes_refs[-1].currentText() != "":
                for i, combobox in enumerate(self.comboboxes_files):
                    # print(combobox.currentText())
                    # print(csts.files[i])
                    csts.files[i] = path_(csts.files[i]).parent.joinpath(combobox.currentText())
                # print(csts.files)
                for i, combobox in enumerate(self.comboboxes_refs):
                    # print(combobox.currentText())
                    if i < len(csts.refs):
                        csts.refs[i] = path_(csts.refs[i]).parent.joinpath(combobox.currentText())
                    elif i >= len(csts.refs):
                        csts.refs.append(path_(csts.refs[0]).parent.joinpath(combobox.currentText()))
                # print(f"refs={(csts.refs)}")

            self.dialog.close()
        except Exception as e:
            print(e)
            self.controler.refreshAll3("Error ")
            return(0)
    
    def verify_last_combo(self):
        # print(f"verifying...")
        if csts.refs and self.comboboxes_refs[-1].currentText()!="":
            self.button_match.setEnabled(True)
        else:
            self.button_match.setEnabled(False)
    
    def show_critical_messagebox(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        
        # ====================== setting message for Message Box ===================== #
        msg.setText("All samples files must be matched to a reference file")
        
        # ===================== setting Message box window title ===================== #
        msg.setWindowTitle("Missing Reference File")
    
        # ===================== declaring buttons on Message Box ===================== #
        # msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        msg.setStandardButtons(QMessageBox.Ok)
    
        # =============================== start the app ============================== #
        retval = msg.exec_()

    def show_warning_messagebox(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)

        # ====================== setting message for Message Box ===================== #
        msg.setText("There is nothing to match !")
    
        # ===================== setting Message box window title ===================== #
        msg.setWindowTitle("No data files were loaded")
    
        # ===================== declaring buttons on Message Box ===================== #
        # msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        msg.setStandardButtons(QMessageBox.Ok)
    
        # =============================== start the app ============================== #
        retval = msg.exec_()
    
    def show_info_noref_messagebox(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        # ====================== setting message for Message Box ===================== #
        msg.setText("There is nothing to match !")
    
        # ===================== setting Message box window title ===================== #
        msg.setWindowTitle("No ref files were loaded")
    
        # ===================== declaring buttons on Message Box ===================== #
        # msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        msg.setStandardButtons(QMessageBox.Ok)
    
        # =============================== start the app ============================== #
        retval = msg.exec_()
    
    
    
class TextBoxWidget(QTextEdit):
    def __init__(self, parent, controler):
        super().__init__(parent)
        self.controler = controler
        self.controler.addClient(self)
        self.controler.addClient3(self)
        self.setReadOnly(True)

        self.setMinimumWidth(560)
        
        references = ['\u00b9 Time should be in ps.\nDatasets in the hdf5 must be named ["timeaxis", "0", "1", ..., "N-1"]\nExample: if we have 1000 time traces, N-1 = 999',
                      '\u00b2 Use to take into account the reference traces in the covariance computation \n(only if the initial data are measures with a sample)\ndon\'t forget to apply the same filters/correction to the reference before',
                      '\u00b3 Sharpness: 100 is almost a step function, 0.1 is really smooth. See graphs in optimization tab.',
                      '\u2074 Plot the noise convolution matrix or the covariance matrix depending on the input \n The matrix must be saved for its computation (Graphical lasso)'
                      ]
        
        references_2=["\u2075 Error options:",
                    "Constant weight: error = \u2016Efit(w)-E(w)\u2016 / \u2016E(w)\u2016",
                    "or error = \u2016Efit(t)-E(t)\u2016 / \u2016E(t)\u2016 in super resolution"]
            
        for reference in references:
            self.append(reference)
            self.append('')
            
    def refresh(self):
        message = self.controler.message
        if type(message)==list:
            for i in message:
                self.append(i)
        else:
            self.append(message)
            



###############################################################################
###############################################################################
#########################   Optimisation tab   ################################
###############################################################################
###############################################################################



class Action_handler(QWidget):
    def __init__(self, parent, controler):
        super().__init__(parent)
        self.parent = parent
        self.controler = controler
        #self.setMaximumWidth(500)

    def refresh(self):
        self.action_widget = Optimization_choices(self,self.controler)

        try:
            deleteLayout(self.layout())
            self.layout().deleteLater()
        except:
            main_layout = QVBoxLayout()
            main_layout.addWidget(self.action_widget)
            self.setLayout(main_layout)

class Optimization_choices(QGroupBox):
    def __init__(self, parent, controler):
        super().__init__(parent)
        self.controler=controler
        self.parent = parent
        #self.setMinimumHeight(150)
        self.emptycellfile = None
        self.setTitle("Optimization")
        self.algo_index = 0
        # Corrective factors
        label_width=200
        action_widget_width=150
        corrective_width_factor=-12
        text_box_height = 22
        # Algorithm choice
        self.label_algo = QLabel("Algorithm - delay/amplitude/dilatation")
        self.label_algo.setMaximumHeight(text_box_height)
        self.options_algo = QComboBox()
        self.options_algo.addItems(['NumPy optimize swarm particle',
                                    '[In Dev]ALPSO without parallelization',
                                    '[In Dev]ALPSO with parallelization',
                                    '[In Dev]SLSQP (pyOpt)',
                                    '[In Dev]SLSQP (pyOpt with parallelization)',
                                    'L-BFGS-B',
                                    'SLSQP (scipy)',
                                    'Dual annealing'])
        self.options_algo.setMaximumHeight(text_box_height)
        self.options_algo.setCurrentIndex(6)
        self.options_algo.currentIndexChanged.connect(self.refresh_param)
        
        
        self.label_algo_ps = QLabel("Algorithm - periodic sampling")
        self.label_algo_ps.setMaximumHeight(text_box_height)
        self.options_algo_ps = QComboBox()
        self.options_algo_ps.addItems(['Dual annealing'])
        self.options_algo_ps.setMaximumHeight(text_box_height)
        self.options_algo_ps.currentIndexChanged.connect(self.refresh_param)
        
        
            # Number of iterations
        self.label_niter = QLabel("\tIterations")
        self.label_niter.setMaximumHeight(text_box_height)
        self.enter_niter = QLineEdit()
        self.enter_niter.setMaximumWidth(action_widget_width +corrective_width_factor)
        self.enter_niter.setMaximumHeight(30)
        self.enter_niter.setMaximumHeight(text_box_height)
        self.enter_niter.setText('1000')
    
            # Number of iterations ps
        self.label_niter_ps = QLabel("Iterations")
        self.label_niter_ps.setMaximumHeight(text_box_height)
        self.enter_niter_ps = QLineEdit()
        self.enter_niter_ps.setMaximumWidth(action_widget_width +corrective_width_factor)
        self.enter_niter_ps.setMaximumHeight(30) 
        self.enter_niter_ps.setMaximumHeight(text_box_height)
        self.enter_niter_ps.setText("1000")
        
                    # SwarmSize
        self.label_swarmsize = QLabel("    Swarmsize")
        self.label_swarmsize.setMaximumHeight(text_box_height)
        self.enter_swarmsize = QLineEdit()
        self.enter_swarmsize.setMaximumWidth(action_widget_width +corrective_width_factor)
        self.enter_swarmsize.setMaximumHeight(30)
        self.enter_swarmsize.setMaximumHeight(text_box_height)
        
        self.refresh_param()
        
            # Button to launch optimization
        self.begin_button = QPushButton("Begin Optimization")
        #TOCHANGE
        #TOVERIFY
        # self.begin_button.clicked.connect(self.begin_optimization)
        self.begin_button.clicked.connect(self.begin_optimizations)
        
        
        # self.begin_button.clicked.connect(self.controler.begin_optimization)
        self.begin_button.pressed.connect(self.pressed_loading)
        self.begin_button.setMaximumHeight(text_box_height)
        self.begin_button.setStyleSheet("background-color : rgb(0, 200, 83)")
        
        #     # Button to stop optimization
        # self.break_button = QPushButton("Stop Optimization (In Dev)")
        # self.break_button.pressed.connect(self.break_loading)
        # self.break_button.clicked.connect(self.controler.stop_optimization)
        # self.break_button.setMaximumHeight(text_box_height)
        # self.break_button.setStyleSheet("background-color : rgb(255, 0, 0)")  # Button to stop opt commented cause it's not working
        
        # #TEST:
        # # ============== add a progress bar to show the overall progress ============= #
        # self.opt_progressbar_label = QLabel(f"Batch progress")
        # self.opt_progressbar = QProgressBar(self,objectName="BlueProgressBar")
        # sub_layoout_pb  = QHBoxLayout()
        # sub_layoout_pb.addWidget(self.opt_progressbar_label)
        # sub_layoout_pb.addWidget(self.opt_progressbar) # Progress bar commented cause it's not working
        
        # COMMENT TO SUPP THE PROGRESS BAR
        
        # # ---------------------------------------------------------------------------- #
        
        
        # Wiget to see how many process are going to be used for the omptimization
        self.label_nb_proc = QLabel("How many process do you want to use?")
        self.enter_nb_proc = QLineEdit()
        self.enter_nb_proc.setText('1')
        self.enter_nb_proc.setMaximumWidth(50)
        self.enter_nb_proc.setMaximumHeight(25)

        # Creation layouts
        sub_layout_h=QHBoxLayout()
        sub_layout_h_2=QHBoxLayout()
        sub_layout_h_3=QHBoxLayout()
        sub_layout_h_7=QHBoxLayout()
        sub_layout_h_8=QHBoxLayout()
        main_layout=QVBoxLayout()
        

        # Organisation layouts
        
        
        # Organisation layouts for optimisation
        sub_layout_h_7.addWidget(self.label_algo,0)
        sub_layout_h_7.addWidget(self.options_algo,0)
        sub_layout_h_3.addWidget(self.label_niter,0)
        sub_layout_h_3.addWidget(self.enter_niter,0)
        sub_layout_h_3.addWidget(self.label_swarmsize,0)
        sub_layout_h_3.addWidget(self.enter_swarmsize,0)
        sub_layout_h_2.addWidget(self.label_algo_ps,0)
        sub_layout_h_2.addWidget(self.options_algo_ps,0)
        sub_layout_h_2.addWidget(self.label_niter_ps,0)
        sub_layout_h_2.addWidget(self.enter_niter_ps,0)
        sub_layout_h_8.addWidget(self.begin_button,0)
        # sub_layout_h_8.addWidget(self.break_button,0) # Button to stop opt cause it's not working

        # Vertical layout   
        main_layout.addLayout(sub_layout_h_7)
        main_layout.addLayout(sub_layout_h_3)
        main_layout.addLayout(sub_layout_h_2)
        main_layout.addLayout(sub_layout_h_8)
        # ====================== progress bar add to main_layout ===================== #
        # main_layout.addLayout(sub_layoout_pb)
        
        # COMMENT TO SUPP THE PROGRESS BAR
        
        # ---------------------------------------------------------------------------- #
        self.setLayout(main_layout)



    def pressed_loading(self):
        self.controler.loading_text3()
    # def break_loading(self):
    #     self.controler.loading_text_break() # Button to stop opt cause it's not working

# ============================================================================ #
#                              Begin optimization                              #
# ============================================================================ #
# ============================================================================ #
#           define new begin_optimization func to loop for all files           #
# ============================================================================ #
# ============================================================================ #
# ======================= be careful with parent.parent ====================== #
# ============================================================================ #

    def begin_optimizations(self):
        # try:
        # self.opt_progressbar.setMaximum(len(csts.files))
        
        
        # COMMENT TO SUPP THE PROGRESS BAR
        
        if self.parent.parent.ref_loaded and len(csts.refs) != len(csts.files):
            raise ValueError
        i=0
        for file in csts.files:
            print(file)
            if len(csts.refs) != 0:
                for ref in csts.refs:
                    self.begin_optimization(file,ref)
            else:
                self.begin_optimization(file)
            i+=1
            self.controler.refreshAll3(f"\nfile {i}/{len(csts.files)} finished optimization")
            self.parent.parent.text_box.repaint()
            if not self.parent.parent.ref_loaded:
                if csts.save_bools["mean"]:
                    self.parent.parent.save_param.save_mean_batch(file)
            elif  self.parent.parent.ref_loaded:
                self.parent.parent.save_param.save_mean_batch(file,ref)
            
            if csts.save_bools["std_time"]:    
                self.parent.parent.save_param.save_std_time_batch(file)
            if csts.save_bools["std_freq"]:    
                self.parent.parent.save_param.save_std_freq_batch(file)
            
            if csts.save_bools["correction_param"]:
                self.parent.parent.save_param.save_param_batch(file)
            
            if csts.save_bools["time_traces"]:
                self.parent.parent.save_param.save_traces_batch(file)
                
            if csts.save_bools["noise_matrix"]:
                self.parent.parent.save_param.save_cov_batch(file)
        
        # except Exception:
            # print(f"you should match references to number of sample files")
            # self.controler.refreshAll3(f"You should match references to number of sample files")
    
        # ============================================================================ #
        #           add all optimized files to list to be able to plot later           #
        # ============================================================================ #
            csts.myinput.append(self.controler.myinput)
            csts.myreferencedata.append(self.controler.myreferencedata)
            csts.ncm.append(self.controler.ncm) 
            csts.ncm_inverse.append(self.controler.ncm_inverse) 
            csts.reference_number.append(self.controler.reference_number) 
            csts.mydatacorrection.append(self.controler.mydatacorrection)
            csts.delay_correction.append(self.controler.delay_correction) 
            csts.dilatation_correction.append(self.controler.dilatation_correction) 
            csts.leftover_correction.append(self.controler.leftover_correction) 
            csts.myglobalparameters.append(self.controler.myglobalparameters) 
            csts.fopt.append(self.controler.fopt) 
            csts.fopt_init.append(self.controler.fopt_init) 
            csts.mode.append(self.controler.mode)
            # self.opt_progressbar.setValue(i)
        self.parent.parent.graph_widget.add_items_to_combobox(csts.files)
        # ---------------------------------------------------------------------------- #

    
    def begin_optimization(self, file,ref_file=None):
        global preview
        # self.print_ui_params()
        # self.method = InitParamWidget(InitParamWidget.parent)
        # InitParamWidget.on_click(InitParamWidget.parent,csts.files[0])
        # InitParamWidget.on_click_param(self)
        # TOVERIFY
        # self.parent.parent.on_click(csts.files[0])
        self.parent.parent.on_click(file,ref_file)
        self.parent.parent.on_click_param()
        
        t1 = time.time()
        global graph_option_2
        submitted = self.submit_algo_param()    #get values from optimisation widget
        if submitted == 1:
            try:
                from mpi4py import MPI
                nb_proc=int(self.enter_nb_proc.text())
            except:
                self.controler.message_log_tab3("You don't have MPI for parallelization, we'll use only 1 process")
                nb_proc=1
            if self.controler.optim.vars_temp_file_5 != None:
                self.controler.begin_optimization(nb_proc)
                graph_option_2='Pulse (E_field)'
            else:
                self.controler.no_temp_file_5()
        preview = 0
        print("\nTime taken by the optimization:")
        self.controler.optim_succeed = 1
        self.controler.refreshAll3("")
        print(time.time()-t1)
    
    def refresh_param(self):
        self.algo_index = self.options_algo.currentIndex()
        if self.algo_index < 3:
            self.label_swarmsize.show()
            self.enter_swarmsize.show()
        else :
            self.label_swarmsize.hide()
            self.enter_swarmsize.hide()
            
            
    def submit_algo_param(self):
        choix_algo=self.algo_index
        swarmsize = 0
        if not self.controler.initialised  or not self.controler.initialised_param:
            self.controler.refreshAll3("Please submit initialization and/or correction parameters first")
            return(0)
        
        if not self.controler.fit_delay  and not self.controler.fit_leftover_noise and not self.controler.fit_periodic_sampling and not self.controler.fit_dilatation:
            self.controler.refreshAll3("Please submit correction parameters first")
            return(0)
        if self.controler.fit_delay or self.controler.fit_leftover_noise or self.controler.fit_dilatation:
            if (choix_algo == 3 or choix_algo == 4):
                try:
                    from pyOpt import SLSQP
                except:
                    self.controler.refreshAll3("SLSQP was not imported successfully and can't be used")
                    return(0)
            if choix_algo <3: #particle swarm
                try:
                    swarmsize=int(self.enter_swarmsize.text())
                    if swarmsize<0:
                        self.controler.invalid_swarmsize()
                        return(0)
                except:
                    self.controler.invalid_swarmsize()
                    return(0)
        try:
            niter = 0
            niter_ps = 0
            if self.controler.fit_delay or self.controler.fit_leftover_noise or self.controler.fit_dilatation:
                niter=int(self.enter_niter.text())
                if niter<0:
                    self.controler.invalid_niter()
                    return(0)
            if self.controler.fit_periodic_sampling :
                niter_ps=int(self.enter_niter_ps.text())
                if niter_ps<0:
                    self.controler.invalid_niter()
                    return(0)
            self.controler.algo_parameters(choix_algo,swarmsize,niter,niter_ps)
        except Exception as e:
            print(e)
            self.controler.invalid_niter()
            return(0)
        return(1)
        
class Saving_parameters(QGroupBox):
    def __init__(self, parent, controler):
        super().__init__(parent)
        self.parent = parent
        self.controler = controler
        text_box_height = 22
        self.setTitle("Saving results")
        
        # Corrective factors
        action_widget_width=150
        corrective_width_factor=-12


        # self.button_save_mean = QPushButton('Mean (.txt)', self)
        # self.button_save_mean.clicked.connect(self.save_mean)
        # self.button_save_mean.setMaximumHeight(text_box_height)    
        
        self.checkbox_save_mean = QCheckBox('Mean', self)
        self.checkbox_save_mean.stateChanged.connect(self.check_save_mean)
        
        self.checkbox_save_std_time = QCheckBox('Std time', self)
        self.checkbox_save_std_time.stateChanged.connect(self.check_save_std_time)
        
        self.checkbox_save_std_freq = QCheckBox('Std freq', self)
        self.checkbox_save_std_freq.stateChanged.connect(self.check_save_std_freq)
        
        self.checkbox_save_traces = QCheckBox('Time traces', self)
        self.checkbox_save_traces.stateChanged.connect(self.check_save_traces)
        
        self.checkbox_save_param = QCheckBox('Correction param', self)
        self.checkbox_save_param.stateChanged.connect(self.check_save_param)
        
        self.checkbox_save_cov = QCheckBox('Noise matrices', self) # "noise matrix" replaced by "Noise matrices"
        self.checkbox_save_cov.stateChanged.connect(self.check_save_cov)
        
        # self.button_save_traces = QPushButton('Time traces (.h5)', self)
        # self.button_save_traces.clicked.connect(self.save_traces)
        # self.button_save_traces.setMaximumHeight(text_box_height) 
        
        # self.button_save_param = QPushButton('Correction param (.txt)', self)
        # self.button_save_param.clicked.connect(self.save_param)
        # self.button_save_param.setMaximumHeight(text_box_height)
        
        # self.button_save_cov = QPushButton('Noise matrix inverse (.h5)', self)
        # self.button_save_cov.clicked.connect(self.save_cov)
        # self.button_save_cov.setMaximumHeight(text_box_height)
        
        # self.button_save_std_time = QPushButton('Std time (.txt)', self)
        # self.button_save_std_time.clicked.connect(self.save_std_time)
        # self.button_save_std_time.setMaximumHeight(text_box_height) 
        
        # self.button_save_std_freq = QPushButton('Std freq (.txt)', self)
        # self.button_save_std_freq.clicked.connect(self.save_std_freq)
        # self.button_save_std_freq.setMaximumHeight(text_box_height) 
        
        # self.button_save_mean.setEnabled(False)
        # self.button_save_traces.setEnabled(False)
        # self.button_save_param.setEnabled(False)
        # self.button_save_cov.setEnabled(False)
        # self.button_save_std_time.setEnabled(False)
        # self.button_save_std_freq.setEnabled(False)3
        
        sub_layout_h1=QHBoxLayout()
        sub_layout_h2=QHBoxLayout()
        sub_layout_h3=QHBoxLayout()            

        # sub_layout_h3.addWidget(self.button_save_mean,0)
        # sub_layout_h3.addWidget(self.button_save_param,0)
        # sub_layout_h3.addWidget(self.button_save_traces,0)
        # sub_layout_h2.addWidget(self.button_save_std_time,0)
        # sub_layout_h2.addWidget(self.button_save_std_freq,0)
        # sub_layout_h2.addWidget(self.button_save_cov,0)
        sub_layout_h3.addWidget(self.checkbox_save_mean,0)
        sub_layout_h3.addWidget(self.checkbox_save_param,0)
        sub_layout_h3.addWidget(self.checkbox_save_traces,0)
        sub_layout_h2.addWidget(self.checkbox_save_std_time,0)
        sub_layout_h2.addWidget(self.checkbox_save_std_freq,0)
        sub_layout_h2.addWidget(self.checkbox_save_cov,0)
        
        self.checkbox_save_mean.setChecked(True)
        self.checkbox_save_traces.setChecked(False)
        self.checkbox_save_param.setChecked(False)
        self.checkbox_save_cov.setChecked(False)
        self.checkbox_save_std_time.setChecked(True)
        self.checkbox_save_std_freq.setChecked(True)

        
        self.main_layout=QVBoxLayout()
        self.main_layout.addLayout(sub_layout_h3)
        self.main_layout.addLayout(sub_layout_h2)
        
        self.setLayout(self.main_layout)
        
    def check_save_mean(self):
        if self.checkbox_save_mean.isChecked():
            csts.save_bools["mean"] = True
            self.controler.refreshAll3(f"\nMean will be saved\n")
        elif not self.checkbox_save_mean.isChecked():
            csts.save_bools["mean"] = False
            self.controler.refreshAll3(f"\nMean will NOT be saved\n")
    
    def check_save_std_time(self):        
        if self.checkbox_save_std_time.isChecked():
            csts.save_bools["std_time"] = True
            self.controler.refreshAll3(f"\nStd time will be saved\n")
        if not self.checkbox_save_std_time.isChecked():
            csts.save_bools["std_time"] = False
            self.controler.refreshAll3(f"\nStd time will NOT be saved\n")
    
    def check_save_std_freq(self):    
        if self.checkbox_save_std_freq.isChecked():
            csts.save_bools["std_freq"] = True
            self.controler.refreshAll3(f"\nStd freq will be saved\n")
        if not self.checkbox_save_std_freq.isChecked():
            csts.save_bools["std_freq"] = False
            self.controler.refreshAll3(f"\nStd freq will NOT be saved\n")
    
    def check_save_traces(self):    
        if self.checkbox_save_traces.isChecked():
            csts.save_bools["time_traces"] = True
            self.controler.refreshAll3(f"\nTime traces will be saved\n")
        elif not self.checkbox_save_traces.isChecked():
            csts.save_bools["time_traces"] = False
            self.controler.refreshAll3(f"\nTime traces will NOT be saved\n")
    
    def check_save_param(self):
        if self.checkbox_save_param.isChecked():
            csts.save_bools["correction_param"] = True
            self.controler.refreshAll3(f"\nCorrection parameters will be saved\n")
        if not self.checkbox_save_param.isChecked():
            csts.save_bools["correction_param"] = False
            self.controler.refreshAll3(f"\nCorrection parameters will NOT be saved\n")
    
    def check_save_cov(self):    
        if self.checkbox_save_cov.isChecked():
            csts.save_bools["noise_matrix"] = True
            self.controler.refreshAll3(f"\nNoise matrices will be saved\n") # "noise matrix" replaced by "Noise matrices"
        if not self.checkbox_save_cov.isChecked():
            csts.save_bools["noise_matrix"] = False
            self.controler.refreshAll3(f"\nNoise matrices will NOT be saved\n") # "noise matrix" replaced by "Noise matrices"
        
        
        # print(csts.save_bools)
        
        # print(f"Will save mean" if csts.save_bools["mean"] else f"Will not save mean")
        # print(f"Will save std time" if csts.save_bools["std_time"] else f"Will not save std time")
        # print(f"Will save std freq" if csts.save_bools["std_freq"] else f"Will not save std freq")
    
    # def save_cov(self):
        # global preview
        # if self.controler.initialised:
            # try:
                # # options = QFileDialog.Options()
                # # options |= QFileDialog.DontUseNativeDialog
                # # fileName, _ = QFileDialog.getSaveFileName(self,"Covariance / Noise convolution matrix","noise.h5","HDF5 (*.h5)", options=options)
                # fileName, _ = QFileDialog.getSaveFileName(self,"Covariance / Noise convolution matrix","noise.h5","HDF5 (*.h5)")
                # try:
                    # name=os.path.basename(fileName)
                    # path = os.path.dirname(fileName)
                    # if name:
                        # saved = self.controler.save_data(name, path, 5)
                        # if saved:
                            # if not self.controler.optim_succeed:
                                # preview = 1
                            # self.controler.refreshAll3(" Saving matrix - Done")
                        # else:
                            # print("Something went wrong")          
                # except:
                    # self.controler.error_message_output_filename()
            # except:
                # self.controler.error_message_output_filename()
                # return(0)
        # else:
            # self.controler.refreshAll3("Please enter initialization data first")
    
    def save_mean_batch(self,filename,ref=None):
        global preview
        if self.controler.initialised:
            name = f"corrected_mean_{path_(filename).stem}.txt"
            # # path = path_(filename).parent.joinpath(f"correct@tds_save_data")
            # print(f"modesuper : {csts.modesuper}")
            # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            
            if csts.modesuper:
                path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
                if not path_(path_init).is_dir():
                    path_(path_init).mkdir()
                    path_(path).mkdir()
            
            else:
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                if not path_(path).is_dir():
                    path_(path).mkdir() 
                # if not path_(path).is_dir():
            # else:
                # pass
            
            if name:
                saved = self.controler.save_data(name, path, 0)
                if saved:
                    if not self.controler.optim_succeed:
                        preview = 1
                    self.controler.refreshAll3(" Saving mean - Done")
                else:
                    print("Something went wrong")          
            # print(path.joinpath(path_(filename).name))
            if sys.platform =="linux" or sys.platform=="darwin":
                subprocess.run(f"cp {filename} {path.joinpath(path_(filename).name)}", shell=True)
            elif sys.platform =="win32" or sys.platform=="cygwin":
                subprocess.run(f"copy {path_(filename)} {path.joinpath(path_(filename).name)}", shell=True)
            if ref != None:
                if sys.platform=="linux" or sys.platform=="darwin":
                    subprocess.run(f"cp {ref} {path.joinpath(path_(ref).name)}", shell=True)
                elif sys.platform=="win32" or sys.platform=="cygwin":
                    subprocess.run(f"copy {path_(ref)} {path.joinpath(path_(ref).name)}", shell=True)
                    
        else:
            self.controler.refreshAll3("Please enter initialization data first")
    
    def save_std_time_batch(self,filename):
        global preview
        if self.controler.initialised:
            name = f"corrected_std_time_{path_(filename).stem}.txt"
            # # path = path_(filename).parent.joinpath(f"correct@tds_save_data")
            # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            
            
            
            if csts.modesuper:
                path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
                if not path_(path_init).is_dir():
                    path_(path_init).mkdir()
                    path_(path).mkdir()
            
            else:
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                if not path_(path).is_dir():
                    path_(path).mkdir() 
            
            # path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if csts.modesuper:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
            # else:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if not path_(path_init).is_dir():
                # path_(path_init).mkdir()
                # path_(path).mkdir()
            # elif path_(path_init).is_dir():
                # if not path_(path).is_dir():
                    # path_(path).mkdir() 
            # else:
                # pass

            if name:
                saved = self.controler.save_data(name, path, 3)
                if saved:
                    if not self.controler.optim_succeed:
                        preview = 1
                    self.controler.refreshAll3(" Saving std in time domain - Done")
                else:
                    print("Something went wrong")          
        else:
            self.controler.refreshAll3("Please enter initialization data first")  
            
    def save_std_freq_batch(self,filename):
        global preview
        if self.controler.initialised:
            name = f"corrected_std_freq_{path_(filename).stem}.txt"
            # # path = path_(filename).parent.joinpath(f"correct@tds_save_data")
            # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            
            if csts.modesuper:
                path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
                if not path_(path_init).is_dir():
                    path_(path_init).mkdir()
                    path_(path).mkdir()
            
            else:
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                if not path_(path).is_dir():
                    path_(path).mkdir() 
            
            # path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if csts.modesuper:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
            # else:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if not path_(path_init).is_dir():
                # path_(path_init).mkdir()
                # path_(path).mkdir()
            # elif path_(path_init).is_dir():
                # if not path_(path).is_dir():
                    # path_(path).mkdir() 
            # else:
                # pass
            
            if name:
                saved = self.controler.save_data(name, path, 4)
                if saved:
                    if not self.controler.optim_succeed:
                        preview = 1
                    self.controler.refreshAll3(" Saving std in frequency domain - Done")
                else:
                    print("Something went wrong")          
        else:
            self.controler.refreshAll3("Please enter initialization data first")

    def save_param_batch(self,filename):
        global preview
        if self.controler.optim_succeed:
            name = f"correction_params_{path_(filename).stem}.txt"
            # # path = path_(filename).parent.joinpath(f"correct@tds_save_data")
            # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            
            if csts.modesuper:
                path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
                if not path_(path_init).is_dir():
                    path_(path_init).mkdir()
                    path_(path).mkdir()
            
            else:
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                if not path_(path).is_dir():
                    path_(path).mkdir() 
            
            
            # path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if csts.modesuper:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
            # else:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if not path_(path_init).is_dir():
                # path_(path_init).mkdir()
                # path_(path).mkdir()
            # elif path_(path_init).is_dir():
                # if not path_(path).is_dir():
                    # path_(path).mkdir() 
            # else:
                # pass
            
            if name:
                saved = self.controler.save_data(name, path, 1)
                if saved:
                    if not self.controler.optim_succeed:
                        preview = 1
                    self.controler.refreshAll3(" Saving parameters - Done")
                else:
                    print("Something went wrong")          
        else:
            self.controler.refreshAll3("Please launch an optimization first")
    
    def save_traces_batch(self,filename):
        global preview
        if self.controler.initialised:
            name = f"corrected_time_traces_{path_(filename).stem}.txt"
            # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            
            if csts.modesuper:
                path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
                if not path_(path_init).is_dir():
                    path_(path_init).mkdir()
                    path_(path).mkdir()
            
            else:
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                if not path_(path).is_dir():
                    path_(path).mkdir() 
            
            # path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if csts.modesuper:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
            # else:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if not path_(path_init).is_dir():
                # path_(path_init).mkdir()
                # path_(path).mkdir()
            # elif path_(path_init).is_dir():
                # if not path_(path).is_dir():
                    # path_(path).mkdir() 
            # else:
                # pass
            
            if name:
                saved = self.controler.save_data(name, path, 2)
                if saved:
                    if not self.controler.optim_succeed:
                        preview = 1
                    self.controler.refreshAll3(" Saving each traces - Done")
                else:
                    print("Something went wrong")          
        else:
            self.controler.refreshAll3("Please enter initialization data first")
    
    def save_cov_batch(self,filename):
        global preview
        if self.controler.initialised:
            name = f"noise_matrix_{path_(filename).stem}.txt"
            # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            
            if csts.modesuper:
                path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
                if not path_(path_init).is_dir():
                    path_(path_init).mkdir()
                    path_(path).mkdir()
            
            else:
                path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
                if not path_(path).is_dir():
                    path_(path).mkdir() 
            
            # path_init = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if csts.modesuper:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}").joinpath(f"superresolution")
            # else:
                # path = path_(filename).parent.joinpath(f"{path_(filename).stem}")
            # if not path_(path_init).is_dir():
                # path_(path_init).mkdir()
                # path_(path).mkdir()
            # elif path_(path_init).is_dir():
                # if not path_(path).is_dir():
                    # path_(path).mkdir() 
            # else:
                # pass
            
            if name:
                saved = self.controler.save_data(name, path, 5)
                if saved:
                    if not self.controler.optim_succeed:
                        preview = 1
                    self.controler.refreshAll3(" Saving matrix - Done")
                else:
                    print("Something went wrong")          
        else:
            self.controler.refreshAll3("Please enter initialization data first")
    
    # def save_mean(self):
        # global preview
        # if self.controler.initialised:
            # try:
                # # options = QFileDialog.Options()
                # # options |= QFileDialog.DontUseNativeDialog
                # # fileName, _ = QFileDialog.getSaveFileName(self,"Mean file","mean.txt","TXT (*.txt)", options=options)
                # fileName, _ = QFileDialog.getSaveFileName(self,"Mean file","mean.txt","TXT (*.txt)")
                # try:
                    # name=os.path.basename(fileName)
                    # path = os.path.dirname(fileName)
                    # if name:
                        # saved = self.controler.save_data(name, path, 0)
                        # if saved:
                            # if not self.controler.optim_succeed:
                                # preview = 1
                            # self.controler.refreshAll3(" Saving mean - Done")
                        # else:
                            # print("Something went wrong")          
                # except:
                    # self.controler.error_message_output_filename()
            # except:
                # self.controler.error_message_output_filename()
                # return(0)
        # else:
            # self.controler.refreshAll3("Please enter initialization data first")
            
    # def save_param(self):
        # global preview
        # if self.controler.optim_succeed:
            # try:
                # # options = QFileDialog.Options()
                # # options |= QFileDialog.DontUseNativeDialog
                # # fileName, _ = QFileDialog.getSaveFileName(self,"Optimization parameters filename","correction_parameters.txt","TXT (*.txt)", options=options)
                # fileName, _ = QFileDialog.getSaveFileName(self,"Optimization parameters filename","correction_parameters.txt","TXT (*.txt)")
                # try:
                    # name=os.path.basename(fileName)
                    # path = os.path.dirname(fileName)
                    # if name:
                        # saved = self.controler.save_data(name, path, 1)
                        # if saved:
                            # if not self.controler.optim_succeed:
                                # preview = 1
                            # self.controler.refreshAll3(" Saving parameters - Done")
                        # else:
                            # print("Something went wrong")          
                # except:
                    # self.controler.error_message_output_filename()
            # except:
                # self.controler.error_message_output_filename()
                # return(0)
        # else:
            # self.controler.refreshAll3("Please launch an optimization first")
        
    # def save_traces(self):
        # global preview
        # if self.controler.initialised:
            # try:
                # # options = QFileDialog.Options()
                # # options |= QFileDialog.DontUseNativeDialog
                # # fileName, _ = QFileDialog.getSaveFileName(self,"Each traces filename","traces.h5","HDF5 (*.h5)", options=options)
                # fileName, _ = QFileDialog.getSaveFileName(self,"Each traces filename","traces.h5","HDF5 (*.h5)")
                # try:
                    # name=os.path.basename(fileName)
                    # path = os.path.dirname(fileName)
                    # if name:
                        # saved = self.controler.save_data(name, path, 2)
                        # if saved:
                            # if not self.controler.optim_succeed:
                                # preview = 1
                            # self.controler.refreshAll3(" Saving each traces - Done")
                        # else:
                            # print("Something went wrong")          
                # except:
                    # self.controler.error_message_output_filename()
            # except:
                # self.controler.error_message_output_filename()
                # return(0)
        # else:
            # self.controler.refreshAll3("Please enter initialization data first")
            
        
    # def save_std_time(self):
        # global preview
        # if self.controler.initialised:
            # try:
                # # options = QFileDialog.Options()
                # # options |= QFileDialog.DontUseNativeDialog
                # # fileName, _ = QFileDialog.getSaveFileName(self,"Std time domain file","std_time.txt","TXT (*.txt)", options=options)
                # fileName, _ = QFileDialog.getSaveFileName(self,"Std time domain file","std_time.txt","TXT (*.txt)")
                # try:
                    # name=os.path.basename(fileName)
                    # path = os.path.dirname(fileName)
                    # if name:
                        # saved = self.controler.save_data(name, path, 3)
                        # if saved:
                            # if not self.controler.optim_succeed:
                                # preview = 1
                            # self.controler.refreshAll3(" Saving std in time domain - Done")
                        # else:
                            # print("Something went wrong")          
                # except:
                    # self.controler.error_message_output_filename()
            # except:
                # self.controler.error_message_output_filename()
                # return(0)
        # else:
            # self.controler.refreshAll3("Please enter initialization data first")

            

    # def save_std_freq(self):
        # global preview
        # if self.controler.initialised:
            # try:
                # # options = QFileDialog.Options()
                # # options |= QFileDialog.DontUseNativeDialog
                # # fileName, _ = QFileDialog.getSaveFileName(self,"Std frequency domain file","std_freq.txt","TXT (*.txt)", options=options)
                # fileName, _ = QFileDialog.getSaveFileName(self,"Std frequency domain file","std_freq.txt","TXT (*.txt)")
                # try:
                    # name=os.path.basename(fileName)
                    # path = os.path.dirname(fileName)
                    # if name:
                        # saved = self.controler.save_data(name, path, 4)
                        # if saved:
                            # if not self.controler.optim_succeed:
                                # preview = 1
                            # self.controler.refreshAll3(" Saving std in frequency domain - Done")
                        # else:
                            # print("Something went wrong")          
                # except:
                    # self.controler.error_message_output_filename()
            # except:
                # self.controler.error_message_output_filename()
                # return(0)
        # else:
            # self.controler.refreshAll3("Please enter initialization data first")
            
            
    def refresh():
        pass
        
    def pressed_loading(self):
        self.controler.loading_text3()
        
    def get_outputdir(self):
        DirectoryName = QFileDialog.getExistingDirectory(self,"Select Directory")
        try:
            self.outputdir=str(DirectoryName)
            name=os.path.basename(str(DirectoryName))
            if name:
                self.button_outputdir.setText(name)
            else:
                self.button_outputdir.setText("browse")
        except:
            self.controler.error_message_path3()        
                

class Graphs_optimisation(QGroupBox):
    def __init__(self, parent, controler):
        super().__init__(parent)
        #TOVERIFY:
        self.parent = parent
        # ---------------------------------------------------------------------------- #
        self.controler = controler
        self.controler.addClient3(self)
        self.setTitle("Graphs")
        # Create objects to plot graphs
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.canvas.draw() 
        

        self.button_E_field_dB = QPushButton('\nE field  [dB]\n', self)
        self.button_E_field_dB.clicked.connect(self.E_field_dB_graph)
        # Pulse (E field)
        self.button_Pulse_E_field = QPushButton('\nPulse E field\n', self)
        self.button_Pulse_E_field.clicked.connect(self.Pulse_E_field_graph)


        self.button_Pulse_E_field_std = QPushButton('\nStd Pulse E field\n', self)
        self.button_Pulse_E_field_std.clicked.connect(self.Pulse_E_field_std_graph)
        # Pulse (E field) [dB] std
        self.button_E_field_dB_std = QPushButton('\nStd E field [dB]\n', self)
        self.button_E_field_dB_std.clicked.connect(self.E_field_std_dB_graph)
        
        #Phase
        self.button_Phase = QPushButton('\nPhase\n', self)
        self.button_Phase.clicked.connect(self.Phase_graph)

        
        #Parameters
        self.button_Correction_param = QPushButton('\nCorrection parameters\n', self)
        self.button_Correction_param.clicked.connect(self.Correction_param_graph)
        
        
        #Covariance
        self.button_Cov_Pulse_E_field= QPushButton('\nNoise matrices \u2074\n', self) # "noise matrix" replaced by "Noise matrices"
        self.button_Cov_Pulse_E_field.clicked.connect (self.Covariance_Pulse)
        
        self.label_window = QLabel("Window")  
        
        # ============================ select file to plot =========================== #
        self.plot_file_label = QLabel(f"choose file to plot") 
        self.plot_file_label.setMaximumSize(170,25)
        self.plot_file = QComboBox(self) 
        self.plot_file.activated.connect(self.plot_batch)

        # Organisation layout
        self.vlayoutmain = QVBoxLayout()
        self.hlayout = QHBoxLayout()
        self.hlayout2 = QHBoxLayout()
        self.hlayout3=QHBoxLayout()
        # ============================== add a combobox ============================== #
        self.hlayout4=QHBoxLayout()
        self.hlayout4.addWidget(self.plot_file_label)
        self.hlayout4.addWidget(self.plot_file)

        self.hlayout.addWidget(self.button_Pulse_E_field)
        self.hlayout.addWidget(self.button_E_field_dB)

        self.hlayout.addWidget(self.button_Pulse_E_field_std)
        self.hlayout.addWidget(self.button_E_field_dB_std)
        self.hlayout.addWidget(self.button_Phase)
        self.hlayout.addWidget(self.button_Correction_param)
        self.hlayout.addWidget(self.button_Cov_Pulse_E_field)
        
        self.label_window.setMaximumHeight(30)
        self.toolbar.setMaximumHeight(30)
        self.hlayout3.addWidget(self.toolbar)
        self.hlayout3.addWidget(self.label_window)
        #window_group = QGroupBox()
        #window_group.setLayout(self.hlayout3)
        #window_group.setMaximumWidth(100)
        #window_group.setMaximumHeight(60)


        self.vlayoutmain.addLayout(self.hlayout3)
        #self.vlayoutmain.addWidget(self.toolbar,1)
        #self.vlayoutmain.addWidget(window_group,Qt.AlignRight)

        self.vlayoutmain.addWidget(self.canvas)
        self.vlayoutmain.addLayout(self.hlayout)
        self.vlayoutmain.addLayout(self.hlayout2)
        # ========================== add combobox to layout ========================== #
        self.vlayoutmain.addLayout(self.hlayout4)
        # ---------------------------------------------------------------------------- #
        self.setLayout(self.vlayoutmain)


    def draw_graph_init(self,myinput, myreferencedata, ncm, ncm_inverse, ref_number, mydatacorrection, delay_correction, dilatation_correction, leftover_correction,
                        myglobalparameters, fopt, fopt_init, mode, preview):
        global graph_option_2
        self.figure.clf()
        
        nsample = len(myinput.pulse[-1])
        windows = np.ones(nsample)
        n_traces = len(myinput.pulse)

        if apply_window:
            windows = signal.tukey(nsample, alpha = 0.05)  #don't forget to modify it in fitc and opt files if it's modify here 

        if graph_option_2=='E_field [dB]':
            self.figure.clf()
            ax1 = self.figure.add_subplot(111)
            
            ax1.set_title('E_field', fontsize=10)
            
            color = 'tab:red'
            ax1.set_xlabel('Frequency [Hz]')
            ax1.set_ylabel('E_field [dB]',color=color)
            ax1.plot(myglobalparameters.freq,20*np.log(abs(TDS.torch_rfft(myinput.moyenne*windows)))/np.log(10), 'b-', label='mean spectre (log)')
            if not preview:
                ax1.plot(myglobalparameters.freq,20*np.log(abs(np.fft.rfft(myreferencedata.Pulseinit*windows)))/np.log(10), 'g-', label='reference spectre (log)')
                ax1.plot(myglobalparameters.freq,20*np.log(abs(np.fft.rfft(mydatacorrection.moyenne*windows)))/np.log(10), 'r-', label='corrected mean spectre (log)')

            if apply_window == 0:
                ax1.plot(myglobalparameters.freq,20*np.log(abs(myinput.freq_std/np.sqrt(n_traces)))/np.log(10), 'b--', label='Standard error (log)')
            else:
                ax1.plot(myglobalparameters.freq,20*np.log(abs(myinput.freq_std_with_window/np.sqrt(n_traces)))/np.log(10), 'b--', label='Standard error (log)')
            if not preview:
                if apply_window == 0:
                    ax1.plot(myglobalparameters.freq,20*np.log(abs(mydatacorrection.freq_std/np.sqrt(n_traces)))/np.log(10), 'r--', label='corrected Standard error (log)')
                else:
                    ax1.plot(myglobalparameters.freq,20*np.log(abs(mydatacorrection.freq_std_with_window/np.sqrt(n_traces)))/np.log(10), 'r--', label='corrected Standard error (log)')
        
            ax1.legend()
            ax1.grid()


        elif graph_option_2=='Pulse (E_field)':
            self.figure.clf()
            ax1 = self.figure.add_subplot(111)
            ax1.set_title('Pulse (E_field)', fontsize=10)
            
            color = 'tab:red'
            ax1.set_xlabel('Time [s]')
            ax1.set_ylabel('Amplitude',color=color)
            ax1.plot(myglobalparameters.t, myinput.moyenne*windows, 'b-', label='mean pulse')
            if not preview:
                ax1.plot(myglobalparameters.t, myreferencedata.Pulseinit*windows, 'g-', label='reference pulse')
                ax1.plot(myglobalparameters.t, mydatacorrection.moyenne*windows, 'r-', label='corrected mean pulse')
            ax1.legend()
            ax1.grid()
           
                
        elif graph_option_2=='Correction parameters':
            self.figure.clf()
            ax1 = self.figure.add_subplot(221)
            ax2 = self.figure.add_subplot(222)
            ax3 = self.figure.add_subplot(223)
            ax4 = self.figure.add_subplot(224)

            ax1.set_title('Delay', fontsize=10)
            ax2.set_title('Coef a - amplitude', fontsize=10)
            
            color = 'tab:red'
            ax1.set_xlabel("Trace index")
            ax1.set_ylabel('Delay [s]',color=color)
            
            ax2.set_xlabel("Trace index")
            ax2.set_ylabel('Coef a',color=color)
            
            ax3.set_title('alpha- dilatation', fontsize=10)
            ax4.set_title('beta - dilatation', fontsize=10)
            
            color = 'tab:red'
            ax3.set_xlabel("Trace index")
            ax3.set_ylabel('alpha',color=color)
            
            ax4.set_xlabel("Trace index")
            ax4.set_ylabel('beta',color=color)
            
            if not preview:
                if delay_correction and len(delay_correction)>1:
                    ax1.plot(delay_correction, "b.-")
                    ax1.plot(ref_number, delay_correction[ref_number], "rv")
                    ax1.annotate(" Reference", (ref_number,delay_correction[ref_number]))
                if leftover_correction and len(leftover_correction)>1:
                    ax2.plot([leftover_correction[i][0] for i in range(len(leftover_correction))], "b.-")
                    ax2.plot(ref_number, leftover_correction[ref_number][0], "rv")
                    ax2.annotate(" Reference", (ref_number,leftover_correction[ref_number][0]))
                if dilatation_correction and len(dilatation_correction)>1:
                    ax3.plot([dilatation_correction[i][0] for i in range(len(dilatation_correction))], "b.-")
                    ax3.plot(ref_number, dilatation_correction[ref_number][0], "rv")
                    ax3.annotate(" Reference", (ref_number,dilatation_correction[ref_number][0]))
                    
                    ax4.plot([dilatation_correction[i][1] for i in range(len(dilatation_correction))], "b.-")
                    ax4.plot(ref_number, dilatation_correction[ref_number][1], "rv")
                    ax4.annotate(" Reference", (ref_number,dilatation_correction[ref_number][1]))
                    
                ax1.grid()
                ax2.grid() 
                ax3.grid()
                ax4.grid()
            
                
        elif graph_option_2=='Std Pulse (E_field)':
            self.figure.clf()
            ax1 = self.figure.add_subplot(111)

            ax1.set_title('Std Pulse (E_field) ', fontsize=10)
            
            color = 'tab:red'
            ax1.set_xlabel('Time [s]')
            ax1.set_ylabel('Standard deviation Pulse (E_field)',color=color)
            ax1.plot(myglobalparameters.t, myinput.time_std*windows, 'b-', label='mean pulse')
            if not preview:
                ax1.plot(myglobalparameters.t, mydatacorrection.time_std*windows, 'r-', label='corrected mean pulse')
            ax1.legend()
            ax1.grid()
            
        elif graph_option_2 == 'Std E_field [dB]':
            self.figure.clf()
            ax1 = self.figure.add_subplot(111)
            ax1.set_title('Standard deviation E_field [dB]', fontsize=10)
            color = 'tab:red'
            ax1.set_xlabel('Frequency [Hz]')
            ax1.set_ylabel('Std E_field [dB]',color=color)
            if apply_window == 0:
                ax1.plot(myglobalparameters.freq,20*np.log(abs(myinput.freq_std))/np.log(10), 'b-', label='std spectre (log)')
            else:
                ax1.plot(myglobalparameters.freq,20*np.log(abs(myinput.freq_std_with_window))/np.log(10), 'b-', label='std spectre (log)')
            if not preview:
                if apply_window == 0:
                    ax1.plot(myglobalparameters.freq,20*np.log(abs(mydatacorrection.freq_std))/np.log(10), 'r-', label='corrected std spectre (log)')
                else:
                    ax1.plot(myglobalparameters.freq,20*np.log(abs(mydatacorrection.freq_std_with_window))/np.log(10), 'r-', label='corrected std spectre (log)')
            ax1.legend()
            ax1.grid()
            
            
        elif graph_option_2 == 'Phase':
            self.figure.clf()
            ax1 = self.figure.add_subplot(111)
            if mode == "basic":
                ax1.set_title('Phase ', fontsize=10)
            else:
                ax1.set_title('Phase', fontsize=10)
            color = 'tab:red'
            ax1.set_xlabel('Frequency [Hz]')
            ax1.set_ylabel('Phase (radians)',color=color)
            ax1.plot(myglobalparameters.freq,np.unwrap(np.angle(np.fft.rfft(myinput.moyenne*windows))), 'b-', label='mean phase')
            if not preview:
                ax1.plot(myglobalparameters.freq,np.unwrap(np.angle(np.fft.rfft(myreferencedata.Pulseinit*windows))), 'g-', label='reference phase')
                ax1.plot(myglobalparameters.freq,np.unwrap(np.angle(np.fft.rfft(mydatacorrection.moyenne*windows))), 'r-', label='corrected mean phase')
            ax1.legend()
            ax1.grid()
        
        elif graph_option_2 == "Noise matrices": # "noise matrix" replaced by "Noise matrices"
            self.figure.clf()
            ax1 = self.figure.add_subplot(121)
            ax2 = self.figure.add_subplot(122)
            
            if ncm is not None:
                ax1.set_title('Noise convolution matrix', fontsize=10)
                
                color = 'tab:red'
                ax1.set_xlabel('Time [s]')
                ax1.set_ylabel('Time [s]',color=color)
                border = max(abs(np.min(ncm)), abs(np.max(ncm)))
                im = ax1.imshow(ncm, vmin=-border, vmax = border, cmap = "seismic",interpolation= "nearest")
                plt.colorbar(im, ax = ax1)    
                
                ax2.set_title('Noise convolution matrix inverse', fontsize=10)
                
                color = 'tab:red'
                ax2.set_xlabel('Time [s]')
                ax2.set_ylabel('Time [s]',color=color)
                border = max(abs(np.min(ncm_inverse)), abs(np.max(ncm_inverse)))
                im = ax2.imshow(ncm_inverse, vmin=-border, vmax=border, cmap = "seismic",interpolation= "nearest")
                plt.colorbar(im, ax = ax2)  
                
            elif mydatacorrection.covariance is not None:
                ax1.set_title('Covariance matrix', fontsize=10)
                
                color = 'tab:red'
                ax1.set_xlabel('Time [s]')
                ax1.set_ylabel('Time [s]',color=color)
                border = max(abs(np.min(mydatacorrection.covariance)), abs(np.max(mydatacorrection.covariance)))
                im = ax1.imshow(mydatacorrection.covariance, vmin=-border, vmax = border, cmap = "seismic",interpolation= "nearest")
                plt.colorbar(im, ax = ax1) 
                
                ax2.set_title('Covariance matrix inverse', fontsize=10)
                
                color = 'tab:red'
                ax2.set_xlabel('Time [s]')
                ax2.set_ylabel('Time [s]',color=color)
                border = max(abs(np.min(mydatacorrection.covariance_inverse)), abs(np.max(mydatacorrection.covariance_inverse)))
                im = ax2.imshow(mydatacorrection.covariance_inverse, vmin=-border, vmax = border, cmap = "seismic",interpolation= "nearest")
                plt.colorbar(im, ax = ax2) 
                    
            elif myinput.covariance is not None:
                ax1.set_title('Covariance matrix', fontsize=10)
                
                color = 'tab:red'
                ax1.set_xlabel('Time [s]')
                ax1.set_ylabel('Time [s]',color=color)
                border = max(abs(np.min(myinput.covariance)), abs(np.max(myinput.covariance)))
                im = ax1.imshow(myinput.covariance, vmin=-border,  vmax = border, cmap = "seismic",interpolation= "nearest")
                plt.colorbar(im, ax = ax1) 

                ax2.set_title('Covariance matrix inverse', fontsize=10)
                
                color = 'tab:red'
                ax2.set_xlabel('Time [s]')
                ax2.set_ylabel('Time [s]',color=color)
                border = max(abs(np.min(myinput.covariance_inverse)), abs(np.max(myinput.covariance_inverse)))
                im = ax2.imshow(myinput.covariance_inverse, vmin=-border, vmax=border, cmap = "seismic",interpolation= "nearest")
                plt.colorbar(im, ax = ax2) 
            
        self.figure.tight_layout()
        self.canvas.draw()



    # def E_field_dB_graph(self):
        # global graph_option_2
        # graph_option_2='E_field [dB]'
        # self.controler.ploting_text3('Ploting E_field [dB]')

    # def Pulse_E_field_graph(self):
        # global graph_option_2
        # graph_option_2='Pulse (E_field)'
        # self.controler.ploting_text3('Ploting pulse E_field')
    
    # def Phase_graph(self):
        # global graph_option_2
        # graph_option_2='Phase'
        # self.controler.ploting_text3('Ploting Phase')
        
    # def Correction_param_graph(self):
        # global graph_option_2
        # graph_option_2='Correction parameters'
        # self.controler.ploting_text3('Ploting Correction parameters')
        
    # def Covariance_Pulse(self):
        # global graph_option_2
        # graph_option_2='Noise matrix'
        # self.controler.ploting_text3('Ploting Covariance\n   Reminder - the matrix must be saved for its computation (Ledoit-Wolf shrinkage)')  
        
    # def Pulse_E_field_std_graph(self):
        # global graph_option_2
        # graph_option_2='Std Pulse (E_field)'
        # self.controler.ploting_text3('Ploting Std pulse E_field')

    # def E_field_std_dB_graph(self):
        # global graph_option_2
        # graph_option_2='Std E_field [dB]'
        # self.controler.ploting_text3('Ploting Std E_field [dB]')
    
    def E_field_dB_graph(self):
        global graph_option_2
        graph_option_2='E_field [dB]'
        # self.controler.ploting_text3('Ploting E_field [dB]')
        self.parent.text_box.append(f'\nPlotting E_field [dB]')
        self.plot_batch()
        
    def Pulse_E_field_graph(self):
        global graph_option_2
        graph_option_2='Pulse (E_field)'
        # self.controler.ploting_text3('Ploting pulse E_field')
        self.parent.text_box.append(f'\nPlotting pulse E_field')
        self.plot_batch()
    
    def Phase_graph(self):
        global graph_option_2
        graph_option_2='Phase'
        # self.controler.ploting_text3('Ploting Phase')
        self.parent.text_box.append(f'\nPlotting Phase')
        self.plot_batch()
    
    def Correction_param_graph(self):
        global graph_option_2
        graph_option_2='Correction parameters'
        # self.controler.ploting_text3('Ploting Correction parameters')
        self.parent.text_box.append(f'\nPlotting Correction parameters')
        self.plot_batch()
    
    def Covariance_Pulse(self):
        global graph_option_2
        graph_option_2='Noise matrices' # "noise matrix" replaced by "Noise matrices"
        # self.controler.ploting_text3('Ploting Covariance\n   Reminder - the matrix must be saved for its computation (Ledoit-Wolf shrinkage)')
        self.parent.text_box.append(f'\nPlotting Covariance\n   Reminder - the matrix must be saved for its computation (Graphical lasso)')  
        self.plot_batch()
    
    def Pulse_E_field_std_graph(self):
        global graph_option_2
        graph_option_2='Std Pulse (E_field)'
        # self.controler.ploting_text3('Ploting Std pulse E_field')
        self.parent.text_box.append(f'\nPlotting Std pulse E_field')
        self.plot_batch()
        
    def E_field_std_dB_graph(self):
        global graph_option_2
        graph_option_2='Std E_field [dB]'
        # self.controler.ploting_text3('Ploting Std E_field [dB]')
        self.parent.text_box.append(f'\nPlotting Std E_field [dB]')
        self.plot_batch()

    # =========================== add itesm to combobox ========================== #
    def add_items_to_combobox(self,items):
        # items = [str(file) for file in csts.files] # items doit etre au format str
        # self.plot_file.addItems(items)
        
        # Ajout de l heure actuelle pour connaitre le file plot
        current_time = datetime.now().strftime("%H:%M:%S")
        items = [f"{str(os.path.basename(file))} ; Added at {current_time}" for file in csts.files]
        self.plot_file.addItems(items)
        
        # Selectionner le dernier element ajoute
        last_index = self.plot_file.count() - 1  # Index du dernier element
        self.plot_file.setCurrentIndex(last_index)
            
        
        
    # ========================= defing file choice signal ======================== #
    def refresh_plot(self):
        self.plot_batch()
    
    # =========== define a new plot function according to selected file ========== #
    def plot_batch(self):
        index = self.plot_file.currentIndex()
        try:
            self.draw_graph_init(csts.myinput[index], csts.myreferencedata[index], csts.ncm[index], csts.ncm_inverse[index], csts.reference_number[index], csts.mydatacorrection[index], csts.delay_correction[index], csts.dilatation_correction[index],csts.leftover_correction[index],csts.myglobalparameters[index],csts.fopt[index], csts.fopt_init[index], csts.mode[index], preview)
        except Exception as e:
            print(e)
    # ---------------------------------------------------------------------------- #
    
    def refresh(self):
        try:
            self.draw_graph_init(self.controler.myinput, self.controler.myreferencedata, self.controler.ncm, self.controler.ncm_inverse, self.controler.reference_number, self.controler.mydatacorrection, self.controler.delay_correction, self.controler.dilatation_correction,self.controler.leftover_correction,self.controler.myglobalparameters,self.controler.fopt, self.controler.fopt_init, self.controler.mode, preview)
                
        except Exception as e:
            pass
        """except:
            print("There is a refresh problem")
            pass"""





###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class MainWindow(QMainWindow):
    def __init__(self, controler):
        super().__init__()
        self.setWindowTitle("Correct@TDS")
        self.mainwidget = MyTableWidget(self,controler)
        self.setCentralWidget(self.mainwidget)
        self.setWindowIcon(QtGui.QIcon('icon.png'))

    def closeEvent(self,event):
        try:
            shutil.rmtree("temp")
        except:
            pass
    
    
def main():

    sys._excepthook = sys.excepthook 
    def exception_hook(exctype, value, traceback):
        print(exctype, value, traceback)
        sys._excepthook(exctype, value, traceback) 
        sys.exit(1) 
    sys.excepthook = exception_hook 
    
    app = QApplication(sys.argv)
    controler = Controler()    
    win = MainWindow(controler)
    qApp.setApplicationName("Correct@TDS")
    #controler.init()
    # win.show()
    win.showMaximized()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
    
