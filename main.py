# %%
"""                                              Import                                                              """
# !/usr/bin/python3
import os
import subprocess
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import scrolledtext
from PIL import Image, ImageTk
import numpy as np
from threading import Thread
import traceback
from inspect import signature
from inspect import getfullargspec
from Source.FittingRoutine import FittingRoutine

import ROOT
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
import imageio

major = 3.5
minor = 1
width = 0.75
plt.rc('text', usetex=True)
plt.rc("axes", labelsize=8)  # 18
plt.rc("xtick", labelsize=6, top=True, direction="in")
plt.rc("ytick", labelsize=6, right=True, direction="in")
plt.rc("axes", titlesize=10)
plt.rc("legend", fontsize=8)
plt.rcParams['font.family'] = "serif"
plt.rcParams['axes.linewidth'] = width
plt.rcParams['xtick.minor.width'] = width
plt.rcParams['xtick.major.width'] = width
plt.rcParams['ytick.minor.width'] = width
plt.rcParams['ytick.major.width'] = width
plt.rcParams['xtick.major.size'] = major
plt.rcParams['xtick.minor.size'] = minor
plt.rcParams['ytick.major.size'] = major
plt.rcParams['ytick.minor.size'] = minor
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=["b", "r", "k", "grey", "magenta", "b", "r", "k", "grey", "magenta"],
                                             ls=["-", "-", "-", "-", "-", "--", "--", "--", "--", "--"])

# %%
"""                                            Description                                                           """

# %%
"""                                             Functions                                                            """


def load_to_root(load_to_root_exe_path, document_name_path):
    output = subprocess.Popen(load_to_root_exe_path.replace(" ", "\\ ") + " " + document_name_path.replace(" ", "\\ "),
                              shell=True, stdout=subprocess.PIPE).communicate()[0]

    if '{}'.format(output)[-4] == '0':
        raise EnvironmentError


def get_all_time_data(get_all_time_data_exe_path, root_file_path, args=None):
    call = get_all_time_data_exe_path.replace(" ", "\\ ") + "/build/./Main" + " " + root_file_path.replace(" ",
                                                                                                           "\\ ") + " "
    if args != None:
        for arg in args:
            call += arg
    os.system(call)


def VMI_filter(VMI_filter_exe_path, root_file_path, args):
    call = VMI_filter_exe_path.replace(" ", "\\ ") + "/build/./Main" + " " + root_file_path.replace(" ", "\\ ") + " "
    for arg in args:
        call += arg
    os.system(call)


def VMI_filter_boot(root, VMI_filter_exe_path, root_file_path, nr_samples, nr_threads='1'):
    call = VMI_filter_exe_path.replace(" ", "\\ ") + "/build/./Main" + " " + root_file_path.replace(" ", "\\ ") + \
           " -b " + nr_samples.replace(" ", "\\ ") + '-n ' + nr_threads.replace(" ", "\\ ")

    p = subprocess.Popen(call, stdout=subprocess.PIPE, shell=True)
    for line in iter(p.stdout.readline, b''):
        root.boot_out.insert(tk.END, line)
        print(line)
    p.stdout.close()
    p.wait()


def validate_int(action, index, value_if_allowed,
                 prior_value, text, validation_type, trigger_type, widget_name):
    # action=1 -> insert
    if action == '1':
        if text in '0123456789':
            try:
                int(value_if_allowed)
                return True
            except ValueError:
                return False
        else:
            return False
    else:
        return True


def validate_float(action, index, value_if_allowed,
                   prior_value, text, validation_type, trigger_type, widget_name):
    # action=1 -> insert
    if action == '1':
        if text in '0123456789.e-':
            try:
                float(value_if_allowed)
                return True
            except ValueError:
                return False
        else:
            return False
    else:
        return True


def search_for_file_path(root, title):
    current_dir = os.getcwd()
    file = filedialog.askopenfilename(parent=root, initialdir=current_dir, title=title)
    return file


def search_for_dir_path(root, title):
    current_dir = os.getcwd()
    file = filedialog.askdirectory(parent=root, initialdir=current_dir, title=title)
    return file


def search_for_dat_file(root, label):
    root.dat_file_path = search_for_file_path(root, 'Please select data file')
    if root.home_dir in root.dat_file_path:
        label.config(text=root.dat_file_path[len(root.home_dir):])
    else:
        label.config(text=root.dat_file_path)


def search_for_load_file(root, label):
    directory = search_for_dir_path(root, 'Please select \'Load data\' directory')
    root.load_file_path = directory
    if root.home_dir in directory:
        label.config(text=directory[len(root.home_dir):])
    else:
        label.config(text=directory)


def search_for_event_file(root, label):
    directory = search_for_dir_path(root, 'Please select \'Event Trig data\' directory')
    root.event_file_path = directory
    if root.home_dir in directory:
        label.config(text=directory[len(root.home_dir):])
    else:
        label.config(text=directory)


def search_for_vmi_file(root, label):
    directory = search_for_dir_path(root, 'Please select \'VMI filter\' directory')
    root.vmi_file_path = directory
    if root.home_dir in directory:
        label.config(text=directory[len(root.home_dir):])
    else:
        label.config(text=directory)


def search_for_gif_file(root, label):
    directory = search_for_dir_path(root, 'Please select \'GIF function\' directory')
    root.gif_file_path = directory
    if root.home_dir in directory:
        label.config(text=directory[len(root.home_dir):])
    else:
        label.config(text=directory)


def search_for_boot_file(root, label):
    directory = search_for_dir_path(root, 'Please select \'Bootstrap function\' directory')
    root.boot_file_path = directory
    if root.home_dir in directory:
        label.config(text=directory[len(root.home_dir):])
    else:
        label.config(text=directory)


def path_check(root, error, error_label):
    if os.path.exists(root.boot_file_path + "/build/Main"):
        if os.path.exists(root.gif_file_path + "/build/Main"):
            if os.path.exists(root.load_file_path + "/build/Main"):
                if os.path.exists(root.event_file_path + "/build/Main"):
                    if os.path.exists(root.vmi_file_path + "/build/Main"):
                        if os.path.exists(root.dat_file_path + ".root"):
                            return True
                        else:
                            try:
                                load_to_root(root.load_file_path + "/build/./Main", root.dat_file_path)
                                get_all_time_data(root.event_file_path, root.dat_file_path + ".root")
                                VMI_filter(root.vmi_file_path, root.dat_file_path + ".root",
                                           ['-r 200', '-o 3', '-v 250', '-c 100'])
                                return True
                            except EnvironmentError:
                                error.insert(tk.END, root.dat_file_path + " : can not be loaded" + "\n")
                    else:
                        error.insert(tk.END, root.vmi_file_path + " : does not include /build/Main" + "\n")
                else:
                    error.insert(tk.END, root.event_file_path + " : does not include /build/Main" + "\n")
            else:
                error.insert(tk.END, root.load_file_path + " : does not include /build/Main" + "\n")
        else:
            error.insert(tk.END, root.gif_file_path + " : does not include /build/Main" + "\n")
    else:
        error.insert(tk.END, root.boot_file_path + " : does not include /build/Main" + "\n")

    error_label.grid(row=8, column=1)
    error.grid(row=9, column=1)
    return False


def confirm_close(root, frame, error, error_label):
    if path_check(root, error, error_label):
        frame.destroy()


# %%
"""                                              Classes                                                             """


class ControlWindow(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self, "Data Analysis Program")
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        geometry = "%ix%i" % (screen_width, screen_height)
        self.geometry(geometry)

        '''             Initial window            '''
        # Initial attempt to set paths
        self.home_dir = os.path.expanduser('~')
        self.current_dir = os.path.dirname(os.path.abspath(__file__))
        self.dat_file_path = self.current_dir
        if os.path.exists(self.current_dir + '/Source/Load data/build/./Main'):
            self.load_file_path = self.current_dir + '/Source/Load data'
        else:
            self.load_file_path = self.current_dir
        if os.path.exists(self.current_dir + '/Source/Event Trig data/build/./Main'):
            self.event_file_path = self.current_dir + '/Source/Event Trig data'
        else:
            self.event_file_path = self.current_dir
        if os.path.exists(self.current_dir + '/Source/VMI filter'):
            self.vmi_file_path = self.current_dir + '/Source/VMI filter'
        else:
            self.vmi_file_path = self.current_dir
        if os.path.exists(self.current_dir + '/Source/VMI GIF'):
            self.gif_file_path = self.current_dir + '/Source/VMI GIF'
        else:
            self.boot_file_path = self.current_dir
        if os.path.exists(self.current_dir + '/Source/VMI Boot'):
            self.boot_file_path = self.current_dir + '/Source/VMI Boot'
        else:
            self.boot_file_path = self.current_dir

        # This code lets you set the needed directories
        self.withdraw()
        load_popup = FilePopup(self)
        self.wait_window(load_popup)
        self.deiconify()

        self.root_file = self.dat_file_path + ".root"

        """             Initial data            """
        self.vmi_resolution = 1
        self.resolution = 200
        self.norder = '4'
        self.abel_invers = False

        self.event_hist = []
        self.event_hist_err = []
        self.event_current_hist = []
        self.event_current_hist_err = []

        self.trig_hist = []
        self.trig_hist_err = []
        self.trig_current_hist = []
        self.trig_current_hist_err = []

        self.peak_sum_hist = []
        self.peak_sum_hist_err = []
        self.peak_sum_current_hist = []
        self.peak_sum_current_hist_err = []

        self.freq_hist = []
        self.freq_hist_err = []

        self.event_time_bin_size = 0.5e-3
        self.trig_time_bin_size = 0.5e-8 / 1e-6  # in micro seconds
        self.peak_sum_bin_size = 20
        self.freq_bin_size = 0.1e-4

        self.background_sub_mode = [1, 0]
        self.event_background_lower = 0
        self.event_background_upper = 1

        self.sys_radial_err = 0.02 * 20

        self.event_time_bins = [0]
        self.trig_time_bins = [0]
        self.peak_sum_bins = [0]
        self.freq_bins = [0]

        self.vmi_laser_on = np.zeros((1, 1))
        self.vmi_laser_off = np.zeros((1, 1))
        self.vmi_fold = np.zeros((1, 1))

        self.radial_projection_values = []
        self.radial_projection_amplitude = []
        self.radial_projection_error = []
        self.radial_projection_amplitude_sum = []
        self.radial_projection_error_sum = []
        self.max_rad = '0'

        self.abel_vmi = np.zeros((1, 1))
        self.combined_vmi = np.zeros((1, 1))
        self.circ_vmi = np.zeros((1, 1))
        self.circularize = False
        self.nr_angular_bin = 100
        self.intensity = []
        self.intensity_err = []

        self.boot_iterations = []
        self.boot_intensity = []
        self.boot_radial = []
        self.intensity_err_lower = []
        self.intensity_err_upper = []
        self.quantiles = []
        self.boot_rad_err = []
        self.boot_rad_sys_err = []
        self.boot_intensity_sys_err = []

        self.center_vals = []

        self.nr_samples = 0
        self.nr_threads = 1
        self.boot_out_str = tk.StringVar()
        self.nr_bootstraps = 1

        self.initialise()

        self.raw_vmi_frames = []
        self.raw_off_vmi_frames = []
        self.raw_circ_vmi_frames = []
        self.inv_vmi_frames = []
        self.intensity_frames = []
        self.radial_frames = []
        self.event_frames = []
        self.circ_frames = False
        self.gif_times = []

        self.current_frame = 0
        self.max_frames = 0
        self.run_gif = False

        self.radius_to_energy_func = {}

        '''             Figures                 '''
        self.event_figure = plt.Figure((2, 1), dpi=200)
        self.event_ax = self.event_figure.subplots(1)

        self.trig_figure = plt.Figure((2, 1), dpi=200)
        self.trig_ax = self.trig_figure.subplots(1)

        self.peak_sum_figure = plt.Figure((2, 1), dpi=200)
        self.peak_sum_ax = self.peak_sum_figure.subplots(1)

        self.laser_on_figure = plt.Figure((2, 1), dpi=200)
        self.laser_on_ax = self.laser_on_figure.subplots(1)
        self.cb_laser_on = self.laser_on_figure.colorbar(self.laser_on_ax.imshow(self.vmi_laser_on))

        self.laser_off_figure = plt.Figure((2, 1), dpi=200)
        self.laser_off_ax = self.laser_off_figure.subplots(1)
        self.cb_laser_off = self.laser_off_figure.colorbar(self.laser_off_ax.imshow(self.vmi_laser_off))

        self.rad_figure = plt.Figure((2, 1), dpi=200)
        self.rad_ax = self.rad_figure.subplots(1)
        self.rad_ax_energy = {}

        self.intensity_figure = plt.Figure((2, 1), dpi=200)
        self.intensity_ax = self.intensity_figure.subplots(1)

        self.boot_intensity_figure = plt.Figure((2, 1), dpi=200)
        self.boot_intensity_ax = self.boot_intensity_figure.subplots(1)

        self.boot_statistic_figure = plt.Figure((2, 1), dpi=200)
        self.boot_statistic_ax = self.boot_statistic_figure.subplots(1)

        self.abel_vmi_figure = plt.Figure((2, 1), dpi=200)
        self.abel_vmi_ax = self.abel_vmi_figure.subplots(1)
        self.cb_abel_vmi = self.abel_vmi_figure.colorbar(self.abel_vmi_ax.imshow(self.abel_vmi))

        self.circ_vmi_figure = plt.Figure((2, 1), dpi=200)
        self.circ_vmi_ax = self.circ_vmi_figure.subplots(1)
        self.cb_circ_vmi = self.circ_vmi_figure.colorbar(self.circ_vmi_ax.imshow(self.circ_vmi))

        self.cali_raw_figure, self.cali_raw_ax = plt.subplots(2, 1, sharex=True, figsize=(6, 6),
                                                              gridspec_kw={'height_ratios': [2, 1]})

        self.laser_hist_figure, self.laser_hist_ax = plt.subplots(2, 1, sharex=True, figsize=(6, 6),
                                                                  gridspec_kw={'height_ratios': [2, 1]})

        self.cali_plot_figure, self.cali_plot_ax = plt.subplots(2, 1, sharex=True, figsize=(6, 6),
                                                                gridspec_kw={'height_ratios': [2, 1]})

        '''                 Frames                '''
        Command = tk.LabelFrame(self, text="Controls", width="270", bg="yellow")
        Command.grid(row=0, column=0, rowspan=2, sticky=tk.W + tk.N + tk.S)
        Command.grid_propagate(False)
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=0)

        tab_bar = ttk.Notebook(self)

        # Raw data page

        Rawdata = ttk.Frame(tab_bar)

        # Tab
        TriggerPeakTab = ttk.Notebook(Rawdata)
        TriggerPeakTab.grid(row=0, column=0, sticky=tk.W + tk.N + tk.S + tk.E)
        Rawdata.grid_rowconfigure(0, weight=1)
        Rawdata.grid_columnconfigure(0, weight=1)
        Trigger = ttk.Frame(TriggerPeakTab)
        PeakSum = ttk.Frame(TriggerPeakTab)

        Event = ttk.Frame(Rawdata)
        Event.grid(row=0, column=1, sticky=tk.W + tk.N + tk.S + tk.E)
        Rawdata.grid_rowconfigure(0, weight=1)
        Rawdata.grid_columnconfigure(1, weight=1)

        VMI_on = ttk.Frame(Rawdata)
        VMI_on.grid(row=1, column=0, sticky=tk.W + tk.N + tk.S + tk.E)
        Rawdata.grid_rowconfigure(1, weight=1)
        Rawdata.grid_columnconfigure(0, weight=1)

        VMI_off = ttk.Frame(Rawdata)
        VMI_off.grid(row=1, column=1, sticky=tk.W + tk.N + tk.S + tk.E)
        Rawdata.grid_rowconfigure(1, weight=1)
        Rawdata.grid_columnconfigure(1, weight=1)

        # Analysed data page

        Analyseddata = ttk.Frame(tab_bar)

        Radial = tk.Frame(Analyseddata)
        Radial.grid(row=0, column=0, sticky=tk.W + tk.N + tk.S + tk.E)
        Analyseddata.grid_rowconfigure(0, weight=1, uniform='group1')
        Analyseddata.grid_columnconfigure(0, weight=1, uniform='group1')

        Intensity = tk.Frame(Analyseddata)
        Intensity.grid(row=0, column=1, sticky=tk.W + tk.N + tk.S + tk.E)
        Analyseddata.grid_rowconfigure(0, weight=1, uniform='group1')
        Analyseddata.grid_columnconfigure(1, weight=1, uniform='group1')

        AbelVMI = tk.Frame(Analyseddata)
        AbelVMI.grid(row=1, column=0, sticky=tk.W + tk.N + tk.S + tk.E)
        Analyseddata.grid_rowconfigure(1, weight=1, uniform='group1')
        Analyseddata.grid_columnconfigure(0, weight=1, uniform='group1')

        Boot = ttk.Notebook(Analyseddata)
        Boot.grid(row=1, column=1, sticky=tk.W + tk.N + tk.S + tk.E)
        Analyseddata.grid_rowconfigure(1, weight=1, uniform='group1')
        Analyseddata.grid_columnconfigure(1, weight=1, uniform='group1')

        Boot_plot = ttk.Frame(Boot)
        Boot_plot.pack(fill='both', expand=True)
        Boot.add(Boot_plot, text='Collected')

        Boot_hist = ttk.Frame(Boot)
        Boot_hist.pack(fill='both', expand=True)
        Boot.add(Boot_hist, text='Statistics')

        # Time resolved data page

        self.Timeresolved = ttk.Frame(tab_bar)

        # Event and controls
        GIFNotebook = tk.Frame(self.Timeresolved)

        GIF_tabs = ttk.Notebook(self.Timeresolved)

        self.GIFEvent = ttk.Frame(GIF_tabs)
        GIF_tabs.add(self.GIFEvent, text='Time Plot')

        self.GIFControls = ttk.Frame(GIF_tabs)
        GIF_tabs.add(self.GIFControls, text='Controls')

        GIF_tabs.grid(row=0, column=0, sticky=tk.W + tk.N + tk.S + tk.E)
        self.Timeresolved.grid_rowconfigure(0, weight=1, uniform='group1')
        self.Timeresolved.grid_columnconfigure(0, weight=1, uniform='group1')

        GIFNotebook.grid(row=0, column=0, sticky=tk.W + tk.N + tk.S + tk.E)
        self.Timeresolved.grid_rowconfigure(0, weight=1, uniform='group1')
        self.Timeresolved.grid_columnconfigure(0, weight=1, uniform='group1')

        # Intensity and radial
        GIFIntensityNotebook = tk.Frame(self.Timeresolved)

        GIFIntensity_tabs = ttk.Notebook(self.Timeresolved)

        self.GIFIntensity = ttk.Frame(GIFIntensity_tabs)
        GIFIntensity_tabs.add(self.GIFIntensity, text='Abel Inverse')

        self.GIFRadial = ttk.Frame(GIFIntensity_tabs)
        GIFIntensity_tabs.add(self.GIFRadial, text='Radial binning')

        GIFIntensity_tabs.grid(row=1, column=0, sticky=tk.W + tk.N + tk.S + tk.E)
        self.Timeresolved.grid_rowconfigure(1, weight=1, uniform='group1')
        self.Timeresolved.grid_columnconfigure(0, weight=1, uniform='group1')

        GIFIntensityNotebook.grid(row=1, column=0, sticky=tk.W + tk.N + tk.S + tk.E)
        self.Timeresolved.grid_rowconfigure(1, weight=1, uniform='group1')
        self.Timeresolved.grid_columnconfigure(0, weight=1, uniform='group1')

        # Inverse abel

        self.GIFAbelVMI = tk.Frame(self.Timeresolved)
        self.GIFAbelVMI.grid(row=0, column=1, sticky=tk.W + tk.N + tk.S + tk.E)
        self.Timeresolved.grid_rowconfigure(0, weight=1, uniform='group1')
        self.Timeresolved.grid_columnconfigure(1, weight=1, uniform='group1')

        # VMI images

        GIFVMINotebook = tk.Frame(self.Timeresolved)

        GIFVMI_tabs = ttk.Notebook(self.Timeresolved)

        self.GIFVMI = ttk.Frame(GIFVMI_tabs)
        GIFVMI_tabs.add(self.GIFVMI, text='VMI On')

        self.GIFVMIOff = ttk.Frame(GIFVMI_tabs)
        GIFVMI_tabs.add(self.GIFVMIOff, text='VMI Off')

        self.GIFVMICirc = ttk.Frame(GIFVMI_tabs)
        GIFVMI_tabs.add(self.GIFVMICirc, text='VMI Circ')

        GIFVMI_tabs.grid(row=1, column=1, sticky=tk.W + tk.N + tk.S + tk.E)
        self.Timeresolved.grid_rowconfigure(1, weight=1, uniform='group1')
        self.Timeresolved.grid_columnconfigure(1, weight=1, uniform='group1')

        GIFVMINotebook.grid(row=1, column=1, sticky=tk.W + tk.N + tk.S + tk.E)
        self.Timeresolved.grid_rowconfigure(1, weight=1, uniform='group1')
        self.Timeresolved.grid_columnconfigure(1, weight=1, uniform='group1')

        # Calibration page

        # Setting up tabs
        Calibration = ttk.Frame(tab_bar)

        CalibrationTabs = ttk.Notebook(Calibration)
        CalibrationTabs.pack(fill='both', expand=True)

        DataFitTab = ttk.Frame(CalibrationTabs)
        LaserFitTab = ttk.Frame(CalibrationTabs)
        CalibrationFitTab = ttk.Frame(CalibrationTabs)

        # Frames
        self.DataFitTabPlot = ttk.Frame(DataFitTab)
        self.DataFitTabPlot.grid(row=0, column=1, sticky=tk.NSEW)
        DataFitTabSettings = ttk.Frame(DataFitTab)
        DataFitTabSettings.grid(row=0, column=0, sticky=tk.NSEW)
        DataFitTab.grid_rowconfigure(0, weight=1, uniform='group1')
        DataFitTab.grid_columnconfigure(0, weight=1, uniform='group1')
        DataFitTab.grid_columnconfigure(1, weight=4, uniform='group1')

        self.LaserFitTabPlot = ttk.Frame(LaserFitTab)
        self.LaserFitTabPlot.grid(row=0, column=1, sticky=tk.NSEW)
        LaserFitTabSettings = ttk.Frame(LaserFitTab)
        LaserFitTabSettings.grid(row=0, column=0, sticky=tk.NSEW)
        LaserFitTab.grid_rowconfigure(0, weight=1, uniform='group1')
        LaserFitTab.grid_columnconfigure(0, weight=1, uniform='group1')
        LaserFitTab.grid_columnconfigure(1, weight=4, uniform='group1')

        self.CalibrationFitTabPlot = ttk.Frame(CalibrationFitTab)
        self.CalibrationFitTabPlot.grid(row=0, column=1, sticky=tk.NSEW)
        CalibrationFitTabSettings = ttk.Frame(CalibrationFitTab)
        CalibrationFitTabSettings.grid(row=0, column=0, sticky=tk.NSEW)
        CalibrationFitTab.grid_rowconfigure(0, weight=1, uniform='group1')
        CalibrationFitTab.grid_columnconfigure(0, weight=1, uniform='group1')
        CalibrationFitTab.grid_columnconfigure(1, weight=4, uniform='group1')

        # Adding tabs
        CalibrationTabs.add(DataFitTab, text='Current Data')
        CalibrationTabs.add(LaserFitTab, text='Laser Data')
        CalibrationTabs.add(CalibrationFitTab, text='Calibration Fit')

        # Adding the tabs to the interface

        tab_bar.add(Rawdata, text="Raw Data")
        tab_bar.add(Analyseddata, text="Analysis")
        tab_bar.add(self.Timeresolved, text="Time Resolved")
        tab_bar.add(Calibration, text="Calibration")

        tab_bar.grid(row=0, column=1, columnspan=2, rowspan=2, sticky=tk.W + tk.N + tk.S + tk.E)
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

        self.frames = {}

        '''           Control frame           '''
        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        vcmdi = (self.register(validate_int),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        # Text field to write in the resolution
        resolution_label = tk.Label(Command, text="Resolution", bg="yellow")
        self.resolution_str_var = tk.StringVar()
        self.resolution_str_var.set(str(self.resolution))
        self.resolution_entry = tk.Entry(Command, validate='key', validatecommand=vcmdi,
                                         width=5, textvariable=self.resolution_str_var)
        self.resolution = self.resolution_entry.get()
        resolution_max_label = tk.Label(Command, text='Max Radial', bg="yellow")
        self.resolution_max_str_var = tk.StringVar()
        self.resolution_max_str_var.set(self.max_rad)
        self.resolution_max_entry = tk.Entry(Command, validate='key', validatecommand=vcmdf,
                                             width=5, textvariable=self.resolution_max_str_var)

        # Text field to write in the angular order
        order_label = tk.Label(Command, text="Angular order", bg="yellow")
        self.order_str_var = tk.StringVar()
        self.order_str_var.set(self.norder)
        self.norder_entry = tk.Entry(Command, validate='key', validatecommand=vcmdi,
                                     width=5, textvariable=self.order_str_var)
        self.norder = self.norder_entry.get()

        # Text field to write in the vmi resolution
        vmi_resolution_label = tk.Label(Command, text="VMI Resolution", bg="yellow")
        self.vmi_resolution_str_var = tk.StringVar()
        self.vmi_resolution_str_var.set(str(self.vmi_resolution))
        self.vmi_resolution_entry = tk.Entry(Command, validate='key', validatecommand=vcmdi,
                                             width=5, textvariable=self.vmi_resolution_str_var)
        self.vmi_resolution = self.vmi_resolution_entry.get()

        # Button to accept settings
        settings_button = tk.Button(Command, text="Set", command=lambda: self.set_settings())

        # Text fields to write trigger time
        trig_time_label = tk.Label(Command, text="Trigger time interval", bg="yellow")
        self.trig_lower_str_var = tk.StringVar()
        self.trig_upper_str_var = tk.StringVar()
        self.trig_lower_str_var.set(str("%.3f" % self.trig_lower))
        self.trig_upper_str_var.set(str("%.3f" % self.trig_upper))
        self.trig_lower_entry = tk.Entry(Command, validate='key', validatecommand=vcmdf, width=9,
                                         textvariable=self.trig_lower_str_var)
        self.trig_upper_entry = tk.Entry(Command, validate='key', validatecommand=vcmdf,
                                         width=11, textvariable=self.trig_upper_str_var)
        self.trig_upper = self.trig_upper_entry.get()
        self.trig_lower = self.trig_lower_entry.get()

        # Text fields to write event time
        event_time_label = tk.Label(Command, text="Event time interval", bg="yellow")
        self.event_lower_str_var = tk.StringVar()
        self.event_upper_str_var = tk.StringVar()
        self.event_lower_str_var.set(str("%.3f" % self.event_lower))
        self.event_upper_str_var.set(str("%.3f" % self.event_upper))
        self.event_lower_entry = tk.Entry(Command, validate='key', validatecommand=vcmdf, width=9,
                                          textvariable=self.event_lower_str_var)
        self.event_upper_entry = tk.Entry(Command, validate='key', validatecommand=vcmdf,
                                          width=11, textvariable=self.event_upper_str_var)
        self.event_upper = self.event_upper_entry.get()
        self.event_lower = self.event_lower_entry.get()

        # Text fields to write peak sum
        peak_sum_label = tk.Label(Command, text="Peak sum interval", bg="yellow")
        self.peak_sum_lower_str_var = tk.StringVar()
        self.peak_sum_upper_str_var = tk.StringVar()
        self.peak_sum_lower_str_var.set(str("%.3f" % self.peak_sum_lower))
        self.peak_sum_upper_str_var.set(str("%.3f" % self.peak_sum_upper))
        self.peak_sum_lower_entry = tk.Entry(Command, validate='key', validatecommand=vcmdf, width=9,
                                             textvariable=self.peak_sum_lower_str_var)
        self.peak_sum_upper_entry = tk.Entry(Command, validate='key', validatecommand=vcmdf,
                                             width=11, textvariable=self.peak_sum_upper_str_var)
        self.peak_sum_upper = self.peak_sum_upper_entry.get()
        self.peak_sum_lower = self.peak_sum_lower_entry.get()

        # Button to accept times
        times_button = tk.Button(Command, text="Set\nTimes", command=lambda: self.set_times(), width=7)
        # Button to check times
        times_check_button = tk.Button(Command, text="Check\nTimes", command=lambda: self.check_times(), width=7)

        # Text field to write in nr of angular bins in circulation
        circular_label = tk.Label(Command, text="Circularization", bg="yellow")
        self.circular_str_var = tk.StringVar()
        self.circular_str_var.set(str(self.nr_angular_bin))
        self.circular_entry = tk.Entry(Command, validate='key', validatecommand=vcmdi,
                                       width=5, textvariable=self.circular_str_var)
        self.nr_angular_bin = self.circular_entry.get()
        # Button to run circularization
        circular_set_button = tk.Button(Command, text='Set', command=lambda: self.set_circular())
        # Button to change bool of circularization
        self.circular_str_truth = tk.StringVar()
        if self.circularize:
            self.circular_str_truth.set('True')
        else:
            self.circular_str_truth.set('False')
        circular_button = tk.Button(Command, textvariable=self.circular_str_truth,
                                    command=lambda: self.change_circular())

        # Button to open GIF window
        boot_button = tk.Button(Command, text='Bootstrap settings', command=lambda: self.calc_boot())

        # Button to open GIF window
        gif_button = tk.Button(Command, text='GIF settings', command=lambda: self.calc_gif())

        # Imoprt/Export of circularization
        CircImportExport = tk.LabelFrame(Command, text='Imoprt/Export Circularization', bg='yellow')
        circ_import_button = tk.Button(CircImportExport, text='Import', command=lambda: self.import_circ())
        circ_export_button = tk.Button(CircImportExport, text='Export', command=lambda: self.export_circ())
        self.circ_import_label_var = tk.StringVar()
        self.circ_import_file = ''
        circ_import_label = tk.Label(CircImportExport, textvariable=self.circ_import_label_var, bg='yellow')
        circ_use_import_label = tk.Label(CircImportExport, text='Use import ', bg='yellow')
        self.circ_use_import_button_var = tk.StringVar()
        self.circ_use_import_button_var.set('False')
        circ_use_import_button = tk.Button(CircImportExport, textvariable=self.circ_use_import_button_var,
                                           command=lambda: self.use_import())
        circ_import_button.grid(row=0, column=0)
        circ_export_button.grid(row=0, column=1)
        circ_import_label.grid(row=2, column=0, columnspan=2)
        circ_use_import_label.grid(row=1, column=0)
        circ_use_import_button.grid(row=1, column=1)

        # Imoprt/Export of calibration
        CaliImportExport = tk.LabelFrame(Command, text='Imoprt Energy Calibration', bg='yellow')
        cali_import_button = tk.Button(CaliImportExport, text='Import', command=lambda: self.import_cali_fit())
        self.cali_import_label_var = tk.StringVar()
        self.cali_import_file = ''
        self.imp_cali_param = []
        self.imp_cali_err = []
        cali_import_label = tk.Label(CaliImportExport, textvariable=self.cali_import_label_var, bg='yellow')
        cali_use_import_label = tk.Label(CaliImportExport, text='Use import ', bg='yellow')
        self.cali_use_import_button_var = tk.StringVar()
        self.cali_use_import_button_var.set('False')
        cali_use_import_button = tk.Button(CaliImportExport, textvariable=self.cali_use_import_button_var,
                                           command=lambda: self.use_calibration_import())
        cali_converter_poppup = tk.Button(CaliImportExport, text='Energy Converter',
                                          command=lambda: self.cali_converter_popup())
        cali_import_button.grid(row=0, column=0)
        cali_import_label.grid(row=2, column=0, columnspan=2)
        cali_use_import_label.grid(row=1, column=0)
        cali_use_import_button.grid(row=1, column=1)
        cali_converter_poppup.grid(row=3, column=0, columnspan=2)

        # Drawing widgets
        resolution_label.grid(row=0, column=0, columnspan=3)
        self.resolution_entry.grid(row=1, column=1, sticky=tk.W + tk.E)
        resolution_max_label.grid(row=0, column=0, sticky=tk.W)
        self.resolution_max_entry.grid(row=1, column=0, sticky=tk.W)

        order_label.grid(row=2, column=0, columnspan=3)
        self.norder_entry.grid(row=3, column=1, sticky=tk.W + tk.E)

        vmi_resolution_label.grid(row=4, column=0, columnspan=3)
        self.vmi_resolution_entry.grid(row=5, column=1, sticky=tk.W + tk.E)

        circular_label.grid(row=6, column=0, columnspan=3)
        self.circular_entry.grid(row=7, column=1, sticky=tk.E + tk.W)
        circular_button.grid(row=7, column=0, sticky=tk.W)

        space1 = tk.Label(Command, text="\n", bg="yellow")
        space1.grid(row=8, column=0, columnspan=3, rowspan=2)

        settings_button.grid(row=0, column=2, rowspan=8, sticky=tk.E + tk.N + tk.S)

        trig_time_label.grid(row=10, column=0, columnspan=2)
        self.trig_lower_entry.grid(row=11, column=0, rowspan=2, sticky=tk.W)
        self.trig_upper_entry.grid(row=11, column=1, rowspan=2, sticky=tk.E)

        event_time_label.grid(row=13, column=0, columnspan=2)
        self.event_lower_entry.grid(row=14, column=0, rowspan=2, sticky=tk.W)
        self.event_upper_entry.grid(row=14, column=1, rowspan=2, sticky=tk.E)

        peak_sum_label.grid(row=16, column=0, columnspan=2)
        self.peak_sum_lower_entry.grid(row=17, column=0, rowspan=2, sticky=tk.W)
        self.peak_sum_upper_entry.grid(row=17, column=1, rowspan=2, sticky=tk.E)

        times_button.grid(row=14, column=2, rowspan=4, sticky=tk.E)
        times_check_button.grid(row=10, column=2, rowspan=4, sticky=tk.E)

        space2 = tk.Label(Command, text="\n", bg="yellow")
        space2.grid(row=19, column=0, columnspan=3, rowspan=2)

        boot_button.grid(row=21, column=0, columnspan=3, rowspan=2, sticky=tk.W + tk.E + tk.N + tk.S, pady=20)

        gif_button.grid(row=23, column=0, columnspan=3, rowspan=2, sticky=tk.W + tk.E + tk.N + tk.S, pady=20)

        CircImportExport.grid(row=25, column=0, columnspan=3, sticky=tk.W + tk.E + tk.N + tk.S, pady=(0, 20))

        CaliImportExport.grid(row=28, column=0, columnspan=3, sticky=tk.W + tk.E + tk.N + tk.S, pady=(0, 20))

        # Menu
        menu_bar = tk.Menu()
        self.config(menu=menu_bar)
        menu_bar.add_command(label='  Open file/path  ', command=lambda: self.menu_open_file())
        menu_bar.add_command(label='  Save data  ', command=lambda: self.save_data())
        menu_bar.add_command(label='  Settings  ', command=lambda: self.general_settings())
        menu_bar.add_command(label='  Exit  ', command=lambda: self.destroy())

        '''         Trigger frame        '''
        self.trigger_figure_canvas = FigureCanvasTkAgg(self.trig_figure, Trigger)
        self.trigger_figure_canvas.draw()
        trigger_toolbar = NavigationToolbar2Tk(self.trigger_figure_canvas, Trigger)
        trigger_toolbar.update()
        self.trigger_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)
        TriggerPeakTab.add(Trigger, text='Trigger times')

        '''         Peak sum frame        '''
        self.peak_sum_figure_canvas = FigureCanvasTkAgg(self.peak_sum_figure, PeakSum)
        self.peak_sum_figure_canvas.draw()
        peak_sum_toolbar = NavigationToolbar2Tk(self.peak_sum_figure_canvas, PeakSum)
        peak_sum_toolbar.update()
        self.peak_sum_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)
        TriggerPeakTab.add(PeakSum, text='Peak sum')

        '''         Event frame        '''
        self.event_figure_canvas = FigureCanvasTkAgg(self.event_figure, Event)
        self.event_figure_canvas.draw()
        event_toolbar = NavigationToolbar2Tk(self.event_figure_canvas, Event)
        event_toolbar.update()
        self.event_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)

        '''         VMI laser on frame        '''
        self.laser_on_figure_canvas = FigureCanvasTkAgg(self.laser_on_figure, VMI_on)
        self.laser_on_figure_canvas.draw()
        laser_on_toolbar = NavigationToolbar2Tk(self.laser_on_figure_canvas, VMI_on)
        laser_on_toolbar.update()
        self.laser_on_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)

        '''         VMI laser off frame        '''
        self.laser_off_figure_canvas = FigureCanvasTkAgg(self.laser_off_figure, VMI_off)
        self.laser_off_figure_canvas.draw()
        laser_off_toolbar = NavigationToolbar2Tk(self.laser_off_figure_canvas, VMI_off)
        laser_off_toolbar.update()
        self.laser_off_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=(20, 0))

        self.laser_off_slider_val = tk.IntVar()
        self.laser_off_slider_val.set(0)
        self.laser_off_slider = ttk.Scale(VMI_off, from_=0, to=1, orient='horizontal',
                                          variable=self.laser_off_slider_val,
                                          command=lambda event: self.plot_vmi_off(event))
        self.laser_off_slider.pack(fill='x', side='top')
        self.laser_off_slider['state'] = 'disabled'
        if self.circularize:
            self.laser_off_slider_val.set(1)
            self.laser_off_slider['state'] = 'normal'

        '''         Radial binning frame        '''
        self.radial_figure_canvas = FigureCanvasTkAgg(self.rad_figure, Radial)
        self.radial_figure_canvas.draw()
        radial_toolbar = NavigationToolbar2Tk(self.radial_figure_canvas, Radial)
        radial_toolbar.update()
        self.radial_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)

        '''         Intensity frame        '''
        self.intensity_figure_canvas = FigureCanvasTkAgg(self.intensity_figure, Intensity)
        self.intensity_figure_canvas.draw()
        intensity_toolbar = NavigationToolbar2Tk(self.intensity_figure_canvas, Intensity)
        intensity_toolbar.update()
        self.intensity_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)

        '''         Bootstrap Intensity frame        '''
        self.boot_intensity_figure_canvas = FigureCanvasTkAgg(self.boot_intensity_figure, Boot_plot)
        self.boot_intensity_figure_canvas.draw()
        boot_intensity_toolbar = NavigationToolbar2Tk(self.boot_intensity_figure_canvas, Boot_plot)
        boot_intensity_toolbar.update()
        self.boot_intensity_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)

        self.boot_statistic_slider_val = tk.IntVar()
        self.boot_statistic_slider_val.set(1)
        self.boot_statistic_slider = ttk.Scale(Boot_hist, from_=0, to=self.nr_bootstraps, orient='horizontal',
                                               variable=self.boot_statistic_slider_val,
                                               command=lambda event: self.boot_slider_changed(event))
        self.boot_statistic_slider.pack(fill='x', side='top')
        self.boot_statistic_slider['state'] = 'disabled'

        self.boot_statistic_figure_canvas = FigureCanvasTkAgg(self.boot_statistic_figure, Boot_hist)
        self.boot_statistic_figure_canvas.draw()
        boot_statistic_toolbar = NavigationToolbar2Tk(self.boot_statistic_figure_canvas, Boot_hist)
        boot_statistic_toolbar.update()
        self.boot_statistic_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=(20, 0))

        '''         Abel VMI frame        '''
        self.abel_vmi_figure_canvas = FigureCanvasTkAgg(self.abel_vmi_figure, AbelVMI)
        self.abel_vmi_figure_canvas.draw()
        abel_vmi_toolbar = NavigationToolbar2Tk(self.abel_vmi_figure_canvas, AbelVMI)
        abel_vmi_toolbar.update()
        self.abel_vmi_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)

        '''         Timeresolved frame        '''
        self.event_gif = ttk.Label(self.GIFEvent)

        self.raw_vmi_gif = ttk.Label(self.GIFVMI)
        self.raw_vmi_gif.pack(fill='both', expand=True, padx=20, pady=20)
        self.raw_off_vmi_gif = ttk.Label(self.GIFVMIOff)
        self.raw_off_vmi_gif.pack(fill='both', expand=True, padx=20, pady=20)
        self.raw_circ_vmi_gif = ttk.Label(self.GIFVMICirc)
        self.raw_circ_vmi_gif.pack(fill='both', expand=True, padx=20, pady=20)
        self.inv_vmi_gif = ttk.Label(self.GIFAbelVMI)
        self.inv_vmi_gif.pack(fill='both', expand=True, padx=20, pady=20)
        self.intensity_gif = ttk.Label(self.GIFIntensity)
        self.intensity_gif.pack(fill='both', expand=True, padx=20, pady=20)
        self.radial_gif = ttk.Label(self.GIFRadial)
        self.radial_gif.pack(fill='both', expand=True, padx=20, pady=20)

        self.gif_frame_slider_val = tk.IntVar()
        self.gif_frame_slider_val.set(self.current_frame + 1)
        self.gif_frame_slider = ttk.Scale(self.GIFEvent, from_=1, to=self.max_frames, orient='horizontal',
                                          variable=self.gif_frame_slider_val,
                                          command=lambda event: self.slider_changed(event))
        self.gif_frame_slider['state'] = 'disabled'
        self.gif_frame_slider.pack(fill='x', expand=True, padx=20, pady=(10, 0))
        self.event_gif.pack(fill='both', expand=True, padx=20, pady=(0, 20))

        # Control frame
        self.pause_play_button_str = tk.StringVar()
        self.pause_play_button_str.set('Pause')
        pause_play_button = tk.Button(self.GIFControls, textvariable=self.pause_play_button_str,
                                      command=lambda: self.play_pause_gif())
        pause_play_button.grid(row=0, column=0)

        wait_label = tk.Label(self.GIFControls, text='Animation speed (ms)')
        wait_str = tk.StringVar()
        wait_str.set('100')
        self.wait_entry = tk.Entry(self.GIFControls, validate='key', validatecommand=vcmdi,
                                   width=5, textvariable=wait_str)
        wait_label.grid(row=1, column=0)
        self.wait_entry.grid(row=1, column=1)

        self.frame_slider_val = tk.IntVar()
        self.frame_slider_val.set(self.current_frame + 1)
        self.frame_slider = ttk.Scale(self.GIFControls, from_=1, to=self.max_frames, orient='horizontal',
                                      variable=self.frame_slider_val, command=lambda event: self.slider_changed(event))
        self.frame_slider['state'] = 'disabled'
        self.frame_slider_str = tk.StringVar()
        self.frame_slider_str.set(f'{self.current_frame + 1}')
        self.frame_slider_label = tk.Label(self.GIFControls, textvariable=self.frame_slider_str)
        self.current_frame_label = tk.Label(self.GIFControls, text='Current frame ')
        self.current_frame_label.grid(row=2, column=0)
        self.frame_slider_label.grid(row=2, column=1)
        self.frame_slider.grid(row=3, column=0, columnspan=3)

        self.frame_time_label = tk.Label(self.GIFControls, text='Event time interval (s)')
        self.frame_time_label_lower_var = tk.StringVar()
        self.frame_time_label_lower = tk.Label(self.GIFControls, textvariable=self.frame_time_label_lower_var)
        self.frame_time_label_upper_var = tk.StringVar()
        self.frame_time_label_upper = tk.Label(self.GIFControls, textvariable=self.frame_time_label_upper_var)
        self.frame_time_label.grid(row=4, column=0, columnspan=3)
        self.frame_time_label_lower.grid(row=5, column=0)
        self.frame_time_label_upper.grid(row=5, column=1)

        '''         Calibration Frames      '''
        '''Raw data frame'''
        self.cali_raw_figure_canvas = FigureCanvasTkAgg(self.cali_raw_figure, self.DataFitTabPlot)
        self.cali_raw_figure_canvas.draw()
        self.cali_raw_toolbar = NavigationToolbar2Tk(self.cali_raw_figure_canvas, self.DataFitTabPlot)
        self.cali_raw_toolbar.update()
        self.cali_raw_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)

        # Preset distributions
        self.t_dist = 'from FittingRoutine import t_dist, t_dist_jac, t_dist_hess\n\ndef f(x,a,mu,sigma,nu,b):\n\treturn t_dist(x,a,mu,sigma,nu,b)\n\ndef jac(x,a,mu,sigma,nu,b):\n\treturn t_dist_jac(x,a,mu,sigma,nu,b)\n\ndef hess(x,a,mu,sigma,nu,b):\n\treturn t_dist_hess(x,a,mu,sigma,nu,b)'
        self.gaus_dist = 'from FittingRoutine import gauss_dist, gauss_dist_jac, gauss_dist_hess\n\ndef f(x,a,mu,sigma,b):\n\treturn gauss_dist(x,a,mu,sigma,b)\n\ndef jac(x,a,mu,sigma,b):\n\treturn gauss_dist_jac(x,a,mu,sigma,b)\n\ndef hess(x,a,mu,sigma,b):\n\treturn gauss_dist_hess(x,a,mu,sigma,b)'
        self.t_dist_double = 'import numpy as np\nfrom FittingRoutine import t_dist, t_dist_jac, t_dist_hess\n\ndef f(x,a1,mu1,sigma1,nu1,b1,a2,mu2,sigma2,nu2,b2):\n\treturn t_dist(x,a1,mu1,sigma1,nu1,b1)+\\\n\t\tt_dist(x,a2,mu2,sigma2,nu2,b2)\n\ndef jac(x,a1,mu1,sigma1,nu1,b1,a2,mu2,sigma2,nu2,b2):\n\tif isinstance(x,np.ndarray):\n\t\tjac = np.array([np.append(t_dist_jac(x[0],a1,mu1,sigma1,nu1,b1),\n\t\tt_dist_jac(x[0],a2,mu2,sigma2,nu2,b2))])\n\t\tfor i in range(1,len(x)):\n\t\t\tjac = np.append(jac,[np.append(t_dist_jac(x[i],a1,mu1,sigma1,nu1,b1),\n\t\tt_dist_jac(x[i],a2,mu2,sigma2,nu2,b2))],axis=0)\n\t\treturn jac.T\n\telse:\n\t\treturn np.append(t_dist_jac(x,a1,mu1,sigma1,nu1,b1),\n\t\tt_dist_jac(x,a2,mu2,sigma2,nu2,b2))\n\ndef hess(x,a1,mu1,sigma1,nu1,b1,a2,mu2,sigma2,nu2,b2):\n\thess_1 = t_dist_hess(x,a1,mu1,sigma1,nu1,b1)\n\thess_2 = t_dist_hess(x,a2,mu2,sigma2,nu2,b2)\n\tsh = hess_1.shape\n\tif len(sh) == 3:\n\t\thess = np.zeros((2*sh[0],2*sh[1],sh[2]))\n\t\thess[:sh[0],:sh[0],:] = hess_1\n\t\thess[sh[0]:2*sh[0],sh[0]:2*sh[0],:] = hess_2\n\telse:\n\t\thess = np.zeros((2*sh[0],2*sh[1]))\n\t\thess[:sh[0],:sh[0]] = hess_1\n\t\thess[sh[0]:2*sh[0],sh[0]:2*sh[0]] = hess_2\n\treturn hess'

        # Function
        self.peak_fit_text = scrolledtext.ScrolledText(DataFitTabSettings, wrap=tk.WORD, height=10)
        self.peak_fit_text.insert(tk.END, self.t_dist)
        self.peak_fit_text.configure(state='disabled')

        self.peak_fit_explanation_txt = tk.StringVar()
        if len(self.boot_intensity) != 0 and len(self.boot_intensity) == len(self.boot_radial):
            self.peak_fit_explanation_txt.set(
                f'Function must have name f, and form\ndef f(x,*params):\nJacobian and Hessian returning \n(len(params),n) and (len(params), len(params, n))\n can be given with jac and hess functions.\n\n Lower and upper bounds are only used for\n global minimization with method \'diff_evol\'\n\n For method \'diff_evol\' errors are only reliable if jac and hess are given\n\nCurrent raidal errors,\nSystematic: {self.sys_radial_err}\nStatistical: {self.boot_radial[-1] - self.boot_radial[-2]}')
        else:
            self.peak_fit_explanation_txt.set(
                f'Function must have name f, and form\ndef f(x,*params):\nJacobian and Hessian returning \n(len(params),n) and (len(params), len(params, n))\n can be given with jac and hess functions.\n\n Lower and upper bounds are only used for\n global minimization with method \'diff_evol\'\n\n For method \'diff_evol\' errors are only reliable if jac and hess are given\n\nCurrent raidal errors,\nSystematic: {self.sys_radial_err}\nStatistical: not calculated')
        peak_fit_explanation = ttk.Label(DataFitTabSettings, textvariable=self.peak_fit_explanation_txt)

        # Presets
        presets_label = ttk.Label(DataFitTabSettings, text='Choose preset')
        preset_dropdown_options = ['Students-t', 'Students-t', 'Two students-t', 'Gauss', 'Custom']
        self.preset_dropdown_var = tk.StringVar()
        self.preset_dropdown_var.trace("w", self.set_peak_fit_preset)
        self.preset_dropdown = ttk.OptionMenu(DataFitTabSettings, self.preset_dropdown_var, *preset_dropdown_options)
        self.preset_dropdown_var.set('Students-t')

        # Lock function
        self.peak_fit_function = {}
        self.peak_fit_function_nr_params = 0
        self.function_lock_var = tk.StringVar()
        self.function_lock_var.set('Function Unlocked')
        self.function_lock_button = tk.Button(DataFitTabSettings, textvariable=self.function_lock_var,
                                              command=self.function_lock)
        self.function_lock_error = scrolledtext.ScrolledText(DataFitTabSettings, wrap=tk.WORD, height=15)
        self.function_lock_error.configure(state='disabled')

        # Parameters table
        self.parameters_table = ttk.Frame(DataFitTabSettings)
        self.parameters_table_head = ['Param', 'init', 'lower', 'upper']
        self.param_names = []
        self.param_init = []
        self.param_lower = []
        self.param_upper = []

        # Fitting
        self.peak_fit_button = ttk.Button(DataFitTabSettings, text='Fit', command=self.run_peak_fit)
        self.peak_fit = {}
        self.method = tk.StringVar()
        self.method.set('BFGS')
        self.kwargs = tk.StringVar()

        self.systematic_errors_var = tk.IntVar()
        self.statistical_errors_var = tk.IntVar()
        self.systematic_errors_var.set(0)
        self.statistical_errors_var.set(1)
        self.radius_low_var = tk.StringVar()
        self.radius_high_var = tk.StringVar()
        self.peak_fit_settings_button = ttk.Button(DataFitTabSettings, text='Settings', command=self.peak_fit_settings)

        # Fitting results
        self.parameters_fit_table = ttk.LabelFrame(DataFitTabSettings, text='Results')
        self.parameters_fit_table_head = ['Param', 'value', 'error', 'select']
        self.param_fit_val = []
        self.param_fit_error = []
        self.param_fit_check = []
        self.peak_fit_export = tk.Button(self.parameters_fit_table, text='Move selected to calibration',
                                         command=self.export_peak_fit)

        # Setting widgets
        peak_fit_explanation.grid(row=0, sticky=tk.W, pady=15)
        self.peak_fit_text.grid(row=1, sticky=tk.EW, pady=5)
        presets_label.grid(row=2, sticky=tk.W)
        self.preset_dropdown.grid(row=2, sticky=tk.E)
        self.function_lock_button.grid(row=3, sticky=tk.EW)
        self.parameters_table.grid(row=5, sticky=tk.EW)
        DataFitTabSettings.grid_columnconfigure(0, weight=1)
        self.parameters_fit_table.grid(row=8, sticky=tk.EW)

        '''Calibration fit frame'''
        self.laser_figure_canvas = FigureCanvasTkAgg(self.laser_hist_figure, self.LaserFitTabPlot)
        self.laser_figure_canvas.draw()
        self.laser_toolbar = NavigationToolbar2Tk(self.laser_figure_canvas, self.LaserFitTabPlot)
        self.laser_toolbar.update()
        self.laser_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)

        # Function
        self.laser_fit_text = scrolledtext.ScrolledText(LaserFitTabSettings, wrap=tk.WORD, height=10)
        self.laser_fit_text.insert(tk.END, self.t_dist)
        self.laser_fit_text.configure(state='disabled')

        self.laser_explanation_txt = tk.StringVar()
        if np.sum(self.freq_hist) != 0:
            self.laser_explanation_txt.set(
                f'Function must have name f, and form\ndef f(x,*params):\nJacobian and Hessian returning \n(len(params),n) and (len(params), len(params, n))\n can be given with jac and hess functions.\n\n Lower and upper bounds are only used for\n global minimization with method \'diff_evol\'\n\n For method \'diff_evol\' errors are only reliable if jac and hess are given\n\nCurrent statistics,\nMean: {np.average(self.freq_bins, weights=self.freq_hist)}\nStd.: {np.sqrt(np.cov(self.freq_bins, aweights=self.freq_hist))}')
        else:
            self.laser_explanation_txt.set(
                f'Function must have name f, and form\ndef f(x,*params):\nJacobian and Hessian returning \n(len(params),n) and (len(params), len(params, n))\n can be given with jac and hess functions.\n\n Lower and upper bounds are only used for\n global minimization with method \'diff_evol\'\n\n For method \'diff_evol\' errors are only reliable if jac and hess are given\n\nCurrent statistics,\nMean: No laser\nStd.: No laser')
        laser_fit_explanation = ttk.Label(LaserFitTabSettings, textvariable=self.laser_explanation_txt)

        # Presets
        laser_presets_label = ttk.Label(LaserFitTabSettings, text='Choose preset')
        laser_preset_dropdown_options = ['Gauss', 'Students-t', 'Two students-t', 'Gauss', 'Custom']
        self.laser_preset_dropdown_var = tk.StringVar()
        self.laser_preset_dropdown_var.trace("w", self.set_laser_fit_preset)
        self.laser_preset_dropdown = ttk.OptionMenu(LaserFitTabSettings, self.laser_preset_dropdown_var,
                                                    *laser_preset_dropdown_options)
        self.laser_preset_dropdown_var.set('Gauss')

        # Lock function
        self.laser_fit_function = {}
        self.laser_fit_function_nr_params = 0
        self.laser_function_lock_var = tk.StringVar()
        self.laser_function_lock_var.set('Function Unlocked')
        self.laser_function_lock_button = tk.Button(LaserFitTabSettings, textvariable=self.laser_function_lock_var,
                                                    command=self.laser_function_lock)
        self.laser_function_lock_error = scrolledtext.ScrolledText(LaserFitTabSettings, wrap=tk.WORD, height=15)
        self.laser_function_lock_error.configure(state='disabled')

        # Parameters table
        self.laser_parameters_table = ttk.Frame(LaserFitTabSettings)
        self.laser_parameters_table_head = ['Param', 'init', 'lower', 'upper']
        self.laser_param_names = []
        self.laser_param_init = []
        self.laser_param_lower = []
        self.laser_param_upper = []

        # Fitting
        self.laser_fit_button = ttk.Button(LaserFitTabSettings, text='Fit', command=self.run_laser_fit)
        self.laser_fit = {}
        self.laser_method = tk.StringVar()
        self.laser_method.set('BFGS')
        self.laser_kwargs = tk.StringVar()

        self.laser_statistical_errors_var = tk.IntVar()
        self.laser_statistical_errors_var.set(1)
        self.laser_energy_low_var = tk.StringVar()
        self.laser_energy_high_var = tk.StringVar()
        self.laser_fit_settings_button = ttk.Button(LaserFitTabSettings, text='Settings',
                                                    command=self.laser_fit_settings)

        # Fitting results
        self.laser_parameters_fit_table = ttk.LabelFrame(LaserFitTabSettings, text='Results')
        self.laser_parameters_fit_table_head = ['Param', 'value', 'error']
        self.laser_param_fit_val = []
        self.laser_param_fit_error = []

        # Setting widgets
        laser_fit_explanation.grid(row=0, sticky=tk.W, pady=15)
        self.laser_fit_text.grid(row=1, sticky=tk.EW, pady=5)
        laser_presets_label.grid(row=2, sticky=tk.W)
        self.laser_preset_dropdown.grid(row=2, sticky=tk.E)
        self.laser_function_lock_button.grid(row=3, sticky=tk.EW)
        self.laser_parameters_table.grid(row=5, sticky=tk.EW)
        LaserFitTabSettings.grid_columnconfigure(0, weight=1)
        self.laser_parameters_fit_table.grid(row=8, sticky=tk.EW)

        '''Calibration fit frame'''
        self.cali_plot_figure_canvas = FigureCanvasTkAgg(self.cali_plot_figure, self.CalibrationFitTabPlot)
        self.cali_plot_figure_canvas.draw()
        self.cali_plot_toolbar = NavigationToolbar2Tk(self.cali_plot_figure_canvas, self.CalibrationFitTabPlot)
        self.cali_plot_toolbar.update()
        self.cali_plot_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=20, pady=20)

        # Preset distributions
        self.linear = 'import numpy as np\n\n# Inverse Linear\ndef f(x,a,b):\n\treturn (x - b)/a\n\n\ndef jac(x,a,b):\n\tif isinstance(x,np.ndarray):\n\t\tb_dev=-np.ones(len(x))/a\n\telse:\n\t\tb_dev=-1/a\n\treturn np.array([(b-x)/a**2,b_dev])\n\ndef hess(x,a,b):\n\treturn np.array([[2*(x-b)/a**3,1/a**2],[1/a**2,x*0]])'
        self.prop = 'import numpy as np\n\n# Inverse Proportional\ndef f(x,a):\n\treturn x/a\n\n\ndef jac(x,a):\n\treturn np.array([-x/a**2])\n\ndef hess(x,a):\n\treturn np.array([2*x/a**3])'
        self.pow = 'import numpy as np\n\n# Inverse Power\ndef f(x,a,b):\n\treturn (x/a)**(1/b)\n\n\ndef jac(x,a,b):\n\treturn np.array([-(x/a)**(1/b)/a/b,-(x/a)**(1/b)*np.log(x/a)/b**2])\n\ndef hess(x,a,b):\n\th0=f(x,a,b)\n\th11 = (b+1)*h0/(a*b)**2\n\th12 = h0*(np.log(x/a)+b)/a/b**3\n\th22 = h0*np.log(x/a)*(np.log(x/a)+2*b)/b**4\n\treturn np.array([[h11,h12],[h12,h22]])'

        # Function
        self.cali_fit_text = scrolledtext.ScrolledText(CalibrationFitTabSettings, wrap=tk.WORD, height=10)
        self.cali_fit_text.insert(tk.END, self.pow)
        self.cali_fit_text.configure(state='disabled')

        cali_fit_explanation_txt = 'Choose fitting function:'
        cali_fit_explanation = ttk.Label(CalibrationFitTabSettings, text=cali_fit_explanation_txt)

        # Presets
        cali_presets_label = ttk.Label(CalibrationFitTabSettings, text='Choose preset')
        cali_preset_dropdown_options = ['Power', 'Power', 'Proportional', 'Linear']
        self.cali_preset_dropdown_var = tk.StringVar()
        self.cali_preset_dropdown = ttk.OptionMenu(CalibrationFitTabSettings, self.cali_preset_dropdown_var,
                                                   *cali_preset_dropdown_options)

        # Lock function
        self.cali_fit_function = {}
        self.cali_fit_function_nr_params = 0
        self.cali_function_lock_error = scrolledtext.ScrolledText(CalibrationFitTabSettings, wrap=tk.WORD, height=15)
        self.cali_function_lock_error.configure(state='disabled')

        # Parameters table
        self.cali_parameters_table = ttk.Frame(CalibrationFitTabSettings)
        self.cali_parameters_table_head = ['Param', 'init', 'lower', 'upper']
        self.cali_param_names = []
        self.cali_param_init = []
        self.cali_param_lower = []
        self.cali_param_upper = []

        # Fitting
        self.cali_fit_button = ttk.Button(CalibrationFitTabSettings, text='Fit', command=self.run_cali_fit)
        self.cali_fit = {}
        self.cali_method = tk.StringVar()
        self.cali_method.set('BFGS')
        self.cali_kwargs = tk.StringVar()

        self.cali_fit_settings_button = ttk.Button(CalibrationFitTabSettings, text='Settings',
                                                   command=self.cali_fit_settings)

        # Fitting results
        self.cali_parameters_fit_table = ttk.LabelFrame(CalibrationFitTabSettings, text='Results')
        self.cali_parameters_fit_table_head = ['Param', 'value', 'error']
        self.cali_param_fit_val = []
        self.cali_param_fit_error = []
        self.cali_fit_export = tk.Button(self.cali_parameters_fit_table, text='Export calibration',
                                         command=self.export_cali_fit)

        # Points frame
        self.cali_points = tk.LabelFrame(CalibrationFitTabSettings, text='Calibration points', height=150)
        self.cali_points_canvas = tk.Canvas(self.cali_points)
        self.cali_points_scroll = tk.Scrollbar(self.cali_points, orient="vertical",
                                               command=self.cali_points_canvas.yview)

        self.cali_points_canvas.configure(yscrollcommand=self.cali_points_scroll.set)

        self.cali_points_scroll.pack(side="right", fill="y")
        self.cali_points_canvas.pack(side='left', fill='both', expand=True)

        self.cali_points_table = tk.Frame(self.cali_points_canvas)
        self.cali_points_table.bind('<Configure>', lambda e: self.cali_points_canvas.configure(
            scrollregion=self.cali_points_canvas.bbox('all')))
        self.cali_points_canvas.create_window((0, 0), window=self.cali_points_table, anchor='nw')
        self.cali_points_canvas.bind_all("<MouseWheel>",
                                         lambda e: self.cali_points_canvas.yview_scroll(-1 * e.delta, "units"))

        # Points table
        self.cali_points_table_head = ['Radius', 'Error', 'Energy', 'Error']
        self.cali_points_exp = []
        self.cali_points_exp_er = []
        self.cali_points_lit = []
        self.cali_points_lit_err = []

        for i in range(len(self.cali_points_table_head)):
            label = ttk.Label(self.cali_points_table, text=self.cali_points_table_head[i])
            label.grid(row=0, column=i)

        size = int((self.cali_fit_text.winfo_width()) / 4 / tk.font.Font(font='TkFixedFont').actual()['size'])

        self.cali_points_button_plus = ttk.Button(self.cali_points_table, text='+',
                                                  command=self.cali_points_plus, width=size)
        self.cali_points_button_minus = ttk.Button(self.cali_points_table, text='-',
                                                   command=self.cali_points_minus, width=size)
        self.cali_points_button_plus.grid(row=1, column=2)
        self.cali_points_button_minus.grid(row=1, column=3)

        # Nice to know
        energy_lookup_table = ttk.Button(CalibrationFitTabSettings, text='Energy Lookup Table',
                                         command=self.open_lookup_tabel)

        # Setting widgets
        energy_lookup_table.grid(row=0, sticky=tk.EW, pady=5)
        self.cali_points.grid(row=1, sticky=tk.EW, pady=10)
        cali_fit_explanation.grid(row=2, sticky=tk.W, pady=15)
        self.cali_fit_text.grid(row=3, sticky=tk.EW, pady=5)
        cali_presets_label.grid(row=4, sticky=tk.W)
        self.cali_preset_dropdown.grid(row=4, sticky=tk.E)
        self.cali_parameters_table.grid(row=7, sticky=tk.EW)
        self.cali_parameters_fit_table.grid(row=10, sticky=tk.EW)
        CalibrationFitTabSettings.grid_columnconfigure(0, weight=1)
        self.cali_fit_button.grid(row=8, column=0, sticky=tk.W)
        self.cali_fit_settings_button.grid(row=8, column=0, sticky=tk.E)
        self.cali_preset_dropdown_var.trace("w", self.set_cali_fit_preset)
        self.cali_preset_dropdown_var.set('Power')
        self.set_cali_fit_preset()

        '''         Initial plot        '''
        self.plot_trig()
        self.plot_event()
        self.plot_peak_sum()
        self.plot_vmi_on()
        self.plot_vmi_off()
        self.plot_rad_bin()
        self.plot_intensity()
        self.plot_abel_vmi()
        self.initialise_gif()
        self.initialise_boot()
        self.update_laser_hist()

    def initialise(self):
        """             Initial data            """
        # Load the file to a root dataframe
        if os.path.exists(self.root_file):
            pass
        else:
            load_to_root(self.load_file_path, self.dat_file_path)
            get_all_time_data(self.event_file_path, self.root_file)
            VMI_filter(self.vmi_file_path, self.root_file, ['-r 200', '-o 3', '-v 250', '-c 100'])

        '''             Read data               '''

        # Point to the dataframe
        root_file = ROOT.TFile(self.root_file, "read")
        trig_hist_data = root_file.Get("TRIG_TOT")
        event_hist_data = root_file.Get("EVENT_TOT")
        peak_hist_data = root_file.Get("PEAK_TOT")
        trig_hist_current_data = root_file.Get("TRIG_NOW")
        event_hist_current_data = root_file.Get("EVENT_NOW")
        peak_hist_current_data = root_file.Get("PEAK_NOW")
        freq_hist_data = root_file.Get("FREQ_TOT")
        freq_bin_data = root_file.Get("FREQ_BIN")
        vmi_laser_on_data = root_file.Get("VMI_LASER_ON")
        vmi_laser_off_data = root_file.Get("VMI_LASER_OFF")
        vmi_fold_data = root_file.Get("VMI_FOLD")
        radial_bin_val_data = root_file.Get("radial_bin_values")
        radial_bin_amp_data = root_file.Get("radial_bin_TOTAL")
        radial_bin_err_data = root_file.Get("radial_bin_error")
        vmi_abel_data = root_file.Get("VMI_ABEL")
        intensity_data = root_file.Get("INTENSITY")
        intensity_err_data = root_file.Get("INTENSITY_ERROR")

        settings_val_data = root_file.Get("SETTINGS_VALUES")
        settings_times_data = root_file.Get("SETTINGS_TIMES")
        settings_bin_data = root_file.Get("SETTINGS_BINS")
        settings_mode_data = root_file.Get("SETTINGS_MODE")
        center_data = root_file.Get("CENTER")

        self.center_vals = []
        for entry in center_data:
            self.center_vals.append(entry.center)

        settings_vals = []
        for entry in settings_val_data:
            settings_vals.append(entry.settings_val)

        settings_times = []
        for entry in settings_times_data:
            settings_times.append(entry.settings_time)

        settings_bins = []
        for entry in settings_bin_data:
            settings_bins.append(entry.settings_bin)

        settings_mode = []
        for entry in settings_mode_data:
            settings_mode.append(entry.settings_mode)

        self.resolution = settings_vals[0]
        self.norder = settings_vals[1] - 1
        self.vmi_resolution = settings_vals[2]
        self.circularize = settings_vals[3]
        self.nr_angular_bin = settings_vals[4]

        self.trig_lower = settings_times[0] / 1e-6
        self.trig_upper = settings_times[1] / 1e-6
        self.event_lower = settings_times[2]
        self.event_upper = settings_times[3]
        self.peak_sum_lower = settings_times[5]
        self.peak_sum_upper = settings_times[6]
        self.max_rad = str(settings_times[4])

        self.event_time_bin_size = settings_bins[0]
        self.trig_time_bin_size = settings_bins[1] / 1e-6  # in micro seconds
        self.peak_sum_bin_size = settings_bins[2]
        self.freq_bin_size = settings_bins[3]

        self.background_sub_mode = [settings_mode[0], settings_mode[1]]
        self.event_background_lower = settings_mode[2]
        self.event_background_upper = settings_mode[3]

        self.event_hist = []
        self.event_hist_err = []
        self.event_current_hist = []
        self.event_current_hist_err = []

        self.trig_hist = []
        self.trig_hist_err = []
        self.trig_current_hist = []
        self.trig_current_hist_err = []

        self.peak_sum_hist = []
        self.peak_sum_hist_err = []
        self.peak_sum_current_hist = []
        self.peak_sum_current_hist_err = []

        self.event_time_bins = [0]
        self.trig_time_bins = [0]
        self.peak_sum_bins = [0]
        self.freq_bins = []

        self.vmi_laser_on = np.zeros((self.vmi_resolution, self.vmi_resolution))
        self.vmi_laser_off = np.zeros((self.vmi_resolution, self.vmi_resolution))
        self.vmi_fold = np.zeros((self.vmi_resolution, self.vmi_resolution))
        self.abel_vmi = np.zeros((self.vmi_resolution, self.vmi_resolution))
        self.circ_vmi = np.zeros((self.vmi_resolution, self.vmi_resolution))

        # Get the event and trigger data from file
        for entry in trig_hist_data:
            self.trig_hist.append(entry.trig_hist)
            self.trig_hist_err.append(np.sqrt(entry.trig_hist))
            self.trig_time_bins.append(self.trig_time_bins[-1] + self.trig_time_bin_size)
        for entry in trig_hist_current_data:
            self.trig_current_hist.append(entry.trig_hist)
            self.trig_current_hist_err.append(np.sqrt(entry.trig_hist))
        for entry in event_hist_data:
            self.event_hist.append(entry.event_hist)
            self.event_hist_err.append(np.sqrt(entry.event_hist))
            self.event_time_bins.append(self.event_time_bins[-1] + self.event_time_bin_size)
        for entry in event_hist_current_data:
            self.event_current_hist.append(entry.event_hist)
            self.event_current_hist_err.append(np.sqrt(entry.event_hist))
        for entry in peak_hist_data:
            self.peak_sum_hist.append(entry.peak_hist)
            self.peak_sum_hist_err.append(np.sqrt(entry.peak_hist))
            self.peak_sum_bins.append(self.peak_sum_bins[-1] + self.peak_sum_bin_size)
        for entry in peak_hist_current_data:
            self.peak_sum_current_hist.append(entry.peak_hist)
            self.peak_sum_current_hist_err.append(np.sqrt(entry.peak_hist))
        for entry in freq_hist_data:
            self.freq_hist.append(entry.freq_hist)
            self.freq_hist_err.append(np.sqrt(entry.freq_hist))
        for entry in freq_bin_data:
            self.freq_bins.append(entry.freq_bin * 4.136 / 1000)  # Convert from THz to eV
        self.event_time_bins.pop(-1)
        self.trig_time_bins.pop(-1)
        self.peak_sum_bins.pop(-1)

        # Get the VMI image data from the file
        count = 0
        for entry in vmi_laser_on_data:
            self.vmi_laser_on[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
            count += 1
        count = 0
        for entry in vmi_laser_off_data:
            self.vmi_laser_off[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
            count += 1
        count = 0
        for entry in vmi_fold_data:
            self.vmi_fold[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
            count += 1
        if self.circularize:
            vmi_circ_data = root_file.Get("VMI_CIRC")
            count = 0
            for entry in vmi_circ_data:
                self.circ_vmi[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
                count += 1

        # Get the radial binning of data form file
        self.radial_projection_values.clear()
        self.radial_projection_amplitude.clear()
        self.radial_projection_error.clear()
        self.radial_projection_amplitude_sum = np.zeros(int(self.resolution))
        self.radial_projection_error_sum = np.zeros(int(self.resolution))
        count_amp = 0
        count_err = 0
        for entry in radial_bin_val_data:
            self.radial_projection_values.append(entry.radial_bin_val)

        for entry in radial_bin_amp_data:
            self.radial_projection_amplitude.append(entry.radial_bin)
            self.radial_projection_amplitude_sum[count_amp] += entry.radial_bin
            count_amp += 1
            if count_amp == int(self.resolution):
                count_amp = 0
        for entry in radial_bin_err_data:
            self.radial_projection_error.append(entry.radial_bin_err)
            self.radial_projection_error_sum[count_err] += entry.radial_bin_err
            count_err += 1
            if count_err == int(self.resolution):
                count_err = 0
        self.radial_projection_error_sum = np.sqrt(np.array(self.radial_projection_error_sum))

        # Get abel inverse image from file
        count = 0
        for entry in vmi_abel_data:
            self.abel_vmi[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
            count += 1
        self.abel_vmi[np.array(self.abel_vmi) < 0] = 0

        # Combine the VMI
        center_col = int(int(self.vmi_resolution) / 2)
        denom = np.max(self.abel_vmi[:, :center_col])
        self.combined_vmi = np.zeros((int(self.vmi_resolution), int(self.vmi_resolution)))
        self.combined_vmi[:, center_col:] = self.vmi_fold[:, center_col:]
        if denom != 0:
            self.combined_vmi[:, :center_col] = self.abel_vmi[:, :center_col] / np.max(self.abel_vmi[:, :center_col]) * \
                                                np.max(self.vmi_fold[:, center_col:])
        else:
            self.combined_vmi[:, :center_col] = self.abel_vmi[:, :center_col]

        # Get abel inversed intensities
        self.intensity.clear()
        for entry in intensity_data:
            self.intensity.append(entry.intensity)

        self.intensity_err.clear()
        for entry in intensity_err_data:
            self.intensity_err.append(entry.error)

        root_file.Close()

    def initialise_gif(self):
        '''      GIF exist       '''
        save_dir = self.root_file
        save_dir = save_dir[:save_dir.rfind('/') + 1]
        self.circ_frames = False

        if os.path.exists(save_dir + 'event.gif') and os.path.exists(save_dir + 'raw_vmi.gif') \
                and os.path.exists(save_dir + 'inv_vmi.gif') and os.path.exists(save_dir + 'intensity.gif') \
                and os.path.exists(save_dir + 'gif_times.txt' and os.path.exists(save_dir + 'radial.gif') \
                                   and os.path.exists(save_dir + 'raw_off_vmi.gif')):
            if os.path.exists(save_dir + 'raw_circ_vmi.gif'):
                self.circ_frames = True

            # Load GIF
            self.event_frames = Image.open(save_dir + 'event.gif')
            self.raw_vmi_frames = Image.open(save_dir + 'raw_vmi.gif')
            self.raw_off_vmi_frames = Image.open(save_dir + 'raw_off_vmi.gif')
            self.inv_vmi_frames = Image.open(save_dir + 'inv_vmi.gif')
            self.intensity_frames = Image.open(save_dir + 'intensity.gif')
            self.radial_frames = Image.open(save_dir + 'radial.gif')
            self.gif_times = np.loadtxt(save_dir + 'gif_times.txt', unpack=True, delimiter=",")

            if self.circ_frames:
                self.raw_circ_vmi_frames = Image.open(save_dir + 'raw_circ_vmi.gif')

            self.max_frames = self.event_frames.n_frames
            self.frame_slider.configure(to=self.max_frames)
            self.gif_frame_slider.configure(to=self.max_frames)

            # Plot the results
            self.run_gif = True
            self.update_gif()
            self.run_gif = False
            self.pause_play_button_str.set('Play')
            self.frame_slider['state'] = 'normal'
            self.gif_frame_slider['state'] = 'normal'

    def initialise_boot(self):
        self.boot_statistic_slider['state'] = 'disabled'
        '''      Boot exist       '''
        if os.path.exists(self.root_file[:-5] + '_boot.root'):
            # Read the file
            root_file = ROOT.TFile(self.root_file[:-5] + '_boot.root', "read")
            intensity_data = root_file.Get("INTENSITY")
            intensity_err_data = root_file.Get("INTENSITY_ERROR")
            radial_data = root_file.Get("RADIUS")
            error_upper_data = root_file.Get("INTENSITY_ERROR_UPPER")
            error_lower_data = root_file.Get("INTENSITY_ERROR_LOWER")
            iterations_data = root_file.Get("BOOTSTRAP_ITERATIONS")

            bins = 0

            self.boot_intensity = []
            self.boot_radial = []
            self.intensity_err_upper = []
            self.intensity_err_lower = []

            for entry in intensity_data:
                self.boot_intensity.append(entry.intensity)
                bins += 1
            for entry in intensity_err_data:
                self.boot_intensity_sys_err.append(entry.error)
            for entry in radial_data:
                self.boot_radial.append(entry.radius)

            for entry in error_upper_data:
                self.intensity_err_upper.append(entry.intensity_err)
            for entry in error_lower_data:
                self.intensity_err_lower.append(entry.intensity_err)

            self.nr_bootstraps = int(iterations_data.GetEntries() / bins)

            self.boot_iterations = np.zeros((bins, self.nr_bootstraps))

            count = 0
            i = 0
            for entry in iterations_data:
                self.boot_iterations[i, count] = entry.iteration
                i += 1
                if i == bins:
                    i = 0
                    count += 1

            root_file.Close()

            self.boot_rad_err = np.ones(len(self.boot_intensity)) * (self.boot_radial[2] - self.boot_radial[1]) / 2
            self.boot_rad_sys_err = np.ones(len(self.boot_intensity)) * 0.02 * 20

            self.boot_statistic_slider.config(to=bins)
            self.boot_statistic_slider['state'] = 'normal'
            self.quantiles = np.array([np.array(self.boot_intensity) - np.array(self.intensity_err_lower),
                                       np.array(self.intensity_err_upper) - np.array(self.boot_intensity)])

            self.plot_bootstrap()
            self.peak_fit_explanation_txt.set(
                f'Function must have name f, and form\ndef f(x,*params):\nJacobian and Hessian returning \n(len(params),n) and (len(params), len(params, n))\n can be given with jac and hess functions.\n\n Lower and upper bounds are only used for\n global minimization with method \'diff_evol\'\n\n For method \'diff_evol\' errors are only reliable if jac and hess are given\n\nCurrent raidal errors,\nSystematic: {self.sys_radial_err}\nStatistical: {self.boot_radial[-1] - self.boot_radial[-2]}')

    def update_laser_hist(self):
        if np.sum(self.freq_hist) != 0:
            self.laser_function_lock()
            self.run_laser_fit()

    def update_settings(self):
        root_file = ROOT.TFile(self.root_file, "read")

        settings_val_data = root_file.Get("SETTINGS_VALUES")
        settings_times_data = root_file.Get("SETTINGS_TIMES")
        center_data = root_file.Get("CENTER")

        center_vals = []
        for entry in center_data:
            center_vals.append(entry.center)

        settings_vals = []
        for entry in settings_val_data:
            settings_vals.append(entry.settings_val)

        settings_times = []
        for entry in settings_times_data:
            settings_times.append(entry.settings_time)

        self.resolution = settings_vals[0]
        self.resolution_str_var.set(str(self.resolution))
        self.norder = settings_vals[1] - 1
        self.order_str_var.set(str(self.norder))
        self.vmi_resolution = settings_vals[2]
        self.vmi_resolution_str_var.set(str(self.vmi_resolution))
        self.circularize = settings_vals[3]
        if self.circularize:
            self.circular_str_truth.set('True')
        else:
            self.circular_str_truth.set('False')
        self.nr_angular_bin = settings_vals[4]
        self.circular_str_var.set(str(self.nr_angular_bin))

        self.trig_lower = settings_times[0] / 1e-6
        self.trig_lower_str_var.set(str("%.3f" % self.trig_lower))
        self.trig_upper = settings_times[1] / 1e-6
        self.trig_upper_str_var.set(str("%.3f" % self.trig_upper))
        self.event_lower = settings_times[2]
        self.event_lower_str_var.set(str("%.3f" % self.event_lower))
        self.event_upper = settings_times[3]
        self.event_upper_str_var.set(str("%.3f" % self.event_upper))
        self.peak_sum_lower = settings_times[5]
        self.peak_sum_lower_str_var.set(str("%.1f" % self.peak_sum_lower))
        self.peak_sum_upper = settings_times[6]
        self.peak_sum_upper_str_var.set(str("%.1f" % self.peak_sum_upper))
        self.max_rad = settings_times[4]
        self.resolution_max_str_var.set(str("%.2f" % self.max_rad))
        self.max_rad = str(self.max_rad)

        root_file.Close()

    def menu_open_file(self):
        load_popup = FilePopup(self)
        self.wait_window(load_popup)
        self.root_file = self.dat_file_path + ".root"
        self.initialise()
        self.update_settings()
        self.laser_off_slider_val.set(0)
        self.laser_off_slider['state'] = 'disabled'
        if self.circularize:
            self.laser_off_slider_val.set(1)
            self.laser_off_slider['state'] = 'normal'
        self.plot_trig()
        self.plot_event()
        self.plot_vmi_on()
        self.plot_vmi_off()
        self.plot_rad_bin()
        self.plot_intensity()
        self.plot_abel_vmi()
        self.initialise_gif()
        self.initialise_boot()

    def plot_event(self):
        self.event_ax.clear()
        self.event_ax.axhline(0, color='k', lw=0.5)
        self.event_ax.axvline(float(self.event_lower), ls='--', color="b", lw=0.5)
        self.event_ax.axvline(float(self.event_upper), ls='--', color="b", lw=0.5)
        if self.background_sub_mode[1] == 1:
            self.event_ax.axvline(self.event_background_lower, ls='--', color="k", lw=0.5)
            self.event_ax.axvline(self.event_background_upper, ls='--', color="k", lw=0.5)
        self.event_ax.errorbar(self.event_time_bins, self.event_hist, self.event_hist_err,
                               fmt='ro-', capsize=2, ms=2, lw=1,
                               errorevery=3, elinewidth=0.5, label='Total')
        self.event_ax.errorbar(self.event_time_bins, self.event_current_hist, self.event_current_hist_err,
                               fmt='bo-', capsize=2, ms=2, lw=1,
                               errorevery=3, elinewidth=0.5, label='Current')
        self.event_ax.set_title('Histogram of event times')
        self.event_ax.set_xlabel('Time (s)')
        self.event_ax.set_ylabel('Counts')
        self.event_ax.set_xlim(self.event_time_bins[0], self.event_time_bins[-1])
        self.event_ax.legend()
        self.event_figure_canvas.draw()

    def plot_peak_sum(self):
        self.peak_sum_ax.clear()
        self.peak_sum_ax.axhline(0, color='k', lw=0.5)
        self.peak_sum_ax.axvline(float(self.peak_sum_lower), ls='--', color="b", lw=0.5)
        self.peak_sum_ax.axvline(float(self.peak_sum_upper), ls='--', color="b", lw=0.5)
        self.peak_sum_ax.errorbar(self.peak_sum_bins, self.peak_sum_hist, self.peak_sum_hist_err,
                                  fmt='ro-', capsize=2, ms=2, lw=1,
                                  errorevery=3, elinewidth=0.5, label='Total')
        self.peak_sum_ax.errorbar(self.peak_sum_bins, self.peak_sum_current_hist, self.peak_sum_current_hist_err,
                                  fmt='bo-', capsize=2, ms=2, lw=1,
                                  errorevery=3, elinewidth=0.5, label='Current')
        self.peak_sum_ax.set_title('Histogram of VMI peak sums')
        self.peak_sum_ax.set_xlabel(r'Amplitude')
        self.peak_sum_ax.set_ylabel('Counts')
        self.peak_sum_ax.set_xlim(self.peak_sum_bins[0], self.peak_sum_bins[-1])
        self.peak_sum_ax.legend()
        self.peak_sum_figure_canvas.draw()

    def plot_trig(self):
        self.trig_ax.clear()
        self.trig_ax.axhline(0, color='k', lw=0.5)
        self.trig_ax.axvline(float(self.trig_lower), ls='--', color="b", lw=0.5)
        self.trig_ax.axvline(float(self.trig_upper), ls='--', color="b", lw=0.5)
        self.trig_ax.errorbar(self.trig_time_bins, self.trig_hist, self.trig_hist_err,
                              fmt='ro-', capsize=2, ms=2, lw=1,
                              errorevery=3, elinewidth=0.5, label='Total')
        self.trig_ax.errorbar(self.trig_time_bins, self.trig_current_hist, self.trig_current_hist_err,
                              fmt='bo-', capsize=2, ms=2, lw=1,
                              errorevery=3, elinewidth=0.5, label='Current')
        self.trig_ax.set_title('Histogram of trigger times')
        self.trig_ax.set_xlabel(r'Time ($\mu$s)')
        self.trig_ax.set_ylabel('Counts')
        self.trig_ax.set_xlim(self.trig_time_bins[0], self.trig_time_bins[-1])
        self.trig_ax.legend()
        self.trigger_figure_canvas.draw()

    def plot_vmi_on(self):
        self.cb_laser_on.remove()
        self.laser_on_ax.clear()
        self.cb_laser_on = self.laser_on_figure.colorbar(self.laser_on_ax.imshow(self.vmi_laser_on))
        self.laser_on_ax.imshow(self.vmi_laser_on, extent=[self.center_vals[0] - self.center_vals[2],
                                                           self.center_vals[0] + self.center_vals[2],
                                                           self.center_vals[1] - self.center_vals[2],
                                                           self.center_vals[1] + self.center_vals[2]])
        self.laser_on_ax.set_xlabel('x (mm)')
        self.laser_on_ax.set_ylabel('y (mm)')
        self.laser_on_ax.set_title('VMI coincidence laser on')
        self.laser_on_figure_canvas.draw()

    def plot_vmi_off(self, event=None):
        self.cb_laser_off.remove()
        self.laser_off_ax.clear()
        if self.laser_off_slider_val.get() == 1:
            self.cb_laser_off = self.laser_off_figure.colorbar(self.laser_off_ax.imshow(self.circ_vmi))
            self.laser_off_ax.imshow(self.circ_vmi, extent=[self.center_vals[0] - self.center_vals[2],
                                                            self.center_vals[0] + self.center_vals[2],
                                                            self.center_vals[1] - self.center_vals[2],
                                                            self.center_vals[1] + self.center_vals[2]])
            self.laser_off_ax.set_xlabel('x (mm)')
            self.laser_off_ax.set_ylabel('y (mm)')
            self.laser_off_ax.set_title('VMI Circularized')
        else:
            self.cb_laser_off = self.laser_off_figure.colorbar(self.laser_off_ax.imshow(self.vmi_laser_off))
            self.laser_off_ax.imshow(self.vmi_laser_off, extent=[self.center_vals[0] - self.center_vals[2],
                                                                 self.center_vals[0] + self.center_vals[2],
                                                                 self.center_vals[1] - self.center_vals[2],
                                                                 self.center_vals[1] + self.center_vals[2]])
            self.laser_off_ax.set_xlabel('x (mm)')
            self.laser_off_ax.set_ylabel('y (mm)')
            self.laser_off_ax.set_title('VMI coincidence laser off')
        self.laser_off_figure_canvas.draw()

    def plot_rad_bin(self):
        self.rad_ax.clear()
        self.rad_ax_energy.clear()
        self.rad_ax.axhline(0, color='k', lw=0.5)
        self.rad_ax.plot(self.radial_projection_values, self.radial_projection_amplitude_sum, label='Total')
        for n in range(int(self.norder) + 1):
            self.rad_ax.plot(self.radial_projection_values,
                             self.radial_projection_amplitude[n * int(self.resolution):(n + 1) * int(self.resolution)],
                             label=f'Harmonic {n}')
        if self.cali_use_import_button_var.get() == 'True':
            self.rad_ax_energy.clear()
            self.rad_ax_energy = self.rad_ax.secondary_xaxis('top', functions=(self.radius_to_energy_func['converter'],
                                                                               self.radius_to_energy_func['inverter']))
            self.rad_ax_energy.set_xlabel('Energy (eV)')
        self.rad_ax.set_title('Radially binned data')
        self.rad_ax.set_xlabel('Radius (mm)')
        self.rad_ax.set_ylabel('Integration of angular orders')
        self.rad_ax.legend()
        self.radial_figure_canvas.draw()

    def plot_intensity(self):
        self.intensity_ax.clear()
        if self.cali_use_import_button_var.get() == 'True':
            self.intensity_ax.axhline(0, color='k', lw=0.5)
            rad_err = np.ones(len(self.intensity)) * self.sys_radial_err
            print(self.radius_to_energy_func['converter'](np.array(self.radial_projection_values)))
            print(self.radius_to_energy_func['error'](np.array(self.radial_projection_values), rad_err))
            self.intensity_ax.errorbar(self.radius_to_energy_func['converter'](np.array(self.radial_projection_values)),
                                       self.intensity, self.intensity_err,
                                       self.radius_to_energy_func['error'](np.array(self.radial_projection_values),
                                                                           rad_err),
                                       fmt='ro', ms=1, lw=0.5, capsize=1.2,
                                       errorevery=1, elinewidth=0.25)
            self.intensity_ax.set_title('Abel Inversed Intensity')
            self.intensity_ax.set_xlabel('Energy (eV)')
            self.intensity_ax.set_ylabel(r'Integration over $\theta^\prime$')
            self.intensity_figure_canvas.draw()
        else:
            self.intensity_ax.axhline(0, color='k', lw=0.5)
            rad_err = np.ones(len(self.intensity)) * self.sys_radial_err
            self.intensity_ax.errorbar(self.radial_projection_values, self.intensity, self.intensity_err, rad_err,
                                       fmt='ro', ms=1, lw=0.5, capsize=1.2,
                                       errorevery=1, elinewidth=0.25)

            self.intensity_ax.set_title('Abel Inversed Intensity')
            self.intensity_ax.set_xlabel('3D Radius (mm)')
            self.intensity_ax.set_ylabel(r'Integration over $\theta^\prime$')
            self.intensity_figure_canvas.draw()

    def plot_bootstrap(self, dot=None):
        if self.cali_use_import_button_var.get() == 'True':
            self.boot_intensity_ax.clear()

            self.boot_intensity_ax.axhline(0, color='k', lw=0.5)
            self.boot_intensity_ax.errorbar(self.radius_to_energy_func['converter'](np.array(self.boot_radial)),
                                            self.boot_intensity, self.quantiles,
                                            self.radius_to_energy_func['error'](np.array(self.boot_radial),
                                                                                np.array(self.boot_rad_err)),
                                            fmt='ro', ms=1, lw=0.5, capsize=1.2,
                                            errorevery=1, elinewidth=0.25)
            if dot is not None:
                self.boot_intensity_ax.axvline(self.radius_to_energy_func['error'](self.boot_radial[dot]), color='k',
                                               lw=0.5)

            self.boot_intensity_ax.set_title('Bootstrapped Intensity')
            self.boot_intensity_ax.set_xlabel('Energy (eV)')
            self.boot_intensity_ax.set_ylabel(r'Integration over $\theta^\prime$')
            self.boot_intensity_figure_canvas.draw()
        else:
            self.boot_intensity_ax.clear()

            self.boot_intensity_ax.axhline(0, color='k', lw=0.5)
            self.boot_intensity_ax.errorbar(self.boot_radial, self.boot_intensity, self.quantiles, self.boot_rad_err,
                                            fmt='ro', ms=1, lw=0.5, capsize=1.2,
                                            errorevery=1, elinewidth=0.25)

            if dot is not None:
                self.boot_intensity_ax.axvline(self.boot_radial[dot], color='k', lw=0.5)

            self.boot_intensity_ax.set_title('Bootstrapped Intensity')
            self.boot_intensity_ax.set_xlabel('3D Radius (mm)')
            self.boot_intensity_ax.set_ylabel(r'Integration over $\theta^\prime$')
            self.boot_intensity_figure_canvas.draw()

    def plot_bootstrap_hist(self, stat):
        self.boot_statistic_ax.clear()

        self.boot_statistic_ax.hist(self.boot_iterations[stat], bins=100, color='b')

        self.boot_statistic_ax.axvline(self.boot_intensity[stat], color='k')
        self.boot_statistic_ax.axvline(self.intensity_err_lower[stat], color='r')
        self.boot_statistic_ax.axvline(self.intensity_err_upper[stat], color='r')

        self.boot_statistic_ax.set_title('Bootstrap\'ed statistic at %.3f mm' % self.boot_radial[stat])
        self.boot_statistic_ax.set_ylabel('Counts')
        self.boot_statistic_ax.set_xlabel(r'Integration over $\theta^\prime$')
        self.boot_statistic_figure_canvas.draw()

    def plot_abel_vmi(self):
        self.cb_abel_vmi.remove()
        self.abel_vmi_ax.clear()
        self.cb_abel_vmi = self.abel_vmi_figure.colorbar(self.abel_vmi_ax.imshow(self.combined_vmi))
        self.abel_vmi_ax.imshow(self.combined_vmi, extent=[self.center_vals[0] - self.center_vals[2],
                                                           self.center_vals[0] + self.center_vals[2],
                                                           self.center_vals[1] - self.center_vals[2],
                                                           self.center_vals[1] + self.center_vals[2]])
        self.abel_vmi_ax.set_xlabel('x (mm)')
        self.abel_vmi_ax.set_ylabel('y (mm)')
        self.abel_vmi_ax.set_title('Abel Inversed VMI')
        self.abel_vmi_figure_canvas.draw()

    def update_vmi(self):
        root_file = ROOT.TFile(self.root_file, "read")
        vmi_laser_on_data = root_file.Get("VMI_LASER_ON")
        vmi_laser_off_data = root_file.Get("VMI_LASER_OFF")
        self.vmi_laser_on = np.zeros((int(self.vmi_resolution), int(self.vmi_resolution)))
        self.vmi_laser_off = np.zeros((int(self.vmi_resolution), int(self.vmi_resolution)))
        self.circ_vmi = np.zeros((int(self.vmi_resolution), int(self.vmi_resolution)))
        count = 0
        for entry in vmi_laser_on_data:
            self.vmi_laser_on[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
            count += 1
        count = 0
        for entry in vmi_laser_off_data:
            self.vmi_laser_off[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
            count += 1
        if self.circularize:
            vmi_circ_data = root_file.Get("VMI_CIRC")
            count = 0
            for entry in vmi_circ_data:
                self.circ_vmi[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
                count += 1
        root_file.Close()
        self.plot_vmi_on()
        self.plot_vmi_off()

    def update_analysed_data(self):
        root_file = ROOT.TFile(self.root_file, "read")
        radial_bin_val_data = root_file.Get("radial_bin_values")
        radial_bin_amp_data = root_file.Get("radial_bin_TOTAL")
        radial_bin_err_data = root_file.Get("radial_bin_error")
        vmi_abel_data = root_file.Get("VMI_ABEL")
        vmi_fold_data = root_file.Get("VMI_FOLD")
        intensity_data = root_file.Get("INTENSITY")
        intensity_err_data = root_file.Get("INTENSITY_ERROR")

        center_data = root_file.Get("CENTER")

        self.center_vals = []
        for entry in center_data:
            self.center_vals.append(entry.center)

        self.abel_vmi = np.zeros((int(self.vmi_resolution), int(self.vmi_resolution)))
        self.vmi_fold = np.zeros((int(self.vmi_resolution), int(self.vmi_resolution)))
        self.circ_vmi = np.zeros((int(self.vmi_resolution), int(self.vmi_resolution)))

        # Get the radial binning of data form file
        self.radial_projection_values.clear()
        self.radial_projection_amplitude.clear()
        self.radial_projection_error.clear()
        self.radial_projection_amplitude_sum = np.zeros(int(self.resolution))
        self.radial_projection_error_sum = np.zeros(int(self.resolution))
        count_amp = 0
        count_err = 0
        for entry in radial_bin_val_data:
            self.radial_projection_values.append(entry.radial_bin_val)
        for entry in radial_bin_amp_data:
            self.radial_projection_amplitude.append(entry.radial_bin)
            self.radial_projection_amplitude_sum[count_amp] += entry.radial_bin
            count_amp += 1
            if count_amp == int(self.resolution):
                count_amp = 0
        for entry in radial_bin_err_data:
            self.radial_projection_error.append(entry.radial_bin_err)
            self.radial_projection_error_sum[count_err] += entry.radial_bin_err
            count_err += 1
            if count_err == int(self.resolution):
                count_err = 0
        self.radial_projection_error_sum = np.sqrt(np.array(self.radial_projection_error_sum))

        # Get abel inverse image from file
        count = 0
        for entry in vmi_abel_data:
            self.abel_vmi[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
            count += 1
        self.abel_vmi[np.array(self.abel_vmi) < 0] = 0
        count = 0
        for entry in vmi_fold_data:
            self.vmi_fold[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
            count += 1

        if self.circularize:
            self.laser_off_slider['state'] = 'normal'
            vmi_circ_data = root_file.Get("VMI_CIRC")
            count = 0
            for entry in vmi_circ_data:
                self.circ_vmi[count // int(self.vmi_resolution)][count % int(self.vmi_resolution)] = entry.image
                count += 1

        # Combine the VMI
        center_col = int(int(self.vmi_resolution) / 2)
        denom = np.max(self.abel_vmi[:, :center_col])
        self.combined_vmi = np.zeros((int(self.vmi_resolution), int(self.vmi_resolution)))
        self.combined_vmi[:, center_col:] = self.vmi_fold[:, center_col:]
        if denom != 0:
            self.combined_vmi[:, :center_col] = self.abel_vmi[:, :center_col] / denom * \
                                                np.max(self.vmi_fold[:, center_col:])
        else:
            self.combined_vmi[:, :center_col] = self.abel_vmi[:, :center_col]

        # Get abel inversed intensities
        self.intensity.clear()
        for entry in intensity_data:
            self.intensity.append(entry.intensity)

        self.intensity_err.clear()
        for entry in intensity_err_data:
            self.intensity_err.append(entry.error)

        root_file.Close()

        self.plot_rad_bin()
        self.plot_abel_vmi()
        self.plot_intensity()

    def update_time_hist(self):
        '''             Read data               '''

        # Point to the dataframe
        root_file = ROOT.TFile(self.root_file, "read")
        trig_hist_data = root_file.Get("TRIG_TOT")
        event_hist_data = root_file.Get("EVENT_TOT")
        peak_hist_data = root_file.Get("PEAK_TOT")
        trig_hist_current_data = root_file.Get("TRIG_NOW")
        event_hist_current_data = root_file.Get("EVENT_NOW")
        peak_hist_current_data = root_file.Get("PEAK_NOW")
        freq_hist_data = root_file.Get("FREQ_TOT")
        freq_bin_data = root_file.Get("FREQ_BIN")

        settings_bin_data = root_file.Get("SETTINGS_BINS")

        settings_bins = []
        for entry in settings_bin_data:
            settings_bins.append(entry.settings_bin)

        self.event_time_bin_size = settings_bins[0]
        self.trig_time_bin_size = settings_bins[1] / 1e-6  # in micro seconds
        self.peak_sum_bin_size = settings_bins[2]
        self.freq_bin_size = settings_bins[3]

        self.event_hist = []
        self.event_hist_err = []
        self.event_current_hist = []
        self.event_current_hist_err = []

        self.trig_hist = []
        self.trig_hist_err = []
        self.trig_current_hist = []
        self.trig_current_hist_err = []

        self.peak_sum_hist = []
        self.peak_sum_hist_err = []
        self.peak_sum_current_hist = []
        self.peak_sum_current_hist_err = []

        self.event_time_bins = [0]
        self.trig_time_bins = [0]
        self.peak_sum_bins = [0]
        self.freq_bins = []

        # Get the event and trigger data from file
        for entry in trig_hist_data:
            self.trig_hist.append(entry.trig_hist)
            self.trig_hist_err.append(np.sqrt(entry.trig_hist))
            self.trig_time_bins.append(self.trig_time_bins[-1] + self.trig_time_bin_size)
        for entry in trig_hist_current_data:
            self.trig_current_hist.append(entry.trig_hist)
            self.trig_current_hist_err.append(np.sqrt(entry.trig_hist))
        for entry in event_hist_data:
            self.event_hist.append(entry.event_hist)
            self.event_hist_err.append(np.sqrt(entry.event_hist))
            self.event_time_bins.append(self.event_time_bins[-1] + self.event_time_bin_size)
        for entry in event_hist_current_data:
            self.event_current_hist.append(entry.event_hist)
            self.event_current_hist_err.append(np.sqrt(entry.event_hist))
        for entry in peak_hist_data:
            self.peak_sum_hist.append(entry.peak_hist)
            self.peak_sum_hist_err.append(np.sqrt(entry.peak_hist))
            self.peak_sum_bins.append(self.peak_sum_bins[-1] + self.peak_sum_bin_size)
        for entry in peak_hist_current_data:
            self.peak_sum_current_hist.append(entry.peak_hist)
            self.peak_sum_current_hist_err.append(np.sqrt(entry.peak_hist))
        for entry in freq_hist_data:
            self.freq_hist.append(entry.freq_hist)
            self.freq_hist_err.append(np.sqrt(entry.freq_hist))
        for entry in freq_bin_data:
            self.freq_bins.append(entry.freq_bin * 4.136 / 1000)  # Convert from THz to eV
        self.event_time_bins.pop(-1)
        self.trig_time_bins.pop(-1)
        self.peak_sum_bins.pop(-1)

        root_file.Close()

    def set_times(self):
        self.trig_upper = self.trig_upper_entry.get()
        self.trig_lower = self.trig_lower_entry.get()
        trig_low = str(float(self.trig_lower) * 1e-6)
        trig_up = str(float(self.trig_upper) * 1e-6)
        self.event_upper = self.event_upper_entry.get()
        self.event_lower = self.event_lower_entry.get()
        self.peak_sum_upper = self.peak_sum_upper_entry.get()
        self.peak_sum_lower = self.peak_sum_lower_entry.get()
        if self.circularize:
            if len(self.circ_import_file) != 0:
                VMI_filter(self.vmi_file_path, self.root_file, ['-t' + trig_low, ' -T' + trig_up +
                                                                ' -e' + self.event_lower + ' -E' + self.event_upper +
                                                                ' -p' + self.peak_sum_lower + ' -P' + self.peak_sum_upper +
                                                                '-C' + self.circ_import_file])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-t' + trig_low, ' -T' + trig_up +
                                                                ' -e' + self.event_lower + ' -E' + self.event_upper +
                                                                ' -p' + self.peak_sum_lower + ' -P' + self.peak_sum_upper +
                                                                ' -c' + str(self.nr_angular_bin)])
        else:
            VMI_filter(self.vmi_file_path, self.root_file, [' -t' + trig_low, ' -T' + trig_up +
                                                            ' -e' + self.event_lower + ' -E' + self.event_upper +
                                                            ' -p' + self.peak_sum_lower + ' -P' + self.peak_sum_upper])
        get_all_time_data(self.event_file_path, self.root_file, [' -t' + trig_low, ' -T' + trig_up +
                                                                 ' -e' + self.event_lower + ' -E' + self.event_upper +
                                                                 ' -p' + self.peak_sum_lower + ' -P' + self.peak_sum_upper])
        self.update_vmi()
        self.update_analysed_data()
        self.update_time_hist()
        self.plot_trig()
        self.plot_event()
        self.plot_peak_sum()

    def check_times(self):
        self.trig_upper = self.trig_upper_entry.get()
        self.trig_lower = self.trig_lower_entry.get()
        self.event_upper = self.event_upper_entry.get()
        self.event_lower = self.event_lower_entry.get()
        self.peak_sum_upper = self.peak_sum_upper_entry.get()
        self.peak_sum_lower = self.peak_sum_lower_entry.get()
        self.plot_trig()
        self.plot_event()
        self.plot_peak_sum()

    def calc_boot(self):
        boot_popup = BootPopup(self, self.boot_file_path, self.root_file)
        self.wait_window(boot_popup)

    def change_circular(self):
        self.nr_angular_bin = self.circular_entry.get()
        if self.circular_str_truth.get() == 'True':
            self.circular_str_truth.set('False')
            self.circularize = False
            VMI_filter(self.vmi_file_path, self.root_file, ['-c -1'])
        else:
            self.circular_str_truth.set('True')
            self.circularize = True
            if len(self.circ_import_file) != 0:
                VMI_filter(self.vmi_file_path, self.root_file, [' -C' + self.circ_import_file])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-c ' + self.nr_angular_bin])
        self.update_analysed_data()
        self.plot_vmi_off()

    def set_circular(self):
        self.nr_angular_bin = self.circular_entry.get()
        if self.circularize:
            VMI_filter(self.vmi_file_path, self.root_file, ['-c ' + self.nr_angular_bin])
            self.update_analysed_data()
            self.plot_vmi_off()
        else:
            pass

    def set_settings(self):
        max_changed = False
        if self.max_rad != self.resolution_max_entry.get():
            self.max_rad = self.resolution_max_entry.get()
            max_changed = True

        if self.resolution != self.resolution_entry.get() and \
                self.norder != self.norder_entry.get() and \
                self.vmi_resolution != self.vmi_resolution_entry.get():
            self.resolution = self.resolution_entry.get()
            self.norder = self.norder_entry.get()
            self.vmi_resolution = self.vmi_resolution_entry.get()
            if self.circularize:
                if len(self.circ_import_file) != 0:
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-o' + self.norder + " " +
                                                                    '-v' + self.vmi_resolution + " " +
                                                                    '-C' + self.circ_import_file + " " +
                                                                    '-m' + self.max_rad])
                else:
                    self.nr_angular_bin = self.circular_entry.get()
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-o' + self.norder + " " +
                                                                    '-v' + self.vmi_resolution + " " +
                                                                    '-c' + self.nr_angular_bin + " " +
                                                                    '-m' + self.max_rad])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-r' + self.resolution + " " +
                                                                '-o' + self.norder + " " +
                                                                '-v' + self.vmi_resolution + " " +
                                                                '-m' + self.max_rad])

        elif self.norder != self.norder_entry.get() and \
                self.vmi_resolution != self.vmi_resolution_entry.get():
            self.norder = self.norder_entry.get()
            self.vmi_resolution = self.vmi_resolution_entry.get()
            if self.circularize:
                if len(self.circ_import_file) != 0:
                    VMI_filter(self.vmi_file_path, self.root_file, ['-o' + self.norder + " " +
                                                                    '-v' + self.vmi_resolution + " " +
                                                                    '-C' + self.circ_import_file + " " +
                                                                    '-m' + self.max_rad])
                else:
                    self.nr_angular_bin = self.circular_entry.get()
                    VMI_filter(self.vmi_file_path, self.root_file, ['-o' + self.norder + " " +
                                                                    '-v' + self.vmi_resolution + " " +
                                                                    '-c' + self.nr_angular_bin + " " +
                                                                    '-m' + self.max_rad])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-o' + self.norder + " "
                                                                                     '-v' + self.vmi_resolution + " "
                                                                                                                  '-m' + self.max_rad])

        elif self.resolution != self.resolution_entry.get() and \
                self.vmi_resolution != self.vmi_resolution_entry.get():
            self.resolution = self.resolution_entry.get()
            self.vmi_resolution = self.vmi_resolution_entry.get()
            if self.circularize:
                if len(self.circ_import_file) != 0:
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-v' + self.vmi_resolution + " " +
                                                                    '-C' + self.circ_import_file + " " +
                                                                    '-m' + self.max_rad])
                else:
                    self.nr_angular_bin = self.circular_entry.get()
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-v' + self.vmi_resolution + " " +
                                                                    '-c' + self.nr_angular_bin + " " +
                                                                    '-m' + self.max_rad])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-r' + self.resolution + " "
                                                                                         '-v' + self.vmi_resolution + " "
                                                                                                                      '-m' + self.max_rad])

        elif self.resolution != self.resolution_entry.get() and \
                self.norder != self.norder_entry.get():
            self.resolution = self.resolution_entry.get()
            self.norder = self.norder_entry.get()
            if self.circularize:
                if len(self.circ_import_file) != 0:
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-o' + self.norder + " " +
                                                                    '-C' + self.circ_import_file + " " +
                                                                    '-m' + self.max_rad])
                else:
                    self.nr_angular_bin = self.circular_entry.get()
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-o' + self.norder + " " +
                                                                    '-c' + self.nr_angular_bin + " " +
                                                                    '-m' + self.max_rad])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-r' + self.resolution + " "
                                                                                         '-o' + self.norder + " "
                                                                                                              '-m' + self.max_rad])

        elif self.vmi_resolution != self.vmi_resolution_entry.get():
            self.vmi_resolution = self.vmi_resolution_entry.get()
            if self.circularize:
                if len(self.circ_import_file) != 0:
                    VMI_filter(self.vmi_file_path, self.root_file, ['-v' + self.vmi_resolution + " " +
                                                                    '-C' + self.circ_import_file + " " +
                                                                    '-m' + self.max_rad])
                else:
                    self.nr_angular_bin = self.circular_entry.get()
                    VMI_filter(self.vmi_file_path, self.root_file, ['-v' + self.vmi_resolution + " " +
                                                                    '-c' + self.nr_angular_bin + " " +
                                                                    '-m' + self.max_rad])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-v' + self.vmi_resolution + " "
                                                                                             '-m' + self.max_rad])

        elif self.resolution != self.resolution_entry.get():
            self.resolution = self.resolution_entry.get()
            if self.circularize:
                if len(self.circ_import_file) != 0:
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-C' + self.circ_import_file + " " +
                                                                    '-m' + self.max_rad])
                else:
                    self.nr_angular_bin = self.circular_entry.get()
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-c' + self.nr_angular_bin + " " +
                                                                    '-m' + self.max_rad])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-r' + self.resolution + " "
                                                                                         '-m' + self.max_rad])

        elif self.norder != self.norder_entry.get():
            self.norder = self.norder_entry.get()
            if self.circularize:
                if len(self.circ_import_file) != 0:
                    VMI_filter(self.vmi_file_path, self.root_file, ['-o' + self.norder + " " +
                                                                    '-C' + self.circ_import_file + " " +
                                                                    '-m' + self.max_rad])
                else:
                    self.nr_angular_bin = self.circular_entry.get()
                    VMI_filter(self.vmi_file_path, self.root_file, ['-o' + self.norder + " " +
                                                                    '-c' + self.nr_angular_bin + " " +
                                                                    '-m' + self.max_rad])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-o ' + self.norder + " "
                                                                                      '-m ' + self.max_rad])

        elif max_changed:
            self.resolution = self.resolution_entry.get()
            self.norder = self.norder_entry.get()
            self.vmi_resolution = self.vmi_resolution_entry.get()
            if self.circularize:
                if len(self.circ_import_file) != 0:
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-o' + self.norder + " " +
                                                                    '-v' + self.vmi_resolution + " " +
                                                                    '-C' + self.circ_import_file + " " +
                                                                    '-m' + self.max_rad])
                else:
                    self.nr_angular_bin = self.circular_entry.get()
                    VMI_filter(self.vmi_file_path, self.root_file, [' -r' + self.resolution + " " +
                                                                    '-o' + self.norder + " " +
                                                                    '-v' + self.vmi_resolution + " " +
                                                                    '-c' + self.nr_angular_bin + " " +
                                                                    '-m' + self.max_rad])
            else:
                VMI_filter(self.vmi_file_path, self.root_file, ['-r' + self.resolution + " "
                                                                                         '-o' + self.norder + " "
                                                                                                              '-v' + self.vmi_resolution + " "
                                                                                                                                           '-m' + self.max_rad])

        self.update_analysed_data()
        self.update_vmi()

    def calc_gif(self):
        gif_popup = GifPopup(self, self.gif_file_path, self.root_file)
        self.wait_window(gif_popup)

    def update_gif(self):
        if self.run_gif:
            try:
                self.frame_time_label_upper_var.set(f'{self.gif_times[1][self.current_frame]}')
                self.frame_time_label_lower_var.set(f'{self.gif_times[0][self.current_frame]} - ')

                self.raw_vmi_frames.seek(self.current_frame)
                self.raw_off_vmi_frames.seek(self.current_frame)
                if self.circ_frames:
                    self.raw_circ_vmi_frames.seek(self.current_frame)
                self.inv_vmi_frames.seek(self.current_frame)
                self.intensity_frames.seek(self.current_frame)
                self.radial_frames.seek(self.current_frame)
                self.event_frames.seek(self.current_frame)

                raw_vmi_image = self.raw_vmi_frames.copy()
                raw_off_vmi_image = self.raw_off_vmi_frames.copy()
                if self.circ_frames:
                    raw_circ_vmi_image = self.raw_circ_vmi_frames.copy()
                inv_vmi_image = self.inv_vmi_frames.copy()
                intensity_image = self.intensity_frames.copy()
                radial_image = self.radial_frames.copy()
                event_image = self.event_frames.copy()

                raw_vmi_image = raw_vmi_image.resize((self.GIFVMI.winfo_width(), self.GIFVMI.winfo_height()),
                                                     Image.LANCZOS)
                raw_off_vmi_image = raw_off_vmi_image.resize(
                    (self.GIFVMIOff.winfo_width(), self.GIFVMIOff.winfo_height()),
                    Image.LANCZOS)
                if self.circ_frames:
                    raw_circ_vmi_image = raw_circ_vmi_image.resize(
                        (self.GIFVMICirc.winfo_width(), self.GIFVMICirc.winfo_height()),
                        Image.LANCZOS)
                inv_vmi_image = inv_vmi_image.resize((self.GIFAbelVMI.winfo_width(), self.GIFAbelVMI.winfo_height()),
                                                     Image.LANCZOS)
                intensity_image = intensity_image.resize(
                    (self.GIFIntensity.winfo_width(), self.GIFIntensity.winfo_height()),
                    Image.LANCZOS)
                radial_image = radial_image.resize(
                    (self.GIFRadial.winfo_width(), self.GIFRadial.winfo_height()),
                    Image.LANCZOS)
                event_image = event_image.resize((self.GIFEvent.winfo_width(), self.GIFEvent.winfo_height()),
                                                 Image.LANCZOS)

                raw_vmi_image = ImageTk.PhotoImage(raw_vmi_image)
                raw_off_vmi_image = ImageTk.PhotoImage(raw_off_vmi_image)
                if self.circ_frames:
                    raw_circ_vmi_image = ImageTk.PhotoImage(raw_circ_vmi_image)
                inv_vmi_image = ImageTk.PhotoImage(inv_vmi_image)
                intensity_image = ImageTk.PhotoImage(intensity_image)
                radial_image = ImageTk.PhotoImage(radial_image)
                event_image = ImageTk.PhotoImage(event_image)

                self.raw_vmi_gif.configure(image=raw_vmi_image)
                self.raw_vmi_gif.image = raw_vmi_image
                self.raw_off_vmi_gif.configure(image=raw_off_vmi_image)
                self.raw_off_vmi_gif.image = raw_off_vmi_image
                if self.circ_frames:
                    self.raw_circ_vmi_gif.configure(image=raw_circ_vmi_image)
                    self.raw_circ_vmi_gif.image = raw_circ_vmi_image
                self.inv_vmi_gif.configure(image=inv_vmi_image)
                self.inv_vmi_gif.image = inv_vmi_image
                self.intensity_gif.configure(image=intensity_image)
                self.intensity_gif.image = intensity_image
                self.radial_gif.configure(image=radial_image)
                self.radial_gif.image = radial_image
                self.event_gif.configure(image=event_image)
                self.event_gif.image = event_image

                self.frame_slider_val.set(self.current_frame + 1)
                self.frame_slider_str.set(f'{self.current_frame + 1}')
                self.current_frame += 1
            except:
                self.current_frame = 0
                self.frame_slider_val.set(self.current_frame + 1)
                self.frame_slider_str.set(f'{self.current_frame + 1}')

                self.raw_vmi_frames.seek(self.current_frame)
                self.raw_off_vmi_frames.seek(self.current_frame)
                if self.circ_frames:
                    self.raw_circ_vmi_frames.seek(self.current_frame)
                self.inv_vmi_frames.seek(self.current_frame)
                self.intensity_frames.seek(self.current_frame)
                self.radial_frames.seek(self.current_frame)
                self.event_frames.seek(self.current_frame)

                raw_vmi_image = self.raw_vmi_frames.copy()
                raw_off_vmi_image = self.raw_off_vmi_frames.copy()
                if self.circ_frames:
                    raw_circ_vmi_image = self.raw_circ_vmi_frames.copy()
                inv_vmi_image = self.inv_vmi_frames.copy()
                intensity_image = self.intensity_frames.copy()
                radial_image = self.radial_frames.copy()
                event_image = self.event_frames.copy()

                raw_vmi_image = raw_vmi_image.resize((self.GIFVMI.winfo_width(), self.GIFVMI.winfo_height()),
                                                     Image.LANCZOS)
                raw_off_vmi_image = raw_off_vmi_image.resize(
                    (self.GIFVMIOff.winfo_width(), self.GIFVMIOff.winfo_height()),
                    Image.LANCZOS)
                if self.circ_frames:
                    raw_circ_vmi_image = raw_circ_vmi_image.resize(
                        (self.GIFVMICirc.winfo_width(), self.GIFVMICirc.winfo_height()),
                        Image.LANCZOS)
                inv_vmi_image = inv_vmi_image.resize((self.GIFAbelVMI.winfo_width(), self.GIFAbelVMI.winfo_height()),
                                                     Image.LANCZOS)
                intensity_image = intensity_image.resize(
                    (self.GIFIntensity.winfo_width(), self.GIFIntensity.winfo_height()),
                    Image.LANCZOS)
                radial_image = radial_image.resize(
                    (self.GIFRadial.winfo_width(), self.GIFRadial.winfo_height()),
                    Image.LANCZOS)
                event_image = event_image.resize((self.GIFEvent.winfo_width(), self.GIFEvent.winfo_height()),
                                                 Image.LANCZOS)

                raw_vmi_image = ImageTk.PhotoImage(raw_vmi_image)
                raw_off_vmi_image = ImageTk.PhotoImage(raw_off_vmi_image)
                if self.circ_frames:
                    raw_circ_vmi_image = ImageTk.PhotoImage(raw_circ_vmi_image)
                inv_vmi_image = ImageTk.PhotoImage(inv_vmi_image)
                intensity_image = ImageTk.PhotoImage(intensity_image)
                radial_image = ImageTk.PhotoImage(radial_image)
                event_image = ImageTk.PhotoImage(event_image)

                self.raw_vmi_gif.configure(image=raw_vmi_image)
                self.raw_vmi_gif.image = raw_vmi_image
                self.raw_off_vmi_gif.configure(image=raw_off_vmi_image)
                self.raw_off_vmi_gif.image = raw_off_vmi_image
                if self.circ_frames:
                    self.raw_circ_vmi_gif.configure(image=raw_circ_vmi_image)
                    self.raw_circ_vmi_gif.image = raw_circ_vmi_image
                self.inv_vmi_gif.configure(image=inv_vmi_image)
                self.inv_vmi_gif.image = inv_vmi_image
                self.intensity_gif.configure(image=intensity_image)
                self.intensity_gif.image = intensity_image
                self.radial_gif.configure(image=radial_image)
                self.radial_gif.image = radial_image
                self.event_gif.configure(image=event_image)
                self.event_gif.image = event_image

            # Call again
            try:
                if int(self.wait_entry.get()) > 1:
                    self.Timeresolved.after(int(self.wait_entry.get()), self.update_gif)
                else:
                    self.Timeresolved.after(250, self.update_gif)
            except IOError:
                self.Timeresolved.after(250, self.update_gif)

    def play_pause_gif(self):
        if self.pause_play_button_str.get() == 'Play':
            self.run_gif = True
            self.update_gif()
            self.pause_play_button_str.set('Pause')
            self.frame_slider['state'] = 'disabled'
            self.gif_frame_slider['state'] = 'disabled'
        else:
            self.run_gif = False
            self.pause_play_button_str.set('Play')
            self.frame_slider['state'] = 'normal'
            self.gif_frame_slider['state'] = 'normal'

    def slider_changed(self, event):
        if self.current_frame != self.frame_slider_val.get() - 1:
            self.current_frame = self.frame_slider_val.get() - 1
            self.gif_frame_slider_val.set(self.frame_slider_val.get())
        else:
            self.current_frame = self.gif_frame_slider_val.get() - 1
            self.frame_slider_val.set(self.gif_frame_slider_val.get())

        self.frame_slider_str.set(f'{self.current_frame + 1}')

        self.frame_time_label_upper_var.set(f'{self.gif_times[1][self.current_frame]}')
        self.frame_time_label_lower_var.set(f'{self.gif_times[0][self.current_frame]} - ')

        self.raw_vmi_frames.seek(self.current_frame)
        self.raw_off_vmi_frames.seek(self.current_frame)
        if self.circ_frames:
            self.raw_circ_vmi_frames.seek(self.current_frame)
        self.inv_vmi_frames.seek(self.current_frame)
        self.intensity_frames.seek(self.current_frame)
        self.radial_frames.seek(self.current_frame)
        self.event_frames.seek(self.current_frame)

        raw_vmi_image = self.raw_vmi_frames.copy()
        raw_off_vmi_image = self.raw_off_vmi_frames.copy()
        if self.circ_frames:
            raw_circ_vmi_image = self.raw_circ_vmi_frames.copy()
        inv_vmi_image = self.inv_vmi_frames.copy()
        intensity_image = self.intensity_frames.copy()
        radial_image = self.radial_frames.copy()
        event_image = self.event_frames.copy()

        raw_vmi_image = raw_vmi_image.resize((self.GIFVMI.winfo_width(), self.GIFVMI.winfo_height()),
                                             Image.LANCZOS)
        raw_off_vmi_image = raw_off_vmi_image.resize((self.GIFVMIOff.winfo_width(), self.GIFVMIOff.winfo_height()),
                                                     Image.LANCZOS)
        if self.circ_frames:
            raw_circ_vmi_image = raw_circ_vmi_image.resize(
                (self.GIFVMICirc.winfo_width(), self.GIFVMICirc.winfo_height()),
                Image.LANCZOS)
        inv_vmi_image = inv_vmi_image.resize((self.GIFAbelVMI.winfo_width(), self.GIFAbelVMI.winfo_height()),
                                             Image.LANCZOS)
        intensity_image = intensity_image.resize(
            (self.GIFIntensity.winfo_width(), self.GIFIntensity.winfo_height()),
            Image.LANCZOS)
        radial_image = radial_image.resize(
            (self.GIFRadial.winfo_width(), self.GIFRadial.winfo_height()),
            Image.LANCZOS)
        event_image = event_image.resize((self.GIFEvent.winfo_width(), self.GIFEvent.winfo_height()),
                                         Image.LANCZOS)

        raw_vmi_image = ImageTk.PhotoImage(raw_vmi_image)
        raw_off_vmi_image = ImageTk.PhotoImage(raw_off_vmi_image)
        if self.circ_frames:
            raw_circ_vmi_image = ImageTk.PhotoImage(raw_circ_vmi_image)
        inv_vmi_image = ImageTk.PhotoImage(inv_vmi_image)
        intensity_image = ImageTk.PhotoImage(intensity_image)
        radial_image = ImageTk.PhotoImage(radial_image)
        event_image = ImageTk.PhotoImage(event_image)

        self.raw_vmi_gif.configure(image=raw_vmi_image)
        self.raw_vmi_gif.image = raw_vmi_image
        self.raw_off_vmi_gif.configure(image=raw_off_vmi_image)
        self.raw_off_vmi_gif.image = raw_off_vmi_image
        if self.circ_frames:
            self.raw_circ_vmi_gif.configure(image=raw_circ_vmi_image)
            self.raw_circ_vmi_gif.image = raw_circ_vmi_image
        self.inv_vmi_gif.configure(image=inv_vmi_image)
        self.inv_vmi_gif.image = inv_vmi_image
        self.intensity_gif.configure(image=intensity_image)
        self.intensity_gif.image = intensity_image
        self.radial_gif.configure(image=radial_image)
        self.radial_gif.image = radial_image
        self.event_gif.configure(image=event_image)
        self.event_gif.image = event_image

    def boot_slider_changed(self, event):
        if self.boot_statistic_slider_val.get() == 0:
            self.plot_bootstrap()
        else:
            self.plot_bootstrap(self.boot_statistic_slider_val.get() - 1)
            self.plot_bootstrap_hist(self.boot_statistic_slider_val.get() - 1)

    def save_data(self):
        save_data_popup = SaveDataPopup(self)
        self.wait_window(save_data_popup)

    def import_circ(self):
        if len(self.circ_import_file) != 0:
            file = filedialog.askopenfilename(parent=self, title='Choose file')
        else:
            file = filedialog.askopenfilename(parent=self, initialfile=self.circ_import_file, title='Choose file')
        if len(file) != 0:
            self.import_file_check_valid(file)
            self.circ_import_file = file.replace(" ", "\\ ")
            self.circ_import_label_var.set('File chosen')

    def export_circ(self):
        if self.circularize:
            file = open(self.root_file[:-5] + '_radial_corrections.txt', 'w')
            root_file = ROOT.TFile(self.root_file, "read")
            rad_corrections_tree = root_file.Get('RAD_CORR')
            rad_corrections_branch = rad_corrections_tree.GetBranch("rad_corrections")

            ang_bins = int(rad_corrections_tree.GetEntries() / 2)
            file.write(f'{ang_bins}\n')

            rad_corrections = np.empty((2, ang_bins))

            for i in range(rad_corrections_branch.GetEntries()):
                rad_corrections_branch.GetEntry(i)
                entry = getattr(rad_corrections_tree, 'rad_corrections/rad_corrections')
                if i < ang_bins:
                    rad_corrections[0, i] = entry
                else:
                    rad_corrections[1, i - ang_bins] = entry
            root_file.Close()
            for i in range(int(self.nr_angular_bin)):
                file.write(f'{rad_corrections[0, i]},{rad_corrections[1, i]}\n')
            file.close()

    def use_import(self):
        if len(self.circ_import_file) != 0:
            if self.circ_use_import_button_var.get() == 'True':
                self.circ_use_import_button_var.set('False')
            else:
                self.circ_use_import_button_var.set('True')
                if self.circularize:
                    VMI_filter(self.vmi_file_path, self.root_file, [' -C' + self.circ_import_file])
                    self.update_analysed_data()
                    self.plot_vmi_off()
        else:
            self.circ_use_import_button_var.set('False')

    def import_file_check_valid(self, file):
        try:
            f = open(file, 'r')
            nr_angular_bins = int(f.readline())
            f.close()
            r, a = np.loadtxt(file, delimiter=',', skiprows=1, unpack=True)
            if len(r) == nr_angular_bins and len(a) == nr_angular_bins:
                return True
            else:
                raise FileExistsError('File not valid format')
        except:
            self.import_circ()

    def set_peak_fit_preset(self, *args):
        if self.preset_dropdown_var.get() == 'Custom':
            self.peak_fit_text.configure(state='normal')
            self.peak_fit_text.delete("1.0", "end")
            self.peak_fit_text.insert(tk.END, 'import numpy as np\n\ndef f(x,*params):\n\treturn')
        else:
            self.peak_fit_text.configure(state='normal')
            if self.preset_dropdown_var.get() == 'Students-t':
                self.peak_fit_text.delete("1.0", "end")
                self.peak_fit_text.insert(tk.END, self.t_dist)
                self.peak_fit_text.configure(state='disabled')
            elif self.preset_dropdown_var.get() == 'Gauss':
                self.peak_fit_text.delete("1.0", "end")
                self.peak_fit_text.insert(tk.END, self.gaus_dist)
                self.peak_fit_text.configure(state='disabled')
            elif self.preset_dropdown_var.get() == 'Two students-t':
                self.peak_fit_text.delete("1.0", "end")
                self.peak_fit_text.insert(tk.END, self.t_dist_double)
                self.peak_fit_text.configure(state='disabled')

    def function_lock(self):
        if self.function_lock_var.get() == 'Function Locked':
            self.function_lock_var.set('Function Unlocked')
            if self.preset_dropdown_var.get() == 'Custom':
                self.peak_fit_text.configure(state='normal')
            self.preset_dropdown.configure(state='normal')
            for widgets in self.parameters_table.winfo_children():
                widgets.destroy()
            for widgets in self.parameters_fit_table.winfo_children():
                if not isinstance(widgets, tk.Button):
                    widgets.destroy()
                else:
                    widgets.grid_forget()
            self.parameters_fit_table.grid_forget()
            self.peak_fit_button.grid_forget()
            self.peak_fit_settings_button.grid_forget()
            self.peak_fit_function = {}
        else:
            try:
                # Checking function
                self.peak_fit_function = {}
                exec(self.peak_fit_text.get("1.0", tk.END), self.peak_fit_function)
                self.param_names = getfullargspec(self.peak_fit_function['f'])[0]
                self.peak_fit_function_nr_params = len(self.param_names)
                ones = np.ones(self.peak_fit_function_nr_params)
                self.peak_fit_function['f'](*ones)

                # Checking jacobian function
                if 'jac' in self.peak_fit_function:
                    self.peak_fit_function['jac'](*ones)
                    self.peak_fit_function['jac'](ones, *ones[1:])
                else:
                    self.peak_fit_function['jac'] = None

                # Checking hessian function
                if 'hess' in self.peak_fit_function:
                    pass  # self.peak_fit_function['hess'](*ones)
                    # self.peak_fit_function['hess'](ones, *ones[1:])
                else:
                    self.peak_fit_function['hess'] = None

                # Removing error
                self.function_lock_error.configure(state='normal')
                self.function_lock_error.delete("1.0", tk.END)
                self.function_lock_error.grid_forget()

                # Parameters table
                self.generate_params_table()

                # Fitting tools
                self.radius_low_var.set(str(np.min(self.boot_radial)))
                self.radius_high_var.set(str(np.max(self.boot_radial)))
                self.peak_fit_button.grid(row=6, column=0, sticky=tk.W)
                self.peak_fit_settings_button.grid(row=6, column=0, sticky=tk.E)

                # Freezing function input
                self.peak_fit_text.configure(state='disabled')
                self.function_lock_var.set('Function Locked')
                self.preset_dropdown.configure(state='disabled')
            except Exception as e:
                print(traceback.format_exc())
                self.function_lock_error.configure(state='normal')
                self.function_lock_error.delete("1.0", tk.END)
                self.function_lock_error.insert(tk.END, str(e))
                self.function_lock_error.configure(state='disabled')
                self.function_lock_error.grid(row=4)
                self.peak_fit_button.grid_forget()
                self.peak_fit_settings_button.grid_forget()

    def generate_params_table(self):
        self.param_init.clear()
        self.param_lower.clear()
        self.param_upper.clear()
        for widgets in self.parameters_table.winfo_children():
            widgets.destroy()

        size = int(self.function_lock_button.winfo_width() / 4 / tk.font.Font(font='TkFixedFont').actual()['size'])

        for i in range(len(self.parameters_table_head)):
            label = ttk.Label(self.parameters_table, text=self.parameters_table_head[i], width=size)
            label.grid(row=0, column=i)

        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        for i in range(1, self.peak_fit_function_nr_params):
            label = ttk.Label(self.parameters_table, text=self.param_names[i], width=size)
            label.grid(row=i, column=0, sticky=tk.W)

        for i in range(1, self.peak_fit_function_nr_params):
            self.param_init.append(tk.StringVar())
            self.param_init[i - 1].set('1')
            text = tk.Entry(self.parameters_table, validate='key', validatecommand=vcmdf,
                            textvariable=self.param_init[i - 1], width=size)
            text.grid(row=i, column=1)

            self.param_lower.append(tk.StringVar())
            self.param_lower[i - 1].set('0')
            text1 = tk.Entry(self.parameters_table, validate='key', validatecommand=vcmdf,
                             textvariable=self.param_lower[i - 1], width=size)
            text1.grid(row=i, column=2)

            self.param_upper.append(tk.StringVar())
            self.param_upper[i - 1].set('1')
            text2 = tk.Entry(self.parameters_table, validate='key', validatecommand=vcmdf,
                             textvariable=self.param_upper[i - 1], width=size)
            text2.grid(row=i, column=3)

    def run_peak_fit(self):
        try:
            mask = np.array([float(self.radius_low_var.get()) <= np.array(self.boot_radial)])[0] * \
                   np.array([np.array(self.boot_radial) <= float(self.radius_high_var.get())])[0] * \
                   np.array([np.array(self.boot_rad_err) != 0])[0] * \
                   np.array([np.array(self.quantiles)[0] != 0])[0] * \
                   np.array([np.array(self.quantiles)[1] != 0])[0]

            if self.systematic_errors_var.get():
                cov_y = np.array([np.diag(np.array(self.quantiles)[0]) ** 2, np.diag(np.array(self.quantiles)[1]) ** 2])
                cov_y += np.einsum('i,j->ij', np.array(self.boot_intensity_sys_err),
                                   np.array(self.boot_intensity_sys_err))
            else:
                cov_y = np.array([np.diag(np.array(self.quantiles)[0]) ** 2, np.diag(np.array(self.quantiles)[1]) ** 2])

            if self.method.get() == 'diff_evol':
                # Call fit
                p0 = np.array([[float(i.get()), float(j.get())] for i, j in zip(self.param_lower, self.param_upper)])
                self.peak_fit = FittingRoutine(self.peak_fit_function['f'], np.array(self.boot_radial),
                                               np.array(self.boot_intensity), covy=cov_y, P0=p0,
                                               method=self.method.get(),
                                               jac=self.peak_fit_function['jac'],
                                               hess=self.peak_fit_function['hess'], mask=mask)
            else:
                # Call fit
                p0 = np.array([float(i.get()) for i in self.param_init])
                self.peak_fit = FittingRoutine(self.peak_fit_function['f'], np.array(self.boot_radial),
                                               np.array(self.boot_intensity), covy=cov_y, P0=p0,
                                               method=self.method.get(),
                                               jac=self.peak_fit_function['jac'],
                                               hess=self.peak_fit_function['hess'], mask=mask)

            # Plot
            self.cali_raw_figure, self.cali_raw_ax = self.peak_fit.Plot(xlabel='Radius (mm)', ylabel='Intensity',
                                                                        init=True, do_return=True, figsize=(6, 6),
                                                                        ms=0.5,
                                                                        capsize=1, elinewidth=0.4, lw=0.7)
            self.cali_raw_figure.set_dpi(200)
            self.cali_raw_figure_canvas.get_tk_widget().destroy()
            self.cali_raw_toolbar.destroy()
            self.cali_raw_figure_canvas = FigureCanvasTkAgg(self.cali_raw_figure, master=self.DataFitTabPlot)
            self.cali_raw_toolbar = NavigationToolbar2Tk(self.cali_raw_figure_canvas, self.DataFitTabPlot)
            self.cali_raw_toolbar.update()
            self.cali_raw_figure_canvas.draw()
            self.cali_raw_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=10, pady=10)

            print(self.peak_fit.params, self.peak_fit.NormRes, self.peak_fit.x, self.peak_fit.y)
            mask = [i in self.peak_fit.x for i in self.peak_fit.x_plot]
            print(self.peak_fit.y_plot[mask], np.sqrt(np.diag(self.peak_fit.cov_y_plot[0]))[mask],
                  np.sqrt(np.diag(self.peak_fit.cov_y_plot[1]))[mask])

            # Parameters table
            self.generate_fit_params_table()

            # Removing error
            self.function_lock_error.configure(state='normal')
            self.function_lock_error.delete("1.0", tk.END)
            self.function_lock_error.grid_forget()
        except Exception as e:
            print(traceback.format_exc())
            self.function_lock_error.configure(state='normal')
            self.function_lock_error.delete("1.0", tk.END)
            self.function_lock_error.insert(tk.END, str(e))
            self.function_lock_error.configure(state='disabled')
            self.function_lock_error.grid(row=7)

    def peak_fit_settings(self):
        peak_fit_settings_popup = PeakSettings(self)
        self.wait_window(peak_fit_settings_popup)

    def generate_fit_params_table(self):
        self.param_fit_val.clear()
        self.param_fit_error.clear()
        self.param_fit_check.clear()
        for widgets in self.parameters_fit_table.winfo_children():
            if not isinstance(widgets, tk.Button):
                widgets.destroy()
            else:
                widgets.grid_forget()

        self.parameters_fit_table.grid(row=8, sticky=tk.EW)

        size = int((self.cali_fit_text.winfo_width()) / 4 / tk.font.Font(font='TkFixedFont').actual()['size'])

        for i in range(len(self.parameters_fit_table_head)):
            label = ttk.Label(self.parameters_fit_table, text=self.parameters_fit_table_head[i], width=size)
            label.grid(row=0, column=i)

        for i in range(1, self.peak_fit_function_nr_params):
            label = ttk.Label(self.parameters_fit_table, text=self.param_names[i], width=size)
            label.grid(row=i, column=0, sticky=tk.W)

        for i in range(1, self.peak_fit_function_nr_params):
            self.param_fit_val.append(tk.StringVar())
            self.param_fit_val[i - 1].set('%1.3e' % self.peak_fit.params[i - 1])
            text = tk.Label(self.parameters_fit_table, textvariable=self.param_fit_val[i - 1], width=size)
            text.grid(row=i, column=1)

            self.param_fit_error.append(tk.StringVar())
            self.param_fit_error[i - 1].set('%1.3e' % self.peak_fit.Error[i - 1])
            text1 = tk.Label(self.parameters_fit_table, textvariable=self.param_fit_error[i - 1], width=size)
            text1.grid(row=i, column=2)

            self.param_fit_check.append(tk.IntVar())
            self.param_fit_check[i - 1].set(0)
            text2 = tk.Checkbutton(self.parameters_fit_table, onvalue=1, offvalue=0,
                                   variable=self.param_fit_check[i - 1], width=size)
            text2.grid(row=i, column=3)

        text = ttk.Label(self.parameters_fit_table, text=('P-Value : %1.2f' % self.peak_fit.Pval), width=4 * size)
        text.grid(row=self.peak_fit_function_nr_params, column=0, columnspan=4, sticky=tk.W, pady=(10, 5))

        text = ttk.Label(self.parameters_fit_table,
                         text=('Chi2, df : %1.2f, %i' % (self.peak_fit.Chi2, self.peak_fit.df)), width=4 * size)
        text.grid(row=self.peak_fit_function_nr_params + 1, column=0, columnspan=4, sticky=tk.W, pady=(0, 10))

        self.peak_fit_export.config(width=4 * size)
        self.peak_fit_export.grid(row=self.peak_fit_function_nr_params + 2, column=0, columnspan=4, sticky=tk.EW,
                                  pady=10)

    def export_peak_fit(self):
        for i in range(len(self.param_fit_check)):
            if self.param_fit_check[i].get() == 1:
                self.cali_points_push(self.param_fit_val[i].get(), self.param_fit_error[i].get())

    def cali_points_push(self, exp, exp_err):
        grid_size = self.cali_points_table.grid_size()
        size = int((self.cali_fit_text.winfo_width()) / 4 / tk.font.Font(font='TkFixedFont').actual()['size'])
        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        plus = self.cali_points_table.grid_slaves(row=grid_size[1] - 1, column=2)[0]
        minus = self.cali_points_table.grid_slaves(row=grid_size[1] - 1, column=3)[0]
        plus.grid_forget()
        minus.grid_forget()
        plus.grid(row=grid_size[1], column=2)
        minus.grid(row=grid_size[1], column=3)

        # Move all slaves one row down
        for row in range(1, grid_size[1] - 1):
            for col in range(grid_size[0]):
                slave = self.cali_points_table.grid_slaves(row=row, column=col)[0]
                slave.grid_forget()
                slave.grid(row=row + 1, column=col)

        self.cali_points_exp.insert(0, tk.StringVar())
        self.cali_points_exp[0].set(exp)
        self.cali_points_exp_er.insert(0, tk.StringVar())
        self.cali_points_exp_er[0].set(exp_err)
        self.cali_points_lit.insert(0, tk.StringVar())
        self.cali_points_lit[0].set('0')
        self.cali_points_lit_err.insert(0, tk.StringVar())
        self.cali_points_lit_err[0].set('1')

        ent1 = tk.Entry(self.cali_points_table, textvariable=self.cali_points_exp[0], width=size, validate='key',
                        validatecommand=vcmdf)
        ent1.grid(row=1, column=0)
        ent2 = tk.Entry(self.cali_points_table, textvariable=self.cali_points_exp_er[0], width=size, validate='key',
                        validatecommand=vcmdf)
        ent2.grid(row=1, column=1)
        ent3 = tk.Entry(self.cali_points_table, textvariable=self.cali_points_lit[0], width=size, validate='key',
                        validatecommand=vcmdf)
        ent3.grid(row=1, column=2)
        ent4 = tk.Entry(self.cali_points_table, textvariable=self.cali_points_lit_err[0], width=size, validate='key',
                        validatecommand=vcmdf)
        ent4.grid(row=1, column=3)

    def set_cali_fit_preset(self, *args):
        self.cali_fit_text.configure(state='normal')
        if self.cali_preset_dropdown_var.get() == 'Linear':
            self.cali_fit_text.delete("1.0", "end")
            self.cali_fit_text.insert(tk.END, self.linear)
            self.cali_fit_text.configure(state='disabled')
        elif self.cali_preset_dropdown_var.get() == 'Proportional':
            self.cali_fit_text.delete("1.0", "end")
            self.cali_fit_text.insert(tk.END, self.prop)
            self.cali_fit_text.configure(state='disabled')
        elif self.cali_preset_dropdown_var.get() == 'Power':
            self.cali_fit_text.delete("1.0", "end")
            self.cali_fit_text.insert(tk.END, self.pow)
            self.cali_fit_text.configure(state='disabled')

        # Set funciton
        self.cali_fit_function = {}
        exec(self.cali_fit_text.get("1.0", tk.END), self.cali_fit_function)
        self.cali_param_names = getfullargspec(self.cali_fit_function['f'])[0]
        self.cali_fit_function_nr_params = len(self.cali_param_names)

        # Create params table
        self.cali_generate_params_table()

    def cali_generate_params_table(self):
        self.cali_param_init.clear()
        self.cali_param_lower.clear()
        self.cali_param_upper.clear()
        for widgets in self.cali_parameters_table.winfo_children():
            widgets.destroy()

        size = int(self.function_lock_button.winfo_width() / 4 / tk.font.Font(font='TkFixedFont').actual()['size'])

        for i in range(len(self.cali_parameters_table_head)):
            label = ttk.Label(self.cali_parameters_table, text=self.cali_parameters_table_head[i], width=size)
            label.grid(row=0, column=i)

        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        for i in range(1, self.cali_fit_function_nr_params):
            label = ttk.Label(self.cali_parameters_table, text=self.cali_param_names[i], width=size)
            label.grid(row=i, column=0, sticky=tk.W)

        for i in range(1, self.cali_fit_function_nr_params):
            self.cali_param_init.append(tk.StringVar())
            self.cali_param_init[i - 1].set('1')
            text = tk.Entry(self.cali_parameters_table, validate='key', validatecommand=vcmdf,
                            textvariable=self.cali_param_init[i - 1], width=size)
            text.grid(row=i, column=1)

            self.cali_param_lower.append(tk.StringVar())
            self.cali_param_lower[i - 1].set('0')
            text1 = tk.Entry(self.cali_parameters_table, validate='key', validatecommand=vcmdf,
                             textvariable=self.cali_param_lower[i - 1], width=size)
            text1.grid(row=i, column=2)

            self.cali_param_upper.append(tk.StringVar())
            self.cali_param_upper[i - 1].set('1')
            text2 = tk.Entry(self.cali_parameters_table, validate='key', validatecommand=vcmdf,
                             textvariable=self.cali_param_upper[i - 1], width=size)
            text2.grid(row=i, column=3)

    def cali_generate_fit_params_table(self):
        self.cali_param_fit_val.clear()
        self.cali_param_fit_error.clear()
        for widgets in self.cali_parameters_fit_table.winfo_children():
            if not isinstance(widgets, tk.Button):
                widgets.destroy()
            else:
                widgets.grid_forget()

        self.cali_parameters_fit_table.grid(row=10, sticky=tk.EW)

        size = int((self.cali_fit_text.winfo_width()) / 3 / tk.font.Font(font='TkFixedFont').actual()['size'])

        for i in range(len(self.cali_parameters_fit_table_head)):
            label = ttk.Label(self.cali_parameters_fit_table, text=self.cali_parameters_fit_table_head[i], width=size)
            label.grid(row=0, column=i)

        for i in range(1, self.cali_fit_function_nr_params):
            label = ttk.Label(self.cali_parameters_fit_table, text=self.cali_param_names[i], width=size)
            label.grid(row=i, column=0, sticky=tk.W)

        for i in range(1, self.cali_fit_function_nr_params):
            self.cali_param_fit_val.append(tk.StringVar())
            self.cali_param_fit_val[i - 1].set('%1.3e' % self.cali_fit.params[i - 1])
            text = tk.Label(self.cali_parameters_fit_table, textvariable=self.cali_param_fit_val[i - 1], width=size)
            text.grid(row=i, column=1)

            self.cali_param_fit_error.append(tk.StringVar())
            self.cali_param_fit_error[i - 1].set('%1.3e' % self.cali_fit.Error[i - 1])
            text1 = tk.Label(self.cali_parameters_fit_table, textvariable=self.cali_param_fit_error[i - 1], width=size)
            text1.grid(row=i, column=2)

        text = ttk.Label(self.cali_parameters_fit_table, text=('P-Value : %1.2f' % self.cali_fit.Pval), width=3 * size)
        text.grid(row=self.cali_fit_function_nr_params, column=0, columnspan=3, sticky=tk.W, pady=(10, 5))

        text = ttk.Label(self.cali_parameters_fit_table,
                         text=('Chi2, df : %1.2f, %i' % (self.cali_fit.Chi2, self.cali_fit.df)), width=3 * size)
        text.grid(row=self.cali_fit_function_nr_params + 1, column=0, columnspan=3, sticky=tk.W, pady=(0, 10))

        self.cali_fit_export.config(width=3 * size)
        self.cali_fit_export.grid(row=self.cali_fit_function_nr_params + 2, column=0, columnspan=3, sticky=tk.EW,
                                  pady=10)

    def run_cali_fit(self):
        try:
            x = np.array([float(i.get()) for i in self.cali_points_lit])
            x_err = np.array([float(i.get()) for i in self.cali_points_lit_err])
            y = np.array([float(i.get()) for i in self.cali_points_exp])
            y_err = np.array([float(i.get()) for i in self.cali_points_exp_er])

            if self.method.get() == 'diff_evol':
                # Call fit
                p0 = np.array(
                    [[float(i.get()), float(j.get())] for i, j in zip(self.cali_param_lower, self.cali_param_upper)])
                self.cali_fit = FittingRoutine(self.cali_fit_function['f'], x, y, y_err, x_err, P0=p0,
                                               method=self.method.get(),
                                               jac=self.cali_fit_function['jac'], hess=self.cali_fit_function['hess'])
            else:
                # Call fit
                p0 = np.array([float(i.get()) for i in self.cali_param_init])
                self.cali_fit = FittingRoutine(self.cali_fit_function['f'], x, y, y_err, x_err, P0=p0,
                                               method=self.method.get(),
                                               jac=self.cali_fit_function['jac'])

            # Plot
            self.cali_plot_figure, self.cali_plot_ax = self.cali_fit.Plot(ylabel='Radius (mm)', xlabel='Energy (eV)',
                                                                          init=True, do_return=True, figsize=(6, 6),
                                                                          ms=0.5,
                                                                          capsize=1, elinewidth=0.4, lw=0.7,
                                                                          legend_loc=0)
            self.cali_plot_figure.set_dpi(200)
            self.cali_plot_figure_canvas.get_tk_widget().destroy()
            self.cali_plot_toolbar.destroy()
            self.cali_plot_figure_canvas = FigureCanvasTkAgg(self.cali_plot_figure, master=self.CalibrationFitTabPlot)
            self.cali_plot_toolbar = NavigationToolbar2Tk(self.cali_plot_figure_canvas, self.CalibrationFitTabPlot)
            self.cali_plot_toolbar.update()
            self.cali_plot_figure_canvas.draw()
            self.cali_plot_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=10, pady=10)

            print(self.cali_fit.params, self.cali_fit.NormRes, self.cali_fit.x, self.cali_fit.y)

            # Parameters table
            self.cali_generate_fit_params_table()

            # Removing error
            self.cali_function_lock_error.configure(state='normal')
            self.cali_function_lock_error.delete("1.0", tk.END)
            self.cali_function_lock_error.grid_forget()
        except Exception as e:
            print(traceback.format_exc())
            self.cali_function_lock_error.configure(state='normal')
            self.cali_function_lock_error.delete("1.0", tk.END)
            self.cali_function_lock_error.insert(tk.END, str(e))
            self.cali_function_lock_error.configure(state='disabled')
            self.cali_function_lock_error.grid(row=9)

    def cali_fit_settings(self):
        cali_fit_settings_popup = CaliSettings(self)
        self.wait_window(cali_fit_settings_popup)

    def cali_points_plus(self):
        grid_size = self.cali_points_table.grid_size()
        size = int((self.cali_fit_text.winfo_width()) / 4 / tk.font.Font(font='TkFixedFont').actual()['size'])
        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        # Remove the buttons to make room for new row
        plus = self.cali_points_table.grid_slaves(row=grid_size[1] - 1, column=grid_size[0] - 1)[0]
        minus = self.cali_points_table.grid_slaves(row=grid_size[1] - 1, column=grid_size[0] - 2)[0]

        plus.grid_forget()
        minus.grid_forget()

        self.cali_points_exp.append(tk.StringVar())
        self.cali_points_exp[-1].set('0')
        self.cali_points_exp_er.append(tk.StringVar())
        self.cali_points_exp_er[-1].set('1')
        self.cali_points_lit.append(tk.StringVar())
        self.cali_points_lit[-1].set('0')
        self.cali_points_lit_err.append(tk.StringVar())
        self.cali_points_lit_err[-1].set('1')

        ent1 = tk.Entry(self.cali_points_table, textvariable=self.cali_points_exp[-1], width=size, validate='key',
                        validatecommand=vcmdf)
        ent1.grid(row=grid_size[1] - 1, column=0)
        ent2 = tk.Entry(self.cali_points_table, textvariable=self.cali_points_exp_er[-1], width=size, validate='key',
                        validatecommand=vcmdf)
        ent2.grid(row=grid_size[1] - 1, column=1)
        ent3 = tk.Entry(self.cali_points_table, textvariable=self.cali_points_lit[-1], width=size, validate='key',
                        validatecommand=vcmdf)
        ent3.grid(row=grid_size[1] - 1, column=2)
        ent4 = tk.Entry(self.cali_points_table, textvariable=self.cali_points_lit_err[-1], width=size, validate='key',
                        validatecommand=vcmdf)
        ent4.grid(row=grid_size[1] - 1, column=3)

        plus.config(width=size)
        minus.config(width=size)
        plus.grid(row=grid_size[1], column=grid_size[0] - 1)
        minus.grid(row=grid_size[1], column=grid_size[0] - 2)

    def cali_points_minus(self):
        grid_size = self.cali_points_table.grid_size()

        if grid_size[1] > 2:
            # Remove the buttons to make room for new row
            plus = self.cali_points_table.grid_slaves(row=grid_size[1] - 1, column=grid_size[0] - 1)[0]
            minus = self.cali_points_table.grid_slaves(row=grid_size[1] - 1, column=grid_size[0] - 2)[0]

            plus.grid_forget()
            minus.grid_forget()

            for col in range(grid_size[0]):
                self.cali_points_table.grid_slaves(row=grid_size[1] - 2, column=col)[0].destroy()

            self.cali_points_exp.pop(grid_size[1] - 3)
            self.cali_points_exp_er.pop(grid_size[1] - 3)
            self.cali_points_lit.pop(grid_size[1] - 3)
            self.cali_points_lit_err.pop(grid_size[1] - 3)

            plus.grid(row=grid_size[1] - 2, column=grid_size[0] - 1)
            minus.grid(row=grid_size[1] - 2, column=grid_size[0] - 2)
        else:
            pass

    def export_cali_fit(self):
        file = open(self.root_file[:-5] + '_' + self.cali_preset_dropdown_var.get() + '_calibration.txt', 'w')
        file.write(self.cali_preset_dropdown_var.get() + '\n')
        for i in range(1, self.cali_fit_function_nr_params):
            file.write(f'{self.cali_param_fit_val[i - 1].get()},{self.cali_param_fit_error[i - 1].get()}\n')
        file.close()

    def open_lookup_tabel(self):
        energy_lookup_tabel = EnergyTable()
        self.wait_window(energy_lookup_tabel)

    def import_file_check_valid_cali(self, file):
        try:
            f = open(file, 'r')
            function = f.readline()
            f.close()
            if 'Linear' in function:
                imp_cali_func = lambda x, a, b: a * x + b
            elif 'Proportional' in function:
                imp_cali_func = lambda x, a: a * x
            elif 'Power' in function:
                imp_cali_func = lambda x, a, b: a * x ** b

            cali_param, _ = np.loadtxt(file, delimiter=",", skiprows=1, unpack=True)

            try:
                imp_cali_func(np.ones(10), *cali_param)
            except:
                raise FileExistsError('File not valid format')
        except:
            self.import_cali_fit()

    def import_cali_fit(self):
        if len(self.cali_import_file) != 0:
            file = filedialog.askopenfilename(parent=self, title='Choose file')
        else:
            file = filedialog.askopenfilename(parent=self, initialfile=self.cali_import_file, title='Choose file')
        if len(file) != 0:
            self.import_file_check_valid_cali(file)
            self.cali_import_file = file
            self.cali_import_label_var.set('File chosen')

            self.imp_cali_param, self.imp_cali_err = np.loadtxt(file, delimiter=",", skiprows=1, unpack=True)
            f = open(self.cali_import_file, 'r')
            function_set = f.readline()
            f.close
            if 'Power' in function_set:
                exec(
                    f'import numpy as np\ndef converter(x):\n\treturn np.nan_to_num({self.imp_cali_param[0]}*x**{self.imp_cali_param[1]})',
                    self.radius_to_energy_func)
                exec(
                    f'import numpy as np\ndef inverter(x):\n\treturn np.nan_to_num((x/{self.imp_cali_param[0]})**{1 / self.imp_cali_param[1]})',
                    self.radius_to_energy_func)
                exec(
                    f'import numpy as np\ndef error(x,err):\n\treturn np.nan_to_num(np.array([(({self.imp_cali_err[0]}*x**{self.imp_cali_param[1]})**2+np.nan_to_num({self.imp_cali_err[1]}*{self.imp_cali_err[0]}*np.log(x)*x**{self.imp_cali_param[1] - 1})**2+({self.imp_cali_param[0]}*x**{self.imp_cali_param[1]}-{self.imp_cali_param[0]}*(x-err)**{self.imp_cali_param[1]})**2)**(1/2),(({self.imp_cali_err[0]}*x**{self.imp_cali_param[1]})**2+np.nan_to_num({self.imp_cali_err[0]}*{self.imp_cali_err[1]}*np.log(x)*x**{self.imp_cali_param[1] - 1})**2+({self.imp_cali_param[0]}*x**{self.imp_cali_param[1]}-{self.imp_cali_param[0]}*(x+err)**{self.imp_cali_param[1]})**2)**(1/2)]))',
                    self.radius_to_energy_func)
            elif 'Proportional' in function_set:
                exec(f'def converter(x):\n\treturn {self.imp_cali_param[0]}*x', self.radius_to_energy_func)
                exec(f'def inverter(x):\n\treturn x/{self.imp_cali_param[0]}',
                     self.radius_to_energy_func)
                exec(f'def error(x,err):\n\treturn (({self.imp_cali_err[0]}*x)**2+err**2)**(1/2)',
                     self.radius_to_energy_func)
            elif 'Linear' in function_set:
                exec(f'def converter(x):\n\treturn {self.imp_cali_param[0]}*x+{self.imp_cali_param[1]}',
                     self.radius_to_energy_func)
                exec(f'def inverter(x):\n\treturn (x-{self.imp_cali_param[1]})/{self.imp_cali_param[0]}',
                     self.radius_to_energy_func)
                exec(
                    f'def error(x,err):\n\treturn np.nan_to_num(np.array([[(({self.imp_cali_err[0]}*x)**2+{self.imp_cali_err[1]}**2+({self.imp_cali_param[0]}*x-{self.imp_cali_param[0]}*(x-err))**2)**(1/2)],[(({self.imp_cali_err[0]}*x)**2+{self.imp_cali_err[1]}**2+({self.imp_cali_param[0]}*x-{self.imp_cali_param[0]}*(x+err))**2)**(1/2)]]))',
                    self.radius_to_energy_func)

    def use_calibration_import(self):
        if len(self.cali_import_file) != 0:
            if self.cali_use_import_button_var.get() == 'True':
                self.cali_use_import_button_var.set('False')
                self.plot_rad_bin()
                self.plot_intensity()
                if len(self.boot_intensity) != 0 and len(self.boot_intensity) == len(self.boot_radial):
                    self.plot_bootstrap()
            else:
                self.cali_use_import_button_var.set('True')
                self.plot_rad_bin()
                self.plot_intensity()
                if len(self.boot_intensity) != 0 and len(self.boot_intensity) == len(self.boot_radial):
                    self.plot_bootstrap()
        else:
            self.cali_use_import_button_var.set('False')

    def cali_converter_popup(self):
        if len(self.cali_import_file) != 0:
            converter_popup = EnergyConverter(self)
            self.wait_window(converter_popup)

    def set_laser_fit_preset(self, *args):
        if self.laser_preset_dropdown_var.get() == 'Custom':
            self.laser_fit_text.configure(state='normal')
            self.laser_fit_text.delete("1.0", "end")
            self.laser_fit_text.insert(tk.END, 'import numpy as np\n\ndef f(x,*params):\n\treturn')
        else:
            self.laser_fit_text.configure(state='normal')
            if self.laser_preset_dropdown_var.get() == 'Students-t':
                self.laser_fit_text.delete("1.0", "end")
                self.laser_fit_text.insert(tk.END, self.t_dist)
                self.laser_fit_text.configure(state='disabled')
            elif self.laser_preset_dropdown_var.get() == 'Gauss':
                self.laser_fit_text.delete("1.0", "end")
                self.laser_fit_text.insert(tk.END, self.gaus_dist)
                self.laser_fit_text.configure(state='disabled')
            elif self.laser_preset_dropdown_var.get() == 'Two students-t':
                self.laser_fit_text.delete("1.0", "end")
                self.laser_fit_text.insert(tk.END, self.t_dist_double)
                self.laser_fit_text.configure(state='disabled')

    def laser_function_lock(self):
        if self.laser_function_lock_var.get() == 'Function Locked':
            self.laser_function_lock_var.set('Function Unlocked')
            if self.laser_preset_dropdown_var.get() == 'Custom':
                self.laser_fit_text.configure(state='normal')
            self.laser_preset_dropdown.configure(state='normal')
            for widgets in self.laser_parameters_table.winfo_children():
                widgets.destroy()
            for widgets in self.laser_parameters_fit_table.winfo_children():
                if not isinstance(widgets, tk.Button):
                    widgets.destroy()
                else:
                    widgets.grid_forget()
            self.laser_parameters_fit_table.grid_forget()
            self.laser_fit_button.grid_forget()
            self.laser_fit_settings_button.grid_forget()
            self.laser_fit_function = {}
        else:
            try:
                # Checking function
                self.laser_fit_function = {}
                exec(self.laser_fit_text.get("1.0", tk.END), self.laser_fit_function)
                self.laser_param_names = getfullargspec(self.laser_fit_function['f'])[0]
                self.laser_fit_function_nr_params = len(self.laser_param_names)
                ones = np.ones(self.laser_fit_function_nr_params)
                self.laser_fit_function['f'](*ones)

                # Checking jacobian function
                if 'jac' in self.laser_fit_function:
                    self.laser_fit_function['jac'](*ones)
                    self.laser_fit_function['jac'](ones, *ones[1:])
                else:
                    self.laser_fit_function['jac'] = None

                # Checking hessian function
                if 'hess' in self.laser_fit_function:
                    pass  # self.peak_fit_function['hess'](*ones)
                    # self.peak_fit_function['hess'](ones, *ones[1:])
                else:
                    self.laser_fit_function['hess'] = None

                # Removing error
                self.laser_function_lock_error.configure(state='normal')
                self.laser_function_lock_error.delete("1.0", tk.END)
                self.laser_function_lock_error.grid_forget()

                # Parameters table
                self.laser_generate_params_table()

                # Fitting tools
                self.laser_energy_low_var.set(str(np.min(self.freq_bins)))
                self.laser_energy_high_var.set(str(np.max(self.freq_bins)))
                self.laser_fit_button.grid(row=6, column=0, sticky=tk.W)
                self.laser_fit_settings_button.grid(row=6, column=0, sticky=tk.E)

                # Freezing function input
                self.laser_fit_text.configure(state='disabled')
                self.laser_function_lock_var.set('Function Locked')
                self.laser_preset_dropdown.configure(state='disabled')
            except Exception as e:
                print(traceback.format_exc())
                self.laser_function_lock_error.configure(state='normal')
                self.laser_function_lock_error.delete("1.0", tk.END)
                self.laser_function_lock_error.insert(tk.END, str(e))
                self.laser_function_lock_error.configure(state='disabled')
                self.laser_function_lock_error.grid(row=4)
                self.laser_fit_button.grid_forget()
                self.laser_fit_settings_button.grid_forget()

    def run_laser_fit(self):
        try:
            mask = np.array([float(self.laser_energy_low_var.get()) <= np.array(self.freq_bins)])[0] * \
                   np.array([np.array(self.freq_bins) <= float(self.laser_energy_high_var.get())])[0] * \
                   np.array([np.array(self.freq_hist_err) != 0])[0]

            if self.laser_method.get() == 'diff_evol':
                # Call fit
                p0 = np.array(
                    [[float(i.get()), float(j.get())] for i, j in zip(self.laser_param_lower, self.laser_param_upper)])
                self.laser_fit = FittingRoutine(self.laser_fit_function['f'], np.array(self.freq_bins),
                                                np.array(self.freq_hist), error_y=self.freq_hist_err, P0=p0,
                                                method=self.laser_method.get(),
                                                jac=self.laser_fit_function['jac'],
                                                hess=self.laser_fit_function['hess'], mask=mask)
            else:
                # Call fit
                p0 = np.array([float(i.get()) for i in self.param_init])
                self.laser_fit = FittingRoutine(self.laser_fit_function['f'], np.array(self.freq_bins),
                                                np.array(self.freq_hist), error_y=self.freq_hist_err, P0=p0,
                                                method=self.laser_method.get(),
                                                jac=self.laser_fit_function['jac'],
                                                hess=self.laser_fit_function['hess'], mask=mask)

            # Plot
            self.laser_hist_figure, self.laser_hist_ax = self.laser_fit.Plot(xlabel='Energy (eV)', ylabel='Counts',
                                                                             init=True, do_return=True, figsize=(6, 6),
                                                                             ms=0.5,
                                                                             capsize=1, elinewidth=0.4, lw=0.7)
            self.laser_hist_figure.set_dpi(200)
            self.laser_figure_canvas.get_tk_widget().destroy()
            self.laser_toolbar.destroy()
            self.laser_figure_canvas = FigureCanvasTkAgg(self.laser_hist_figure, master=self.LaserFitTabPlot)
            self.laser_toolbar = NavigationToolbar2Tk(self.laser_figure_canvas, self.LaserFitTabPlot)
            self.laser_toolbar.update()
            self.laser_figure_canvas.draw()
            self.laser_figure_canvas.get_tk_widget().pack(fill='both', expand=True, padx=10, pady=10)

            # Parameters table
            self.laser_generate_fit_params_table()

            # Removing error
            self.laser_function_lock_error.configure(state='normal')
            self.laser_function_lock_error.delete("1.0", tk.END)
            self.laser_function_lock_error.grid_forget()
        except Exception as e:
            print(traceback.format_exc())
            self.laser_function_lock_error.configure(state='normal')
            self.laser_function_lock_error.delete("1.0", tk.END)
            self.laser_function_lock_error.insert(tk.END, str(e))
            self.laser_function_lock_error.configure(state='disabled')
            self.laser_function_lock_error.grid(row=7)

    def laser_fit_settings(self):
        laser_fit_settings_popup = LaserSettings(self)
        self.wait_window(laser_fit_settings_popup)

    def laser_generate_params_table(self):
        self.laser_param_init.clear()
        self.laser_param_lower.clear()
        self.laser_param_upper.clear()
        for widgets in self.laser_parameters_table.winfo_children():
            widgets.destroy()

        size = int(
            self.laser_function_lock_button.winfo_width() / 4 / tk.font.Font(font='TkFixedFont').actual()['size'])

        for i in range(len(self.laser_parameters_table_head)):
            label = ttk.Label(self.laser_parameters_table, text=self.laser_parameters_table_head[i], width=size)
            label.grid(row=0, column=i)

        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        for i in range(1, self.laser_fit_function_nr_params):
            label = ttk.Label(self.laser_parameters_table, text=self.laser_param_names[i], width=size)
            label.grid(row=i, column=0, sticky=tk.W)

        for i in range(1, self.laser_fit_function_nr_params):
            self.laser_param_init.append(tk.StringVar())
            self.laser_param_init[i - 1].set('1')
            text = tk.Entry(self.laser_parameters_table, validate='key', validatecommand=vcmdf,
                            textvariable=self.laser_param_init[i - 1], width=size)
            text.grid(row=i, column=1)

            self.laser_param_lower.append(tk.StringVar())
            self.laser_param_lower[i - 1].set('0')
            text1 = tk.Entry(self.laser_parameters_table, validate='key', validatecommand=vcmdf,
                             textvariable=self.laser_param_lower[i - 1], width=size)
            text1.grid(row=i, column=2)

            self.laser_param_upper.append(tk.StringVar())
            self.laser_param_upper[i - 1].set('1')
            text2 = tk.Entry(self.laser_parameters_table, validate='key', validatecommand=vcmdf,
                             textvariable=self.laser_param_upper[i - 1], width=size)
            text2.grid(row=i, column=3)

    def laser_generate_fit_params_table(self):
        self.laser_param_fit_val.clear()
        self.laser_param_fit_error.clear()
        for widgets in self.laser_parameters_fit_table.winfo_children():
            if not isinstance(widgets, tk.Button):
                widgets.destroy()
            else:
                widgets.grid_forget()

        self.laser_parameters_fit_table.grid(row=8, sticky=tk.EW)

        size = int((self.laser_fit_text.winfo_width()) / 3 / tk.font.Font(font='TkFixedFont').actual()['size'])

        for i in range(len(self.laser_parameters_fit_table_head)):
            label = ttk.Label(self.laser_parameters_fit_table, text=self.laser_parameters_fit_table_head[i], width=size)
            label.grid(row=0, column=i)

        for i in range(1, self.laser_fit_function_nr_params):
            label = ttk.Label(self.laser_parameters_fit_table, text=self.laser_param_names[i], width=size)
            label.grid(row=i, column=0, sticky=tk.W)

        for i in range(1, self.laser_fit_function_nr_params):
            self.laser_param_fit_val.append(tk.StringVar())
            self.laser_param_fit_val[i - 1].set('%1.3e' % self.laser_fit.params[i - 1])
            text = tk.Label(self.laser_parameters_fit_table, textvariable=self.laser_param_fit_val[i - 1], width=size)
            text.grid(row=i, column=1)

            self.laser_param_fit_error.append(tk.StringVar())
            self.laser_param_fit_error[i - 1].set('%1.3e' % self.laser_fit.Error[i - 1])
            text1 = tk.Label(self.laser_parameters_fit_table, textvariable=self.laser_param_fit_error[i - 1],
                             width=size)
            text1.grid(row=i, column=2)

        text = ttk.Label(self.laser_parameters_fit_table, text=('P-Value : %1.2f' % self.laser_fit.Pval),
                         width=3 * size)
        text.grid(row=self.laser_fit_function_nr_params, column=0, columnspan=3, sticky=tk.W, pady=(10, 5))

        text = ttk.Label(self.laser_parameters_fit_table,
                         text=('Chi2, df : %1.2f, %i' % (self.laser_fit.Chi2, self.laser_fit.df)), width=3 * size)
        text.grid(row=self.laser_fit_function_nr_params + 1, column=0, columnspan=3, sticky=tk.W, pady=(0, 10))

    def general_settings(self):
        settings_popup = GeneralSettings(self)
        self.wait_window(settings_popup)


class FilePopup(tk.Toplevel):
    def __init__(self, root):
        tk.Toplevel.__init__(self, root)

        data_txt_label = tk.Label(self, text='Data file: ')
        data_txt_label.grid(row=0, column=0, sticky=tk.E)
        data_label = tk.Label(self, text=root.dat_file_path)
        data_label.grid(row=0, column=1, sticky=tk.W)
        data_button = tk.Button(self, text="Choose file", command=lambda: search_for_dat_file(root, data_label))
        data_button.grid(row=0, column=2, sticky=tk.E)

        load_txt_label = tk.Label(self, text='Load data function directory: ')
        load_txt_label.grid(row=1, column=0, sticky=tk.E)
        if root.home_dir in root.load_file_path:
            load_label = tk.Label(self, text=root.load_file_path[len(root.home_dir):])
        else:
            load_label = tk.Label(self, text=root.load_file_path)
        load_label.grid(row=1, column=1, sticky=tk.W)
        load_button = tk.Button(self, text="Choose directory", command=lambda: search_for_load_file(root, load_label))
        load_button.grid(row=1, column=2, sticky=tk.E)

        event_txt_label = tk.Label(self, text='Event/Trig data function directory: ')
        event_txt_label.grid(row=2, column=0, sticky=tk.E)
        if root.home_dir in root.load_file_path:
            event_label = tk.Label(self, text=root.event_file_path[len(root.home_dir):])
        else:
            event_label = tk.Label(self, text=root.event_file_path)
        event_label.grid(row=2, column=1, sticky=tk.W)
        event_button = tk.Button(self, text="Choose directory",
                                 command=lambda: search_for_event_file(root, event_label))
        event_button.grid(row=2, column=2, sticky=tk.E)

        vmi_txt_label = tk.Label(self, text='VMI filter function directory: ')
        vmi_txt_label.grid(row=3, column=0, sticky=tk.E)
        if root.home_dir in root.load_file_path:
            vmi_label = tk.Label(self, text=root.vmi_file_path[len(root.home_dir):])
        else:
            vmi_label = tk.Label(self, text=root.vmi_file_path)
        vmi_label.grid(row=3, column=1, sticky=tk.W)
        vmi_button = tk.Button(self, text="Choose directory", command=lambda: search_for_vmi_file(root, vmi_label))
        vmi_button.grid(row=3, column=2, sticky=tk.E)

        gif_txt_label = tk.Label(self, text='GIF function directory: ')
        gif_txt_label.grid(row=4, column=0, sticky=tk.E)
        if root.home_dir in root.load_file_path:
            gif_label = tk.Label(self, text=root.gif_file_path[len(root.home_dir):])
        else:
            gif_label = tk.Label(self, text=root.gif_file_path)
        gif_label.grid(row=4, column=1, sticky=tk.W)
        gif_button = tk.Button(self, text="Choose directory", command=lambda: search_for_gif_file(root, gif_label))
        gif_button.grid(row=4, column=2, sticky=tk.E)

        boot_txt_label = tk.Label(self, text='Bootstrap function directory: ')
        boot_txt_label.grid(row=5, column=0, sticky=tk.E)
        if root.home_dir in root.load_file_path:
            boot_label = tk.Label(self, text=root.boot_file_path[len(root.home_dir):])
        else:
            boot_label = tk.Label(self, text=root.boot_file_path)
        boot_label.grid(row=5, column=1, sticky=tk.W)
        boot_button = tk.Button(self, text="Choose directory", command=lambda: search_for_boot_file(root, boot_label))
        boot_button.grid(row=5, column=2, sticky=tk.E)

        error_label = tk.Label(self, text='Error')
        error = tk.Text(self)
        button = tk.Button(self, text="Click to continue",
                           command=lambda: confirm_close(root, self, error, error_label))
        button.grid(row=6, column=1)


class BootCalcPopup(tk.Toplevel):
    def __init__(self, root, boot_file_path, root_file, args, circ):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "Go get some coffee")

        self.geometry('300x100')
        self.args = args
        self.root = root
        self.circ = circ
        self.boot_file_path = boot_file_path
        self.root_file = root_file

        boot_label = tk.Label(self, text="Calculating Bootstrap")
        self.boot_str = tk.StringVar()
        self.boot_out = tk.Label(self, textvariable=self.boot_str)
        boot_label.pack(fill='x', expand=True)
        self.boot_out.pack(fill='x', expand=True)

        t = Thread(target=lambda: self.sudo())
        t.start()

    def sudo(self):

        call = self.boot_file_path.replace(" ", "\\ ") + "/build/./Main" + " " + self.root_file.replace(" ",
                                                                                                        "\\ ") + " "

        for arg in self.args:
            call += arg

        print(call)

        p = subprocess.Popen(call, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1,
                             universal_newlines=True, shell=True)
        p.poll()

        while True:
            line = p.stdout.readline().strip("\n")
            p_line = self.boot_str.get().split("\n")
            if len(p_line) >= 2:
                p_line = p_line[1]
                self.boot_str.set(p_line + "\n" + line)
            else:
                self.boot_str.set(line + "\n" + line)
            if not line and p.poll is not None: break

        while True:
            err = p.stderr.readline().strip("\n")
            p_line = self.boot_str.get().split("\n")
            if len(p_line) >= 2:
                p_line = p_line[1]
                self.boot_str.set(p_line + "\n" + err)
            else:
                self.boot_str.set(err + "\n" + err)
            if not err and p.poll is not None: break

        self.root.initialise_boot()

        self.destroy()
        self.update()


class BootPopup(tk.Toplevel):
    def __init__(self, root, boot_file_path, root_file):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "Bootstrap settings")

        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        vcmdi = (self.register(validate_int),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        # Initial parameters needed
        self.geometry('650x375')
        self.root = root
        self.boot_file_path = boot_file_path
        self.root_file = root_file

        ''' labels and widgets '''
        radial_resolution_label = tk.Label(self, text='Radial resolution')
        self.radial_resolution_str = tk.StringVar()
        self.radial_resolution_str.set(root.resolution_entry.get())
        self.radial_resolution_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                                width=5, textvariable=self.radial_resolution_str)

        radial_max_label = tk.Label(self, text='Maximal radius (mm)')
        self.radial_max_str = tk.StringVar()
        self.radial_max_str.set(root.resolution_max_entry.get())
        self.radial_max_entry = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                         width=5, textvariable=self.radial_max_str)

        angular_order_label = tk.Label(self, text='Angular orders')
        self.angular_order_str = tk.StringVar()
        self.angular_order_str.set(root.norder_entry.get())
        self.angular_order_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                            width=5, textvariable=self.angular_order_str)

        vmi_resolution_label = tk.Label(self, text='VMI resolution')
        self.vmi_resolution_str = tk.StringVar()
        self.vmi_resolution_str.set(root.vmi_resolution_entry.get())
        self.vmi_resolution_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                             width=5, textvariable=self.vmi_resolution_str)

        circularize_order_label = tk.Label(self, text='Circularization order')
        self.circularize_order_str = tk.StringVar()
        self.circularize_order_str.set(root.circular_entry.get())
        self.circularize_order_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                                width=5, textvariable=self.circularize_order_str)
        self.circularize_button_str = tk.StringVar()
        if root.circularize:
            self.circularize_button_str.set('On')
        else:
            self.circularize_button_str.set('Off')
        self.circularize_button = tk.Button(self, textvariable=self.circularize_button_str, width=4,
                                            command=lambda: self.circularize_on_off())

        trigger_time_label = tk.Label(self, text='Trigger time interval (mu s)')
        trigger_time_label_lower = tk.Label(self, text='Lower')
        self.trigger_time_str_lower = tk.StringVar()
        self.trigger_time_str_lower.set(root.trig_lower_entry.get())
        self.trigger_time_entry_lower = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                                 width=7, textvariable=self.trigger_time_str_lower)
        trigger_time_label_upper = tk.Label(self, text='Upper')
        self.trigger_time_str_upper = tk.StringVar()
        self.trigger_time_str_upper.set(root.trig_upper_entry.get())
        self.trigger_time_entry_upper = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                                 width=7, textvariable=self.trigger_time_str_upper)

        event_time_label = tk.Label(self, text='Event time interval (s)')
        event_time_label_lower = tk.Label(self, text='Lower')
        self.event_time_str_lower = tk.StringVar()
        self.event_time_str_lower.set(root.event_lower_entry.get())
        self.event_time_entry_lower = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                               width=7, textvariable=self.event_time_str_lower)
        event_time_label_upper = tk.Label(self, text='Upper')
        self.event_time_str_upper = tk.StringVar()
        self.event_time_str_upper.set(root.event_upper_entry.get())
        self.event_time_entry_upper = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                               width=7, textvariable=self.event_time_str_upper)

        peak_sum_label = tk.Label(self, text='Peak sum interval')
        peak_sum_label_lower = tk.Label(self, text='Lower')
        self.peak_sum_str_lower = tk.StringVar()
        self.peak_sum_str_lower.set(root.peak_sum_lower_entry.get())
        self.peak_sum_entry_lower = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                             width=7, textvariable=self.peak_sum_str_lower)
        peak_sum_label_upper = tk.Label(self, text='Upper')
        self.peak_sum_str_upper = tk.StringVar()
        self.peak_sum_str_upper.set(root.peak_sum_upper_entry.get())
        self.peak_sum_entry_upper = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                             width=7, textvariable=self.peak_sum_str_upper)

        samples_label = tk.Label(self, text='Number of bootstrap samples')
        self.samples_str = tk.StringVar()
        self.samples_str.set('100')
        self.samples_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                      width=8, textvariable=self.samples_str)

        threads_label = tk.Label(self, text='Number of threads')
        self.threads_str = tk.StringVar()
        self.threads_str.set('1')
        self.threads_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                      width=8, textvariable=self.threads_str)

        self.boot_button = tk.Button(self, text='Confirm and calculate', command=lambda: self.calc_boot())

        ''' Set the widgets'''
        radial_resolution_label.grid(row=0, column=0, sticky=tk.W, pady=15, padx=10)
        self.radial_resolution_entry.grid(row=0, column=1, pady=15)

        radial_max_label.grid(row=1, column=0, sticky=tk.W, pady=15, padx=10)
        self.radial_max_entry.grid(row=1, column=1, pady=15)

        angular_order_label.grid(row=2, column=0, sticky=tk.W, pady=15, padx=10)
        self.angular_order_entry.grid(row=2, column=1, pady=15)

        vmi_resolution_label.grid(row=3, column=0, sticky=tk.W, pady=15, padx=10)
        self.vmi_resolution_entry.grid(row=3, column=1, pady=15)

        circularize_order_label.grid(row=3, column=0, sticky=tk.W, pady=15, padx=10)
        self.circularize_order_entry.grid(row=3, column=1, pady=15)
        self.circularize_button.grid(row=3, column=2, pady=15)

        space = tk.Label(self, text="                           ")
        space.grid(row=0, column=3, sticky=tk.W + tk.E)

        trigger_time_label.grid(row=0, column=4, columnspan=2, pady=10, padx=10)
        trigger_time_label_lower.grid(row=1, column=4, padx=10, sticky=tk.E + tk.N)
        trigger_time_label_upper.grid(row=1, column=4, padx=10, sticky=tk.E + tk.S)
        self.trigger_time_entry_lower.grid(row=1, column=5, padx=10, sticky=tk.E + tk.N)
        self.trigger_time_entry_upper.grid(row=1, column=5, padx=10, sticky=tk.E + tk.S)

        event_time_label.grid(row=2, column=4, columnspan=2, pady=10, padx=10)
        event_time_label_lower.grid(row=3, column=4, padx=10, sticky=tk.E + tk.N)
        event_time_label_upper.grid(row=3, column=4, padx=10, sticky=tk.E + tk.S)
        self.event_time_entry_lower.grid(row=3, column=5, padx=10, sticky=tk.E + tk.N)
        self.event_time_entry_upper.grid(row=3, column=5, padx=10, sticky=tk.E + tk.S)

        peak_sum_label.grid(row=4, column=4, columnspan=2, pady=10, padx=10)
        peak_sum_label_lower.grid(row=5, column=4, padx=10, sticky=tk.E + tk.N)
        peak_sum_label_upper.grid(row=5, column=4, padx=10, sticky=tk.E + tk.S)
        self.peak_sum_entry_lower.grid(row=5, column=5, padx=10, sticky=tk.E + tk.N)
        self.peak_sum_entry_upper.grid(row=5, column=5, padx=10, sticky=tk.E + tk.S)

        samples_label.grid(row=4, column=0, columnspan=3, pady=15, sticky=tk.S + tk.W)
        self.samples_entry.grid(row=4, column=2, columnspan=2)

        threads_label.grid(row=5, column=0, columnspan=3, pady=15, sticky=tk.S + tk.W)
        self.threads_entry.grid(row=5, column=2, columnspan=2)

        self.boot_button.grid(row=6, column=2, columnspan=3, pady=15, sticky=tk.S + tk.W + tk.E)

    def circularize_on_off(self):
        if self.circularize_button_str.get() == 'On':
            self.circularize_button_str.set('Off')
        else:
            self.circularize_button_str.set('On')

    def calc_boot(self):
        # Call to pass
        call = '-r' + self.radial_resolution_entry.get() + " "
        call += '-v' + self.vmi_resolution_entry.get() + ' '
        call += '-o' + self.angular_order_entry.get() + ' '
        call += '-t' + str(float(self.trigger_time_entry_lower.get()) * 1.e-6) + ' '
        call += '-T' + str(float(self.trigger_time_entry_upper.get()) * 1.e-6) + ' '
        call += '-e' + self.event_time_entry_lower.get() + ' '
        call += '-E' + self.event_time_entry_upper.get() + ' '
        call += '-m' + self.radial_max_entry.get() + ' '
        call += '-b' + self.samples_entry.get() + ' '
        call += '-n' + self.threads_entry.get() + ' '
        call += '-p' + self.peak_sum_entry_lower.get() + ' '
        call += '-P' + self.peak_sum_entry_upper.get() + ' '

        if self.root.background_sub_mode[0] == 1:
            call += '-L '
        else:
            call += '-l '
        if self.root.background_sub_mode[1] == 1:
            call += '-S' + str(self.root.background_sub_mode[3]) + ' '
            call += '-s' + str(self.root.background_sub_mode[2]) + ' '
        else:
            call += '-s-1 '

        if len(self.root.circ_import_file) != 0:
            if self.root.circ_use_import_button_var.get() == 'True':
                call += '-C' + self.root.circ_import_file + ' '
        if self.circularize_button_str.get() == 'On':
            call += '-c' + self.circularize_order_entry.get() + ' '

        # Do the calculation
        calculation_popup = BootCalcPopup(self.root, self.boot_file_path, self.root_file, call,
                                          (self.circularize_button_str.get() == 'On'))
        self.wait_window(calculation_popup)


class GifCalcPopup(tk.Toplevel):
    def __init__(self, root, gif_file_path, root_file, args, circ, boot):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "Go get some coffee")

        self.geometry('300x100')
        self.args = args
        self.root = root
        self.circ = circ
        self.boot = boot
        self.gif_file_path = gif_file_path
        self.root_file = root_file

        boot_label = tk.Label(self, text="Calculating GIF")
        self.boot_str = tk.StringVar()
        self.boot_out = tk.Label(self, textvariable=self.boot_str)
        boot_label.pack(fill='x', expand=True)
        self.boot_out.pack(fill='x', expand=True)

        t = Thread(target=lambda: self.sudo())
        t.start()

    def sudo(self):

        call = self.gif_file_path.replace(" ", "\\ ") + "/build/./Main" + " " + self.root_file.replace(" ", "\\ ") + " "

        for arg in self.args:
            call += arg

        p = subprocess.Popen(call, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1,
                             universal_newlines=True, shell=True)
        p.poll()

        while True:
            line = p.stdout.readline().strip("\n")
            p_line = self.boot_str.get().split("\n")
            if len(p_line) >= 2:
                p_line = p_line[1]
                self.boot_str.set(p_line + "\n" + line)
            else:
                self.boot_str.set(line + "\n" + line)
            if not line and p.poll is not None: break

        while True:
            err = p.stderr.readline().strip("\n")
            p_line = self.boot_str.get().split("\n")
            if len(p_line) >= 2:
                p_line = p_line[1]
                self.boot_str.set(p_line + "\n" + err)
            else:
                self.boot_str.set(err + "\n" + err)
            if not err and p.poll is not None: break
        self.boot_str.set('Calculation done\nCreating GIF files')

        # Empty figures
        intensity_figure = plt.Figure((6, 3), dpi=200)
        intensity_ax = intensity_figure.subplots(1)
        radius_figure = plt.Figure((6, 3), dpi=200)
        radius_ax = radius_figure.subplots(1)
        event_figure = plt.Figure((6, 3), dpi=200)
        event_ax = event_figure.subplots(1)
        raw_figure = plt.Figure((6, 4), dpi=200)
        raw_ax = raw_figure.subplots(1)
        raw_off_figure = plt.Figure((6, 4), dpi=200)
        raw_off_ax = raw_off_figure.subplots(1)
        inv_vmi_figure = plt.Figure((6, 4), dpi=200)
        inv_vmi_ax = inv_vmi_figure.subplots(1)

        if self.circ:
            raw_circ_figure = plt.Figure((6, 4), dpi=200)
            raw_circ_ax = raw_circ_figure.subplots(1)
            raw_circ_vmi_files = []

        # create the gifs!
        gif_file = ROOT.TFile(self.root_file[:-5] + "_gif.root", 'READ')
        raw_vmi_files = []
        raw_off_vmi_files = []
        inv_vmi_files = []
        intensity_files = []
        radius_files = []
        event_files = []
        times = []

        # Values
        event_min = np.min(np.array(self.root.event_hist))
        event_max = np.max(np.array(self.root.event_hist))
        time_min = self.root.event_time_bins[0]
        time_max = self.root.event_time_bins[-1]
        frame_tree = gif_file.Get("SETTINGS")
        settings_branch = frame_tree.GetBranch("settings")

        settings_branch.GetEntry(0)
        intensity_min = getattr(frame_tree, "settings/settings")
        settings_branch.GetEntry(1)
        intensity_max = getattr(frame_tree, "settings/settings")
        settings_branch.GetEntry(2)
        radius_max = getattr(frame_tree, "settings/settings")
        settings_branch.GetEntry(4)
        radial_tot_min = getattr(frame_tree, "settings/settings")
        settings_branch.GetEntry(5)
        radial_tot_max = getattr(frame_tree, "settings/settings")

        file_nr = 0

        for key in gif_file.GetListOfKeys():
            if key.GetName() != "SETTINGS":
                self.boot_str.set(f'Creating frame {file_nr + 1}')
                # Keys and branches
                frame_tree = gif_file.Get(key.GetName())
                circ_image_branch = frame_tree.GetBranch("circ_image")
                inv_image_branch = frame_tree.GetBranch("inv_image")
                raw_image_branch = frame_tree.GetBranch("raw_image")
                raw_off_image_branch = frame_tree.GetBranch("raw_off_image")
                intensity_branch = frame_tree.GetBranch("intensity")
                intensity_err_branch = frame_tree.GetBranch("intensity_err")
                rad_tot_branch = frame_tree.GetBranch("radial_bin_TOTAL")
                rad_val_branch = frame_tree.GetBranch("radial_bin_values")
                rad_err_branch = frame_tree.GetBranch("radial_bin_error")
                time_branch = frame_tree.GetBranch("time")
                settings_branch = frame_tree.GetBranch("settings")
                if self.boot:
                    I_err_l_branch = frame_tree.GetBranch("intensity_err_lower")
                    I_err_u_branch = frame_tree.GetBranch("intensity_err_upper")

                time_branch.GetEntry(0)
                time = [getattr(frame_tree, "time/time")]
                time_branch.GetEntry(1)
                time.append(getattr(frame_tree, "time/time"))

                times.append([time[0], time[1]])

                settings_branch.GetEntry(0)
                resolution = getattr(frame_tree, "settings/settings")
                settings_branch.GetEntry(1)
                max_radial = getattr(frame_tree, "settings/settings")
                settings_branch.GetEntry(2)
                vmi_resolution = getattr(frame_tree, "settings/settings")

                # Arrays to append to
                circ_image = np.zeros((int(vmi_resolution), int(vmi_resolution)))
                inv_image = np.zeros((int(vmi_resolution), int(vmi_resolution)))
                image_val = np.zeros((int(vmi_resolution), int(vmi_resolution)))
                image_off_val = np.zeros((int(vmi_resolution), int(vmi_resolution)))
                intensity = np.zeros(int(resolution))
                intensity_err = np.zeros(int(resolution))
                radial_tot = np.zeros(int(resolution))
                radial_err = np.zeros(int(resolution))
                radial = np.zeros(int(resolution))
                rad_val_err = np.ones(int(resolution)) * self.root.sys_radial_err
                if self.boot:
                    I_err_l = np.zeros(int(resolution))
                    I_err_u = np.zeros(int(resolution))

                # Reading the vmi images
                count = 0
                for i in range(inv_image_branch.GetEntries()):
                    if self.circ:
                        circ_image_branch.GetEntry(i)
                    inv_image_branch.GetEntry(i)
                    raw_image_branch.GetEntry(i)
                    raw_off_image_branch.GetEntry(i)

                    if self.circ:
                        circ_image[count // int(vmi_resolution)][count % int(vmi_resolution)] += \
                            getattr(frame_tree, 'circ_image/circ_image')
                    inv_image[count // int(vmi_resolution)][count % int(vmi_resolution)] += \
                        getattr(frame_tree, 'inv_image/inv_image')
                    image_val[count // int(vmi_resolution)][count % int(vmi_resolution)] += \
                        getattr(frame_tree, 'raw_image/raw_image')
                    image_off_val[count // int(vmi_resolution)][count % int(vmi_resolution)] += \
                        getattr(frame_tree, 'raw_off_image/raw_off_image')
                    count += 1

                combined_vmi = np.zeros((int(vmi_resolution), int(vmi_resolution)))
                if self.circ:
                    # Combine circ and inverse
                    center_col = int(int(vmi_resolution) / 2)
                    if np.max(inv_image[:, :center_col]) != 0:
                        combined_vmi[:, center_col:] = circ_image[:, center_col:]
                        combined_vmi[:, :center_col] = inv_image[:, :center_col] / np.max(inv_image[:, :center_col]) * \
                                                       np.max(circ_image[:, center_col:])
                    else:
                        combined_vmi[:, center_col:] = circ_image[:, center_col:]
                        combined_vmi[:, :center_col] = inv_image[:, :center_col]

                # Reading the intensity
                for i in range(intensity_branch.GetEntries()):
                    intensity_branch.GetEntry(i)
                    intensity_err_branch.GetEntry(i)
                    rad_tot_branch.GetEntry(i)
                    rad_val_branch.GetEntry(i)
                    rad_err_branch.GetEntry(i)
                    if self.boot:
                        I_err_l_branch.GetEntry(i)
                        I_err_u_branch.GetEntry(i)
                        I_err_l[i] = getattr(frame_tree, 'intensity_err_lower/intensity_err_lower')
                        I_err_u[i] = getattr(frame_tree, 'intensity_err_upper/intensity_err_upper')

                    intensity[i] = getattr(frame_tree, 'intensity/intensity')
                    intensity_err[i] = getattr(frame_tree, 'intensity_err/intensity_err')
                    radial_tot[i] = getattr(frame_tree, 'radial_bin_TOTAL/radial_bin_TOTAL')
                    radial_err[i] = getattr(frame_tree, 'radial_bin_error/radial_bin_error')
                    radial[i] = getattr(frame_tree, 'radial_bin_values/radial_bin_values')
                rad_val_err_stat = np.ones(int(resolution)) * (radial[-2] - radial[-3])
                if self.boot:
                    quantiles = np.array([intensity - I_err_l,
                                          I_err_u - intensity])

                # Make the plots
                event_ax.clear()
                event_ax.axhline(0, color='k', lw=0.5)
                event_ax.axvline(time[0], ls='--', color="b", lw=0.5)
                event_ax.axvline(time[1], ls='--', color="b", lw=0.5)
                event_ax.errorbar(self.root.event_time_bins, self.root.event_hist, self.root.event_hist_err,
                                  fmt='ro-', capsize=2, ms=2, lw=1,
                                  errorevery=3, elinewidth=0.5)
                event_ax.set_title('Histogram of event times')
                event_ax.set_xlabel('Time (s)')
                event_ax.set_ylabel('Counts')
                event_ax.set_xlim(time_min, time_max)
                event_ax.set_ylim(event_min, event_max)

                intensity_ax.clear()
                intensity_ax.axhline(0, color='k', lw=0.5)
                if self.boot:
                    if self.root.cali_use_import_button_var.get() == 'True':
                        intensity_ax.errorbar(
                            self.root.radius_to_energy_func['converter'](radial),
                            intensity, intensity_err,
                            self.root.radius_to_energy_func['error'](radial,
                                                                     rad_val_err),
                            fmt='ro', ms=0, lw=0.5, capsize=1.2,
                            errorevery=1, elinewidth=0.25)
                        intensity_ax.errorbar(
                            self.root.radius_to_energy_func['converter'](radial),
                            intensity, quantiles,
                            self.root.radius_to_energy_func['error'](radial,
                                                                     rad_val_err_stat),
                            fmt='bo', ms=1, lw=0.5, capsize=1.2,
                            errorevery=1, elinewidth=0.25)
                        intensity_ax.set_xlim(0, self.root.radius_to_energy_func['converter'](radius_max))
                        intensity_ax.set_xlabel('Energy (eV)')
                    else:
                        intensity_ax.errorbar(radial,
                                              intensity, intensity_err, rad_val_err,
                                              fmt='ro', ms=0, lw=0.5, capsize=1.2,
                                              errorevery=1, elinewidth=0.25)
                        intensity_ax.errorbar(radial,
                                              intensity, quantiles, rad_val_err_stat,
                                              fmt='bo', ms=1, lw=0.5, capsize=1.2,
                                              errorevery=1, elinewidth=0.25)
                        intensity_ax.set_xlim(0, radius_max)
                        intensity_ax.set_xlabel('3D Radius (mm)')
                else:
                    if self.root.cali_use_import_button_var.get() == 'True':
                        intensity_ax.errorbar(
                            self.root.radius_to_energy_func['converter'](radial),
                            intensity, intensity_err,
                            self.root.radius_to_energy_func['error'](radial,
                                                                     rad_val_err),
                            fmt='ro', ms=1, lw=0.5, capsize=1.2,
                            errorevery=1, elinewidth=0.25)
                        intensity_ax.set_xlim(0, self.root.radius_to_energy_func['converter'](radius_max))
                        intensity_ax.set_xlabel('Energy (eV)')
                    else:
                        intensity_ax.errorbar(radial,
                                              intensity, intensity_err, rad_val_err,
                                              fmt='ro', ms=1, lw=0.5, capsize=1.2,
                                              errorevery=1, elinewidth=0.25)
                        intensity_ax.set_xlim(0, radius_max)
                        intensity_ax.set_xlabel('3D Radius (mm)')
                intensity_ax.set_title('Abel Inversed Intensity')
                intensity_ax.set_ylabel(r'Integration over $\theta^\prime$')
                intensity_ax.set_ylim(intensity_min, intensity_max)

                radius_ax.clear()
                radius_ax.axhline(0, color='k', lw=0.5)
                if self.root.cali_use_import_button_var.get() == 'True':
                    radius_ax.errorbar(self.root.radius_to_energy_func['converter'](radial), radial_tot,
                                       np.sqrt(radial),
                                       self.root.radius_to_energy_func['error'](np.array(radial),
                                                                                rad_val_err_stat),
                                       fmt='ro', ms=1, lw=0.5, capsize=1.2,
                                       errorevery=1, elinewidth=0.25)
                    radius_ax.set_xlabel('Energy (eV)')
                    radius_ax.set_xlim(0, self.root.radius_to_energy_func['converter'](radius_max))
                else:
                    radius_ax.errorbar(radial, radial_tot, np.sqrt(radial), rad_val_err_stat, fmt='ro', ms=1, lw=0.5,
                                       capsize=1.2,
                                       errorevery=1, elinewidth=0.25)
                    radius_ax.set_xlabel('2D Radius (mm)')
                    radius_ax.set_xlim(0, radius_max)
                radius_ax.set_title('Radial projection')
                radius_ax.set_ylabel(r'Integration over $\theta^\prime$')
                radius_ax.set_ylim(radial_tot_min, radial_tot_max)

                raw_ax.clear()
                raw_ax.imshow(image_val)
                raw_ax.set_title('VMI Laser On')

                raw_off_ax.clear()
                raw_off_ax.imshow(image_off_val)
                raw_off_ax.set_title('VMI Laser Off')

                if self.circ:
                    raw_circ_ax.clear()
                    raw_circ_ax.imshow(circ_image)
                    raw_circ_ax.set_title('VMI Circularized')

                inv_vmi_ax.clear()
                if self.circ:
                    inv_vmi_ax.imshow(combined_vmi)
                    inv_vmi_ax.set_title('Abel inverse / Circularized raw image')
                else:
                    inv_vmi_ax.imshow(inv_image)
                    inv_vmi_ax.set_title('Abel inverse image')

                # Save the plots
                raw_filename = key.GetName() + '_raw.png'
                raw_vmi_files.append(raw_filename)
                raw_off_filename = key.GetName() + '_raw_off.png'
                raw_off_vmi_files.append(raw_off_filename)
                inv_filename = key.GetName() + '_inv.png'
                inv_vmi_files.append(inv_filename)
                intensity_filename = key.GetName() + '_int.png'
                intensity_files.append(intensity_filename)
                radius_filename = key.GetName() + '_rad.png'
                radius_files.append(radius_filename)
                event_filename = key.GetName() + '_evt.png'
                event_files.append(event_filename)

                raw_figure.savefig(raw_filename)
                raw_off_figure.savefig(raw_off_filename)
                inv_vmi_figure.savefig(inv_filename)
                intensity_figure.savefig(intensity_filename)
                radius_figure.savefig(radius_filename)
                event_figure.savefig(event_filename)

                if self.circ:
                    raw_circ_filename = key.GetName() + '_raw_circ.png'
                    raw_circ_vmi_files.append(raw_circ_filename)
                    raw_circ_figure.savefig(raw_circ_filename)

                file_nr += 1

        # build gif
        save_dir = self.root_file
        save_dir = save_dir[:save_dir.rfind('/') + 1]

        self.boot_str.set(f'Building GIF\'s\nVMI laser on image')
        with imageio.get_writer(save_dir + 'raw_vmi.gif', mode='I') as writer:
            for filename in raw_vmi_files:
                image = imageio.v3.imread(filename)
                writer.append_data(image)

        self.boot_str.set(f'Building GIF\'s\nVMI laser off image')
        with imageio.get_writer(save_dir + 'raw_off_vmi.gif', mode='I') as writer:
            for filename in raw_off_vmi_files:
                image = imageio.v3.imread(filename)
                writer.append_data(image)

        if self.circ:
            self.boot_str.set(f'Building GIF\'s\nVMI circularized image')
            with imageio.get_writer(save_dir + 'raw_circ_vmi.gif', mode='I') as writer:
                for filename in raw_circ_vmi_files:
                    image = imageio.v3.imread(filename)
                    writer.append_data(image)

        self.boot_str.set(f'Building GIF\'s\nInverse/Circular VMI image')
        with imageio.get_writer(save_dir + 'inv_vmi.gif', mode='I') as writer:
            for filename in inv_vmi_files:
                image = imageio.v3.imread(filename)
                writer.append_data(image)

        self.boot_str.set(f'Building GIF\'s\nIntensity plot')
        with imageio.get_writer(save_dir + 'intensity.gif', mode='I') as writer:
            for filename in intensity_files:
                image = imageio.v3.imread(filename)
                writer.append_data(image)

        self.boot_str.set(f'Building GIF\'s\nRadial binning plot')
        with imageio.get_writer(save_dir + 'radial.gif', mode='I') as writer:
            for filename in radius_files:
                image = imageio.v3.imread(filename)
                writer.append_data(image)

        self.boot_str.set(f'Building GIF\'s\nEvent times plot')
        with imageio.get_writer(save_dir + 'event.gif', mode='I') as writer:
            for filename in event_files:
                image = imageio.v3.imread(filename)
                writer.append_data(image)

        # Write times to file
        self.root.gif_times = np.array(times).T
        file = open(save_dir + 'gif_times.txt', 'w')
        for i, j in zip(self.root.gif_times[0], self.root.gif_times[1]):
            file.write(f'{i},{j}\n')
        file.close()

        self.boot_str.set(f'Removing files')
        # Remove files
        for filename in set(raw_vmi_files):
            os.remove(filename)
        for filename in set(raw_off_vmi_files):
            os.remove(filename)
        for filename in set(inv_vmi_files):
            os.remove(filename)
        for filename in set(intensity_files):
            os.remove(filename)
        for filename in set(radius_files):
            os.remove(filename)
        for filename in set(event_files):
            os.remove(filename)
        if self.circ:
            for filename in set(raw_circ_vmi_files):
                os.remove(filename)

        self.boot_str.set(f'Loading GIF')

        self.root.event_frames = Image.open(save_dir + 'event.gif')
        self.root.raw_vmi_frames = Image.open(save_dir + 'raw_vmi.gif')
        self.root.raw_off_vmi_frames = Image.open(save_dir + 'raw_off_vmi.gif')
        self.root.inv_vmi_frames = Image.open(save_dir + 'inv_vmi.gif')
        self.root.intensity_frames = Image.open(save_dir + 'intensity.gif')
        self.root.radial_frames = Image.open(save_dir + 'radial.gif')

        if self.circ:
            self.root.circ_frames = True
            self.root.raw_circ_vmi_frames = Image.open(save_dir + 'raw_circ_vmi.gif')

        self.root.max_frames = file_nr
        self.root.frame_slider.configure(to=file_nr)
        self.root.gif_frame_slider.configure(to=file_nr)

        self.destroy()
        self.update()


class GifPopup(tk.Toplevel):
    def __init__(self, root, gif_file_path, root_file):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "Video time")

        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        vcmdi = (self.register(validate_int),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        # Initial parameters needed
        self.geometry('650x375')
        self.root = root
        self.gif_file_path = gif_file_path
        self.root_file = root_file

        ''' labels and widgets '''
        radial_resolution_label = tk.Label(self, text='Radial resolution')
        self.radial_resolution_str = tk.StringVar()
        self.radial_resolution_str.set(root.resolution_entry.get())
        self.radial_resolution_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                                width=5, textvariable=self.radial_resolution_str)

        radial_max_label = tk.Label(self, text='Maximal radius (mm)')
        self.radial_max_str = tk.StringVar()
        self.radial_max_str.set(root.resolution_max_entry.get())
        self.radial_max_entry = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                         width=5, textvariable=self.radial_max_str)

        angular_order_label = tk.Label(self, text='Angular orders')
        self.angular_order_str = tk.StringVar()
        self.angular_order_str.set(root.norder_entry.get())
        self.angular_order_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                            width=5, textvariable=self.angular_order_str)

        vmi_resolution_label = tk.Label(self, text='VMI resolution')
        self.vmi_resolution_str = tk.StringVar()
        self.vmi_resolution_str.set(root.vmi_resolution_entry.get())
        self.vmi_resolution_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                             width=5, textvariable=self.vmi_resolution_str)

        circularize_order_label = tk.Label(self, text='Circularization order')
        self.circularize_order_str = tk.StringVar()
        self.circularize_order_str.set(root.circular_entry.get())
        self.circularize_order_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                                width=5, textvariable=self.circularize_order_str)
        self.circularize_button_str = tk.StringVar()
        if root.circularize:
            self.circularize_button_str.set('On')
        else:
            self.circularize_button_str.set('Off')
        self.circularize_button = tk.Button(self, textvariable=self.circularize_button_str, width=4,
                                            command=lambda: self.circularize_on_off())

        bootstrap_order_label = tk.Label(self, text='Nr. bootstrap samples')
        self.bootstrap_order_str = tk.StringVar()
        self.bootstrap_order_str.set('100')
        self.bootstrap_order_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                              width=5, textvariable=self.bootstrap_order_str)
        self.bootstrap_button_str = tk.StringVar()
        self.bootstrap_button_str.set('Off')
        self.bootstrap_button = tk.Button(self, textvariable=self.bootstrap_button_str, width=4,
                                          command=lambda: self.bootstrap_on_off())

        trigger_time_label = tk.Label(self, text='Trigger time interval (mu s)')
        trigger_time_label_lower = tk.Label(self, text='Lower')
        self.trigger_time_str_lower = tk.StringVar()
        self.trigger_time_str_lower.set(root.trig_lower_entry.get())
        self.trigger_time_entry_lower = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                                 width=7, textvariable=self.trigger_time_str_lower)
        trigger_time_label_upper = tk.Label(self, text='Upper')
        self.trigger_time_str_upper = tk.StringVar()
        self.trigger_time_str_upper.set(root.trig_upper_entry.get())
        self.trigger_time_entry_upper = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                                 width=7, textvariable=self.trigger_time_str_upper)

        event_time_label = tk.Label(self, text='Event time interval (s)')
        event_time_label_lower = tk.Label(self, text='Lower')
        self.event_time_str_lower = tk.StringVar()
        self.event_time_str_lower.set(root.event_lower_entry.get())
        self.event_time_entry_lower = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                               width=7, textvariable=self.event_time_str_lower)
        event_time_label_upper = tk.Label(self, text='Upper')
        self.event_time_str_upper = tk.StringVar()
        self.event_time_str_upper.set(root.event_upper_entry.get())
        self.event_time_entry_upper = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                               width=7, textvariable=self.event_time_str_upper)

        peak_sum_label = tk.Label(self, text='Peak sum interval')
        peak_sum_label_lower = tk.Label(self, text='Lower')
        self.peak_sum_str_lower = tk.StringVar()
        self.peak_sum_str_lower.set(root.peak_sum_lower_entry.get())
        self.peak_sum_entry_lower = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                             width=7, textvariable=self.peak_sum_str_lower)
        peak_sum_label_upper = tk.Label(self, text='Upper')
        self.peak_sum_str_upper = tk.StringVar()
        self.peak_sum_str_upper.set(root.peak_sum_upper_entry.get())
        self.peak_sum_entry_upper = tk.Entry(self, validate='key', validatecommand=vcmdf,
                                             width=7, textvariable=self.peak_sum_str_upper)

        frames_label = tk.Label(self, text='Number of frames')
        self.frames_str = tk.StringVar()
        self.frames_str.set('10')
        self.frames_entry = tk.Entry(self, validate='key', validatecommand=vcmdi,
                                     width=8, textvariable=self.frames_str)

        self.gif_button = tk.Button(self, text='Confirm and calculate frames', command=lambda: self.calc_gif())

        ''' Set the widgets'''
        radial_resolution_label.grid(row=0, column=0, sticky=tk.W, pady=15, padx=10)
        self.radial_resolution_entry.grid(row=0, column=1, pady=15)

        radial_max_label.grid(row=1, column=0, sticky=tk.W, pady=15, padx=10)
        self.radial_max_entry.grid(row=1, column=1, pady=15)

        angular_order_label.grid(row=2, column=0, sticky=tk.W, pady=15, padx=10)
        self.angular_order_entry.grid(row=2, column=1, pady=15)

        vmi_resolution_label.grid(row=3, column=0, sticky=tk.W, pady=15, padx=10)
        self.vmi_resolution_entry.grid(row=3, column=1, pady=15)

        circularize_order_label.grid(row=3, column=0, sticky=tk.W, pady=15, padx=10)
        self.circularize_order_entry.grid(row=3, column=1, pady=15)
        self.circularize_button.grid(row=3, column=2, pady=15)

        bootstrap_order_label.grid(row=4, column=0, sticky=tk.W, pady=15, padx=10)
        self.bootstrap_order_entry.grid(row=4, column=1, pady=15)
        self.bootstrap_button.grid(row=4, column=2, pady=15)

        space = tk.Label(self, text="                           ")
        space.grid(row=0, column=3, sticky=tk.W + tk.E)

        trigger_time_label.grid(row=0, column=4, columnspan=2, pady=10, padx=10)
        trigger_time_label_lower.grid(row=1, column=4, padx=10, sticky=tk.E + tk.N)
        trigger_time_label_upper.grid(row=1, column=4, padx=10, sticky=tk.E + tk.S)
        self.trigger_time_entry_lower.grid(row=1, column=5, padx=10, sticky=tk.E + tk.N)
        self.trigger_time_entry_upper.grid(row=1, column=5, padx=10, sticky=tk.E + tk.S)

        event_time_label.grid(row=2, column=4, columnspan=2, pady=10, padx=10)
        event_time_label_lower.grid(row=3, column=4, padx=10, sticky=tk.E + tk.N)
        event_time_label_upper.grid(row=3, column=4, padx=10, sticky=tk.E + tk.S)
        self.event_time_entry_lower.grid(row=3, column=5, padx=10, sticky=tk.E + tk.N)
        self.event_time_entry_upper.grid(row=3, column=5, padx=10, sticky=tk.E + tk.S)

        peak_sum_label.grid(row=4, column=4, columnspan=2, pady=10, padx=10)
        peak_sum_label_lower.grid(row=5, column=4, padx=10, sticky=tk.E + tk.N)
        peak_sum_label_upper.grid(row=5, column=4, padx=10, sticky=tk.E + tk.S)
        self.peak_sum_entry_lower.grid(row=5, column=5, padx=10, sticky=tk.E + tk.N)
        self.peak_sum_entry_upper.grid(row=5, column=5, padx=10, sticky=tk.E + tk.S)

        frames_label.grid(row=5, column=0, pady=15, sticky=tk.S + tk.W + tk.E)
        self.frames_entry.grid(row=5, column=1)
        self.gif_button.grid(row=6, column=2, columnspan=3, pady=15, sticky=tk.S + tk.W + tk.E)

    def circularize_on_off(self):
        if self.circularize_button_str.get() == 'On':
            self.circularize_button_str.set('Off')
        else:
            self.circularize_button_str.set('On')

    def bootstrap_on_off(self):
        if self.bootstrap_button_str.get() == 'On':
            self.bootstrap_button_str.set('Off')
        else:
            self.bootstrap_button_str.set('On')

    def calc_gif(self):
        # Call to pass
        call = '-r' + self.radial_resolution_entry.get() + " "
        call += '-v' + self.vmi_resolution_entry.get() + ' '
        call += '-o' + self.angular_order_entry.get() + ' '
        call += '-t' + str(float(self.trigger_time_entry_lower.get()) * 1.e-6) + ' '
        call += '-T' + str(float(self.trigger_time_entry_upper.get()) * 1.e-6) + ' '
        call += '-e' + self.event_time_entry_lower.get() + ' '
        call += '-E' + self.event_time_entry_upper.get() + ' '
        call += '-m' + self.radial_max_entry.get() + ' '
        call += '-g' + self.frames_entry.get() + ' '
        call += '-p' + self.peak_sum_entry_lower.get() + ' '
        call += '-P' + self.peak_sum_entry_upper.get() + ' '

        if self.root.background_sub_mode[0] == 1:
            call += '-L '
        else:
            call += '-l '
        if self.root.background_sub_mode[1] == 1:
            call += '-S' + str(self.root.background_sub_mode[3]) + ' '
            call += '-s' + str(self.root.background_sub_mode[2]) + ' '
        else:
            call += '-s-1 '

        if len(self.root.circ_import_file) != 0:
            if self.root.circ_use_import_button_var.get() == 'True':
                call += '-C' + self.root.circ_import_file + ' '
        if self.circularize_button_str.get() == 'On':
            call += '-c' + self.circularize_order_entry.get() + ' '

        if self.bootstrap_button_str.get() == 'On':
            call += '-b' + self.bootstrap_order_entry.get() + ' '

        # Do the calculation
        calculation_popup = GifCalcPopup(self.root, self.gif_file_path, self.root_file, call,
                                         (self.circularize_button_str.get() == 'On') or
                                         ((len(self.root.circ_import_file) != 0) and
                                          (self.root.circ_use_import_button_var.get() == 'True')),
                                         (self.bootstrap_button_str.get() == 'On'))
        self.wait_window(calculation_popup)

        # Plot the results
        self.root.run_gif = True
        self.root.update_gif()


class SaveDataPopup(tk.Toplevel):
    def __init__(self, root):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "Save Data")

        # Initial parameters needed
        self.root = root
        self.filepath = os.path.expanduser("~")

        ''' Check boxes '''
        check_boxes_frame = tk.Frame(self)

        self.time_hist_check_val = tk.IntVar()
        self.time_hist_check_val.set(0)
        time_hist_check = tk.Checkbutton(check_boxes_frame, text="Time Histograms", variable=self.time_hist_check_val,
                                         onvalue=1,
                                         offvalue=0)
        time_hist_check.pack(expand=True, anchor='w')

        self.vmi_images_check_val = tk.IntVar()
        self.vmi_images_check_val.set(0)
        vmi_images_check = tk.Checkbutton(check_boxes_frame, text="VMI Images", variable=self.vmi_images_check_val,
                                          onvalue=1,
                                          offvalue=0)
        vmi_images_check.pack(expand=True, anchor='w')

        self.radial_bin_check_val = tk.IntVar()
        self.radial_bin_check_val.set(0)
        radial_bin_check = tk.Checkbutton(check_boxes_frame, text="Radial Binned Data",
                                          variable=self.radial_bin_check_val, onvalue=1,
                                          offvalue=0)
        radial_bin_check.pack(expand=True, anchor='w')

        self.abel_inv_check_val = tk.IntVar()
        self.abel_inv_check_val.set(0)
        abel_inv_check = tk.Checkbutton(check_boxes_frame, text="Abel Inversed Intensity",
                                        variable=self.abel_inv_check_val, onvalue=1,
                                        offvalue=0)
        abel_inv_check.pack(expand=True, anchor='w')

        self.boot_check_val = tk.IntVar()
        self.boot_check_val.set(0)
        boot_check = tk.Checkbutton(check_boxes_frame, text="Bootstrapped Data", variable=self.boot_check_val,
                                    onvalue=1,
                                    offvalue=0)
        boot_check.pack(expand=True, anchor='w')

        self.gif_check_val = tk.IntVar()
        self.gif_check_val.set(0)
        gif_check = tk.Checkbutton(check_boxes_frame, text="Time Resolved Data", variable=self.gif_check_val, onvalue=1,
                                   offvalue=0)
        gif_check.pack(expand=True, anchor='w')

        check_boxes_frame.grid(row=0, column=0, sticky=tk.NSEW, pady=10, padx=10)

        ''' Save parameters frame '''
        save_param_frame = tk.Frame(self)

        delimiter_label = tk.Label(save_param_frame, text='Delimiter', width=12)
        self.delimiter_var = tk.StringVar()
        self.delimiter_var.set(',')
        delimiter_val = tk.Entry(save_param_frame, width=5, textvariable=self.delimiter_var)

        self.path_label = tk.Label(save_param_frame, text=self.filepath)
        self.choose_path = tk.Button(save_param_frame, text='Select Path', command=lambda: self.search_for_path(),
                                     width=18)

        delimiter_label.grid(row=0, column=0, sticky=tk.W, pady=10)
        delimiter_val.grid(row=0, column=1, sticky=tk.W, pady=10)

        self.choose_path.grid(row=1, column=0, columnspan=2)
        self.path_label.grid(row=2, column=0, columnspan=2)

        save_param_frame.grid(row=0, column=1, sticky=tk.NSEW, pady=10, padx=10)

        ''' Save button '''
        save_button_frame = tk.Frame(self)

        self.name_var = tk.StringVar()
        self.name_var.set('*.txt')
        name = tk.Entry(save_button_frame, textvariable=self.name_var)
        save_button = tk.Button(save_button_frame, text='Save', command=lambda: self.save_function())
        name.grid(row=0, column=0, columnspan=3, sticky=tk.EW)
        save_button.grid(row=0, column=3, sticky=tk.E)

        save_button_frame.grid(row=1, column=1, sticky=tk.NS + tk.E)

    def search_for_path(self):
        directory = self.search_for_dir_path('Please select save path')
        self.filepath = directory
        if os.path.expanduser("~") in directory:
            self.path_label.config(text=os.path.expanduser("~") + '\n' + directory[len(os.path.expanduser("~")):])
        else:
            self.path_label.config(text=directory)

    def search_for_dir_path(self, title):
        current_dir = os.path.expanduser("~")
        file = filedialog.askdirectory(parent=self, initialdir=current_dir, title=title)
        return file

    def save_function(self):

        delim = self.delimiter_var.get()
        extension = ''
        save_file = self.filepath + '/' + self.name_var.get()[:self.name_var.get().find('.')]
        if '.' in self.name_var.get():
            extension = self.name_var.get()[self.name_var.get().find('.'):]

        if self.time_hist_check_val.get():
            # Trigger times
            f = open(save_file + '_trigger_times' + extension, "a")
            f.write('time bins(mu s)' + delim + 'total counts' + delim + 'counts error' +
                    delim + 'current counts' + delim + 'counts error' + '\n')
            for i in range(len(self.root.trig_time_bins)):
                f.write(f'{self.root.trig_time_bins[i]}' + delim + f'{self.root.trig_hist[i]}' + delim +
                        f'{self.root.trig_hist_err[i]}' + delim + f'{self.root.trig_current_hist[i]}' +
                        delim + f'{self.root.trig_current_hist_err[i]}' + '\n')
            f.close()

            # Event times
            f = open(save_file + '_event_times' + extension, "a")
            f.write('time bins(s)' + delim + 'total counts' + delim + 'counts error' +
                    delim + 'current counts' + delim + 'counts error' + '\n')
            for i in range(len(self.root.event_time_bins)):
                f.write(f'{self.root.event_time_bins[i]}' + delim + f'{self.root.event_hist[i]}' + delim +
                        f'{self.root.event_hist_err[i]}' + delim + f'{self.root.event_current_hist[i]}' + delim +
                        f'{self.root.event_current_hist_err[i]}' + '\n')
            f.close()

            # Peak sum
            f = open(save_file + '_peak_sums' + extension, "a")
            f.write('sum bins(s)' + delim + 'total counts' + delim + 'counts error' + delim +
                    'current counts' + delim + 'counts error' + '\n')
            for i in range(len(self.root.peak_sum_bins)):
                f.write(f'{self.root.peak_sum_bins[i]}' + delim + f'{self.root.peak_sum_hist[i]}' + delim +
                        f'{self.root.peak_sum_hist_err[i]}' + delim + f'{self.root.peak_sum_current_hist[i]}' + delim +
                        f'{self.root.peak_sum_current_hist_err[i]}' + '\n')
            f.close()

        if self.vmi_images_check_val.get():
            # Raw laser on
            f = open(save_file + '_VMI_laser_on' + extension, "a")
            for i in self.root.vmi_laser_on:
                for j in range(len(i) - 1):
                    f.write(f'{i[j]}' + delim)
                f.write(f'{i[len(i) - 1]}' + '\n')
            f.close()

            # Raw laser off
            f = open(save_file + '_VMI_laser_off' + extension, "a")
            for i in self.root.vmi_laser_off:
                for j in range(len(i) - 1):
                    f.write(f'{i[j]}' + delim)
                f.write(f'{i[len(i) - 1]}' + '\n')
            f.close()

            # Abel vmi
            f = open(save_file + '_VMI_abel_inv' + extension, "a")
            for i in self.root.abel_vmi:
                for j in range(len(i) - 1):
                    f.write(f'{i[j]}' + delim)
                f.write(f'{i[len(i) - 1]}' + '\n')
            f.close()

            if self.root.circularize:
                # Circular vmi
                f = open(save_file + '_VMI_circular' + extension, "a")
                for i in self.root.circ_vmi:
                    for j in range(len(i) - 1):
                        f.write(f'{i[j]}' + delim)
                    f.write(f'{i[len(i) - 1]}' + '\n')
                f.close()

        if self.radial_bin_check_val.get():
            # Radial bins
            f = open(save_file + '_radial_bins' + extension, "a")
            string = 'radius(mm)' + delim + 'total' + delim
            for n in range(int(self.root.norder)):
                string += f'order-{n}' + delim
            string += f'order-{int(self.root.norder)}\n'
            f.write(string)
            for i in range(len(self.root.radial_projection_values)):
                string = f'{self.root.radial_projection_values[i]}' + delim + f'{self.root.radial_projection_amplitude_sum[i]}'
                for n in range(int(self.root.norder) + 1):
                    string += delim + f'{self.root.radial_projection_amplitude[n * int(self.root.resolution) + i]}'
                f.write(string + '\n')
            f.close()

        if self.abel_inv_check_val.get():
            # Abel inv bins
            rad_err = 0.02 * 20
            f = open(save_file + '_abel_inverse' + extension, "a")
            f.write(
                'radius(mm)' + delim + 'intensity' + delim + 'systematic intensity error' +
                delim + 'systematic radius error\n')
            for i in range(len(self.root.intensity)):
                f.write(f'{self.root.radial_projection_values[i]}' + delim + f'{self.root.intensity[i]}' + delim +
                        f'{self.root.intensity_err[i]}' + delim + f'{rad_err}' + '\n')
            f.close()

        if self.boot_check_val.get():
            # Abel inv bins
            rad_err = (self.root.boot_radial[2] - self.root.boot_radial[1]) / 2
            f = open(save_file + '_bootstrap_abel_inverse' + extension, "a")
            f.write(
                'radius(mm)' + delim + 'intensity' + delim + 'intensity error lower' + 'intensity error upper' +
                delim + 'radius error\n')
            for i in range(len(self.root.boot_radial)):
                f.write(f'{self.root.boot_radial[i]}' + delim + f'{self.root.boot_intensity[i]}' + delim +
                        f'{self.root.quantiles[0][i]}' + delim + f'{self.root.quantiles[1][i]}' + delim +
                        f'{rad_err}' + '\n')
            f.close()

            # Iterations
            f = open(save_file + '_bootstrap_iterations' + extension, "a")
            string = f'iteration-{1}'
            for i in range(1, len(self.root.boot_iterations[0])):
                string += delim + f'iteration-{i + 1}'
            f.write(string + '\n')
            for i in range(len(self.root.boot_iterations)):
                string = f'{self.root.boot_iterations[i][0]}'
                for j in range(1, len(self.root.boot_iterations[0])):
                    string += delim + f'{self.root.boot_iterations[i][j]}'
                f.write(string + '\n')
            f.close()

        if self.gif_check_val.get():
            if os.path.exists(self.root.root_file[:-5] + "_gif.root"):
                gif_file = ROOT.TFile(self.root.root_file[:-5] + "_gif.root", 'READ')
                frame_tree = gif_file.Get("SETTINGS")
                settings_branch = frame_tree.GetBranch("settings")

                #Loading settings

                settings_branch.GetEntry(3)
                nr_frames = int(getattr(frame_tree, "settings/settings"))

                boot = False
                frame_tree = gif_file.Get('FRAME0')
                for branch in frame_tree.GetListOfBranches():
                    if branch.GetName() == "intensity_err_lower":
                        boot = True
                settings_branch = frame_tree.GetBranch("settings")
                settings_branch.GetEntry(0)
                resolution = int(getattr(frame_tree, "settings/settings"))
                settings_branch.GetEntry(2)
                vmi_resolution = (getattr(frame_tree, "settings/settings"))

                #Setting up data structures

                circ_image = np.zeros((nr_frames, vmi_resolution, vmi_resolution))
                inv_image = np.zeros((nr_frames, vmi_resolution, vmi_resolution))
                image_val = np.zeros((nr_frames, vmi_resolution, vmi_resolution))
                image_off_val = np.zeros((nr_frames, vmi_resolution, vmi_resolution))
                intensity = np.zeros((nr_frames, resolution))
                intensity_err = np.zeros((nr_frames, resolution))
                radial_tot = np.zeros((nr_frames, resolution))
                radial_err = np.zeros((nr_frames, resolution))
                radial = np.zeros((nr_frames, resolution))
                rad_val_err = np.ones(resolution) * self.root.sys_radial_err
                rad_val_err_stat = np.ones((nr_frames,resolution))
                if boot:
                    I_err_l = np.zeros((nr_frames, resolution))
                    I_err_u = np.zeros((nr_frames, resolution))

                frame_nr = 0

                # Loading from file

                for key in gif_file.GetListOfKeys():
                    if key.GetName() != "SETTINGS":
                        # Keys and branches
                        frame_tree = gif_file.Get(key.GetName())
                        circ_image_branch = frame_tree.GetBranch("circ_image")
                        inv_image_branch = frame_tree.GetBranch("inv_image")
                        raw_image_branch = frame_tree.GetBranch("raw_image")
                        raw_off_image_branch = frame_tree.GetBranch("raw_off_image")
                        intensity_branch = frame_tree.GetBranch("intensity")
                        intensity_err_branch = frame_tree.GetBranch("intensity_err")
                        rad_tot_branch = frame_tree.GetBranch("radial_bin_TOTAL")
                        rad_val_branch = frame_tree.GetBranch("radial_bin_values")
                        rad_err_branch = frame_tree.GetBranch("radial_bin_error")
                        if boot:
                            I_err_l_branch = frame_tree.GetBranch("intensity_err_lower")
                            I_err_u_branch = frame_tree.GetBranch("intensity_err_upper")

                        # Reading the vmi images
                        count = 0
                        for i in range(inv_image_branch.GetEntries()):
                            circ_image_branch.GetEntry(i)
                            inv_image_branch.GetEntry(i)
                            raw_image_branch.GetEntry(i)
                            raw_off_image_branch.GetEntry(i)

                            circ_image[frame_nr][count // int(vmi_resolution)][count % int(vmi_resolution)] += \
                                getattr(frame_tree, 'circ_image/circ_image')
                            inv_image[frame_nr][count // int(vmi_resolution)][count % int(vmi_resolution)] += \
                                getattr(frame_tree, 'inv_image/inv_image')
                            image_val[frame_nr][count // int(vmi_resolution)][count % int(vmi_resolution)] += \
                                getattr(frame_tree, 'raw_image/raw_image')
                            image_off_val[frame_nr][count // int(vmi_resolution)][count % int(vmi_resolution)] += \
                                getattr(frame_tree, 'raw_off_image/raw_off_image')
                            count += 1

                        # Reading the intensity
                        for i in range(intensity_branch.GetEntries()):
                            intensity_branch.GetEntry(i)
                            intensity_err_branch.GetEntry(i)
                            rad_tot_branch.GetEntry(i)
                            rad_val_branch.GetEntry(i)
                            rad_err_branch.GetEntry(i)
                            if boot:
                                I_err_l_branch.GetEntry(i)
                                I_err_u_branch.GetEntry(i)
                                I_err_l[frame_nr][i] = getattr(frame_tree, 'intensity_err_lower/intensity_err_lower')
                                I_err_u[frame_nr][i] = getattr(frame_tree, 'intensity_err_upper/intensity_err_upper')

                            intensity[frame_nr][i] = getattr(frame_tree, 'intensity/intensity')
                            intensity_err[frame_nr][i] = getattr(frame_tree, 'intensity_err/intensity_err')
                            radial_tot[frame_nr][i] = getattr(frame_tree, 'radial_bin_TOTAL/radial_bin_TOTAL')
                            radial_err[frame_nr][i] = getattr(frame_tree, 'radial_bin_error/radial_bin_error')
                            radial[frame_nr][i] = getattr(frame_tree, 'radial_bin_values/radial_bin_values')
                        rad_val_err_stat[frame_nr][:] = np.ones(int(resolution)) * (radial[frame_nr][-2] - radial[frame_nr][-3])

                        frame_nr += 1

                #Saving to files
                f = open(save_file + '_gif_vmi_on' + extension, "a")
                for fr in range(nr_frames):
                    f.write(f'Frame {fr}\n')
                    for i in image_val[fr]:
                        for j in range(len(i) - 1):
                            f.write(f'{i[j]}' + delim)
                        f.write(f'{i[len(i) - 1]}' + '\n')
                f.close()

                f = open(save_file + '_gif_vmi_off' + extension, "a")
                for fr in range(nr_frames):
                    f.write(f'Frame {fr}\n')
                    for i in image_off_val[fr]:
                        for j in range(len(i) - 1):
                            f.write(f'{i[j]}' + delim)
                        f.write(f'{i[len(i) - 1]}' + '\n')
                f.close()

                f = open(save_file + '_gif_vmi_circ' + extension, "a")
                for fr in range(nr_frames):
                    f.write(f'Frame {fr}\n')
                    for i in circ_image[fr]:
                        for j in range(len(i) - 1):
                            f.write(f'{i[j]}' + delim)
                        f.write(f'{i[len(i) - 1]}' + '\n')
                f.close()

                f = open(save_file + '_gif_inv_vmi' + extension, "a")
                for fr in range(nr_frames):
                    f.write(f'Frame {fr}\n')
                    for i in inv_image[fr]:
                        for j in range(len(i) - 1):
                            f.write(f'{i[j]}' + delim)
                        f.write(f'{i[len(i) - 1]}' + '\n')
                f.close()

                f = open(save_file + '_gif_abel_inverse' + extension, "a")
                for fr in range(nr_frames):
                    f.write(f'{fr} '+ 'radius(mm)' + delim + f'{fr} '+'intensity' + delim + f'{fr} '+'intensity error lower' + f'{fr} '+'intensity error upper' +
                    delim + f'{fr} '+'radius error')
                f.write('\n')
                for i in range(resolution):
                    for fr in range(nr_frames):
                        f.write(f'{radial[fr][i]}' + delim + f'{intensity[fr][i]}' + delim +
                                f'{I_err_l[fr][i]}' + delim + f'{I_err_u[fr][i]}' + delim +
                                f'{rad_val_err_stat[fr][i]}')
                    f.write('\n')
                f.close()

                f = open(save_file + '_gif_radial' + extension, "a")
                for fr in range(nr_frames):
                    f.write(f'{fr} '+ 'radius(mm)' + delim + f'{fr} '+'integrated val' + delim + f'{fr} '+'integrated error' + f'{fr} '+'radius error')
                f.write('\n')
                for i in range(resolution):
                    for fr in range(nr_frames):
                        f.write(f'{radial[fr][i]}' + delim + f'{radial_tot[fr][i]}' + delim +
                                f'{radial_err[fr][i]}' + delim + f'{rad_val_err[i]}')
                    f.write('\n')
                f.close()


class PeakSettings(tk.Toplevel):
    def __init__(self, root):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "Raw Data Fit Settings")

        self.geometry('600x250')
        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        methods = [root.method.get(), 'diff_evol', 'BFGS']

        self.root = root

        self.systematic_errors = tk.Checkbutton(self, text="Systematic Error", variable=root.systematic_errors_var,
                                                onvalue=1,
                                                offvalue=0, command=self.sys_err_check)
        self.statistical_errors = tk.Checkbutton(self, text="Statistical Error", variable=root.statistical_errors_var,
                                                 onvalue=1,
                                                 offvalue=0, command=self.stat_err_check)

        self.methods_dropdown = ttk.OptionMenu(self, root.method, *methods)

        self.radius_label = ttk.Label(self, text='Radius interval')
        self.radius_low_entry = ttk.Entry(self, textvariable=root.radius_low_var, validate='key', validatecommand=vcmdf)
        self.radius_high_entry = ttk.Entry(self, textvariable=root.radius_high_var, validate='key',
                                           validatecommand=vcmdf)

        self.kwargs_text = tk.Entry(self, textvariable=root.kwargs)
        kwargs_label = tk.Label(self, text='kwargs')

        # Setting widgets
        self.systematic_errors.grid(row=0, sticky=tk.W, pady=5)
        self.statistical_errors.grid(row=1, sticky=tk.W, pady=5)
        self.methods_dropdown.grid(row=2, sticky=tk.W, pady=5)
        self.radius_label.grid(row=3, sticky=tk.EW, pady=(10, 5))
        self.radius_low_entry.grid(row=4, sticky=tk.W, pady=5)
        self.radius_high_entry.grid(row=4, sticky=tk.E, pady=5)
        kwargs_label.grid(row=5, sticky=tk.W, pady=(5, 0))
        self.kwargs_text.grid(row=6, sticky=tk.EW, pady=(0, 5))
        self.grid_columnconfigure(0, weight=1, uniform='group1')

    def sys_err_check(self):
        if self.root.statistical_errors_var.get() == 0 \
                and self.root.systematic_errors_var.get() == 0:
            self.root.statistical_errors_var.set(1)

    def stat_err_check(self):
        if self.root.statistical_errors_var.get() == 0 \
                and self.root.systematic_errors_var.get() == 0:
            self.root.systematic_errors_var.set(1)


class CaliSettings(tk.Toplevel):
    def __init__(self, root):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "Calibration Fit Settings")

        self.geometry('600x250')

        methods = [root.method.get(), 'diff_evol', 'BFGS']

        self.root = root

        self.methods_dropdown = ttk.OptionMenu(self, root.method, *methods)

        self.kwargs_text = tk.Entry(self, textvariable=root.kwargs)
        kwargs_label = tk.Label(self, text='kwargs')

        # Setting widgets
        self.methods_dropdown.grid(row=2, sticky=tk.W, pady=5)
        kwargs_label.grid(row=5, sticky=tk.W, pady=(5, 0))
        self.kwargs_text.grid(row=6, sticky=tk.EW, pady=(0, 5))
        self.grid_columnconfigure(0, weight=1, uniform='group1')


class LaserSettings(tk.Toplevel):
    def __init__(self, root):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "Laser Energy Fit Settings")

        self.geometry('600x250')
        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        methods = [root.laser_method.get(), 'diff_evol', 'BFGS']

        self.root = root

        self.statistical_errors = tk.Checkbutton(self, text="Statistical Error",
                                                 variable=root.laser_statistical_errors_var,
                                                 onvalue=1,
                                                 offvalue=0, command=self.stat_err_check)

        self.methods_dropdown = ttk.OptionMenu(self, root.laser_method, *methods)

        self.radius_label = ttk.Label(self, text='Radius interval')
        self.radius_low_entry = ttk.Entry(self, textvariable=root.laser_energy_low_var, validate='key',
                                          validatecommand=vcmdf)
        self.radius_high_entry = ttk.Entry(self, textvariable=root.laser_energy_high_var, validate='key',
                                           validatecommand=vcmdf)

        self.kwargs_text = tk.Entry(self, textvariable=root.laser_kwargs)
        kwargs_label = tk.Label(self, text='kwargs')

        # Setting widgets
        self.statistical_errors.grid(row=1, sticky=tk.W, pady=5)
        self.methods_dropdown.grid(row=2, sticky=tk.W, pady=5)
        self.radius_label.grid(row=3, sticky=tk.EW, pady=(10, 5))
        self.radius_low_entry.grid(row=4, sticky=tk.W, pady=5)
        self.radius_high_entry.grid(row=4, sticky=tk.E, pady=5)
        kwargs_label.grid(row=5, sticky=tk.W, pady=(5, 0))
        self.kwargs_text.grid(row=6, sticky=tk.EW, pady=(0, 5))
        self.grid_columnconfigure(0, weight=1, uniform='group1')

    def stat_err_check(self):
        if self.root.laser_statistical_errors_var.get() == 0:
            self.root.laser_statistical_errors_var.set(1)


class EnergyTable(tk.Toplevel):
    def __init__(self):
        tk.Toplevel.__init__(self)
        tk.Toplevel.wm_title(self, "Energy Lookup Table")
        self.geometry('400x150')

        Cs = tk.LabelFrame(self, text='Cs- --> Cs + e-')
        Cs_lvl1 = tk.Label(Cs, text='S0 --> S1/2\t\t0.471593(36)eV')
        Cs_lvl2 = tk.Label(Cs, text='S0 --> P1/2\t\t1.857521(36)eV')
        Cs_lvl3 = tk.Label(Cs, text='S0 --> P3/2\t\t1.926214(36)eV')
        Cs_lvl1.pack(fill='x', expand=True)
        Cs_lvl2.pack(fill='x', expand=True)
        Cs_lvl3.pack(fill='x', expand=True)
        Cs.pack(fill='both', expand=True)


class EnergyConverter(tk.Toplevel):
    def __init__(self, root):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "Energy Converter")
        self.geometry('450x250')
        self.root = root

        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        radius_label = tk.Label(self, text='Radius')
        energy_label = tk.Label(self, text='Energy')
        radius_err_label = tk.Label(self, text='Error')
        energy_err_label = tk.Label(self, text='Error')

        self.energy_val = tk.StringVar()
        self.energy_val.set('0')
        self.energy_err_val = tk.StringVar()
        self.energy_err_val.set('0')
        energy_val_label = tk.Label(self, textvariable=self.energy_val)
        energy_err_val_label = tk.Label(self, textvariable=self.energy_err_val)

        self.radius_val = tk.StringVar()
        self.radius_val.set('0')
        self.radius_err_val = tk.StringVar()
        self.radius_err_val.set('0')
        radius_entry = tk.Entry(self, textvariable=self.radius_val, validate='key', validatecommand=vcmdf)
        radius_err_entry = tk.Entry(self, textvariable=self.radius_err_val, validate='key', validatecommand=vcmdf)

        calculate_button = tk.Button(self, text='Calculate', command=lambda: self.calculate())

        radius_label.grid(row=0, column=0, pady=5, padx=5)
        energy_label.grid(row=0, column=1, pady=5, padx=5)
        radius_entry.grid(row=1, column=0, pady=5, padx=5)
        energy_val_label.grid(row=1, column=1, pady=5, padx=5)
        radius_err_label.grid(row=2, column=0, pady=5, padx=5)
        energy_err_label.grid(row=2, column=1, pady=5, padx=5)
        radius_err_entry.grid(row=3, column=0, pady=5, padx=5)
        energy_err_val_label.grid(row=3, column=1, pady=5, padx=5)
        calculate_button.grid(row=4, column=0, columnspan=2, pady=5, padx=5)

    def calculate(self):
        self.energy_val.set(str(self.root.radius_to_energy_func['converter'](float(self.radius_val.get()))))
        self.energy_err_val.set(str(
            self.root.radius_to_energy_func['error'](float(self.radius_val.get()), float(self.radius_err_val.get()))))


class GeneralSettings(tk.Toplevel):
    def __init__(self, root):
        tk.Toplevel.__init__(self, root)
        tk.Toplevel.wm_title(self, "General Settings")
        self.root = root

        vcmdf = (self.register(validate_float),
                 '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

        # Frames
        BinSize = tk.LabelFrame(self, text='Bin size control')
        BackgroundMode = tk.LabelFrame(self, text='Background subtraction mode')

        # Bin Size frame
        trigger_label = tk.Label(BinSize, text='Trigger time')
        self.trigger_val = tk.StringVar()
        self.trigger_val.set(str(root.trig_time_bin_size))
        trigger_entry = tk.Entry(BinSize, textvariable=self.trigger_val, validate='key', validatecommand=vcmdf)
        trigger_unit = tk.Label(BinSize, text=r'\mu s')

        event_label = tk.Label(BinSize, text='Event time')
        self.event_val = tk.StringVar()
        self.event_val.set(str(root.event_time_bin_size))
        event_entry = tk.Entry(BinSize, textvariable=self.event_val, validate='key', validatecommand=vcmdf)
        event_unit = tk.Label(BinSize, text=r's')

        peak_label = tk.Label(BinSize, text='Peak sum')
        self.peak_val = tk.StringVar()
        self.peak_val.set(str(root.peak_sum_bin_size))
        peak_entry = tk.Entry(BinSize, textvariable=self.peak_val, validate='key', validatecommand=vcmdf)

        freq_label = tk.Label(BinSize, text='Peak sum')
        self.freq_val = tk.StringVar()
        self.freq_val.set(str(root.freq_bin_size))
        freq_entry = tk.Entry(BinSize, textvariable=self.freq_val, validate='key', validatecommand=vcmdf)

        bin_set = tk.Button(BinSize, text='Set', command=self.set_bin_size)

        # Bin Size draw
        trigger_label.grid(row=0, column=0, sticky=tk.EW, pady=(10, 0))
        trigger_entry.grid(row=1, column=0, sticky=tk.EW)
        trigger_unit.grid(row=1, column=1, sticky=tk.E)

        event_label.grid(row=2, column=0, sticky=tk.EW, pady=(10, 0))
        event_entry.grid(row=3, column=0, sticky=tk.EW)
        event_unit.grid(row=3, column=1, sticky=tk.E)

        peak_label.grid(row=4, column=0, sticky=tk.EW, pady=(10, 0))
        peak_entry.grid(row=5, column=0, sticky=tk.EW)

        freq_label.grid(row=6, column=0, sticky=tk.EW, pady=(10, 0))
        freq_entry.grid(row=7, column=0, sticky=tk.EW)

        bin_set.grid(row=0, column=2, rowspan=8, sticky=tk.NS)

        # Background mode frame
        mode_label = tk.Label(BackgroundMode, text='Laser mode')
        mode_dropdown_options = ['Laser On - Laser Off', 'Laser On - Laser Off', 'Laser Off']
        self.mode_dropdown_var = tk.StringVar()
        self.mode_dropdown = ttk.OptionMenu(BackgroundMode, self.mode_dropdown_var, *mode_dropdown_options)
        self.mode_dropdown_var.set('Laser On - Laser Off')
        if root.background_sub_mode[0] == 0:
            self.mode_dropdown_var.set('Laser Off')
        self.mode_dropdown_var.trace("w", self.set_laser_mode)

        event_time_label = tk.Label(BackgroundMode, text='Event time background window')
        self.event_time_lower_var = tk.StringVar()
        self.event_time_lower_var.set(root.event_background_lower)
        event_time_lower = tk.Entry(BackgroundMode, textvariable=self.event_time_lower_var, validate='key',
                                    validatecommand=vcmdf)
        self.event_time_upper_var = tk.StringVar()
        self.event_time_upper_var.set(root.event_background_upper)
        event_time_upper = tk.Entry(BackgroundMode, textvariable=self.event_time_upper_var, validate='key',
                                    validatecommand=vcmdf)

        self.use_background_window = tk.StringVar()
        if root.background_sub_mode[1] == 0:
            self.use_background_window.set('Not in use')
        else:
            self.use_background_window.set('In use')
        use_background_window = tk.Button(BackgroundMode, textvariable=self.use_background_window,
                                          command=self.use_window)

        # Draw Background mode frame
        mode_label.grid(row=0, column=0, columnspan=2, sticky=tk.EW, pady=(10, 5))
        self.mode_dropdown.grid(row=1, column=0, columnspan=2, sticky=tk.EW, pady=(0, 10))

        event_time_label.grid(row=2, column=0, columnspan=2, sticky=tk.EW, pady=(10, 0))
        event_time_lower.grid(row=3, column=0, sticky=tk.EW)
        event_time_upper.grid(row=3, column=1, sticky=tk.EW)

        use_background_window.grid(row=4, column=0, columnspan=2, sticky=tk.EW)

        # Set frames
        BinSize.grid(row=0, column=0, sticky=tk.NSEW, padx=(0, 5))
        BackgroundMode.grid(row=0, column=1, sticky=tk.NSEW, padx=(5, 0))
        self.grid_columnconfigure(0, weight=1, uniform='group1')
        self.grid_columnconfigure(1, weight=1, uniform='group1')

    def set_bin_size(self):
        self.root.trig_upper = self.root.trig_upper_entry.get()
        self.root.trig_lower = self.root.trig_lower_entry.get()
        trig_low = str(float(self.root.trig_lower) * 1e-6)
        trig_up = str(float(self.root.trig_upper) * 1e-6)
        self.root.event_upper = self.root.event_upper_entry.get()
        self.root.event_lower = self.root.event_lower_entry.get()
        self.root.peak_sum_upper = self.root.peak_sum_upper_entry.get()
        self.root.peak_sum_lower = self.root.peak_sum_lower_entry.get()

        trig_bin_size = str(float(self.trigger_val.get()) * 1e-6)
        self.root.event_time_bin_size = float(self.event_val.get())
        self.root.trig_time_bin_size = float(self.trigger_val.get()) * 1e-6 / 1e-6  # in micro seconds
        self.root.peak_sum_bin_size = float(self.peak_val.get())
        self.root.freq_bin_size = float(self.freq_val.get())

        get_all_time_data(self.root.event_file_path, self.root.root_file, [' -t' + trig_low, ' -T' + trig_up +
                                                                           ' -e' + self.root.event_lower + ' -E' + self.root.event_upper +
                                                                           ' -p' + self.root.peak_sum_lower + ' -P' + self.root.peak_sum_upper +
                                                                           ' -f' + self.freq_val.get() +
                                                                           ' -g' + self.event_val.get() +
                                                                           ' -h' + trig_bin_size +
                                                                           ' -j' + self.peak_val.get()])
        self.root.update_time_hist()
        self.root.plot_trig()
        self.root.plot_event()
        self.root.plot_peak_sum()

    def set_laser_mode(self, *args):
        if self.mode_dropdown_var.get == 'Laser Off':
            self.root.background_sub_mode[0] = 0
            if self.root.circularize and len(self.root.circ_import_file) != 0:
                VMI_filter(self.root.vmi_file_path, self.root.root_file,
                           ['-C' + self.root.circ_import_file + " " + '-l'])
            else:
                VMI_filter(self.root.vmi_file_path, self.root.root_file, [' -l'])
        else:
            self.root.background_sub_mode[0] = 1
            if self.root.circularize and len(self.root.circ_import_file) != 0:
                VMI_filter(self.root.vmi_file_path, self.root.root_file,
                           ['-C' + self.root.circ_import_file + " " + '-L'])
            else:
                VMI_filter(self.root.vmi_file_path, self.root.root_file, [' -L'])
        self.root.update_analysed_data()
        self.root.update_vmi()

    def use_window(self):
        if self.use_background_window.get == 'Not in use':
            self.use_background_window.set('In use')
            self.root.background_sub_mode[1] = 1
            self.root.background_sub_mode[2] = float(self.event_time_lower_var.get())
            self.root.background_sub_mode[3] = float(self.event_time_upper_var.get())
            if self.root.circularize and len(self.root.circ_import_file) != 0:
                VMI_filter(self.root.vmi_file_path, self.root.root_file,
                           ['-C' + self.root.circ_import_file + ' -s' + self.event_time_lower_var.get() +
                            ' -S' + self.event_time_upper_var.get()])
            else:
                VMI_filter(self.root.vmi_file_path, self.root.root_file, [' -s' + self.event_time_lower_var.get() +
                                                                          ' -S' + self.event_time_upper_var.get()])
            self.root.update_analysed_data()
            self.root.update_vmi()
            self.root.plot_event()
        else:
            self.use_background_window.set('Not in use')
            self.root.background_sub_mode[1] = 0
            if self.root.circularize and len(self.root.circ_import_file) != 0:
                VMI_filter(self.root.vmi_file_path, self.root.root_file,
                           ['-C' + self.root.circ_import_file + " " + '-s-1'])
            else:
                VMI_filter(self.root.vmi_file_path, self.root.root_file, [' -s-1'])
            self.root.update_analysed_data()
            self.root.update_vmi()
            self.root.plot_event()


# %%
"""                                              Main                                                                """
if __name__ == '__main__':
    '''              GUI                  '''

    control_app = ControlWindow()
    control_app.mainloop()
