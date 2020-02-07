# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:37:00 2019

@author: lzinkiewicz
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import scrolledtext, ttk
import physics as p
from PIL import ImageTk,Image


class MainWindow:
    
    
    
    def __init__(self, master):
        # create tkinter window
        self.window = master
        
         # define variables
        self.optics_model = tk.StringVar()
        self.is_curvature_graph_visible = tk.BooleanVar()
        self.is_legend_visible = tk.BooleanVar()
        self.initial_waist = tk.DoubleVar()
        self.target_waist = tk.DoubleVar()
        self.max_distance = tk.DoubleVar()
        self.wavelength = tk.DoubleVar()
        self.rayleigh_range = tk.DoubleVar()
        self.initial_angle = tk.DoubleVar()
        self.initial_NA = tk.DoubleVar()
        self.lenses = [25.0, 35.0, 50.0, 75.0, 100.0, 125.0, 150.0, 200.0, 250.0, 500.0, 750.0, 1000.0]  # list of available lenses
        self.points = 1000
        self.colors = {True : 'b', False : 'orange'}
        self.is_updating_fields = False
        
        # initialize variables 
        self.optics_model.set("gaussian")
        self.is_curvature_graph_visible.set(False)
        self.is_legend_visible.set(True)
        self.initial_waist.set(0.1)
        self.target_waist.set(2.0)
        self.max_distance.set(500.0)
        self.wavelength.set(532.0)
        self.initial_angle.set(0.05)
        self.initial_angle_changed()
        self.initial_waist_changed()
        self.optics_model.trace('w', self.optics_model_changed)
        self.initial_waist.trace('w', self.initial_waist_changed)
        self.wavelength.trace('w', self.wavelength_changed)
        self.rayleigh_range.trace('w', self.rayleigh_range_changed)
        self.initial_angle.trace('w', self.initial_angle_changed)
        self.initial_NA.trace('w', self.initial_NA_changed)
        
        # other functions preparing layout
        self.vcmd_positive = master.register(self.validate_positive)
        self.vcmd_angle = master.register(self.validate_angle)
        self.vcmd_NA = master.register(self.validate_NA)
        self.create_and_configure_gui_objects()
        self.print_lenses()
        self.prepare_figures()
        
    def create_and_configure_gui_objects(self):
        
        
        # create GUI objects
        self.label_waist_initial = tk.Label(self.window, text="Initial waist size (mm):")
        self.input_waist_initial = tk.Entry(self.window, width=10, textvariable=self.initial_waist, validate='key', validatecommand=(self.vcmd_positive, "%P"))
        self.label_waist_target = tk.Label(self.window, text="Target waist size (mm):")
        self.input_waist_target = tk.Entry(self.window, width=10, textvariable=self.target_waist, validate='key', validatecommand=(self.vcmd_positive, "%P"))
        self.label_size = tk.Label(self.window, text="Maximum distance (mm):")
        self.input_size = tk.Entry(self.window, width=10, textvariable=self.max_distance, validate='key', validatecommand=(self.vcmd_positive, "%P"))
        self.label_wavelength = tk.Label(self.window, text="Wavelength (nm):")
        self.input_wavelength = tk.Entry(self.window, width=10, textvariable=self.wavelength, validate='key', validatecommand=(self.vcmd_positive, "%P"))
        self.label_rayleigh = tk.Label(self.window, text="Rayleigh range (mm):")
        self.input_rayleigh = tk.Entry(self.window, width=10, textvariable=self.rayleigh_range, validate='key', validatecommand=(self.vcmd_positive, "%P"))
        self.label_angle = tk.Label(self.window, text="Initial angle (rad):")
        self.input_angle = tk.Entry(self.window, width=10, textvariable=self.initial_angle, validate='key', validatecommand=(self.vcmd_angle, "%P"))
        self.label_NA = tk.Label(self.window, text="Initial numerical aperture (NA):")
        self.input_NA = tk.Entry(self.window, width=10, textvariable=self.initial_NA, validate='key', validatecommand=(self.vcmd_NA, "%P"))
        self.button_run = tk.Button(self.window, text="CALCULATE", command=self.run_algorithm)
        self.radiobutton_gauss = tk.Radiobutton(self.window, variable = self.optics_model, value = "gaussian", text = "Gaussian beams", command=self.run_algorithm) 
        self.radiobutton_ray = tk.Radiobutton(self.window, variable = self.optics_model, value = "ray", text = "Ray optics", command=self.run_algorithm)
        self.checkbutton_show_curvature = tk.Checkbutton(self.window, text="Show curvature/angle", variable=self.is_curvature_graph_visible, command=self.draw_beam_evolution_graph)
        self.checkbutton_show_legend = tk.Checkbutton(self.window, text="Show legend", variable=self.is_legend_visible, command=self.show_legend)
        self.label_focal_lengths = tk.Label(self.window, text="Available lenses f = (mm):")
        self.scrolledtext_focal_lengths = scrolledtext.ScrolledText(self.window, width=20, height=14)
        self.button_update_lenses = tk.Button(self.window, text="UPDATE", command=self.update_lenses)
        self.label_optimal_parameters = tk.Label(self.window, text="Optimal parameters (mm):")
        self.scrolledtext_optimal_parameters = scrolledtext.ScrolledText(self.window, width=20, height=10)
        self.figure1 = plt.Figure(figsize=(10,4.5))
        self.figure2 = plt.Figure(figsize=(10,4.5))
        self.image1 = FigureCanvasTkAgg(self.figure1, self.window)
        self.image2 = FigureCanvasTkAgg(self.figure2, self.window)
        self.graphics = Image.open("scheme_1.jpg").resize((720,350), Image.ANTIALIAS)
        self.image3 = ImageTk.PhotoImage(image=self.graphics) 
        self.label_scheme = tk.Label(self.window, image=self.image3) 
        self.label_scheme.image = self.image3
        
        # initial configuration of GUI objects
        self.input_angle.config(state='disabled')
        self.input_NA.config(state='disabled')
        self.checkbutton_show_curvature.config(state='disabled')
        
        # position GUI objects
        self.label_waist_initial.grid(column=0, row=0, sticky=tk.W)
        self.input_waist_initial.grid(column=1, row=0)
        self.label_waist_target.grid(column=0, row=1, sticky=tk.W)
        self.input_waist_target.grid(column=1, row=1)
        self.label_size.grid(column=0, row=2, sticky=tk.W)
        self.input_size.grid(column=1, row=2)
        self.label_wavelength.grid(column=2, row=1, sticky=tk.E)
        self.input_wavelength.grid(column=3, row=1)
        self.label_rayleigh.grid(column=2, row=2, sticky=tk.E)
        self.input_rayleigh.grid(column=3, row=2)
        self.label_angle.grid(column=4, row=1, sticky=tk.E)
        self.input_angle.grid(column=5, row=1)
        self.label_NA.grid(column=4, row=2, sticky=tk.E)
        self.input_NA.grid(column=5, row=2)
        self.radiobutton_gauss.grid(column=2, columnspan=2, row=0)
        self.radiobutton_ray.grid(column=4, columnspan=2, row=0)
        self.button_run.grid(column=7, row=0, rowspan=3)
        self.label_focal_lengths.grid(column=7, row=3, sticky=tk.S)
        self.scrolledtext_focal_lengths.grid(column=7, row=4)
        self.button_update_lenses.grid(column=7, row=5, sticky=tk.N)
        self.label_optimal_parameters.grid(column=7, row=6, sticky=tk.S)
        self.scrolledtext_optimal_parameters.grid(column=7, row=7, sticky=tk.N)
        self.checkbutton_show_curvature.grid(column=7, row=8)
        self.checkbutton_show_legend.grid(column=7, row=9)
        self.image1.get_tk_widget().grid(columnspan=6, rowspan=3, row=3)
        self.label_scheme.grid(columnspan=6, rowspan=4, row=6) 
        self.image2.get_tk_widget().grid(columnspan=6, rowspan=4, row=6) 
        self.image2.get_tk_widget().grid_remove()
        tk.ttk.Separator(self.window, orient=tk.VERTICAL).grid(column=6, row=0, rowspan=10, sticky='ns')
        
        
    def prepare_figures(self):
        
        # create axes
        self.ax1 = self.figure1.add_subplot(111)
        self.ax2 = self.figure2.add_subplot(111)
        self.ax3 = self.ax2.twinx()

        # configure axes layout
        self.ax1.set_xlabel("d1 (mm)")
        self.ax1.set_ylabel("d2 (mm)")
        self.ax1.set_xlim(0, self.max_distance.get())
        self.ax1.set_ylim(0, self.max_distance.get())
        self.ax2.set_xlim(0, self.max_distance.get())
        self.ax2.set_ylim(0, self.max_distance.get())
        self.ax1.set_facecolor("white")
        self.ax1.grid(b=False, which='major', axis='both')
        self.ax2.grid(b=False, which='major', axis='both')
        self.ax3.grid(b=False, which='major', axis='both')
        self.ax1.spines['bottom'].set_color('black')
        self.ax1.spines['top'].set_color('black')
        self.ax1.spines['left'].set_color('black')
        self.ax1.spines['right'].set_color('black')
        self.ax1.tick_params('both', reset=True)
        self.ax1.tick_params('both', direction='in')
        self.ax2.set_facecolor("white")
        self.ax3.set_facecolor("white")
        self.ax2.grid(False)
        self.ax3.grid(False)
        self.ax2.spines['bottom'].set_color('black')
        self.ax2.spines['top'].set_color('black')
        self.ax2.spines['left'].set_color('black')
        self.ax2.spines['right'].set_color('black') 
        self.ax3.spines['bottom'].set_color('black')
        self.ax3.spines['top'].set_color('black')
        self.ax3.spines['left'].set_color('black')
        self.ax3.spines['right'].set_color('black') 
        self.ax2.tick_params('both', reset=True)
        self.ax2.tick_params('both', direction='in', colors="white")
        self.ax3.tick_params('both', reset=True)
        self.ax3.tick_params('both', direction='in', colors="white")
    
    def update_lenses(self):
        lines = self.scrolledtext_focal_lengths.get('1.0', tk.END).splitlines()
        lenses_set = set()
        for text in lines:
            text_split = text.split(' ')
            if len(text_split) > 0:
                for f in text_split:
                    if f != ' ' and f != '':
                        try:
                            lenses_set.add(float(f))
                        except:
                           print("Wrong input: ", f, '$')                  
        self.lenses = list(lenses_set)
        self.lenses.sort()
        self.print_lenses()
    
    def print_lenses(self):
        self.scrolledtext_focal_lengths.delete(1.0, tk.END)
        str_to_dsp = ''
        for lens in self.lenses:
            str_to_dsp += f"{lens:.1f}" + '\n'
        self.scrolledtext_focal_lengths.insert(tk.INSERT, str_to_dsp)
    
    def run_algorithm(self):
        
        self.new_problem = p.Problem(self.optics_model.get(), self.lenses, self.initial_waist.get(), self.target_waist.get(), self.initial_angle.get(), self.max_distance.get(), self.points)
        self.solutions = self.new_problem.solve()
        self.draw_solution_map()        
        print(len(self.solutions), "solution(s) were found.")
          
        self.figure1.canvas.mpl_connect('pick_event', self.onpick)
        self.figure1.canvas.mpl_connect('motion_notify_event', self.hover)
       
        
    def draw_solution_map(self):
        self.ax1.cla()
        self.ax2.cla()
        self.ax3.cla()
        self.ax1.set_xlabel("d1 (mm)")
        self.ax1.set_ylabel("d2 (mm)")
        self.annot = self.ax1.annotate('', xy=(0,0))
        self.annot.set_visible(False)
        self.plotted_pts = self.ax1.scatter([s.d1 for s in self.solutions], [s.d2 for s in self.solutions], color = [self.colors.get(s.is_picked) for s in self.solutions], picker=3)
        self.ax1.set_xlim(0, self.max_distance.get())
        self.ax1.set_ylim(0, self.max_distance.get()) 
        self.figure1.canvas.draw()
        
    def draw_beam_evolution_graph(self):
        self.label_scheme.grid_remove()
        self.is_legend_visible.set(False)
        self.image2.get_tk_widget().grid()
        self.checkbutton_show_curvature.config(state='normal')
        self.ax2.cla()
        self.ax3.cla()
         
        positions, waists, curvatures = p.trace_beam(self.optics_model.get(), self.chosen_solution.d1, self.chosen_solution.d2, self.chosen_solution.f1, self.chosen_solution.f2, self.initial_waist.get(), self.wavelength.get()*10**(-6), self.initial_angle.get(), self.max_distance.get()/1000)
        print(len(positions), len(waists), len(curvatures))
        print()
        self.ax2.plot(positions, np.abs(waists), color='black')
        self.ax2.tick_params('both', colors="black", right=False)
        self.ax2.set_xlim(0, np.amax(positions))
        self.ax2.set_xlabel('Position (mm)')
        self.ax2.axvline(x=self.chosen_solution.d1, color='blue', linestyle='--')
        self.ax2.axvline(x=self.chosen_solution.d1+self.chosen_solution.d2, color='blue', linestyle='--')
        self.ax2.annotate('f1 = '+f"{self.chosen_solution.f1:.1f}"+' mm', xy=(self.chosen_solution.d1/positions[-1]-0.16, 1.05), xycoords='axes fraction',
                        bbox=dict(boxstyle="round", fc="orange"))
        self.ax2.annotate('f2 = '+f"{self.chosen_solution.f2:.1f}"+' mm', xy=((self.chosen_solution.d1+self.chosen_solution.d2)/positions[-1], 1.05), xycoords='axes fraction',
                        bbox=dict(boxstyle="round", fc="orange"))
        
        if self.is_curvature_graph_visible.get():
            self.ax3.plot(positions, curvatures, color='orange')
            self.ax3.tick_params('both', colors="black", left=False, right=True, labelleft=False, labelright=True, grid_alpha=0)
            self.ax3.tick_params('y', colors='orange')
            self.ax3.axhline(y=0, color='grey', linestyle='dashdot')
            self.ax3.set_ylim(1.0*np.amin(curvatures), 1.0*np.amax(curvatures))
            if self.optics_model.get() == 'gaussian':
                self.ax3.set_ylabel('1/R (1/mm)', color='orange')
            else:
                self.ax3.set_ylabel('Angle to axis (rad)', color='orange')
        else:
           self.ax3.tick_params('both', colors="black", left=False, right=False, labelleft=False, labelright=False)
           self.ax3.grid(False) 
           
        if self.optics_model.get() == 'gaussian':
            self.ax2.set_ylabel('Waist size (mm)')
        else:
            self.ax2.set_ylabel('Beam size (mm)')

        self.figure2.canvas.draw()        
        
        
        
    def onpick(self, event):
        click_pos_x = event.mouseevent.xdata
        click_pos_y = event.mouseevent.ydata 
        self.chosen_solution = self.point_to_solution(click_pos_x, click_pos_y)
        print("Solution = ", self.chosen_solution.d1, self.chosen_solution.d2)
        self.draw_solution_map()
        self.draw_beam_evolution_graph()
        self.print_optimal_parameters()  
     
    def point_to_solution(self, x_clk, y_clk):
        distance = self.max_distance.get()
        idx = 0
        for i in range(len(self.solutions)):
            self.solutions[i].is_picked = False
            d = dst(x_clk, y_clk, self.solutions[i].d1, self.solutions[i].d2)
            if d < distance:
                distance = d
                idx = i
        self.solutions[idx].is_picked = True
        return self.solutions[idx]

    def print_optimal_parameters(self):
        str_to_dsp = 'd1 = ' + f"{self.chosen_solution.d1:.1f}" + '\nd2 = ' + f"{self.chosen_solution.d2:.1f}" + '\nf1 = ' + f"{self.chosen_solution.f1:.1f}" + '\nf2 = ' + f"{self.chosen_solution.f2:.1f}" + '\nD tot. = ' + f"{self.chosen_solution.d1+self.chosen_solution.d2:.1f}"
        self.scrolledtext_optimal_parameters.delete(1.0, tk.END)
        self.scrolledtext_optimal_parameters.insert(tk.INSERT, str_to_dsp)
        
    def hover(self, event):
        if event.inaxes == self.ax1:
            cont, ind = self.plotted_pts.contains(event)
            if cont:
                self.annot.set_visible(False)
                pos = self.plotted_pts.get_offsets()[ind["ind"][0]]
                for s in self.solutions:
                    if pos[0] == s.d1 and pos[1] == s.d2:
                        text = 'f1 = ' + f"{s.f1:.1f}" + ' mm  f2 = ' + f"{s.f2:.1f}" + ' mm\nD = ' + f"{s.d1+s.d2:.1f}" + ' mm'
                        if s.d1 < 0.3*self.max_distance.get():
                            xx = 0
                        else:
                            xx = -150
                        if s.d2 < 0.1*self.max_distance.get():
                            yy = 50
                        else:
                            yy = -50
                        break
                self.annot = self.ax1.annotate(text, xy=pos, xytext=(xx, yy),textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="orange"),
                            arrowprops=dict(arrowstyle="-", color='black'))
                self.annot.set_visible(True)
            else:
                self.annot.set_visible(False)
        
            self.figure1.canvas.draw_idle()

    def show_legend(self):
        if self.is_legend_visible.get() == True:
            self.label_scheme.grid()
            self.image2.get_tk_widget().grid_remove()
        else:
            self.label_scheme.grid_remove()
            self.image2.get_tk_widget().grid()
    
    def optics_model_changed(self, *args):
        if self.optics_model.get() == 'gaussian':
            self.input_angle.config(state='disabled')
            self.input_NA.config(state='disabled')
            self.input_wavelength.config(state='normal')
            self.input_rayleigh.config(state='normal')
            self.initial_waist_changed()
        else:
            self.input_angle.config(state='normal')
            self.input_NA.config(state='normal')
            self.input_wavelength.config(state='disabled')
            self.input_rayleigh.config(state='disabled')
        
    def wavelength_changed(self, *args):
        if self.is_updating_fields: return
        new_rayleigh_range = np.pi*self.initial_waist.get()**2/(self.wavelength.get()*10**(-6))
        self.is_updating_fields = True
        self.rayleigh_range.set(f"{new_rayleigh_range:.3f}")
        self.is_updating_fields = False
        
    def rayleigh_range_changed(self, *args):
        if self.is_updating_fields: return
        new_wavelength = (10**6)*np.pi*self.initial_waist.get()**2/self.rayleigh_range.get()
        self.is_updating_fields = True
        self.wavelength.set(f"{new_wavelength:.1f}")
        self.is_updating_fields = False
        
    def initial_waist_changed(self, *args):
        if self.optics_model.get() == 'gaussian':
            if self.is_updating_fields: return
            new_rayleigh_range = np.pi*self.initial_waist.get()**2/(self.wavelength.get()*10**(-6))
            self.is_updating_fields = True
            self.rayleigh_range.set(f"{new_rayleigh_range:.3f}")
            self.is_updating_fields = False
        
    def initial_angle_changed(self, *args):
        if self.is_updating_fields: return
        self.is_updating_fields = True
        self.initial_NA.set(f"{np.sin(self.initial_angle.get()):.3f}")
        self.is_updating_fields = False
        
    def initial_NA_changed(self, *args):
        if self.is_updating_fields: return
        self.is_updating_fields = True
        self.initial_angle.set(f"{np.arcsin(self.initial_NA.get()):.3f}")
        self.is_updating_fields = False
    
    def validate_positive(self, number):
        try:
            if float(number) >= 0:
                return True
            else:
                return False
        except:
            return False
        
    def validate_angle(self, number):
        try:
            if float(number) > -np.pi/2 and float(number) < np.pi/2:
                return True
            else:
                return False
        except:
            return False
        
    def validate_NA(self, number):
        try:
            if float(number) >= 0 and float(number) < 1:
                return True
            else:
                return False
        except:
            return False
        
        
def dst(x1, y1, x2, y2):
    return np.sqrt((x1-x2)**2+(y1-y2)**2)



