import tkinter as tk
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from scipy import stats
import numpy as np
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from main import *


exitflag = False
root = tk.Tk()

class app(tk.Frame):

    def __init__(self, master=None):

        self.output = None
        self.viewdict = []
        self.dictvars = []
        self.viewopt = None
        self.filename = ''
        master.title("FishRand")
        super().__init__(master)
        self.grid()
        self.add_widg()

    def add_widg(self):
        ###Top Labels###
        tk.Label(self, text="Input:", font=("Times New Roman", 18), height= 2).grid(row=2,column=0, columnspan=1)
        tk.Label(self, text="View Distributions", font=("Times New Roman", 18), height= 2).grid(row=2,column=2, columnspan=1, sticky=tk.W)
        tk.Label(self, text="Save Results", font=("Times New Roman", 18), height= 2).grid(row=2,column=5, columnspan=1, sticky=tk.W)
        ttk.Separator(self, orient="horizontal").grid(row=1, column=0, columnspan=8, sticky= 'ew')
        ttk.Separator(self, orient="horizontal").grid(row=3, column=0, columnspan=8, sticky= 'ew')

        ###getting reference to region animal chemical###

        ident = ["Time:", 0, "Fish:", 0, "Chemical:"]
        #To be replaced
        time = ["None"]
        fishs = ["None"]
        chemicals = ["None"]
        options = [time, 0, fishs, 0,chemicals]

        for i in range (3,9,2):

            tk.Label(self, text=ident[i-3]).grid(row=i+1, column=2, sticky=tk.W)
            optvar = tk.StringVar(self)
            optvar.set("None")
            opt = tk.OptionMenu(self, optvar, *options[i-3])
            opt.grid(row=i+2, column=2, sticky=tk.W)
            self.viewdict.append(opt)
            self.dictvars.append(optvar)
        
            
        ###choosing what graphs to view###
        
        tk.Label(self, text="Graph Options:").grid(column=3,row=4, sticky=tk.W)
        options = [("CDFs with all fits", 0), ("PDFs with all fits", 1),
                   ("Both with all fits", 2), ("CDF and PDF of selected fit:", 3)
                   ]
        count = 4
        self.viewopt = tk.IntVar()
        for text, option in options:
            button = tk.Radiobutton(self, text=text, variable=self.viewopt, value=option)
            button.grid(column=3, row = count, sticky=tk.W)
            count +=1

        self.optvar1 = tk.StringVar(self)
        self.optvar1.set("None")
        tk.OptionMenu(self, self.optvar1, "Normal", 'Lognormal', 'Uniform', 'Gamma').grid(row=8, column=3, sticky=tk.E)
        tk.Button(self, text="Show Distributions", command=self.show_dist).grid(row=9, column=3, sticky=tk.SW)
        tk.Button(self, text="Show Time Graph", command=self.show_time_graph).grid(row=10, column=3, sticky=tk.SW)

        
        ttk.Separator(self,orient=tk.VERTICAL).grid(row=1, column=1 , rowspan=10, sticky='ns')
        ttk.Separator(self, orient=tk.VERTICAL).grid(row=1, column=4 , rowspan=10, sticky='ns')

        ###saving data###
        
        tk.Label(self, text="Timesteps to Save and Display:").grid(column=0,row=6, padx=20)
        self.timeentry = tk.Entry(self)
        self.timeentry.grid(column=0,row=7, padx=20, pady=2)
        tk.Label(self, text="Filename:").grid(column=5, row=4, padx=20, pady=2)
        self.filenameentry = tk.Entry(self)
        self.filenameentry.grid(column=5,row=5, padx=20,pady=2)
        direcbutton = tk.Button(self, text="Choose Directory", command=self.askdirectory, width = 18)
        direcbutton.grid(column=5,row=6,padx=20,pady=2)
        savebutton = tk.Button(self, text="Save", command=self.save_to_excel, width=18)
        savebutton.grid(column=5, row=7, padx=20,pady=2)

        ###inputs###

        inputbutton = tk.Button(self, text="Choose File", command=self.askfile, width=18)
        inputbutton.grid(column=0, row=5)
        tk.Button(self, text="Run", command=self.run_model).grid(column=0, row=8)


         ###Fish Image###                                                                                                                         
        fish_image = tk.PhotoImage(file="fishrand.gif")
        label_image0 = tk.Label(self, image=fish_image)
        label_image0.image = fish_image
        label_image0.grid(row=0,column=0, sticky=tk.W)

    def askfile(self):

        self.filename = filedialog.askopenfilename()
        self.parse_filename()

    def askdirectory(self):
        self.savedirec = filedialog.askdirectory()

    def run_model(self):

        if self.filename == '':
            print('Please choose an input file')
        else:

            # get timesteps

            time_entry = self.timeentry.get().split(",")
            print(time_entry)
            if time_entry[0] == '':
                print('Need at least 1 Timestep.')
            else:
                self.time_entry = [int(i) for i in time_entry]

            # run the code

            self.output = filter_cases(self.filename, self.time_entry, 'output_test.xlsx')

            if self.output[0] == 'YES':

                self.to_write = self.output[1]
                self.stat_check = self.output[2]



            if self.output[0] == 'NO':

                self.to_write = self.output[1]
                self.stat_check = self.output[2]
                self.region_areas = self.output[3]
                self.graph_data = self.output[4]

            ########### TO SET UP INDIVDUAL GRAPHS ###########################

                if self.stat_check == True:
                    # update time steps
                    times = time_entry
                    # get fish names
                    fishs = list(self.to_write[0][1].keys())
                    chemicals = list(list(self.to_write[0][1].values())[0].keys())
                    chemicals =  chemicals
                    reset_list  = [times, fishs, chemicals]


                    ###refresh menus###
                    for i in range (len(self.viewdict)):
                        menu = self.viewdict[i]
                        var = self.dictvars[i]

                        var.set('None')
                        menu['menu'].delete(0, 'end')

                        for entry in reset_list[i]:
                            menu['menu'].add_command(label=entry, command=tk._setit(var, entry))

            ##################################################################


    def parse_filename(self):

        pieces = self.filename.split('/')
        return pieces[len(pieces) - 1]

    def show_dist(self):


        if self.stat_check == True:

            type_index_list = ['Normal', 'Lognormal', 'Uniform', 'Gamma']
            type_index = self.optvar1.get()

            dist_type = self.viewopt.get()


            where  = []
            for opt in self.dictvars:
                where.append(opt.get())

            dist_to_show = self.graph_data[int(where[0])][where[1]][where[2]]
            dist_to_show.display = dist_type

            dist_to_show.show(type_index_list.index(type_index))

        if self.stat_check == False:

            print('Concentration is deterministic. No Distrubtion to Plot.')


    def show_time_graph(self):

        if self.output[0] == 'NO':

            fish = self.dictvars[1].get()
            chemical = self.dictvars[2].get()

            params = []
            types = []
            mean_stds = []
            for dic in self.graph_data:
                params.append(dic[fish][chemical].best_para[1])
                types.append(dic[fish][chemical].index)
                mean_stds.append(dic[fish][chemical].v_mean_std)



            norm = lambda x, param, offset: offset - stats.norm.pdf(x, loc=param[0], scale=param[1])
            lognorm = lambda x, param, offset: offset - stats.lognorm.pdf(x, s=param[0], loc=param[1], scale=param[2])
            uniform = lambda x, param, offset: offset - stats.uniform.pdf(x, loc=param[0], scale=param[1])
            gamma = lambda x, param, offset: offset - stats.gamma.pdf(x, a=param[0], loc=param[1], scale=param[2])

            type_functions = [norm, lognorm, uniform, gamma]
            fig, ax = plt.subplots()
            ax.tick_params(bottom=False, labelbottom=False)

            max_step = 0
            for param, type, mean_std in zip(params, types, mean_stds):

                if type == 1 or type == 3:
                    y_plot = np.linspace(mean_std[0] - (4*mean_std[1]), mean_std[0] + (6*mean_std[1]), num=500)
                else:
                    y_plot = np.linspace(mean_std[0] - (3 * mean_std[1]), mean_std[0] + (3 * mean_std[1]), num=500)
                values = type_functions[type](y_plot, param, 0)

                if -min(values) > max_step:
                    max_step = -min(values)

            step = max_step + (.1 * max_step)

            timesteps = [i*step for i in range (1, len(params)+1)]

            data = zip(params, types, timesteps, mean_stds)
            timelabel = list(range(len(params)))
            count = 1
            for param, type, time, mean_std in data:

                y_plot = np.linspace(mean_std[0] - (4*mean_std[1]), mean_std[0] + (6*mean_std[1]), num=500)
                ax.plot(type_functions[type](y_plot, param, time), y_plot, color='b')

                mean_point_x = type_functions[type](mean_std[0], param, time)
                ax.plot(mean_point_x, mean_std[0], 'o', color='b')

                ax.annotate(str(count), xy=(time, mean_std[0]), xytext=(time, max(y_plot) + (max(y_plot)*.05)))
                count += 1

            blue_line  = matplotlib.lines.Line2D([], [], color='blue', label='Best Fit Distribution of Timestep')
            blue_dot = matplotlib.lines.Line2D([], [], color='blue', marker = 'o', linestyle = 'None', label='Mean of Samples during Timestep')

            ax.legend(handles = [blue_line, blue_dot])
            ax.set_xlabel('timesteps')
            ax.set_ylabel('Concentration')
            plt.show()

        else:
            print('There is no time graph for steady state.')


    def save_to_excel(self):

        save_name = self.filenameentry.get()
        use_name = self.savedirec + '/' +save_name

        if self.output[0] == 'YES':

            steady_state_output(self.to_write, self.stat_check, use_name)

        if self.output[0] == 'NO':

            temporal_output(self.stat_check, self.to_write, use_name, self.time_entry, self.region_areas)






def closing():
    global exitflag
    global root

    if messagebox.askokcancel("Quit", "Do you want to quit?"):
        exitflag = True
        root.destroy()

def main():

    global exitflag

    root.protocol("WM_DELETE_WINDOW", closing)
    appwind = app(master=root)

    appwind.mainloop()

main()

