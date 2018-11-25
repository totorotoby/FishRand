import tkinter as tk
import matplotlib
matplotlib.use("TkAgg")
from scipy import stats
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
import networkx as nx
from main import *
import subprocess
import sys
import os


exitflag = False
root = tk.Tk()

class app(tk.Frame):

    def __init__(self, master=None):

        self.output = None
        self.viewdict = []
        self.dictvars = []
        self.viewopt = None
        self.filename = ''
        self.mypath = os.path.dirname(os.path.abspath(__file__))
        master.title("FishRand")
        super().__init__(master)
        self.grid()
        self.add_widg()

    def add_widg(self):
        ###Top Labels###
        tk.Label(self, text="Input", font=("Times New Roman", 18), height= 2).grid(row=4,column=0, columnspan=1)
        tk.Label(self, text="View Distributions", font=("Times New Roman", 18), height= 2).grid(row=4,column=2, columnspan=1, sticky=tk.W)
        tk.Label(self, text= "Other Info", font=("Times New Roman", 18), height=2).grid(row=0, column=2, columnspan=1, sticky=tk.SW)
        tk.Label(self, text="Save Results", font=("Times New Roman", 18), height= 2).grid(row=4,column=5, columnspan=1, sticky=tk.W)
        ttk.Separator(self, orient="horizontal").grid(row=3, column=0, columnspan=8, sticky= 'ew')
        ttk.Separator(self, orient="horizontal").grid(row=5, column=0, columnspan=8, sticky= 'ew')

        ###getting reference to region animal chemical###

        ident = ["Time:", 0, "Fish:", 0, "Chemical:"]
        #To be replaced
        time = ["None"]
        fishs = ["None"]
        chemicals = ["None"]
        options = [time, 0, fishs, 0,chemicals]

        for i in range (5,11,2):

            tk.Label(self, text=ident[i-5]).grid(row=i+1, column=2, sticky=tk.W)
            optvar = tk.StringVar(self)
            optvar.set("None")
            opt = tk.OptionMenu(self, optvar, *options[i-5])
            opt.grid(row=i+2, column=2, sticky=tk.W)
            self.viewdict.append(opt)
            self.dictvars.append(optvar)
        
            
        ###choosing what graphs to view###
        
        tk.Label(self, text="Graph Options:").grid(column=3,row=6, sticky=tk.W)
        options = [("CDFs with all fits", 0), ("PDFs with all fits", 1),
                   ("Both with all fits", 2), ("CDF and PDF of selected fit:", 3)
                   ]
        count = 6
        self.viewopt = tk.IntVar()
        for text, option in options:
            button = tk.Radiobutton(self, text=text, variable=self.viewopt, value=option)
            button.grid(column=3, row = count, sticky=tk.W)
            count +=1

        self.optvar1 = tk.StringVar(self)
        self.optvar1.set("Normal")
        tk.OptionMenu(self, self.optvar1, "Normal", 'Lognormal', 'Uniform', 'Gamma').grid(row=10, column=3, sticky=tk.E)
        tk.Button(self, text="Show Distributions", command=self.show_dist).grid(row=11, column=3, sticky=tk.SW)
        tk.Button(self, text="Show Time Graph", command=self.show_time_graph).grid(row=12, column=3, sticky=tk.SW)

        
        ttk.Separator(self,orient=tk.VERTICAL).grid(row=3, column=1 , rowspan=10, sticky='ns')
        ttk.Separator(self, orient=tk.VERTICAL).grid(row=3, column=4 , rowspan=10, sticky='ns')


        ###Other info####

        tk.Button(self, text='Show Regions', command=self.show_regions).grid(row=1, column=2, columnspan=1, sticky=tk.NW)
        tk.Button(self, text='Show Foodweb', command=self.show_foodweb).grid(row=2, column=2, columnspan=1, sticky=tk.NW)


        ###saving data###
        self.filebox = tk.Entry(self)
        self.filebox.grid(column=0, row=8)
        tk.Label(self, text="Timesteps to Save and Display:").grid(column=0,row=9)
        self.timeentry = tk.Entry(self)
        self.timeentry.grid(column=0,row=10)


        tk.Label(self, text='Type of Distrubtion to save: ').grid(column=5, row=7)
        self.save_dist_type  = tk.StringVar(self)
        self.save_dist_type.set('Normal')
        tk.OptionMenu(self, self.save_dist_type, "Normal", 'Lognormal', 'Uniform', 'Gamma', "KS Best").grid(column=5, row=8, sticky=tk.E)

        direcbutton = tk.Button(self, text="Save", command=self.askdirectory, width = 18)
        direcbutton.grid(column=5,row=10)

        ###inputs###
        inputbutton = tk.Button(self, text="Choose File", command=self.askfile, width=18)
        inputbutton.grid(column=0, row=7)
        tk.Label(self, text='File: ').grid(column=0, row=6)
        tk.Button(self, text="Run", command=self.run_model).grid(column=0, row=11)


         ###Fish Image###
        fish_image = tk.PhotoImage(file= self.mypath + "/fishrand.gif")
        label_image0 = tk.Label(self, image=fish_image)
        label_image0.image = fish_image
        label_image0.grid(row=0,column=0, sticky=tk.W)

    def askfile(self):
        self.filebox.delete(0, tk.END)
        filename = filedialog.askopenfilename(initialdir=self.mypath + '/sheets/input')
        self.filebox.insert(tk.END, filename)

    def askdirectory(self):

        self.savepath = filedialog.asksaveasfilename(initialdir=self.mypath + '/sheets/output', defaultextension=".xlsx")
        dist_type = self.save_dist_type.get()

        if self.output[0] == 'YES':

            steady_state_output(self.to_write, self.stat_check, self.savepath, dist_type)

        if self.output[0] == 'NO':

            temporal_output(self.stat_check, self.to_write, self.savepath, self.time_entry, self.region_areas, dist_type)

    def run_model(self):
        self.filename = self.filebox.get()

        if self.filename == '':
            print('Please choose an input file')

        else:
            # get timesteps

            time_entry = self.timeentry.get().split(",")
            try:
                self.time_entry = [int(i) for i in time_entry]
            except:
                self.time_entry = [0]

            # run the code
            self.output = filter_cases(self.filename, self.time_entry)

            if self.output[0] == 'YES':

                self.to_write = self.output[1]
                self.stat_check = self.output[2]
                self.foodweb_graph = self.output[3]


            if self.output[0] == 'NO':
                #print(self.output)
                self.to_write = self.output[1]
                self.stat_check = self.output[2]
                self.region_areas = self.output[3]
                self.graph_data = self.output[4]
                self.timescale = self.output[5]
                self.region_info = self.output[6]
                self.foodweb_graph = self.output[7]


            ########### TO SET UP INDIVDUAL GRAPHS ###########################
            # update time steps
            times = time_entry
            count = 0
            # get fish names
            if self.stat_check == True and self.output[0] == 'NO':
                fishs = list(self.to_write[0][1].keys())
                chemicals = list(list(self.to_write[0][1].values())[0].keys())
                count = 1
            elif self.stat_check == False and self.output[0] == 'NO':
                fishs = list(self.to_write[0][2].keys())
                chemicals = list(list(self.to_write[0][1].values())[0].keys())
                count = 1
            elif self.stat_check == True and self.output[0] == 'YES':
                fishs = list(list(self.to_write.values())[0].keys())
                chemicals = list(list(list(self.to_write.values())[0].values())[0].keys())
                count = 1


            if count == 1:

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

            if self.output[0] == 'NO':
                dist_to_show = self.graph_data[int(where[0])][where[1]][where[2]]
            else:
                dist_to_show = list(self.to_write.values())[0][where[1]][where[2]]
            dist_to_show.display = dist_type

            dist_to_show.show(type_index_list.index(type_index))

        if self.stat_check == False:

            print('Concentration is deterministic. No Distrubtion to Plot.')


    def show_time_graph(self):

        if self.output[0] == 'NO':


            fish = self.dictvars[1].get()
            chemical = self.dictvars[2].get()

            if self.stat_check == True:

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
                        y_plot = np.linspace(mean_std[0] - (2*mean_std[1]), mean_std[0] + (2*mean_std[1]), num=500)
                    else:
                        y_plot = np.linspace(mean_std[0] - (2 * mean_std[1]), mean_std[0] + (2 * mean_std[1]), num=500)
                    values = type_functions[type](y_plot, param, 0)

                    if -min(values) > max_step:
                        max_step = -min(values)

                step = max_step + (.1 * max_step)

                timesteps = [i*step for i in range (1, len(params)+1)]

                data = zip(params, types, timesteps, mean_stds)
                timelabel = list(range(len(params)))
                count = 1
                for param, type, time, mean_std in data:

                    y_plot = np.linspace(mean_std[0] - (3*mean_std[1]), mean_std[0] + (3*mean_std[1]), num=500)
                    ax.plot(type_functions[type](y_plot, param, time), y_plot, color='b')

                    mean_point_x = type_functions[type](mean_std[0], param, time)
                    ax.plot(mean_point_x, mean_std[0], 'o', color='b')

                    ax.annotate(str(count), xy=(time, mean_std[0]), xytext=(time, max(y_plot) + (max(y_plot)*.05)))
                    count += 1

                blue_line  = matplotlib.lines.Line2D([], [], color='blue', label='Distrubtion of Best Fit during Timestep')
                blue_dot = matplotlib.lines.Line2D([], [], color='blue', marker = 'o', linestyle = 'None', label='Mean of Samples during Timestep')

                ax.legend(handles = [blue_line, blue_dot])
                ax.set_xlabel('Timesteps (' + self.timescale + ')')
                ax.set_ylabel('Concentration of ' + chemical + ' in ' + fish + ' (ng/g)')
                ax.set_title('Distrubtions of Best Fit over Time')
                fig.tight_layout()
                plt.show()

            else:

                graph_by_time(self.graph_data, fish, chemical, self.timescale)

        else:
            print('There is no time graph for steady state.')


    def show_regions(self):

        if self.output[0] == 'NO':

            b_xy = self.region_info[0].exterior.xy

            reg_xys = []
            rep_point = []
            for reg in self.region_info[1]:
                rep_point.append(reg[1].centroid.coords)
                reg_xys.append(reg[1].exterior.xy)

            hotspot_xys=[]
            hotspot_point = []
            for hotspot in self.region_info[2]:
                hotspot_point.append(hotspot.polygon.centroid.coords)
                hotspot_xys.append(hotspot.polygon.exterior.xy)

            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            ax.plot(b_xy[0], b_xy[1],color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
            count = 0
            for pair in reg_xys:
                ax.plot(pair[0], pair[1],color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
                ax.annotate(self.region_info[1][count][0],  xy=(rep_point[count][0][0],rep_point[count][0][1]), xytext=(rep_point[count][0][0],rep_point[count][0][1]), ha='center')
                count += 1

            count = 0
            for pair in hotspot_xys:
                ax.plot(pair[0], pair[1], color='green', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
                ax.annotate(self.region_info[2][count].name, xy=(hotspot_point[count][0][0], hotspot_point[count][0][1]), xytext=(hotspot_point[count][0][0], hotspot_point[count][0][1]), ha='center')
                count += 1

            fig.tight_layout()
            plt.show()

        else:
            print('Steady state can not be run with regions.')

    def show_foodweb(self):

        pos = nx.spectral_layout(self.foodweb_graph)

        pos_higher = {}
        for k, v in pos.items():
            pos_higher[k] = (v[0] + .15, v[1])

        labels = {}
        for label in self.foodweb_graph.nodes:
            labels[label] = label

        nx.draw_networkx_labels(self.foodweb_graph, pos=pos_higher, labels=labels, font_size=10)
        nx.draw(self.foodweb_graph, pos, arrows=True)

        lowerx, upperx = plt.xlim()
        lowery, uppery = plt.ylim()
        plt.xlim(lowerx - .05, upperx + .08)
        plt.ylim(lowery - .05, uppery + .08)
        plt.show()



def closing():
    global exitflag
    global root

    if messagebox.askokcancel("Quit", "Do you want to quit?"):
        exitflag = True
        root.destroy()


def install(package):

    subprocess.call([sys.executable, "-m", "pip", "install", package])

def main():

    global exitflag

    root.protocol("WM_DELETE_WINDOW", closing)
    appwind = app(master=root)

    appwind.mainloop()

main()

