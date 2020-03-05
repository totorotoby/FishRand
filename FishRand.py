import tkinter as tk
import matplotlib
matplotlib.use("TkAgg")
from scipy import stats
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
import networkx as nx
from main import *
from FR_Input_Output import *
import subprocess
import sys
import os
from threading import Thread

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
        tk.Label(self, text="Save Results", font=("Times New Roman", 18), height= 2).grid(row=4,column=5, columnspan=1, sticky=tk.W)
        ttk.Separator(self, orient="horizontal").grid(row=3, column=0, columnspan=8, sticky= 'ew')
        ttk.Separator(self, orient="horizontal").grid(row=5, column=0, columnspan=8, sticky= 'ew')

        ###getting reference to region animal chemical###

        ident = ["Time:", 0, "Species:", 0, "Chemical:"]
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
        
        
        options = [("Cumulative Distribution Function(s)", 0), ("Probability Density Function(s)", 1),
                   ("Both CDF(s) and PDF(s)", 2), ("CDF and PDF of selected fit:", 3)
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

        tk.Label(self, text='Type of Distribution to save: ').grid(column=5, row=7)
        self.save_dist_type  = tk.StringVar(self)
        self.save_dist_type.set('Normal')
        tk.OptionMenu(self, self.save_dist_type, "Normal", 'Lognormal', 'Uniform', 'Gamma', "KS Best").grid(column=5, row=8, sticky=tk.E)

        direcbutton = tk.Button(self, text="Save", command=self.askdirectory, width = 18)
        direcbutton.grid(column=5,row=10)

        ###inputs###
        inputbutton = tk.Button(self, text="Choose File", command=self.askfile, width=18)
        inputbutton.grid(column=0, row=7)
        tk.Label(self, text='File: ').grid(column=0, row=6)
        self.filebox = tk.Entry(self)
        self.filebox.grid(column=0, row=8)
        tk.Label(self, text="Timesteps to Save and Display:").grid(column=0,row=9)
        self.timeentry = tk.Entry(self)
        self.timeentry.grid(column=0,row=10)
        self.curvefit = tk.IntVar()
        tk.Checkbutton(self, text="Curve Fit Output Samples", variable=self.curvefit).grid(column=0, row=13)
        tk.Button(self, text='View Foodweb', command=self.show_foodweb).grid(column=0, row=11)
        tk.Button(self, text='View Map', command=self.show_map).grid(column=0, row=12)
        tk.Button(self, text="Run", command=self.loading).grid(column=0, row=14)


         ###Fish Image###
        fish_image = tk.PhotoImage(file= self.mypath + "/fishrand.gif")
        label_image0 = tk.Label(self, image=fish_image)
        label_image0.image = fish_image
        label_image0.grid(row=0,column=0, sticky=tk.W)

        ###Loading window###
        
        self.Loading = tk.Toplevel(self)
        tk.Label(self.Loading, text="Running Simulations...").pack()
        self.bar = tk.ttk.Progressbar(self.Loading, orient='horizontal')
        self.bar.pack()
        self.Loading.withdraw()

    def askfile(self):
        self.filebox.delete(0, tk.END)
        filename = filedialog.askopenfilename(initialdir=self.mypath + '/sheets/input')
        self.filebox.insert(tk.END, filename)

    def askdirectory(self):

        self.savepath = filedialog.asksaveasfilename(initialdir=self.mypath + '/sheets/output', defaultextension=".xlsx")
        dist_type = self.save_dist_type.get()

        if self.output[0] == 'YES':

            steady_state_output(self.to_write, self.stat_check, self.savepath, dist_type, self.inputs, self.data[1], self.tofit)

        if self.output[0] == 'NO':

            temporal_output(self.stat_check, self.to_write, self.savepath, self.time_entry, self.region_info, dist_type, self.inputs, self.data[1], self.tofit)

    def loading(self):

        self.run_thread(self.run_model)

    def run_thread(self, func):
        Thread(target=self.run_func, args=[func]).start()

    def run_func(self, func):

        self.Loading.update()
        self.Loading.deiconify()
        self.bar.start()
        self.run_model()
        self.bar.stop()
        self.Loading.withdraw()
        
    def run_model(self):

        self.filename = self.filebox.get()

        if self.filename == '':
            print('Please choose an input file')

        else:
            # get timesteps

            time_entry = self.timeentry.get().split(",")
            self.tofit = self.curvefit.get()
            try:
                self.time_entry = [int(i)-1 for i in time_entry]
            except:
                self.time_entry = [0]

            # run the code
            self.data = convert_to_lists(self.filename, max(self.time_entry))
            self.output, self.inputs = filter_cases(self.data, self.time_entry, self.tofit)

            if self.output[0] == 'YES':

                self.to_write = self.output[1]
                self.stat_check = self.output[2]
                self.foodweb_graph = self.output[3]


            if self.output[0] == 'NO':
                
                self.to_write = self.output[1]
                self.stat_check = self.output[2]
                self.region_areas = self.output[3]
                self.graph_data = self.output[4]
                self.timescale = self.output[5]
                self.region_info = self.output[6]
                self.foodweb_graph = self.output[7]
                self.lower_graph_data = self.output[8]

       

            ########### TO SET UP INDIVDUAL GRAPHS ###########################
            # update time steps
            times = time_entry
            count = 0
            # get fish names
            if self.stat_check == True and self.output[0] == 'NO':
                try:
                    r = list(self.to_write[0][0])[0]
                    lower = list(self.to_write[0][0][r].keys())
                    fishs = list(self.to_write[0][1].keys())
                    chemicals = list(list(self.to_write[0][1].values())[0].keys())
                    total = lower + fishs
                except:
                    print('You asked to display a timestep that is outside the range between 0, and ((End Time - Start Time) / step) - 1)')
                    exit(0)
                count = 1
            elif self.stat_check == False and self.output[0] == 'NO':
                try:
                    lower = list(self.lower_graph_data[0].keys())
                    fishs = list(self.to_write[0][2].keys())
                    total = lower + fishs
                except:
                    print('Not enough data in either tempatures, chemical concenrations, or abundance for the number of timesteps input')
                    exit(0)
                chemicals = list(list(self.to_write[0][1].values())[0].keys())
                count = 1
            elif self.stat_check == True and self.output[0] == 'YES':
                fishs = list(list(self.to_write.values())[0].keys())
                chemicals = list(list(list(self.to_write.values())[0].values())[0].keys())
                count = 1


            if count == 1:

                reset_list  = [times, total, chemicals]


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
                try:
                    dist_to_show = self.graph_data[int(where[0])-1][where[1]][where[2]]
                except:
                    print('You are missing a graphing argument from the time, fish, chemical selection.')
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

              
                mean = []
                std = []
                try:
                    for dic in self.graph_data:
                        mean.append(dic[fish][chemical].v_mean_std[0])
                        std.append(dic[fish][chemical].v_mean_std[1])

                    fig, ax = plt.subplots()
                    ax.set_ylabel('(ng/g) of ' + chemical + ' in ' + fish, size='large')
                    ax.set_xlabel('Timestep')
                    ax.errorbar([i for i in range(len(self.graph_data))], mean, yerr=std, fmt='o', capsize=2, capthick=2)
                    plt.show()
                    
                except:
                    
                    regs = list(self.lower_graph_data[0].keys())
                    for dic in self.lower_graph_data:
                        m = []
                        for reg in regs:
                            m.append(dic[reg][fish][chemical].v_mean_std[0])

                        mean.append(np.average(m, weights=self.region_areas))

                    fig, ax = plt.subplots()
                    times = list(range(1, len(mean)+1))
                    ax.plot(times, mean, 'ro')
                    ax.set_xlabel('Timesteps (' + self.timescale + ')')
                    ax.set_ylabel('Concentration of ' + chemical+ ' in ' + fish + ' (ng/g)')
                    ax.set_xlim(0)
                    ax.set_ylim(0)
                    plt.show()
                       


            else:
                
                if fish in self.graph_data[0].keys():
                    graph_by_time(self.graph_data, fish, chemical, self.timescale)
                else:
                    graph_by_time(self.lower_graph_data, fish, chemical, self.timescale)

        else:
            print('There is no time graph for steady state.')


    def show_map(self):

        self.filename = self.filebox.get()

        if self.filename == '':
            print('Please choose an input file')
            
        else:

            everything = convert_to_lists(self.filename, 0)
            
            if everything[0][8] == 'NO':
                
                spatial_data = everything[4]
                region_info = pre_run_loc_data(spatial_data)
                b_xy = region_info[0].exterior.xy

                reg_xys = []
                rep_point = []
                for reg in region_info[1]:
                    rep_point.append(reg[1].centroid.coords)
                    reg_xys.append(reg[1].exterior.xy)
    
                
                hotspot_xys=[]
                hotspot_point = []
                attractions = []
                
                for hotspot in region_info[2]:
                    hotspot_point.append(hotspot.polygon.centroid.coords)
                    hotspot_xys.append(hotspot.polygon.exterior.xy)
                    attractions.append(hotspot.attraction)
                    
                fig = plt.figure(1)
                ax = fig.add_subplot(111)
                ax.plot(b_xy[0], b_xy[1],color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
                count = 0
                
                for pair in reg_xys:
                    ax.plot(pair[0], pair[1],color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
                    ax.annotate(region_info[1][count][0],  xy=(rep_point[count][0][0],rep_point[count][0][1]), xytext=(rep_point[count][0][0],rep_point[count][0][1]), ha='center')
                    count += 1

                count = 0
                for pair in hotspot_xys:
                    if pair != b_xy:
                        ax.plot(pair[0], pair[1], color='green', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
                        ax.annotate(region_info[2][count].name + ': ' + str(attractions[count]), xy=(hotspot_point[count][0][0], hotspot_point[count][0][1]), xytext=(hotspot_point[count][0][0], hotspot_point[count][0][1]), ha='center')
                    count += 1

                self.print_calc_regions(reg_xys, [reg[0] for reg in region_info[1]], [reg[1].area for reg in region_info[1]])

                fig.tight_layout()
                plt.show()

                


            else:
                print('Steady state can not be run with regions.')


    def print_calc_regions(self, reg_xys, names, areas):
        
        print('\n')
        max_points = max([len(reg[0]) for reg in reg_xys])

        print('Thessian Polygon generated coordinates:')
        print(f'|{" ":^10}', end='')
        for i in range (1, max_points):
            print(f'|{"Coord. " + str(i):^14}', end='')
        print('|')
        for i in range (len(names)):
            print(f'|{names[i]:^10}|', end='')
            for j in range(max_points-1):
                try:
                    x = round(reg_xys[i][0][j], 2)
                    y = round(reg_xys[i][1][j], 2)
                    print(f'{str(x) + ", " + str(y):^14}|', end='')
                except:
                    print(f'{" ":^14}|', end='')

            print()

        print("Polygon Areas:")
        for i in range(len(names)):
            area = round(areas[i], 2)
            print(f'|{names[i]:^10}|' + f'{str(area):^10}|', end='')
            print()



    def show_foodweb(self):

        filename = self.filebox.get()

        if filename == '':
            print('Please choose an input file')
            
        else:

            foodweb_graph = convert_to_lists(filename, 0)[5]
            
            
            pos = nx.spring_layout(foodweb_graph)

            pos_higher = {}
            for k, v in pos.items():
                pos_higher[k] = (v[0] + .15, v[1])

            labels = {}
            for label in foodweb_graph.nodes:
                labels[label] = label

            nx.draw_networkx_labels(foodweb_graph, pos=pos_higher, labels=labels, font_size=10)
            nx.draw(foodweb_graph, pos, arrows=True)

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

