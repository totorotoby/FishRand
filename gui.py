import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from Bioaccum import run_bio


exitflag = False
root = tk.Tk()

class app(tk.Frame):

    def __init__(self, master=None):

        self.r_dictionary = None
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

        ident = ["Region:", "Organism:", "Chemical:"]
        #To be replaced
        regs = ["None"]
        animals = ["None"]
        chemicals = ["None"]
        options = [regs, animals, chemicals]

        for i in range (3,9,2):

            tk.Label(self, text=ident[i//4]).grid(row=i+1, column=2, sticky=tk.W)
            optvar = tk.StringVar(self)
            optvar.set("None")
            opt = tk.OptionMenu(self, optvar, *options[i//4])
            opt.grid(row=i+2, column=2, sticky=tk.W)
            self.viewdict.append(opt)
            self.dictvars.append(optvar)
        
            
        ###choosing what graphs to view###
        
        tk.Label(self, text="Graph Options:").grid(column=3,row=4, sticky=tk.W)
        options = [("CDFs with all fits", 0), ("PDFs with all fits", 1),
                   ("Both with all fits", 2), ("CDF and PDF of specified fit", 3)
                   ]
        count = 4
        self.viewopt = tk.IntVar()
        for text, option in options:
            button = tk.Radiobutton(self, text=text, variable=self.viewopt, value=option)
            button.grid(column=3, row = count, sticky=tk.W)
            count +=1
            #self.viewopt.append(viewnum)
        #command= needs to be edited                                                                  
        getgraphs = tk.Button(self, text="Show Distributions", command=self.show_dist).grid(row=8, column=3, sticky=tk.W)

        
        ttk.Separator(self,orient=tk.VERTICAL).grid(row=1, column=1 , rowspan=10, sticky='ns')
        ttk.Separator(self, orient=tk.VERTICAL).grid(row=1, column=4 , rowspan=10, sticky='ns')

        ###saving data###
        
        tk.Label(self, text="Timestep to save:").grid(column=5,row=4, padx=20)
        timeentry = tk.Entry(self)
        timeentry.grid(column=5,row=5, padx=20, pady=2)
        tk.Label(self, text="Filename:").grid(column=5, row=6, padx=20, pady=2)
        filenameentry = tk.Entry(self)
        filenameentry.grid(column=5,row=7, padx=20,pady=2)
        direcbutton = tk.Button(self, text="Choose Directory", command=filedialog.askdirectory, width = 18)
        direcbutton.grid(column=5,row=8,padx=20,pady=2)
        #Command
        savebutton = tk.Button(self, text="Save", width=18)
        savebutton.grid(column=5, row=9, padx=20,pady=2)

        ###inputs###

        inputbutton = tk.Button(self, text="Choose File", command=self.askfile, width=18)
        inputbutton.grid(column=0, row=6)
        tk.Button(self, text="Run", command=self.run_bio).grid(column=0, row=7)


         ###Fish Image###                                                                                                                         
        fish_image = tk.PhotoImage(file="fishrand.gif")
        label_image0 = tk.Label(self, image=fish_image)
        label_image0.image = fish_image
        label_image0.grid(row=0,column=0, sticky=tk.W)

    def askfile(self):

        self.filename = filedialog.askopenfilename()
        self.parse_filename()

    def run_bio(self):

        if self.filename == '':
            print('need exception here')
        else:
            self.r_dictionary = run_bio(1, self.filename, self.parse_filename())

            regions = list(self.r_dictionary.keys())
            regions = ['None'] + regions

            animalstotal = []

            for value in self.r_dictionary.values():
                keys = list(value.keys())
                chemicals = list(list(value.values())[0].keys())
                for animal in keys:
                    animalstotal.append(animal)

            animalstotal = ['None'] + animalstotal
            chemicals = ['None'] + chemicals
            reset_list  = [regions, animalstotal, chemicals]


            ###refresh menus###
            for i in range (len(self.viewdict)):
                menu = self.viewdict[i]
                var = self.dictvars[i]

                var.set('None')
                menu['menu'].delete(0, 'end')

                for entry in reset_list[i]:
                    menu['menu'].add_command(label=entry, command=tk._setit(var, entry))

    def parse_filename(self):

        pieces = self.filename.split('/')
        return pieces[len(pieces) - 1]

    def show_dist(self):

        print('blah')


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

