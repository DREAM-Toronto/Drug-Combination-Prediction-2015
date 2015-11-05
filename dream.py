from tkinter import *
from tkinter import messagebox
import csv, os, webbrowser

# Directory 
script_dir = os.getcwd() #<-- absolute dir the script is in
mono = "../Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_ch2_monoTherapy"
comb = "../Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_training_combinations"
abs_mono = os.path.join(script_dir, mono)
abs_comb = os.path.join(script_dir, comb)

# Training data file
path = os.path.join(script_dir, "../Drug Synergy Data")
file_lea = open(path+'/ch1_leaderBoard_monoTherapy.csv', "rt")
file_tra = open(path+'/ch1_train_combination_and_monoTherapy.csv', "rt")
data = csv.reader(file_tra, delimiter=' ', quotechar='|')
data_lea = csv.reader(file_lea, delimiter=' ', quotechar='|')

# Store data in dict
train_dict = {}
lea_dict = {}

row_no = 1
for dat_a in data:
    values = dat_a[0].split(',')
    if (values[0][1:-1]!='CELL_LINE'): #and values[0]=='HCC1143'):
        row_no+=1
        train_dict[values[-1][1:-1]+"."+values[0][1:-1]] = [values[0][1:-1],
                                                            values[1][1:-1],values[2][1:-1],values[3][1:-1],
                                                            values[4][1:-1],values[5],values[6],values[7],
                                                            values[8],values[9],values[10],values[11],
                                                            values[12],values[13][1:-1],row_no]

row_no_lea = 1
b = 1
for dat_a_lea in data_lea:
    #print(dat_a_lea)
    if  dat_a_lea != []:
        values_lea = dat_a_lea[0].split(',')
    if (values_lea[0]!='CELL_LINE'): #and values_lea[0]=='HCC1143'):
        row_no_lea+=1
        lea_dict [values_lea[-1]+"."+values_lea[0]] = [values_lea[0],
                                                            values_lea[1],values_lea[2],values_lea[3],
                                                            values_lea[4],values_lea[5],values_lea[6],values_lea[7],
                                                            values_lea[8],values_lea[9],values_lea[10],values_lea[11],
                                                            values_lea[12],values_lea[13],row_no_lea]
        b = -1
    b = 1

print(train_dict)   
print(lea_dict)


class TraDataGUI:
    def __init__(self, root):
        self.root = root
        root.title("DREAM Toronto")
        self.seaSea = StringVar()
        self.inSea = StringVar()
        
        self.scrollbar = Scrollbar(fra, orient="vertical", command=self.yview)
        self.scrollbar.pack(side=RIGHT, fill=Y)
        
        self.sear = Label(fra1, text='Search: "CompoundX" or "CompA.CompB" or ?CompX.CELL_LINE" or "RowNo"').pack(side=TOP)
        self.ent = Entry(fra1,textvariable=self.inSea).pack(side=LEFT,padx=50,pady=10)
        
        self.butS = Button(fra1, text ="Search",
                           command = self.callSearch).pack(side=LEFT,padx=5,pady=5)      
        
        self.col_names = Label(fra, text=" COMPOUND_A.COMPOUND_B.CELL_LINE  | SYNERGY_SCORE | ROW_NO ").pack(side=TOP)
        
        self.comb_names = Listbox(fra, yscrollcommand = self.scrollbar.set, 
                                  selectmode=SINGLE)
        self.row = Listbox(fra, yscrollcommand = self.scrollbar.set, 
                                         selectmode=SINGLE)
        self.act_syn = Listbox(fra, yscrollcommand = self.scrollbar.set, 
                                                 selectmode=SINGLE)        
        self.row.pack(side = RIGHT, fill = Y)
        self.row.configure(width=12, height=30) 
        self.comb_names.pack(side = LEFT, fill = Y)
        self.comb_names.configure(width=35, height=30) 
        self.act_syn.pack(side = RIGHT, fill = Y)
        self.act_syn.configure(width=12, height=30)        
        self.comb_names.bind("<MouseWheel>", self.OnMouseWheel)     
        self.row.bind("<MouseWheel>", self.OnMouseWheel)
        self.act_syn.bind("<MouseWheel>", self.OnMouseWheel)
        root.bind('<Return>', self.callSearch)

        #adding drug combinations to list
        self.j = -1
        for self.dcomb in train_dict:
            self.common_Method()
                                                      
        self.butA = Button(fra2, text ="Combination",command = self.callA)
        self.butA.pack(side=LEFT,padx=5,pady=10)
        
        self.butB = Button(fra2, text ="Monothereapy",command = self.callB)
        self.butB.pack(side=LEFT,padx=5,pady=10)
        
        self.butC = Button(fra2, text ="TrainingData",command = self.callC)
        self.butC.pack(side=LEFT,padx=5,pady=10)  
        
        self.butD = Button(fra2, text ="TestData",command = self.callD)
        self.butD.pack(side=LEFT,padx=5,pady=10)         
                
    def callA(self):
        if self.comb_names.curselection():
            cur_selection = self.comb_names.get(self.comb_names.curselection()) \
            + ".Rep1.csv"
            webbrowser.open_new(abs_comb + "//" + cur_selection)
        else:
            messagebox.showerror("Error", "Please select a Drug Combination First!")
        
    def callB(self):
        if self.comb_names.curselection():
            cur_selection = self.comb_names.get(self.comb_names.curselection()) \
                + ".Rep1.csv"
            webbrowser.open_new(abs_mono + "//" + cur_selection)
        else:
            messagebox.showerror("Error", "Please select a Drug Combination First!")
    
    def callC(self):
        webbrowser.open_new(path+'/ch1_train_combination_and_monoTherapy.csv')   
    def callD(self):
        webbrowser.open_new(path+'/ch1_test_monoTherapy.csv')   
        
    def callSearch(self, event=None):
        self.comb_names.delete(0,END)
        self.row.delete(0,END)
        self.act_syn.delete(0,END)
        self.seaSea = self.inSea.get()
        self.seaSea = self.seaSea.upper()
        if self.seaSea == "":
            self.j = -1
            for self.dcomb in train_dict:
                self.common_Method()
        else:
            self.j = -1
            for self.dcomb in train_dict:
                y = 0
                if self.seaSea[0] == "?":
                    char = self.seaSea.find(".")
                    self.seaSea1 = self.seaSea[1:char]
                    self.seaSea2 = self.seaSea[char+1:]
                    y = 1
                    
                if y == 0:
                    if (train_dict[self.dcomb][0] == self.seaSea or 
                        train_dict[self.dcomb][1] == self.seaSea or 
                        train_dict[self.dcomb][2] == self.seaSea or 
                        train_dict[self.dcomb][-2] == self.seaSea or 
                        str(train_dict[self.dcomb][-1]) == self.seaSea):
                        self.common_Method()
                elif y==1:
                    if ((train_dict[self.dcomb][0] == self.seaSea2 or 
                         train_dict[self.dcomb][1] == self.seaSea2  or 
                         train_dict[self.dcomb][2] == self.seaSea2 or 
                         train_dict[self.dcomb][-2] == self.seaSea2 or 
                         str(train_dict[self.dcomb][-1]) == self.seaSea2) and 
                        ( train_dict[self.dcomb][0] == self.seaSea1 or 
                          train_dict[self.dcomb][1] == self.seaSea1 or 
                          train_dict[self.dcomb][2] == self.seaSea1 or 
                          train_dict[self.dcomb][-2] == self.seaSea1 or 
                          str(train_dict[self.dcomb][-1]) == self.seaSea1)):
                        self.common_Method()
    
    def yview(self, *args):
            self.comb_names.yview(*args)
            self.row.yview(*args) 
            self.act_syn.yview(*args)
    
    def OnMouseWheel(self, *args):
        if self.comb_names.yview() != self.row.yview():
            self.row.yview_moveto(args[0])
        self.scrollbar.set(*args)
           
    def common_Method(self): 
        if train_dict[self.dcomb][-3] == '1':
            self.comb_names.insert(END,  self.dcomb)
            self.row.insert(END, train_dict[ self.dcomb][-1])
            self.act_syn.insert(END, train_dict[ self.dcomb][-4])
            #print (str(train_dict[self.dcomb][-3]) +":!:"+ str(self.j))
            self.j+=1
        else:
            self.j+=1  
            self.comb_names.insert(END,  self.dcomb)
            self.row.insert(END, train_dict[ self.dcomb][-1])
            self.act_syn.insert(END, train_dict[ self.dcomb][-4])              
            self.comb_names.itemconfig(self.j, bg='red')
            self.row.itemconfig(self.j, bg='red')
            self.act_syn.itemconfig(self.j, bg='red')
            #print (str(train_dict[self.dcomb][-3]) +":-:"+ str(self.j))
                    
root = Tk()
fra = Frame(root, relief=SOLID, borderwidth=1,background='black')
fra1 = Frame(root, relief=RIDGE, borderwidth=2,background='black')
fra2 = Frame(root, relief=RIDGE, borderwidth=3,background='black')
fra.pack(fill=BOTH,expand=1)
fra1.pack(fill=BOTH,expand=1)
fra2.pack(fill=BOTH,expand=1)
datagui = TraDataGUI(root)
root.mainloop()


class LeaDataGUI:
    def __init__(self, root):
        self.root = root
        root.title("DREAM Lea Toronto")
        self.seaSea = StringVar()
        self.inSea = StringVar()
        
        self.scrollbar = Scrollbar(fra, orient="vertical", command=self.yview)
        self.scrollbar.pack(side=RIGHT, fill=Y)
        
        self.sear = Label(fra1, text='Search: "CompoundX" or "CompA.CompB" or ?CompX.CELL_LINE" or "RowNo"').pack(side=TOP)
        self.ent = Entry(fra1,textvariable=self.inSea).pack(side=LEFT,padx=50,pady=10)
        
        self.butS = Button(fra1, text ="Search",
                           command = self.callSearch).pack(side=LEFT,padx=5,pady=5)      
        
        self.col_names = Label(fra, text=" COMPOUND_A.COMPOUND_B.CELL_LINE  | SYNERGY_SCORE | ROW_NO ").pack(side=TOP)
        
        self.comb_names = Listbox(fra, yscrollcommand = self.scrollbar.set, 
                                  selectmode=SINGLE)
        self.row = Listbox(fra, yscrollcommand = self.scrollbar.set, 
                                         selectmode=SINGLE)
        self.act_syn = Listbox(fra, yscrollcommand = self.scrollbar.set, 
                                                 selectmode=SINGLE)        
        self.row.pack(side = RIGHT, fill = Y)
        self.row.configure(width=12, height=30) 
        self.comb_names.pack(side = LEFT, fill = Y)
        self.comb_names.configure(width=35, height=30) 
        self.act_syn.pack(side = RIGHT, fill = Y)
        self.act_syn.configure(width=12, height=30)        
        self.comb_names.bind("<MouseWheel>", self.OnMouseWheel)     
        self.row.bind("<MouseWheel>", self.OnMouseWheel)
        self.act_syn.bind("<MouseWheel>", self.OnMouseWheel)
        root.bind('<Return>', self.callSearch)

        #adding drug combinations to list
        self.j = -1
        for self.dcomb in lea_dict:
            self.common_Method()
                                                      
        self.butA = Button(fra2, text ="Combination",command = self.callA)
        self.butA.pack(side=LEFT,padx=5,pady=10)
        
        self.butB = Button(fra2, text ="Monothereapy",command = self.callB)
        self.butB.pack(side=LEFT,padx=5,pady=10)
        
        self.butC = Button(fra2, text ="TrainingData",command = self.callC)
        self.butC.pack(side=LEFT,padx=5,pady=10)  
        
        self.butD = Button(fra2, text ="TestData",command = self.callD)
        self.butD.pack(side=LEFT,padx=5,pady=10)         
                
    def callA(self):
        if self.comb_names.curselection():
            cur_selection = self.comb_names.get(self.comb_names.curselection()) \
            + ".Rep1.csv"
            webbrowser.open_new(abs_comb + "//" + cur_selection)
        else:
            messagebox.showerror("Error", "Please select a Drug Combination First!")
        
    def callB(self):
        if self.comb_names.curselection():
            cur_selection = self.comb_names.get(self.comb_names.curselection()) \
                + ".Rep1.csv"
            webbrowser.open_new(abs_mono + "//" + cur_selection)
        else:
            messagebox.showerror("Error", "Please select a Drug Combination First!")
    
    def callC(self):
        webbrowser.open_new(path+'/ch1_train_combination_and_monoTherapy.csv')   
    def callD(self):
        webbrowser.open_new(path+'/ch1_test_monoTherapy.csv')   
        
    def callSearch(self, event=None):
        self.comb_names.delete(0,END)
        self.row.delete(0,END)
        self.act_syn.delete(0,END)
        self.seaSea = self.inSea.get()
        self.seaSea = self.seaSea.upper()
        if self.seaSea == "":
            self.j = -1
            for self.dcomb in lea_dict:
                self.common_Method()
        else:
            self.j = -1
            for self.dcomb in lea_dict:
                y = 0
                if self.seaSea[0] == "?":
                    char = self.seaSea.find(".")
                    self.seaSea1 = self.seaSea[1:char]
                    self.seaSea2 = self.seaSea[char+1:]
                    y = 1
                    
                if y == 0:
                    if (lea_dict[self.dcomb][0] == self.seaSea or 
                        lea_dict[self.dcomb][1] == self.seaSea or 
                        lea_dict[self.dcomb][2] == self.seaSea or 
                        lea_dict[self.dcomb][-2] == self.seaSea or 
                        str(lea_dict[self.dcomb][-1]) == self.seaSea):
                        self.common_Method()
                elif y==1:
                    if ((lea_dict[self.dcomb][0] == self.seaSea2 or 
                         lea_dict[self.dcomb][1] == self.seaSea2  or 
                         lea_dict[self.dcomb][2] == self.seaSea2 or 
                         lea_dict[self.dcomb][-2] == self.seaSea2 or 
                         str(lea_dict[self.dcomb][-1]) == self.seaSea2) and 
                        ( lea_dict[self.dcomb][0] == self.seaSea1 or 
                          lea_dict[self.dcomb][1] == self.seaSea1 or 
                          lea_dict[self.dcomb][2] == self.seaSea1 or 
                          lea_dict[self.dcomb][-2] == self.seaSea1 or 
                          str(lea_dict[self.dcomb][-1]) == self.seaSea1)):
                        self.common_Method()
    
    def yview(self, *args):
            self.comb_names.yview(*args)
            self.row.yview(*args) 
            self.act_syn.yview(*args)
    
    def OnMouseWheel(self, *args):
        if self.comb_names.yview() != self.row.yview():
            self.row.yview_moveto(args[0])
        self.scrollbar.set(*args)
           
    def common_Method(self): 
        if lea_dict[self.dcomb][-3] == '1':
            self.comb_names.insert(END,  self.dcomb)
            self.row.insert(END, lea_dict[ self.dcomb][-1])
            self.act_syn.insert(END, lea_dict[ self.dcomb][-4])
            #print (str(lea_dict[self.dcomb][-3]) +":!:"+ str(self.j))
            self.j+=1
        else:
            self.j+=1  
            self.comb_names.insert(END,  self.dcomb)
            self.row.insert(END, lea_dict[ self.dcomb][-1])
            self.act_syn.insert(END, lea_dict[ self.dcomb][-4])              
            self.comb_names.itemconfig(self.j, bg='red')
            self.row.itemconfig(self.j, bg='red')
            self.act_syn.itemconfig(self.j, bg='red')
            #print (str(lea_dict[self.dcomb][-3]) +":-:"+ str(self.j))
                    
root = Tk()
fra = Frame(root, relief=SOLID, borderwidth=1,background='black')
fra1 = Frame(root, relief=RIDGE, borderwidth=2,background='black')
fra2 = Frame(root, relief=RIDGE, borderwidth=3,background='black')
fra.pack(fill=BOTH,expand=1)
fra1.pack(fill=BOTH,expand=1)
fra2.pack(fill=BOTH,expand=1)
datagui = LeaDataGUI(root)
root.mainloop()


