#! /usr/bin/env python

import Tkinter as tk
import tkFileDialog
root = tk.Tk()
root.withdraw()
filename = tkFileDialog.askopenfilename(parent=root,title="Select QE MD output file" )


from pwscf_analyzer import PwscfAnalyzer 


def run_statsandplot(Event):
    begin = tk.IntVar()
    begin=e1.get()
    removeskips=removeSkips.get()
    startfromone=startfromOne.get()
    # default statistical analysis
    #   weak attempt at equilibration detection
    #   automatic estimation of autocorrelation
    #print pa.md_stats
    # equilibration time specified, estimate autocorrelation
    print ("Cutting off {0} steps for equilibration".format(begin))
    print("Mean    StdErrorMean   Correlation Time (steps)")
    print ("begin=",begin)
    print (pa.md_statistics(equil=int(begin)))
    # equilibration and autocorrelation time specified
    #print pa.md_statistics(equil=20,autocorr=45)
    #end if
    pa.md_plots(filename,int(begin),removeskips=removeskips,startfromone=startfromone)

    
pa = PwscfAnalyzer(filename,md_only=True,analyze=True)
md = pa.md_data
# widget for cutting off first equil steps
eqwin=tk.Toplevel()
eqwin.title=("analyze_pwscf_md")
begin = tk.IntVar()
begin = 0
tk.Label(eqwin , text="Equilibration cutoff").grid(row=0)
e1 = tk.Entry(eqwin)
e1.grid(row=0, column=1)
e1.insert(0,0)
removeSkips = tk.IntVar()
startfromOne = tk.IntVar()
tk.Checkbutton(eqwin, text="Remove Skips", variable=removeSkips).grid(row=1,column=0)
tk.Checkbutton(eqwin, text="Start from 0  ", variable=startfromOne).grid(row=2,column=0)
tk.Button(eqwin,text="Quit",command=quit).grid(row=3,column=0)
eqwin.focus_set()
e1.bind('<Return>',run_statsandplot)
eqwin.mainloop()
