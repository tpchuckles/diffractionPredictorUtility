import tkinter as tk
import matplotlib,os,ast
from tkinter import *
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
matplotlib.use("Agg") # FOR SOME REASON, WE HAD "TkAgg" BEFORE, AND IT SUDDENLY DECIDED TO START OPENING NEW WINDOWS FOR EACH PLOT? comparing to gui.py > TDTR_fitting.py > niceplot.py, which works correctly (updated plot opens within window), gui.py specifies TkAgg, BUT, niceplot specifies Agg
import matplotlib.pyplot as plt

#matplotlib.rc("svg", **{'fonttype':'none'}) # UNCOMMENT THIS IF SAVING SVG, SO TEXT IS EDITABLE IN INKSCAPE
#matplotlib.rc("font",**{'size':16})

window = tk.Tk()

figsize=(6,6)

livePlot={} # global dist to store active fig elements
# TRIANGLE 110 111 100 PLOT. USER CLICKS ON THE PLOT TO CHOOSE CRYSTAL ORIENTATION. located upper-left of 4-panel GUI
def genFig(x,y):
	fig,ax=plt.subplots(figsize=figsize)			# new empty matplotlib figure/axes objects
	ax.fill_between([0,1], [0,1], [0,0])			# triangle
	ax.plot([x],[y],marker='o',color='r')			# and red dot
	# for some crystals, 100 is the same as 010, for others (e.g. wurtzite) it's not. 
	directions={"ij":("100","110","111"),"ik":("001","101","111"),"ji":("010","110","111"),"jk":("010","011","111"),"ki":("001","101","111"),"kj":("001","011","111")}[ordering]
	xys=[(0,0),(1,0),(1,1)]
	for direc,xy in zip(directions,xys):
		ax.annotate("["+",".join([c for c in direc])+"]",xy=xy,xytext=xy)
	livePlot["fig"],livePlot["ax"]=fig,ax			# store the fig,ax objects to our dict
	if "canvas" in livePlot.keys():				# if this isn't the first time, destroy the old "canvas" object
		livePlot["canvas"].get_tk_widget().destroy()
	canvas = FigureCanvasTkAgg(fig, master=window)		# make a new canvas object to hold the figure
	canvas.get_tk_widget().grid(row=0,column=0,sticky='ew')	# pack it
	livePlot["canvas"]=canvas				# and store off the canvas object too
	canvas.mpl_connect('button_press_event',callback)	# link callback function to clicks

lattices={} # lattices are defined by atom positions, lattice constant(s), atomic masses.
def readStructureFile(filename):
	global lattices
	lines=open(filename).readlines()
	for l in lines:			# input file format, each line should have: # structureName: { "pos":[...],"abc":[...],"atomic_number":[...] }
		l=l.replace(" ","")
		k=l.split(":")[0]
		l=l.replace(k+":","")
		k=k.replace("\"","")
		lattices[k]=ast.literal_eval(l) # a safer way to convert string to dict: BEWARE: division doesn't seem to work "1/2" won't work ".5" will 
		print(k)
if os.path.exists("structures/defaults.txt"):
	readStructureFile("structures/defaults.txt")
else:
	lattices={ "diamond" : { "pos":[[0,0,0],[.25,.25,.25],[.5,.5,0],[.5,0,.5],[0,.5,.5],[.75,.25,.75],[.25,.75,.75],[.75,.75,.25]] , "abc":5 },
		"FCC" : { "pos":[[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]] , "abc":5 },
		"cubic" : { "pos":[[0,0,0]] , "abc":5 }, #} 
		"hexagonal" : { "pos":[[1/3,2/3,0],[2/3,1/3,1/2],[1/3,2/3,3/8],[2/3,1/3,7/8]] , "abc":[3.19,3.19,5.19,90,90,120] , "atomic_number":[14,14,28,28] },
		"BCC" : { "pos":[[0,0,0],[.5,.5,.5]],"abc":5} }

def loadMoreStructures(event):
	# query user for file
	files=tk.filedialog.askopenfilenames(title="select file with structures dict. see structures/defaults.txt as an example")
	# load them into lattices dict
	for f in files:
		readStructureFile(f)
	# regenerate dropdown menu
	global dropDown
	dropDown.pack_forget()
	dropDown = OptionMenu( frame , structureDropdownValue , *list(lattices.keys()) )
	dropDown.pack(after=structureLabel)

# contains various user-enterable options (crystal structure drop-down, etc). located upper-right of 4-panel GUI
def paramsPanel():
	global structureDropdownValue,kmaxField,ijButton,dropDown,frame,structureLabel
	frame = tk.Frame(window)
	frame.grid(row=0,column=1,sticky='ew')

	structureLabel = tk.Label(frame,text="structure") 
	structureLabel.pack()
	structureDropdownValue=tk.StringVar(window) ; structureDropdownValue.set("hexagonal") # tk stringvar holds value of dropdown menu, so other funcs can retreive it
	dropDown = OptionMenu( frame , structureDropdownValue , *list(lattices.keys()) )
	dropDown.pack()

	label = tk.Label(frame,text="k_max") 
	label.pack()
	kmaxField=tk.Entry(frame)
	kmaxField.insert(0,"2.0")
	kmaxField.pack()

	buttonFrame=tk.Frame(frame)
	buttonFrame.pack()
	
	global ijkFuncs ; ijkFuncs=[]
	for text in orderings:
		def makeButtonFunc(text): # why nest? this tricks python(?) into giving us a copy of each func in the list! 
			def buttonFunc(event):
				print(text)
				global ordering ; ordering=text
				for button in ijkButtons:
					button.config(background="grey")
				i=orderings.index(text)
				ijkButtons[i].config(background="white")
			return buttonFunc
		ijkFuncs.append(makeButtonFunc(text))

	global ijkButtons ; ijkButtons=[]
	for i,text in enumerate(orderings):
		button=tk.Button(buttonFrame,text=text)
		button.bind("<Button-1>",ijkFuncs[i])
		button.grid(row=0,column=i)
		ijkButtons.append(button)
	
	global ordering
	ordering="ij" ; ijkButtons[0].config(background="white")

	structuresButton=tk.Button(frame,text="load more structures")
	structuresButton.bind("<Button-1>",loadMoreStructures)
	structuresButton.pack()

orderings=["ij","ik","ji","jk","ki","kj"] 	# these handle non-symmetric crystal structures. 100 is not the same as 001 in
ordering="ij"					# wurtzite, so use ij to view a=1,b=0,c=0 or kj indices to view a=0,b=0,c=1 etc

import py4DSTEM
settings={"structure":"","atomic_number":28,"k_max":0}

# we use py4DSTEM to preview real-space 3D atomic structure (showStructure) and generate diffraction pattern (showDiffraction). this sets up parameters for both 
def regenerate():
	global crystal

	regen=False

	structure=structureDropdownValue.get()
	if structure!=settings["structure"]:
		regen=True
		settings["structure"]=structure
	k_max=float(kmaxField.get())
	if k_max!=settings["k_max"]:
		regen=True
		settings["k_max"]=k_max

	if not regen:
		return

	print("regenerating")

	lattice=lattices[settings["structure"]]
	points=lattice["pos"] ; abc=lattice["abc"]

	#atom_num=28 ; a=5.4

	# generate crystal, and preview it
	atom_num=lattice.get("atomic_number",28)
	crystal = py4DSTEM.process.diffraction.Crystal(points, atom_num, abc)
	# crystal.plot_structure(zone_axis_lattice=[3,2,1])

	# maximum scattering vector. need to look this up
	#k_max = 2.0
	# generate structure factors, and plot
	crystal.calculate_structure_factors(settings["k_max"])
	#crystal.plot_structure_factors(zone_axis_lattice=[3,2,1])

	# accelerating voltage
	crystal.setup_diffraction(300e3)

# Example of generating a diffraction pattern for a given orientation
#zone_axis_test=(0,0,1) # result also depends on k_max
#bragg_peaks = crystal.generate_diffraction_pattern(zone_axis_lattice = zone_axis_test,sigma_excitation_error=0.02)
#py4DSTEM.process.diffraction.plot_diffraction_pattern(bragg_peaks, add_labels=False)

# Create and preview an "orientation plan" (allows easy rotation)
# crystal.orientation_plan()
#crystal.plot_orientation_zones()

# upon clicking on the 100-110-111 plot, we generate an oriented 3D atomic structure located in lower-left of 4-panel GUI
# (note: turning the 3D structure doesn't regenerate the diffraction pattern. re-click the 100-110-111 plot)
def showStructure(i,j,k):
	fig,ax=crystal.plot_structure(zone_axis_lattice=[i,j,k],returnfig=True,figsize=figsize)
	livePlot["fig2"],livePlot["ax2"]=fig,ax			# store the fig,ax objects to our dict
	if "canvas2" in livePlot.keys():			# if this isn't the first time, destroy the old "canvas" object
		livePlot["canvas2"].get_tk_widget().destroy()
	canvas = FigureCanvasTkAgg(fig, master=window)		# make a new canvas object to hold the figure
	canvas.get_tk_widget().grid(row=1,column=0,sticky='ew')	# pack it
	livePlot["canvas2"]=canvas				# and store off the canvas object too

# upon clicking on the 100-110-111 plot, we generate a diffraction pattern for the given structure/orientation. located in lower-right of 4-panel GUI
def showDiffraction(i,j,k):
	bragg_peaks = crystal.generate_diffraction_pattern(zone_axis_lattice = (i,j,k),sigma_excitation_error=0.02)
	fig,ax=py4DSTEM.process.diffraction.plot_diffraction_pattern(bragg_peaks, add_labels=True,returnfig=True,figsize=figsize) # EDIT /home/athena/.local/lib/python3.11/site-packages/py4DSTEM/process/diffraction/crystal_viz.py > add_labels == True font size IF YOU WANT FONT SIZE CHANGED FOR DOT LABELS
	livePlot["fig3"],livePlot["ax3"]=fig,ax			# store the fig,ax objects to our dict
	if "canvas3" in livePlot.keys():			# if this isn't the first time, destroy the old "canvas" object
		livePlot["canvas3"].get_tk_widget().destroy()
	canvas = FigureCanvasTkAgg(fig, master=window)		# make a new canvas object to hold the figure
	canvas.get_tk_widget().grid(row=1,column=1,sticky='ew')	# pack it
	livePlot["canvas3"]=canvas				# and store off the canvas object too
	#plt.savefig("diffractions.svg") # UNCOMMENT THIS TO SAVE DIFFRACTION PATTERN

# This is how we read the crystal orientation from clicks on the 100-110-111 plot.
def callback(event): # (https://stackoverflow.com/questions/27565939/getting-the-location-of-a-mouse-click-in-matplotlib-using-tkinter)
	x,y=event.xdata,event.ydata				# when the user clicks on the plot...
	print("clicked at", x,y)
	if x>1:							# detect out-of-bounds clicks
		x=1
	if x<0:
		x=0
	if y<0:
		y=0
	if y>x:
		y=x
	regenerate()						# update globals for py4DSTEM 3D atom visual and diffraction pattern
	genFig(x,y)						# update the 100-110-111 figure
	i,j,k=1,x,y

	i,j,k={"ij":(i,j,k),"ik":(i,k,j),"ji":(j,i,k),"jk":(k,i,j),"ki":(j,k,i),"kj":(k,j,i)}[ordering]

	print(i,j,k)

	showStructure(i,j,k)					# update the 3D atom structure visualizer
	showDiffraction(i,j,k)					# update the diff pattern
	window.update()						# and update the window


def quit_me():				# weird thing, when we generate a matplotlib plot, tkinter main loop doesn't exit when we click the x.
	#log("Quitting")		# to deal with that, we detect a "delete window" and then use that to quit.
	window.quit()			# https://stackoverflow.com/questions/55201199/the-python-program-is-not-ending-when-tkinter-window-is-closed
	window.destroy()
	#global done
	#done=True
genFig(0,0)				# to initialize, start with just the 100-110-111 fig, and the user options panel
paramsPanel()
window.protocol("WM_DELETE_WINDOW", quit_me)

window.mainloop()




