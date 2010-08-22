#!BPY
# coding: latin-1

"""
Name: 'Optimal Cutting 3D'
Blender: 100
Group: 'Misc'
Tooltip: 'Cut a volume to reduce its genus'
Author: 'Quentin Fortier'
"""

import Blender
import off_import
#import off_export
import BPyMessages
import Blender.Window
import bpy
import threading
import string 

from Blender.BGL import *
from Blender.Draw import *
from Blender.Window import *

import subprocess
import os
import time

PATH_TETGEN = "C:\\tetgen.exe" 
PATH_OC3D = "C:\\oc3d_debug.exe"
PATH_DEFAULT = "C:\\Program Files\\Blender Foundation\\Mesh"
EXT_CUT = "_cut_"
EXT_DUAL = "_dual.off"
EXT_THIN = "_thin"

Num_Cut = Create(0)
Neighbors_Toggle = Create(1)
Continue_Toggle = Create(1)
Details_Toggle = Create(1)

last_cmd = []
numToCut = []	
tet_file = ""
base_file = ""
oc3d = 0
dual_made = False

CMD_PASS = 0
CMD_MAKE_DUAL = CMD_PASS + 1
CMD_OPT = CMD_MAKE_DUAL + 1
CMD_SAVE_CUT = CMD_OPT + 1
CMD_SAVE_ALL_CUTS = CMD_SAVE_CUT + 1
CMD_LOAD_THICK = CMD_SAVE_CUT + 1
CMD_LOAD_TET = CMD_LOAD_THICK + 1

# Events
EVENT_NOEVENT = 0 

EVENT_LOAD_MESH = EVENT_NOEVENT + 1
EVENT_SELECTED_MESH = EVENT_LOAD_MESH + 1
EVENT_MAKE_DUAL = EVENT_SELECTED_MESH + 1

EVENT_NUM_CUT = EVENT_MAKE_DUAL + 1
EVENT_LOAD_CUT = EVENT_NUM_CUT + 1
EVENT_SAVE_CUT = EVENT_LOAD_CUT + 1
EVENT_SAVE_THICK_CUT = EVENT_SAVE_CUT + 1

EVENT_OPTIMIZE_CUT = EVENT_SAVE_THICK_CUT + 1
EVENT_UNDO = EVENT_OPTIMIZE_CUT + 1
EVENT_DETAILS = EVENT_UNDO + 1
EVENT_OPTIMIZE_ALL = EVENT_DETAILS + 1

EVENT_CONTINUE = EVENT_OPTIMIZE_ALL + 1
EVENT_NEIGHBORHOOD = EVENT_CONTINUE + 1
EVENT_LOAD_ALL = EVENT_NEIGHBORHOOD + 1
EVENT_EXIT = EVENT_LOAD_ALL + 1

def file_cut(num):
	return base_file + EXT_CUT + repr(num) + ".cut"

def clean(file):
	try:
		os.remove(file)
	except:
		#print "Warning: can't remove old file " + file
		return False
	return True

def read_all(filename):
	global base_file, numToCut, oc3d, last_cmd, dual_made
	split_file = os.path.splitext(filename)
	base_file = split_file[0]
	set_layer(1)
	read(base_file + ".off")
	tet_file = base_file + ".1"
	set_layer(2)
	read(tet_file + ".off")
	dual_made = False
	oc3d.stdin.write("load " + tet_file + ".ele" + "\n")
	last_cmd.append((CMD_PASS, 0))
	oc3d.stdin.write("make_dual\n")
	last_cmd.append((CMD_MAKE_DUAL, base_file + EXT_DUAL))
	while(dual_made == False): time.sleep(0.1)
	i = 0
	while read_cut(i, file_cut(i)):
		oc3d.stdin.write("load_cut " + str(i) + "\n")
		last_cmd.append((CMD_PASS, 0))
		i = i+1
	if(i > Num_Cut.val):
		select(Num_Cut.val)
		
def write_off(filename):
	file = open(filename, 'wb')
	scn= Blender.Scene.GetCurrent()
	object= scn.objects.active
	
	if not object or object.type != 'Mesh':
		BPyMessages.Error_NoMeshActive()
		return
	
	Blender.Window.WaitCursor(1)
	mesh = object.getData(mesh=1)

	file.write('OFF\n')
	file.write('%d %d %d\n' % (len(mesh.verts), len(mesh.faces), len(mesh.edges)))

	for i, v in enumerate(mesh.verts):
		file.write('%.6f %.6f %.6f\n' % tuple(v.co))

	for i, f in enumerate(mesh.faces):
		file.write('%i' % len(f))
		for v in reversed(f.v):
			file.write(' %d' % v.index)
		file.write('\n')
	
	for i, e in enumerate(mesh.edges):
		file.write("2 %d %d\n" % (e.v1.index, e.v2.index))
		
	file.close()
	Blender.Window.WaitCursor(0)
	message = 'Successfully exported "%s"' % Blender.sys.basename(filename)
	
def read_cut(num, filename):
    global numToCut
    try:
        file = open(filename, "rb")
    except:
        print (str("Error: cannot read cut file ") + filename)
        return False
    lEdges = file.readline()
    E = map(int, lEdges.split())
    while(num >= len(numToCut)):
        numToCut.append([])
    numToCut[num] = []
    mesh = bpy.data.scenes.active.objects.active.getData(mesh=1)
    for l in file:
        try:
            u, v = map(int, l.split())
            e = mesh.findEdges(u, v)
            numToCut[num].append(mesh.edges[e])
        except:
            print ("Error: can't read edge (" + str(u) + ", " + str(v) + ")")
    file.close()
    print ("Read cut " + str(num) + " with " + str(len(numToCut[num])) + "edges in file " +  filename)
    return True
	
def write_cut(num, filename):
	global numToCut
	file = open(filename, 'wb')
	file.write(str(len(numToCut[num])) + "\n")
	for e in numToCut[num]: file.write(str(e.v1.index) + " " + str(e.v2.index) + "\n")
	file.close()
	
def save_cut(num):
	global numToCut
	if num == len(numToCut):
		numToCut.append([])
	in_editmode = Blender.Window.EditMode()
	if in_editmode: Blender.Window.EditMode(0)
	mesh = bpy.data.scenes.active.objects.active.getData(mesh=1)
	numToCut[num] = [e for e in mesh.edges if e.sel]
	if in_editmode: Blender.Window.EditMode(1)

def select(num):	
	global numToCut
	try:
		Blender.Window.EditMode(0)
		unselect()
		for e in numToCut[num]:
			e.sel = 1
		Blender.Window.EditMode(1)
	except:
		print ("Error: Can't select cut %d" % num)
			
def unselect():
	try:
		in_editmode = Blender.Window.EditMode()
		if in_editmode: Blender.Window.EditMode(0)
		mesh = bpy.data.scenes.active.objects.active.getData(mesh=1)
		mesh.sel = False
		mesh.update()
		if in_editmode: Blender.Window.EditMode(1)
	except:
		print ("Warning: can't unselect")
	
def write(file):
	try:
		write_off(file)
	except IOError:
		print ("Error: can't write " + file)
		return False
	return True

def read(file):
	try:
		off_import.read(file)
	except:
		print ("Error: can't read " + file)
		return False
	print ("Read " + file)
	return True
		
def tetgen(file):
	global tet_file, base_file
	split_file = os.path.splitext(file)
	base_file = split_file[0]
	tet_file = base_file + ".1"
	clean(tet_file + ".off")
	clean(tet_file + ".ele")
	clean(tet_file + ".face")
	clean(tet_file + ".node")
	clean(tet_file + ".smesh")
	Blender.Window.WaitCursor(1)
	try:
		p = subprocess.Popen(args=[PATH_TETGEN,"-Oq",file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=False)
	except:
		print ("Error: can't execute " + PATH_TETGEN)
		return False
	outputlines = p.stdout.readlines()
	p.wait()
	for l in outputlines:
		print (l, )
	set_layer(2)
	read(tet_file + ".off")
	Blender.Window.WaitCursor(0)
	return True

def read_output():
	import string
	global oc3d, last_cmd, dual_made, Num_Cut, numToCut
	while(True):
		#time.sleep(0.1)
		l = oc3d.stdout.readline()
		l_split = str.split(l, "\n")
		nEnd = string.count(l, "END")
		if(nEnd == 0):
			print (l_split[0])
		while(nEnd >= 1):
			nEnd = nEnd - 1
			if(last_cmd):
				cmd, info = last_cmd.pop(0)
				if(cmd == CMD_SAVE_CUT):
					read_cut(info, file_cut(info))
					select(info)
				elif(cmd == CMD_SAVE_ALL_CUTS):
					for i in range(len(numToCut)):
						read_cut(i, file_cut(i))
				elif(cmd == CMD_MAKE_DUAL):
					set_layer(3)
					read(info)
					dual_made = True
				elif(cmd == CMD_LOAD_THICK):
					read_cut(info, base_file + EXT_CUT + repr(Num_Cut.val) + EXT_THIN +  ".cut")
					select(info)
				elif(cmd == CMD_PASS):
					pass
				elif(cmd == CMD_LOAD_TET):
					oc3d.stdin.write("make_dual\n")
					last_cmd.append((CMD_MAKE_DUAL, base_file + EXT_DUAL))
					
			else: print ("Error: can't find last command")

def set_layer(num):
	ViewLayers([num])
	SetActiveLayer(1<<(num-1))

def gui():
	global Num_Cut, Neighbors_Toggle, Continue_Toggle, Details_Toggle, oc3d

	glClear(GL_COLOR_BUFFER_BIT)
	
	area = GetAreaSize()
	bottom = area[1] / 20
	left = area[0] / 20
	width = area[0] / 8
	height = min(area[1] / 5, 45)
	col_separation = int(left + 1.5*width)
	# TetGen
	Button("Make dual", EVENT_MAKE_DUAL, left, bottom + 2*height, width, height, "Make dual graph")
	Button("TetGen..", EVENT_LOAD_MESH, left, bottom + height, width, height, "Load and tetgen a mesh")
	Button("TetGen selected", EVENT_SELECTED_MESH, left, bottom, width, height, "Export as .off and tetgen it")
		
	# Cuts
	
	def num_changed(event, num):
		unselect()
		if (num < len(numToCut)): select(num)
	
	Num_Cut = Number("Cut n°", EVENT_NOEVENT, col_separation, bottom, width, height, Num_Cut.val, 0, len(numToCut), "Select the i-th cut", num_changed, 1)
	#Button("Save thick cut", EVENT_SAVE_THICK_CUT, col_separation, bottom + 3*height, width, height, "Save selected thick cut as " + file_cut(Num_Cut.val))
	Button("Save cut", EVENT_SAVE_CUT, col_separation, bottom + 2*height, width, height, "Save selected cut as " + file_cut(Num_Cut.val))
	Button("Load cut", EVENT_LOAD_CUT, col_separation, bottom + height, width, height, "Load " + file_cut(Num_Cut.val))

	# Optimize
	# Button("Undo / redo", EVENT_UNDO, col_separation*2, bottom + 2*height, width, height, "Undo / redo last optimization")
	Details_Toggle = Toggle("Optimization details", EVENT_DETAILS, col_separation*2, bottom + 2*height, width, height, Details_Toggle.val, "Gives more information in the console")
	Button("Optimize cut", EVENT_OPTIMIZE_CUT, col_separation*2, bottom + height, width, height, "Optimize cut " + str(Num_Cut.val))
	Button("Optimize all", EVENT_OPTIMIZE_ALL, col_separation*2, bottom, width, height, "Optimize cut " + str(Num_Cut.val))
	
	# Options
	#Button("Save all", EVENT_SAVE_CONFIG, col_separation*3, bottom + 2*height, width, height)
	Continue_Toggle = Toggle("Step", EVENT_CONTINUE, col_separation*3, bottom + 3*height, width, height, Continue_Toggle.val, str("Enable the variant of neighborhood algorithm searching all possible paths before augmentating "))
	Neighbors_Toggle = Toggle("Neighborhood", EVENT_NEIGHBORHOOD, col_separation*3, bottom + 2*height, width, height, Neighbors_Toggle.val, "Enable the use of neighborhood during max flow")
	Button("Load all", EVENT_LOAD_ALL, col_separation*3, bottom + height, width, height, "Select a .off file and load all related files")
	Button("Exit",EVENT_EXIT , col_separation*3, bottom, width, height)
		
def event(evt, val):
	global oc3d
	if evt == QKEY and not val:
		oc3d.kill()
		Exit()
		
def load_tet(file):
	global base_file
	split_file = os.path.splitext(file)
	base_file = (os.path.splitext(split_file[0]))[0]
	oc3d.stdin.write("load " + file + "\n")
	last_cmd.append((CMD_PASS, 0))
	oc3d.stdin.write("make_dual\n")
	last_cmd.append((CMD_MAKE_DUAL, base_file + EXT_DUAL))
	
def bevent(evt):
	global numToCut, oc3d, tet_file, base_file, dual_made
	
	if evt == EVENT_EXIT:
		if(oc3d):
			oc3d.kill()
		Exit()
	elif evt == EVENT_LOAD_MESH:
		FileSelector(tetgen, "Load a .off mesh", PATH_DEFAULT)
		
	elif evt == EVENT_SELECTED_MESH:
		obj = Blender.Object.GetSelected()
		if(obj):
			def save(file):
				split_file = os.path.splitext(file)
				clean(split_file[0] + ".off")
				if write(split_file[0] + ".off"):
					tetgen(split_file[0] + ".off")
			FileSelector(save, "Save as .off mesh", PATH_DEFAULT + "\\" + bpy.data.scenes.active.objects.active.getName() + ".off")
			
	elif evt == EVENT_MAKE_DUAL:
		dual_made = False	
		oc3d.stdin.write("load " + tet_file + ".ele" + "\n")
		last_cmd.append((CMD_PASS, 0))
		oc3d.stdin.write("make_dual\n")
		last_cmd.append((CMD_MAKE_DUAL, base_file + EXT_DUAL))
		while(dual_made == False): time.sleep(0.1)
			
	elif evt == EVENT_LOAD_CUT:
		read_cut(Num_Cut.val, file_cut(Num_Cut.val))
		select(Num_Cut.val)
		oc3d.stdin.write("load_cut " + str(Num_Cut.val) + "\n")
		last_cmd.append((CMD_PASS, 0))
		
	elif evt == EVENT_SAVE_CUT:
		save_cut(Num_Cut.val)
		write_cut(Num_Cut.val, file_cut(Num_Cut.val))
		oc3d.stdin.write("load_cut " + str(Num_Cut.val) + "\n")
		last_cmd.append((CMD_PASS, 0))
	
	elif evt == EVENT_SAVE_THICK_CUT:
		save_cut(Num_Cut.val)
		write_cut(Num_Cut.val, file_cut(Num_Cut.val))
		oc3d.stdin.write("load_thick_cut " + str(Num_Cut.val) + "\n")
		last_cmd.append((CMD_LOAD_THICK, Num_Cut.val))
		
	elif evt == EVENT_OPTIMIZE_CUT:
		if(Num_Cut.val >= len(numToCut)):
			print ("Error: cut " + str(Num_Cut.val) + " not saved")
			return
		oc3d.stdin.write("init\n")
		last_cmd.append((CMD_PASS, 0))
		oc3d.stdin.write("opt " + str(Num_Cut.val) + "\n")
		last_cmd.append((CMD_PASS, 0))
		oc3d.stdin.write("save_cut " + str(Num_Cut.val) + "\n")
		last_cmd.append((CMD_SAVE_CUT, Num_Cut.val))

	elif evt == EVENT_OPTIMIZE_ALL:
		oc3d.stdin.write("init\n")
		last_cmd.append((CMD_PASS, 0))
		oc3d.stdin.write("opt\n")
		last_cmd.append((CMD_PASS, 0))
		oc3d.stdin.write("save_cut" + "\n")
		last_cmd.append((CMD_SAVE_ALL_CUTS, 0))
		
	elif evt == EVENT_LOAD_ALL:
		FileSelector(read_all, "Load a base .off mesh", PATH_DEFAULT + ".off")
	
	elif evt == EVENT_NEIGHBORHOOD:
		oc3d.stdin.write("neighbors\n")
		last_cmd.append((CMD_PASS, 0))
	
	elif evt == EVENT_CONTINUE:
		oc3d.stdin.write("step\n")
		last_cmd.append((CMD_PASS, 0))
	
	elif evt == EVENT_DETAILS:
		oc3d.stdin.write("details\n")
		last_cmd.append((CMD_PASS, 0))
		
oc3d = subprocess.Popen(args=[PATH_OC3D], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin = subprocess.PIPE, shell=False)	
a = threading.Thread(target=read_output)
a.start()

Register(gui, event, bevent)