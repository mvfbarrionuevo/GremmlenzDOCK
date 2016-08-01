#!/usr/bin/python
#################################################################
# GremmlenzDOCK.py - created by Manoel Barrionuevo - 2016       #
#################################################################
import sys, time, os, shutil
#################################################################

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("Invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            print " "
            print "You have hit [ENTER] thus I'm going to assume you haven't a glue about what is going on. Bye."
            print " "
	    quit()
        elif choice in valid:
            if choice == 'no' or choice == 'n':
		print " "
		print "So you have chosen '%s', thus you may like to take my advise first." % choice
		print " "
		quit()
	    elif choice == 'yes' or choice == 'y':
		print " "
		print "You have chosen '%s'." %choice		
		print " "
		return False 
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

#################################################################
#			SETTING VARIABLES			#
# These variables must to be changed accordingly. Please, before#
#you run this script, make sure you have set up all paths.      #
# This script was firstly made by taking it inside the receptors#
#folder, thus it explains why the variables are set as it is. If#
#you decide to keep unchanged, make sure you have receptors and #
#ligands folder inside MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDock#
#Tools/Utilities24 path, otherwise your script is going to fail.#
#################################################################

#receptor_foler ='./' #path to receptor proteins inside MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/
#ligand_folder='../ligands/' #path to ligands to be docked inside MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24
#utilities='../' #path to MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24
#pyshell='../../../../bin/' #path to MGLTools-1.5.7rc1/bin/ took from MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/receptors/

#################################################################
#   VERIFICATION STEP - ARE ALL SCRIPTS WHERE THEY SHOULD BE?   #
#################################################################

print " "
print "Hi, first thing first. If you didn't get (or didn't known before) the scripts for autodockZN, MGLTools and GremmlenzDOCK you will be able to download them from here:"
print " "
print " 1) http://autodock.scripps.edu/resources/autodockzn-forcefield"
print " 2) http://mgltools.scripps.edu/downloads/mgltools-1-5-7rc1"
print " 3) http://GremmlenzDOCKgithubhere "
print " "
print "Once you get them, please uncompress and install MGLTools, then move all autodockZN and GremmlenzDOCK scripts to the path /path/to/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/ (wherever this path may look like).'

print " "
query_yes_no("So, have you done everything I told you?")
print " "
print "Hnmm, it seems you know what you are doing. Lets take some addresses in order to deal with everything. Please, carefully write down the addresses I will ask you. I will not take any responsability about your typos. In case you agree with the default paths just hit [ENTER]."
print " "
utilities=raw_input("[~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/] What's the Utilities24 path?")
if(not utilities or utilities.isspace()):
    utilities='~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'
print " "
receptor_folder=raw_input("[~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/receptors] What's the path of your receptors?")
if(not receptor_folder or receptor_folder.isspace()):
    receptor_folder='~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/receptors/'
print " "
ligand_folder=raw_input("[~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/ligands/] What's the path of you ligands?")
if(not ligand_folder or ligand_folder.isspace()):
    ligand_folder='~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/ligands/'
print " "
pyshell=raw_input("[~/MGLTools-1.5.7rc1/bin/] What's the pythonsh path?")
if(not pyshell or pyshell.isspace()):
    pyshell='~/MGLTools-1.5.7rc1/bin/'
print " "
print "Now that every path has been set up, I'm going to verify through a very simple way if you really got all scripts where they should be (I'll know if you read me before saying 'yes' to me). If something is missing I'll quit the execution, otherwise your docking will not run smoothly."
print " "

os.chdir(utilities)
#now we're inside ~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24
if receptor_folder=='~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/receptors/':
	os.system('mkdir ./receptors/')
	os.system('cp Lister* '+receptor_folder)
	print " "
	print "Please, take a moment and manually copy all your receptors to the path: "+receptor_folder+"."
	print " "
	os.system('read -s -n 1 -p "I am waiting for you, please, press [ENTER] after you have copied all receptor files..."')
	print
else:
	os.system('cp Lister* '+receptor_folder)

if ligand_folder=='~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/ligands/':
	os.system('mkdir ./ligands/')
	os.system('cp Lister-stp2.py '+ligand_folder)
	os.system('cp Lister-stp3.py '+ligand_folder)
	print " "
	print "Please, take a moment and manually copy all your ligands to the path: "+ligand_folder+"."
	print " "
	os.system('read -s -n 1 -p "I am waiting for you, please, press [ENTER] after you have copied all ligand files..."')
	print
else:
	os.system('cp Lister-stp2.py '+ligand_folder)
	os.system('cp Lister-stp3.py '+ligand_folder)

if os.path.exists('./prepare_gpf4zn.py')==True:
	print "prepare_gpf4zn.py has been found."
else:
	print "Sorry, prepare_gpf4zn.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'."
	print " "
	quit()
if os.path.exists('./prepare_dpf42.py')==True:
	print "prepare_dpf42.py has been found."
else:
	print "Sorry, prepare_dpf42.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'."
	print " "
	quit()
if os.path.exists('./prepare_receptor4.py')==True:
	print "prepare_receptor4.py has been found."
else:
	print "Sorry, prepare_receptor4.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'."
	print " "
	quit()
if os.path.exists('./prepare_ligand4.py')==True:
		print "prepare_gpf4zn.py has been found."
else:
	print "Sorry, prepare_ligand4.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'."
	print " "
	quit()
if os.path.exists('./zinc_pseudo.py')==True:
	print "zinc_pseudo.py has been found."
else:
	print "Sorry, zinc_pseudo.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'."
	print " "
	quit()
if os.path.exists('./autodock4')==True:
	print "autodock4.py has been found."
else:
	print "Sorry, autodock4 isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'."
	print " "
	quit()
if os.path.exists('./autogrid4.2.5.x.20131125')==True:
	print "autogrid4.2.5.x.20131125.py has been found."
else:
	print "Sorry, autogrid4.2.5.x.20131125 isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'."
	print " "
	quit()
if os.path.exists('./AD4Zn.dat')==True:
	print "AD4Zn.dat has been found."
	print "Copying AD4Zn.dat to working directory: "+receptor_folder+"."
	ad4=utilities+'AD4Zn.dat'	
	os.system('cp '+ad4+' '+receptor_folder)
else:
	print "Sorry, AD4Zn.dat isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'."
	print " "
	quit()

os.chdir(pyshell)
#now we're inside ~/MGLTools-1.5.7rc1/bin/
if os.path.exists('./pythonsh')==True:
	print "pythonsh has been found."
else:
	print "Sorry, pythonsh isn't in 'MGLTools-1.5.7rc1/bin/'."
	print " "
	quit()

os.chdir(ligand_folder)
#now we're insde ~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/ligands/
if os.path.exists('./Lister-stp2.py')==True:
	print "Lister-stp2.py has been found."
else:
	print "Sorry, Lister-stp2.py isn't in '"+ligand_folder+"'."
	print " "
	quit()
if os.path.exists('./Lister-stp3.py')==True:
	print "Lister-stp2.py has been found."
else:
	print "Sorry, Lister-stp3.py isn't in '"+ligand_folder+"'."
	print " "
	quit()

os.chdir(receptor_folder)
#now we're inside ~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/receptors/
if os.path.exists('./Lister-stp2.py')==True:
	print "Lister-stp2.py has been found."
else:
	print "Sorry, Lister-stp2.py isn't in '"+receptor_folder+"'."
	print " "
	quit()
if os.path.exists('./Lister-stp3.py')==True:
	print "Lister-stp3.py has been found."
else:
	print "Sorry, Lister-stp3.py isn't in '"+receptor_folder+"'."
	print " "
	quit()
if os.path.exists('./Lister-stp4.py')==True:
	print "Lister-stp4.py has been found."
else:
	print "Sorry, Lister-stp4.py isn't in '"+receptor_folder+"'."
	print " "
	quit()

#################################################################
#     PREPARATION STEP - CREATING RECEPTOR AND LIGAND PDBQT     #
#################################################################

#################################################################
#			   RECEPTOR				#
#################################################################

# Creating a list

lines=[]

# Calling Lister to get a list of pdb files to work with

os.chdir(receptor_folder)

if os.path.exists('./Lister-stp2.py')==True:
	print " "
	print "Hey, I'm going to call a friend of mine to get a list of all pdbs found in receptor directory."
	os.system('python '+receptor_folder+'Lister-stp2.py '+receptor_folder)
	if os.path.exists(receptor_folder+'exit.txt')==True:
	  print "You told to my friend you wouldn't like to get a new list, hope you fix it and come back. See you then."
	  os.remove(receptor_folder+'exit.txt')
	  quit()

else:
	print " "
	print "Sorry, I'm affraid 'Lister-stp2.py' isn't in receptors directory. May you call it to join in here?"
	print " "
	quit()

if os.path.exists(receptor_folder+"listpdb.txt")==True:
	print " "
	print "Yeah, 'Lister-stp2.py' has listed all receptors. Now I'll rock things up!"
	print " "
else:
	print " "
	print "Something still wrong here, I found 'Lister-stp2.py' but he didn't give me back a receptor's pdb list, could you please check if he's working well?"
	print " "
	quit()

with open(receptor_folder+'listpdb.txt') as lst:
	for line in lst:
		lines.append(line.strip())
lst.close()

# Creating pdbqt files

index=0
for index in range(len(lines)):
	os.system(pyshell+'pythonsh '+utilities+'prepare_receptor4.py -r '+lines[index]+'.pdb -o '+lines[index].replace('-zn-f','')+'.pdbqt')

# Ends everything

print " "
print "Great, pdbqt of receptor files are ready! :D"
print " "

del lines=[:]

#################################################################
#			   LIGAND				#
#################################################################

os.chdir(ligand_folder)

# Calling Lister to get a list of pdb files to work with

if os.path.exists(ligand_folder+'Lister-stp2.py')==True:
	print " "
	print "Hey, I'm going to call a friend of mine to get a list of all pdbs found in ligand directory."
	os.system('python Lister-stp2.py '+ligand_folder)
	if os.path.exists(ligand_folder+'exit.txt')==True:
	  print "You told to my friend you wouldn't like to get a new list, hope you fix it and come back. See you then."
	  os.remove(ligand_folder+'exit.txt')
	  quit()

else:
	print " "
	print "Sorry, I'm affraid 'Lister-stp2.py' isn't in ligand directory. May you call it to join in here?"
	print " "
	quit()

if os.path.exists(ligand_folder+"listpdb.txt")==True:
	print " "
	print "Yeah, 'Lister-stp2.py' has listed all ligands. Now I'll rock things up!"
	print " "
else:
	print " "
	print "Something still wrong here, I found 'Lister-stp2.py' but he didn't give me back a ligand's pdb list, could you please check if he's working well?"
	print " "
	quit()

with open(ligand_folder+'listpdb.txt') as lst:
	for line in lst:
		lines.append(line.strip())
lst.close()

# Creating pdbqt files

index=0
for index in range(len(lines)):
	os.system(pyshell+'pythonsh '+utilities+'prepare_receptor4.py -r '+lines[index]+'.pdb -o '+lines[index].replace('-zn-f','')+'.pdbqt')

# Ends everything

print " "
print "Great, pdbqt of ligand files are ready! :D"
print " "

del lines=[:]

#################################################################
#      PREPARATION STEP - CREATING ZINC PSEUDO ATOM CENTRES     #
#################################################################

os.chdir(receptor_folder)

# Calling Lister to get a list of pdbqt files to work with

with open(receptor_folder+'listpdb.txt') as lst:
	for line in lst:
		lines.append(line.strip())
lst.close()

# Creating pdbqt files

index=0
for index in range(len(lines)):
	os.system(pyshell+'pythonsh '+utilities+'zinc_pseudo.py -r '+lines[index].replace('-zn-f','')+'.pdbqt -o '+lines[index]+'_tz.pdbqt')

# Ends everything

print " "
print "Great, pseudo zinc atoms were added to receptor pdbqt files! :D"
print " "

del lines=[:]

#################################################################
#		GENERATING THE GRID PARAMETER FILE		#
#################################################################

# we're still inside receptor's folder

# Creating lists

linesr=[]
linesl=[]

# Getting lists of pdbqt files to work with

with open(receptor_folder+'listpdb.txt') as lstr:
	for line in lstr:
		linesr.append(line.strip())
lstr.close()

with open(ligand_folder+'listpdb.txt') as lstl:
	for line in lstl:
		linesl.append(line.strip())
lstl.close()

# Creating gpf files

indexr=0
for indexr in range(len(linesr)):
	indexl=0	
	for indexl in range(len(linesl)):
		os.system(pyshell+'pythonsh '+utilities+'prepare_gpf4zn.py -l '+ligand_folder+linesl[indexl]+'.pdbqt -r '+receptor_folder+linesr[indexr].replace('-zn-f','')+'.pdbqt -o '+receptor_folder+linesr[indexr].replace('_tz','')+'_%d.gpf -p npts=36,36,36 -p gridcenter=35,40,34 -p parameter_file='+utilities+'AD4Zn.dat' % indexl)

# Ends everything

print " "
print "Great, grid point files (*.gpf) are ready! :D"
print " "

del lines=[:]
del linesr=[:]
del linesl=[:]

#################################################################
#		GENERATING THE MAPS PARAMETER FILE		#
#################################################################

#we're still inside receptor folder

if os.path.exists(receptor_folder+'Lister-stp3.py')==True:
	print " "
	print "Hey, I'm going to call a friend of mine to get a list of all gpf found in here."
	os.system('python '+receptor_folder+'Lister-stp3.py '+receptor_folder)
	if os.path.exists(receptor_folder+'exit.txt')==True:
	  print "You told to my friend you wouldn't like to get a new list, hope you fix it and come back. See you then."
	  os.remove(receptor_folder+'exit.txt')
	  quit()

else:
	print " "
	print "Sorry, I'm affraid 'Lister-stp3.py' isn't in receptors directory. May you call it to join in here?"
	print " "
	quit()

if os.path.exists(receptor_folder+"listgpf.txt")==True:
	print " "
	print "Yeah, 'Lister-stp3.py' has listed all gpf files. Now I'll rock things up!"
	print " "
else:
	print " "
	print "Something still wrong here, I found 'Lister-stp3.py' but he didn't give me back a list of gpfs, could you please check if he's working well?"
	print " "
	quit()

# Getting lists of gpf files to work with

with open(receptor_folder+'listgpf.txt') as lstr:
	for line in lstr:
		linesr.append(line.strip())
lstr.close()

# Creating map files

indexr=0
for indexr in range(len(linesr)):
	os.system(utilities+'/autogrid4.2.5.x.20131125 -p '+linesr[indexr]+'.gpf')

# Ends everything

print " "
print "Great, maps files are ready! :D"
print " "

del linesr=[:]

#################################################################
#		GENERATING DOCKING PARAMETER FILE		#
#################################################################

#we're still inside receptor folder

# Getting lists of pdbqt files to work with

with open(receptor_folder+'listpdb.txt') as lstr:
	for line in lstr:
		linesr.append(line.strip())
lstr.close()

with open(ligand_folder+'listpdb.txt') as lstl:
	for line in lstl:
		linesl.append(line.strip())
lstl.close()

# Creating dpf files

indexr=0
for indexr in range(len(linesr)):
	indexl=0	
	for indexl in range(len(linesl)):
		os.system(pyshell+'pythonsh '+utilities+'prepare_dpf42.py -l '+ligand_folder+linesl[indexl]+'.pdbqt -r '+receptor_folder+linesr[indexr].replace('-zn-f','')+'_tz.pdbqt -o '+receptor_folder+linesr[indexr].replace('-zn-f','')+'_%d.dpf ' % indexl)

# Ends everything

print " "
print "Great, docking parameter files are ready! :D"
print " "

del linesr=[:]
del linesl=[:]

#################################################################
#			     DOCKING				#
#################################################################

#we're still inside receptor folder

if os.path.exists(receptor_folder+'Lister-stp4.py')==True:
	print " "
	print "Hey, I'm going to call a friend of mine to get a list of all dpf found in here."
	os.system('python '+receptor_folder+'Lister-stp4.py '+receptor_folder)
	if os.path.exists(receptor_folder+'exit.txt')==True:
	  print "You told to my friend you wouldn't like to get a new list, hope you fix it and come back. See you then."
	  os.remove(receptor_folder+'exit.txt')
	  quit()

else:
	print " "
	print "Sorry, I'm affraid 'Lister-stp4.py' isn't in receptors directory. May you call it to join in here?"
	print " "
	quit()

if os.path.exists(receptor_folder+"listdpf.txt")==True:
	print " "
	print "Yeah, 'Lister-stp4.py' has listed all dpf files. Now I'll rock things up!"
	print " "
else:
	print " "
	print "Something still wrong here, I found 'Lister-stp4.py' but he didn't give me back a list of dpfs, could you please check if he's working well?"
	print " "
	quit()

# Getting lists of dpf files to work with

with open(receptor_folder+'listdpf.txt') as lstr:
	for line in lstr:
		linesr.append(line.strip())
lstr.close()

os.system('cp '+ligands_folder+'*.pdbqt ./')

indexr=0
for indexr in range(len(linesr)):
	print "Docking "+linesr[indexr]+" "
	os.system(utilities+'autodock4 -p '+receptor_folder+linesr[indexr]+'.dpf -l '+receptor_folder+linesr[indexr]+'.log')
	print "Done."
	print " "

os.system('mkdir ./Log')
os.system('mkdir ./Map')
os.system('mkdir ./PDBQT')
os.system('mkdir ./PDB')
os.system('mkdir ./DPF')
os.system('mkdir ./GPF')
os.system('mv *.log ./Log')
os.system('mv *.map ./Map')
os.system('mv *.xyz ./Map')
os.system('mv *.fld ./Map')
os.system('mv *.pdbqt ./PDBQT')
os.system('mv *.PDB ./PDB')
os.system('mv *.dpf ./DPF')
os.system('mv *.gpf ./GPF')
os.system('mv *.dpf ./DPF')
# Ends everything

print " "
print "Great, everything is done, now it's time to statiscally evaluate your docking results! Good luck and ~take care~! ;D"
print " "
quit()
