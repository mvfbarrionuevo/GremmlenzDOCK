#!/usr/bin/python
#################################################################
# GremmlenzDOCK.py - created by Manoel Barrionuevo - 2016       #
#################################################################
import sys, time, os, shutil, multiprocessing, random
import getopt, string, os.path
from os.path import expanduser
home = expanduser("~")
#################################################################

def usage():
    print "Usage: GremmlenzDOCK.py [option] [file]"
    print "    -i interactive use"
    print "    -F parameter file [needs parameter file]"
    print "    -h help"
    print
    print "Parameter file:"
    print "    <filename>.txt"
    print
    print "Parameter file should contain all input data. Please, see further information inside example.txt."
    print

def listing(p):
	lines=[]
	with open(p+'listpdb.txt') as lst:
		for line in lst:
			lines.append(line.strip())
	lst.close()
	return lines

def listlog(p):
	# Creating txt file to write down a new list
	if not os.path.exists(p+'listlog.txt'):
		print " "
		print "Hey, I'm going do write down a new file name 'listlog.txt' that will be our reference file to all log files herein."
		print " "
		lst = open(p+"listlog.txt","wr+")
	else:
		print " "
		lister_yes_no("Sorry, but you already have a 'listlog.txt' file herein. Would you like me to delete it? ")
		print "I'm going to delete 'listlog.txt' and create a new file."
		print " "
		os.remove(p+'listlog.txt')
		lst = open(p+"listlog.txt","wr+")
	# Creating lista and writing it down to a txt file
	dirs = os.listdir(p)
	listalog = []
	listaint = []
	listafin = []
	for files in dirs:
	   if files.endswith("log"):
	      listalog.append(files.replace('log',''))
	   if files.endswith(".dpf"):
	      listaint.append(files.replace('.dpf',''))
	listalog.sort()
	listaint.sort()
	listafin = [x for x in listaint if x not in listalog]
	lst.write("\n".join(listafin))
	lst.close()
	# Ends everything
	print " "
	print "Great, your backup list is now completed! :D"
	print " "
	return True

def lister2(p):
	# Creating txt file to write down a new list
	if not os.path.exists(p+'listpdb.txt'):
		print " "
		print "Hey, I'm going do write down a new file name 'listpdb.txt' that will be our reference file to all pdbs herein."
		print " "
		lst = open(p+"listpdb.txt","wr+")
	else:
		print " "
		lister_yes_no("Sorry, but you already have a 'listpdb.txt' file herein. Would you like me to delete it? ")
		print "I'm going to delete 'listpdb.txt' and create a new file."
		print " "
		os.remove(p+'listpdb.txt')
		lst = open(p+"listpdb.txt","wr+")
	# Creating lista and writing it down to a txt file
	dirs = os.listdir(p)
	lista = []
	for files in dirs:
	   if files.endswith(".pdb"):
	      lista.append(files.replace('.pdb',''))
	lista.sort()
	lst.write("\n".join(lista))
	lst.close()
	# Ends everything
	print " "
	print "Great, your list is now completed! :D"
	print " "
	return True

def lister3(p):
	# Creating txt file to write down a new list
	if not os.path.exists(p+'listgpf.txt'):
		print " "
		print "Hey, I'm going do write down a new file name 'listgpf.txt' that will be our reference file to all gpf files herein."
		print " "
		lst = open(p+"listgpf.txt","wr+")
	else:
		print " "
		lister_yes_no("Sorry, but you already have a 'listgpf.txt' file herein. Would you like me to delete it? ")
		print "I'm going to delete 'listgpf.txt' and create a new file."
		print " "
		os.remove(p+'listgpf.txt')
		lst = open(p+"listgpf.txt","wr+")
	# Creating lista and writing it down to a txt file
	dirs = os.listdir(p)
	lista = []
	for files in dirs:
	   if files.endswith(".gpf"):
	      lista.append(files.replace('.gpf',''))
	lista.sort()
	lst.write("\n".join(lista))
	lst.close()
	# Ends everything
	print " "
	print "Great, your gpf list is now completed! :D"
	print " "
	return True

def lister4(p):
	# Creating txt file to write down a new list
	if not os.path.exists(p+'listdpf.txt'):
		print " "
		print "Hey, I'm going do write down a new file name 'listdpf.txt' that will be our reference file to all dpf files herein."
		print " "
		lst = open(p+"listdpf.txt","wr+")
	else:
		print " "
		lister_yes_no("Sorry, but you already have a 'listdpf.txt' file herein. Would you like me to delete it? ")
		print "I'm going to delete 'listdpf.txt' and create a new file."
		print " "
		os.remove(p+'listdpf.txt')
		lst = open(p+"listdpf.txt","wr+")
	# Creating lista and writing it down to a txt file
	dirs = os.listdir(p)
	lista = []
	for files in dirs:
	   if files.endswith(".dpf"):
	      lista.append(files.replace('.dpf',''))
	lista.sort()
	lst.write("\n".join(lista))
	lst.close()
	# Ends everything
	print " "
	print "Great, your dpf list is now completed! :D"
	print " "
	return True

def lister_yes_no(question, default="yes"):
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
            print "You have hit [ENTER] thus I'm going to proceed as default."
            print " "
	    return False
        elif choice in valid:
            if choice == 'no' or choice == 'n':
		print " "
		print "So you have chosen '%s', thus you may like to change it by your own." % choice
		print " "
		quit()
	    elif choice == 'yes' or choice == 'y':
		print " "
		print "You have chosen '%s'." %choice		
		print " "
		return False 
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

def isnumber(s):
   try:
     float(s)
     return 'True'
   except ValueError:
     return 'False'

def isint(s,t):
   try:
     int(s) <= t
     return 'True'
   except ValueError:
     print " "
     print "Sorry your computer has %d CPUs, try to type an integer value less or equal than it." % t     
     print " "
     return 'False'


def receptor(lst,pyshell,utilities,receptor_folder):
	os.system('mkdir '+receptor_folder+'OldPDB')
	for index in range(len(lst)):
		print "Preparing receptor "+lst[index]+" "		
		os.system(utilities+'reduce.3.23.130521 -NOFLIP '+receptor_folder+''+lst[index]+'.pdb > '+receptor_folder+''+lst[index]+'_H.pdb')
		os.system('mv '+receptor_folder+''+lst[index]+'.pdb '+receptor_folder+'OldPDB/')		
		os.system(pyshell+'pythonsh '+utilities+'prepare_receptor4.py -r '+receptor_folder+''+lst[index]+'_H.pdb -A hydrogens -o '+receptor_folder+''+lst[index]+'.pdbqt')
		print "Done."
	return

def lig(lst,pyshell,ligand_folder,utilities):
	for index in range(len(lst)):
		print "Preparing ligands "+lst[index]+" "
		os.system(pyshell+'pythonsh '+utilities+'prepare_ligand4.py -l '+ligand_folder+''+lst[index]+'.pdb -A hydrogens -U lps -B amide -o '+ligand_folder+''+lst[index]+'.pdbqt')		
		print "Done."
	return

def zinc(lst,pyshell,utilities):
	for index in range(len(lst)):
		print "Preparing pseudo zinc of receptor "+lst[index]+" "
		os.system(pyshell+'pythonsh '+utilities+'zinc_pseudo.py -r '+lst[index]+'.pdbqt -o '+lst[index]+'_tz.pdbqt')
		print "Done."
	return

def gpf(lst1,lst2,x,y,z,xc,yc,zc,ligand_folder,receptor_folder,pyshell,utilities,spacing):
	indexr=0
	for indexr in range(len(lst1)):
		indexl=0	
		for indexl in range(len(lst2)):
			print "Creating gpf file for "+lst1[indexr]+"_%d" % indexl
			os.system(pyshell+'pythonsh '+utilities+'prepare_gpf4zn.py -l '+ligand_folder+''+lst2[indexl]+'.pdbqt -r '+receptor_folder+''+lst1[indexr]+'_tz.pdbqt -o '+receptor_folder+''+lst1[indexr]+'_%d.gpf -p npts=%s,%s,%s -p gridcenter=%s,%s,%s -p parameter_file=../AD4Zn.dat -p spacing=%s' % (indexl,x,y,z,xc,yc,zc,spacing))
			print "Done."	
	return

def autogrid(lst,utilities):
	for indexr in range(len(lst)):
		print "Preparing autogrid of "+lst[indexr]+" "
		os.system(utilities+'/autogrid4.2.5.x.20131125 -p '+lst[indexr]+'.gpf')
		print "Done."
	return

def dpf(lst1,lst2,pyshell,utilities,ligand_folder,receptor_folder,torsdof, rmstol, extnrg, ga_pop_size, ga_num_evals, ga_num_generations, ga_elitism, ga_mutation_rate, ga_crossover_rate, ga_window_size, ga_cauchy_alpha, ga_cauchy_beta, sw_max_its, sw_max_succ, sw_max_fail, sw_rho, sw_lb_rho, ls_search_freq, ga_run,):
	indexr=0
	for indexr in range(len(lst1)):
		indexl=0	
		for indexl in range(len(lst2)):
			print "Creating dpf files of "+lst1[indexr]+"_%d " % indexl
			os.system(pyshell+"pythonsh "+utilities+"prepare_dpf42.py -l "+ligand_folder+lst2[indexl]+".pdbqt -r "+receptor_folder+lst1[indexr]+"_tz.pdbqt -o "+receptor_folder+lst1[indexr]+"_"+str(indexl)+".dpf -p torsdof="+torsdof+" -p rmstol="+rmstol+" -p extnrg="+extnrg+"-p ga_pop_size="+ga_pop_size+" -p ga_num_evals="+ga_num_evals+" -p ga_num_generations="+ga_num_generations+" -p ga_elitism="+ga_elitism+" -p ga_mutation_rate="+ga_mutation_rate+" -p ga_crossover_rate="+ga_crossover_rate+" -p ga_window_size="+ga_window_size+" -p ga_cauchy_alpha="+ga_cauchy_alpha+" -p ga_cauchy_beta="+ga_cauchy_beta+" -p sw_max_its="+sw_max_its+" -p sw_max_succ="+sw_max_succ+" -p sw_max_fail="+sw_max_fail+" -p sw_rho="+sw_rho+" -p sw_lb_rho="+sw_lb_rho+" -p ls_search_freq="+ls_search_freq+" -p ga_run="+ga_run+"")
			print "Done."	
	return

def dock(lst,utilities,receptor_folder,pyshell):
	for i in range(len(lst)):
		print "Docking "+lst[i]+"_ "
		os.system(utilities+'autodock4 -p '+receptor_folder+''+lst[i]+'.dpf -l '+receptor_folder+''+lst[i]+'log')
		os.system(pyshell+'pythonsh '+utilities+'summarize_docking.py -a -l '+lst[i]+'log')		
		print "Done."
	return

def split_seq(line,ncpus):
	sublst = []
	splitsize = 1.0/ncpus*len(line)
	for i in range(ncpus):
		sublst.append(line[int(round(i*splitsize)):int(round((i+1)*splitsize))])
	return sublst

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
#   VERIFICATION STEP - ARE ALL SCRIPTS WHERE THEY SHOULD BE?   #
#################################################################

def automatic(par_file):
	if not os.path.exists("./"+par_file):
			print " "
			print "Sorry, but you have provided a non-existent parameter file. Please, make sure the parameter file you are providing is in the same folder of the running script."
			print " "
			quit()
	else:
			print " "
			print "The parameter file you have provided is: "+str(par_file)+"."
			p = {}
			with open(par_file) as f:
			    for line in f:
			       (key, val) = line.split()
			       p[str(key)] = val
			print " "
	
	value = False
	while value == False:
		value = check_script(p['utilities'],p['pyshell'],p['receptor_folder'])
	
	value = False
	while value == False:
		value = prep_receptor(p['receptor_folder'],p['ncpus'],p['pyshell'],p['utilities'])

	value = False
	while value == False:
		value = prep_ligand(p['ligand_folder'],p['ncpus'],p['pyshell'],p['utilities'])

	value = False
	while value == False:
		value = prep_zinc(p['receptor_folder'],p['ncpus'],p['pyshell'],p['utilities'])

	value = False
	while value == False:
		sublines=[]
		linesr=[]
		linesl=[]
		linesr = listing(p['receptor_folder'])
		linesl = listing(p['ligand_folder'])
		sublines = split_seq(linesr,int(p['ncpus']))
		value = grid_par(sublines,linesl,p['x'],p['y'],p['z'],p['xc'],p['yc'],p['zc'],p['ligand_folder'],p['receptor_folder'],p['pyshell'],p['utilities'],p['spacing'],p['ncpus'])

	value = False
	while value == False:
		value = map_par(p['receptor_folder'],p['ncpus'],p['utilities'])

	value = False
	while value == False:
		value = dock_par(p['receptor_folder'],p['ligand_folder'],p['ncpus'],p['pyshell'], p['utilities'], p['torsdof'], p['rmstol'], p['extnrg'], p['ga_pop_size'], p['ga_num_evals'], p['ga_num_generations'], p['ga_elitism'], p['ga_mutation_rate'], p['ga_crossover_rate'], p['ga_window_size'], p['ga_cauchy_alpha'], p['ga_cauchy_beta'], p['sw_max_its'], p['sw_max_succ'], p['sw_max_fail'], p['sw_rho'], p['sw_lb_rho'], p['ls_search_freq'], p['ga_run'])

	value = False
	while value == False:
		value = docking(p['ligand_folder'],p['receptor_folder'],p['ncpus'],p['utilities'],p['pyshell'])

	print
	return True

def check_script(utilities,pyshell,receptor_folder):
	os.chdir(utilities)
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
		print "Copying AD4Zn.dat to working directory: "+receptor_folder+" "	
		os.system('cp ./AD4Zn.dat '+receptor_folder)
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
	return True

def interactive():
	print " "
	print "Hi, first thing first. If you didn't get (or didn't known before) the scripts for AutoDockZn, MGLTools and Reduce you will be able to download them from here:"
	print " "
	print " 1) http://autodock.scripps.edu/resources/autodockzn-forcefield"
	print " 2) http://mgltools.scripps.edu/downloads/mgltools-1-5-7rc1"
	print " 3) http://kinemage.biochem.duke.edu/software/reduce.php"
	print " "
	print "Once you get them, please uncompress and install MGLTools, then move all autodockZN and reduce scripts to the path /path/to/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/ (wherever this path may look like according to your previous path specification)."
	print " "
	print "PLEASE DO NOT FORGET TO CITE EACH SOFTWARE YOU HAVE USED THROUGH THIS SCRIPT. YOU ARE ADVISED TO SEE HOW TO CITE EACH SOFTWARE BY ACCESSING THE PROVIDED WEB SITES."
	print " "

	print " "
	query_yes_no("So, have you done everything I told you?")
	print " "
	print "Hnmm, it seems you know what you are doing here. Lets take some addresses in order to deal with everything. Please, carefully write down the addresses I will ask you. I will not take any responsability about your typos. In case you agree with the default paths just hit [ENTER]."
	print " "

	utilities=raw_input("[$HOME/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/] What's the Utilities24 path? \n")
	if(not utilities or utilities.isspace()):
	    utilities=home+'/MGLTools-1.5.7rc1/MGLToolsPckgs/AutoDockTools/Utilities24/'
	print " "

	receptor_folder=raw_input("[$HOME/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/receptors] What's the path of your receptors? \n")
	if(not receptor_folder or receptor_folder.isspace()):
		receptor_folder=home+'/MGLTools-1.5.7rc1/MGLToolsPckgs/AutoDockTools/Utilities24/receptors/'
		print " "
		print "Standard folder created at "+receptor_folder
		print " "
		os.system('mkdir '+receptor_folder)
	else:
		os.system('mkdir '+receptor_folder)
	print " "

	ligand_folder=raw_input("[$HOME/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/ligands/] What's the path of you ligands? \n")
	if(not ligand_folder or ligand_folder.isspace()):
		ligand_folder=home+'/MGLTools-1.5.7rc1/MGLToolsPckgs/AutoDockTools/Utilities24/ligands/'
		print " "
		print "Standard folder created at "+ligand_folder
		print " "
		os.system('mkdir '+ligand_folder)
	else:
		os.system('mkdir '+ligand_folder)
	print " "

	pyshell=raw_input("[$HOME/MGLTools-1.5.7rc1/bin/] What's the pythonsh path? \n")
	if(not pyshell or pyshell.isspace()):
	    pyshell=home+'/MGLTools-1.5.7rc1/bin/'
	print " "

	var = 'False'
	t = multiprocessing.cpu_count()
	while var == 'False':
		ncpus = raw_input('Your computer has %d CPUs, how many jobs would you like to run? \n' % t)
		var = isint(ncpus,t)

	print "Now that every path has been set up, I'm going to verify through a very simple way if you really got all scripts where they should be (I'll know if you read me before saying 'yes' to me). If something is missing I'll quit the execution, otherwise your docking will not run smoothly."
	print " "

	os.chdir(utilities)
	#now we're inside ~/MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24
	if os.path.exists(receptor_folder) == True:
		print " "
		print "Please, take a moment and manually copy all your receptors to the path: "+receptor_folder+"."
		print " "
		raw_input("I am waiting for you, please, press [ENTER] after you have copied all receptor files...")
		print " "
	else:
		print " "
		print "Sorry, something went wrong, I couldn't find the path "+receptor_folder
		print " "

	if os.path.exists(ligand_folder)==True:
		print " "
		print "Please, take a moment and manually copy all your ligands to the path: "+ligand_folder+"."
		print " "
		raw_input("I am waiting for you, please, press [ENTER] after you have copied all ligand files...")
		print " "
	else:
		print " "
		print "Sorry, something went wrong, I couldn't find the path "+ligand_folder
		print " "
	value = False
	while value == False:
		value = check_script(utilities,pyshell,receptor_folder)
	
	value = False
	while value == False:
		value = prep_receptor(receptor_folder,ncpus,pyshell,utilities)

	value = False
	while value == False:
		value = prep_ligand(ligand_folder,ncpus,pyshell,utilities)
	
	value = False
	while value == False:
		value = prep_zinc(receptor_folder,ncpus,pyshell,utilities)
	
	value = False
	while value == False:
		value = grid_ipar(receptor_folder,ligand_folder,ncpus,pyshell,utilities)

	value = False
	while value == False:
		value = map_par(receptor_folder,ncpus,utilities)

	value = False
	while value == False:
		value = dock_ipar(receptor_folder,ligand_folder,ncpus,pyshell,utilities)
	
	value = False
	while value == False:
		value = docking(ligand_folder,receptor_folder,ncpus,utilities,pyshell)
	return True

#################################################################
#     PREPARATION STEP - CREATING RECEPTOR AND LIGAND PDBQT     #
#################################################################

#################################################################
#			   RECEPTOR				#
#################################################################

def prep_receptor(receptor_folder,ncpus,pyshell,utilities):
	# Creating a list
	lines=[]
	sublines=[]

	# Calling lister2 to get a list of pdb files to work with

	print "Hey, I'm going to call a friend of mine to get a list of all pdbs found in receptor directory."
	lister2(receptor_folder)

	if os.path.exists(receptor_folder+"listpdb.txt")==True:
		print " "
		print "Yeah, my friend 'lister2' has listed all receptors. Now I'll rock things up!"
		print " "
	else:
		print " "
		print "Something still wrong here, could you please check if my friend 'lister2' is working well?"
		print " "
		quit()
	
	lines = listing(receptor_folder)
	sublines = split_seq(lines,int(ncpus)) 

	# Creating pdbqt files

	if __name__ == '__main__':
		jobs = []
		for i in range(int(ncpus)):
			p = multiprocessing.Process(target=receptor, args=(sublines[i],pyshell,utilities,receptor_folder,))
			jobs.append(p)
			p.start()
		for job in jobs:
			job.join()

	# Ends everything

	print " "
	print "Great, pdbqt of receptor files are ready! :D"
	print " "
	return True
#################################################################
#			   LIGAND				#
#################################################################

def prep_ligand(ligand_folder,ncpus,pyshell,utilities):
	lines=[]
	sublines=[]

	# Calling lister2 to get a list of pdb files to work with

	lister2(ligand_folder)

	if os.path.exists(ligand_folder+"listpdb.txt")==True:
		print " "
		print "Yeah, 'lister2' has listed all ligands. Now I'll rock things up!"
		print " "
	else:
		print " "
		print "Something still wrong here, could you please check if my friend 'lister2' is working well?"
		print " "
		quit()

	lines = listing(ligand_folder)
	sublines = split_seq(lines,int(ncpus)) 

	# Creating pdbqt files

	if __name__ == '__main__':
		jobs = []
		for i in range(int(ncpus)):
			p = multiprocessing.Process(target=lig, args=(sublines[i],pyshell,ligand_folder,utilities))
			jobs.append(p)
			p.start()
		for job in jobs:
			job.join()

	# Ends everything

	print " "
	print "Great, pdbqt of ligand files are ready! :D"
	print " "
	return True

#################################################################
#      PREPARATION STEP - CREATING ZINC PSEUDO ATOM CENTRES     #
#################################################################

def prep_zinc(receptor_folder,ncpus,pyshell,utilities):
	lines=[]
	sublines=[]	
	os.chdir(receptor_folder)

	lines = listing(receptor_folder)
	sublines = split_seq(lines,int(ncpus)) 

	# Creating _tz.pdbqt files

	if __name__ == '__main__':
		jobs = []
		for i in range(int(ncpus)):
			p = multiprocessing.Process(target=zinc, args=(sublines[i],pyshell,utilities))
			jobs.append(p)
			p.start()
		for job in jobs:
			job.join()

	# Ends everything

	print " "
	print "Great, pseudo zinc atoms were added to receptor pdbqt files! :D"
	print " "
	return True

#################################################################
#		GENERATING THE GRID PARAMETER FILE		#
#################################################################

def grid_ipar(receptor_folder,ligand_folder,ncpus,pyshell,utilities):

	sublines=[]
	linesr=[]
	linesl=[]

	linesr = listing(receptor_folder)
	linesl = listing(ligand_folder)
	sublines = split_seq(linesr,int(ncpus))

	print " "
	print "Now, in order to dock I need some box parameters, please type them carefully!"
	print " "

	var = 'False'
	while var == 'False':	
		x = raw_input("Please type 'x' box size: ")
		var = isnumber(x)
	var = 'False'
	while var == 'False':	
		y = raw_input("Please type 'y' box size: ")
		var = isnumber(y)
	var = 'False'
	while var == 'False':	
		z = raw_input("Please type 'z' box size: ")
		var = isnumber(z)
	var = 'False'
	while var == 'False':	
		xc = raw_input("Please type 'xc' box position: ")
		var = isnumber(xc)
	var = 'False'
	while var == 'False':	
		yc = raw_input("Please type 'yc' box position: ")
		var = isnumber(yc)
	var = 'False'
	while var == 'False':	
		zc = raw_input("Please type 'zc' box position: ")
		var = isnumber(zc)
	var = 'False'
	while var == 'False':
		spacing = raw_input("Please type the box spacing: ")
		var = isnumber(spacing)	
	print " "
	value = False	
	while value == False:
		value = grid_par(sublines,linesl,x,y,z,xc,yc,zc,ligand_folder,receptor_folder,pyshell,utilities,spacing,ncpus)
	return

def grid_par(sublines,linesl,x,y,z,xc,yc,zc,ligand_folder,receptor_folder,pyshell,utilities,spacing,ncpus):
	# Creating gpf files

	if __name__ == '__main__':
		jobs = []
		for i in range(int(ncpus)):
			p = multiprocessing.Process(target=gpf, args=(sublines[i],linesl,x,y,z,xc,yc,zc,ligand_folder,receptor_folder,pyshell,utilities,spacing))
			jobs.append(p)
			p.start()
		for job in jobs:
			job.join()

	# Ends everything

	print " "
	print "Great, grid point files (*.gpf) are ready! :D"
	print " "
	return True
#################################################################
#		GENERATING THE MAPS PARAMETER FILE		#
#################################################################

def map_par(receptor_folder,ncpus,utilities):
	lines=[]
	sublines=[]
	linesr=[]
	os.chdir(receptor_folder)

	lister3(receptor_folder)

	if os.path.exists(receptor_folder+"listgpf.txt")==True:
		print " "
		print "Yeah, 'lister3' has listed all gpf files. Now I'll rock things up!"
		print " "
	else:
		print " "
		print "Something still wrong here, could you please check if my friend 'lister3' is working well?"
		print " "
		quit()

	# Getting lists of gpf files to work with

	with open(receptor_folder+'listgpf.txt') as lstr:
		for line in lstr:
			linesr.append(line.strip())
	lstr.close()

	sublines = split_seq(linesr,int(ncpus)) 

	# Creating map files

	if __name__ == '__main__':
		jobs = []
		for i in range(int(ncpus)):
			p = multiprocessing.Process(target=autogrid, args=(sublines[i],utilities,))
			jobs.append(p)
			p.start()
		for job in jobs:
			job.join()

	# Ends everything

	print " "
	print "Great, maps files are ready! :D"
	print " "
	return True

#################################################################
#		GENERATING DOCKING PARAMETER FILE		#
#################################################################

def dock_ipar(receptor_folder,ligand_folder,ncpus,pyshell,utilities):
	print " "
	torsdof = raw_input("Please, type the torsional degrees of freedom [default 5]: ")
	if(not torsdof or torsdof.isspace()):
	    torsdof='5'
	print " "

	rmstol = raw_input("Please, type the RMS cluster tolerance [default 2.0]: ")
	if(not rmstol or rmstol.isspace()):
	    rmstol='2.0'
	print " "

	extnrg = raw_input("Please, type the external grid energy [default 1000.0]: ")
	if(not extnrg or extnrg.isspace()):
	    extnrg='1000.0'
	print " "

	ga_pop_size = raw_input("Please, type the number of individuals in population [default 150]: ")
	if(not ga_pop_size or ga_pop_size.isspace()):
	    ga_pop_size='150'
	print " "

	ga_num_evals = raw_input("Please, type the maximum number of energy evaluations [default 2500000]: ")
	if(not ga_pop_size or ga_pop_size.isspace()):
	    ga_num_evals='2500000'
	print " "

	ga_num_generations = raw_input("Please, type the maximum number of generations [default 27000]: ")
	if(not ga_num_generations or ga_num_generations.isspace()):
	    ga_num_generations='27000'
	print " "

	ga_elitism = raw_input("Please, type the number of top individuals to survive to next generations [default 1]: ")
	if(not ga_elitism or ga_elitism.isspace()):
	    ga_elitism='1'
	print " "

	ga_mutation_rate = raw_input("Please, type the rate of mutation [default 0.02]: ")
	if(not ga_mutation_rate or ga_mutation_rate.isspace()):
	    ga_mutation_rate='0.02'
	print " "

	ga_crossover_rate = raw_input("Please, type the rate of crossover [default 0.8]: ")
	if(not ga_crossover_rate or ga_crossover_rate.isspace()):
	    ga_crossover_rate='0.8'
	print " "

	ga_window_size = raw_input("Please, type the number of generations for picking the worst individual [default 10]: ")
	if(not ga_window_size or ga_window_size.isspace()):
	    ga_window_size='10'
	print " "

	ga_cauchy_alpha = raw_input("Please, type the alpha parameter of Cauchy distribution [default 0.0]: ")
	if(not ga_cauchy_alpha or ga_cauchy_alpha.isspace()):
	    ga_cauchy_alpha='0.0'
	print " "

	ga_cauchy_beta = raw_input("Please, type the beta parameter of Cauchy distribution [default 1.0]: ")
	if(not ga_cauchy_beta or ga_cauchy_beta.isspace()):
	    ga_cauchy_beta='1.0'
	print " "

	sw_max_its = raw_input("Please, type the maximum iterations of Solis & Wets local search [default 300]: ")
	if(not sw_max_its or sw_max_its.isspace()):
	    sw_max_its='300'
	print " "

	sw_max_succ = raw_input("Please, type the maximum consecutive successes before changing rho of Solis & Wets [default 4]: ")
	if(not sw_max_succ or sw_max_succ.isspace()):
	    sw_max_succ='4'
	print " "

	sw_max_fail = raw_input("Please, type the maximum consecutive failures before changing rho of Solis & Wets [default 4]: ")
	if(not sw_max_fail or sw_max_fail.isspace()):
	    sw_max_fail='4'
	print " "

	sw_rho = raw_input("Please, type the size of local search space to sample of Solis & Wets [default 1.0]: ")
	if(not sw_rho or sw_rho.isspace()):
	    sw_rho='1.0'
	print " "

	sw_lb_rho = raw_input("Please, type the lower bound on rho of Solis & Wets [default 0.01]: ")
	if(not sw_lb_rho or sw_lb_rho.isspace()):
	    sw_lb_rho='0.01'
	print " "

	ls_search_freq = raw_input("Please, type the probability of performing local search on individual [default 0.06]: ")
	if(not ls_search_freq or ls_search_freq.isspace()):
	    ls_search_freq='0.06'
	print " "

	ga_run = raw_input("Please, type the number of hybrid GA-LS runs [default 10]: ")
	if(not ga_run or ga_run.isspace()):
	    ga_run='10'
	print " "

	value = False
	while value == False:
		value = dock_par(receptor_folder,ligand_folder,ncpus,pyshell, utilities, torsdof, rmstol, extnrg, ga_pop_size, ga_num_evals, ga_num_generations, ga_elitism, ga_mutation_rate, ga_crossover_rate, ga_window_size, ga_cauchy_alpha, ga_cauchy_beta, sw_max_its, sw_max_succ, sw_max_fail, sw_rho, sw_lb_rho, ls_search_freq, ga_run)
	return True

def dock_par(receptor_folder,ligand_folder,ncpus,pyshell, utilities, torsdof, rmstol, extnrg, ga_pop_size, ga_num_evals, ga_num_generations, ga_elitism, ga_mutation_rate, ga_crossover_rate, ga_window_size, ga_cauchy_alpha, ga_cauchy_beta, sw_max_its, sw_max_succ, sw_max_fail, sw_rho, sw_lb_rho, ls_search_freq, ga_run):
	sublines=[]
	linesr=[]
	linesl=[]
	
	linesr = listing(receptor_folder)
	linesl = listing(ligand_folder)
	sublines = split_seq(linesr,int(ncpus)) 

	# Creating dpf files

	if __name__ == '__main__':
		jobs = []
		for i in range(int(ncpus)):
			p = multiprocessing.Process(target=dpf, args=(sublines[i],linesl,pyshell,utilities,ligand_folder,receptor_folder,torsdof, rmstol, extnrg, ga_pop_size, ga_num_evals, ga_num_generations, ga_elitism, ga_mutation_rate, ga_crossover_rate, ga_window_size, ga_cauchy_alpha, ga_cauchy_beta, sw_max_its, sw_max_succ, sw_max_fail, sw_rho, sw_lb_rho, ls_search_freq, ga_run,))
			jobs.append(p)
			p.start()
		for job in jobs:
			job.join()

	# Ends everything

	print " "
	print "Great, docking parameter files are ready! :D"
	print " "
	return True


#################################################################
#	GETTING BACK FROM THE LAST ENDING POINT & DOCKING	#
#################################################################

# Do not change these coments below, unless you know what you are doing. Remember you are going to take all responsibility for possible malfunctioning of this script.
#listlog(receptor_folder)

#if os.path.exists(receptor_folder+"listlog.txt")==True:
#	print " "
#	print "There is a backup point, I'll try to handle this."
#	print " "
#else:
#	print " "
#	print "Something went wrong, it was supposed to have a backup list of log files here."
#	print " "
#	quit()

def docking(ligand_folder,receptor_folder,ncpus,utilities,pyshell):
	lines=[]
	sublines=[]
	linesr=[]
	linesl=[]
	os.chdir(receptor_folder)
	os.system('cp '+ligand_folder+'*.pdbqt '+receptor_folder)
	# Getting lists of dpf files to work with
	listlog(receptor_folder)
	with open(receptor_folder+'listlog.txt') as lstr:
		for line in lstr:
			linesr.append(line.strip())
	lstr.close()

	if len(linesr) < 1:
		del linesr[:]
		del sublines[:]
		lister4(receptor_folder)
		if os.path.exists(receptor_folder+"listdpf.txt")==True:
			print " "
			print "Yeah, 'lister4' has listed all dpf files. Now I'll rock things up!"
			print " "
		else:
			print " "
			print "Something still wrong here, could you please check if my friend 'lister4' is working well?"
			print " "
			quit()
		# Getting lists of dpf files to work with
		with open(receptor_folder+'listdpf.txt') as lstr:
			for line in lstr:
				linesr.append(line.strip())
		lstr.close()
		sublines = split_seq(linesr,int(ncpus)) 
		os.system('cp '+ligand_folder+'*.pdbqt ./')
		if __name__ == '__main__':
			jobs = []
			for i in range(int(ncpus)):
				p = multiprocessing.Process(target=dock, args=(sublines[i],utilities,receptor_folder,pyshell,))
				jobs.append(p)
				p.start()
			for job in jobs:
				job.join()
	else:
		sublines = split_seq(linesr,int(ncpus))
		if __name__ == '__main__':
			jobs = []
			for i in range(int(ncpus)):
				p = multiprocessing.Process(target=dock, args=(sublines[i],utilities,receptor_folder,pyshell,))
				jobs.append(p)
				p.start()
			for job in jobs:
				job.join()

	os.system('mkdir ./Log')
	os.system('mkdir ./Map')
	os.system('mkdir ./PDBQT')
	os.system('mkdir ./PDB')
	os.system('mkdir ./DPF')
	os.system('mkdir ./GPF')
	os.system('mv *log ./Log')
	os.system('mv *map ./Map')
	os.system('mv *xyz ./Map')
	os.system('mv *fld ./Map')
	os.system('mv *pdbqt ./PDBQT')
	os.system('mv *pdb ./PDB')
	os.system('mv *dpf ./DPF')
	os.system('mv *gpf ./GPF')
	return True

if __name__ == '__main__':
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'ihF:',["file="])
    except getopt.GetoptError, msg:
        print 'GremmelenzDOCK_pseudo.py: %s' % msg
        usage()
        sys.exit(2)

    for o, a in opt_list:
        if o in ('-h', '--'):
            usage()
            sys.exit()
        if o in ('-i', '--i'):
            print
	    print "Interactive command line has been chosen."
	    value = False
	    while value == False:
	    	value = interactive()
	    print
	if o in ('-F', '--F'):
	    print
	    print "Taking parameter file."
	    print
	    par_file = a
	    value = False
	    while value == False:
		value = automatic(par_file)

# Ends everything

	print " "
	print "Great, everything is done, now it's time to statiscally evaluate your docking results! Good luck and ~take care~! ;D"
	print " "
	quit()
