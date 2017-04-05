#!/usr/bin/python
#################################################################
# GremmlenzDOCK.py - created by Manoel Barrionuevo - 2016       #
# Last update: 04/04/2017					#
#################################################################
import sys, time, os, shutil, multiprocessing, random
import getopt, string, os.path
from os.path import expanduser
home = expanduser("~")
#################################################################

def usage():
    print "\nUsage: GremmlenzDOCK.py [option] [file]"
    print "    Use one of the parameter file options  [needs parameter file]: -f --f -F --F --file --File"
    print "    Print this help section: -h --help -H --Help\n"
    print "\nParameter file:"
    print "    <filename>.txt\n"
    print "\nExample of file option usage:\n"
    print "    ./GremmelenzDOCK.py -f myfile.txt\n"
    print "\nParameter file should contain all input data. Please, see further information inside example.txt.\n"

def listing(p):
        lines=[]
        dirs = os.listdir(p)
        for files in dirs:
                if files.endswith(".pdb"):
                        lines.append(files.replace(".pdb",""))
        lines.sort()
        return lines

def listlog(p):
	# Creating lista and writing it down to a txt file
	dirs = os.listdir(p)
	lstgpf = []
	lstgpf = listgpf(p)
	lstgpf.sort()
	lstdpf = []
	lstdpf = listdpf(p)
	lstdpf.sort()
	lstremainder = []
	lstremainder = [x for x in lstgpf if x not in lstdpf]
	# Ends everything
	print "\nGreat, your backup list is now completed! :D\n"
	return lstremainder

def listgpf(p):
	# Creating list of gpf files
	lstgpf = []
        subs = os.listdir(p)
        os.chdir(p)
        for dirs in subs:
         if os.path.isdir(dirs):
          os.chdir(dirs)
          lst = os.listdir('./')
          for files in lst:
           if files.endswith(".gpf"):
                lstgpf.append(files.replace('.gpf',''))
          os.chdir('../')
        lstgpf.sort()
	# Ends everything
	print "\nGreat, your gpf list is now completed! :D\n"
	return lstgpf 

def listdpf(p):
	# Creating list of dpf files
	lstdpf = []
        subs = os.listdir(p)
        os.chdir(p)
        for dirs in subs:
         if os.path.isdir(dirs):
          os.chdir(dirs)
          lst = os.listdir('./')
          for files in lst:
           if files.endswith(".dpf"):
                lstdpf.append(files.replace('.dpf',''))
          os.chdir('../')
        lstdpf.sort()
	# Ends everything
	print "\nGreat, your dpf list is now completed! :D\n"
	return lstdpf

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
		print "Preparing receptor "+lst[index]+"."		
		os.system(utilities+'reduce.3.23.130521 -NOFLIP '+receptor_folder+''+lst[index]+'.pdb > '+receptor_folder+''+lst[index]+'_H.pdb')
		os.system('mv '+receptor_folder+''+lst[index]+'.pdb '+receptor_folder+'OldPDB/')		
		os.system(pyshell+'pythonsh '+utilities+'prepare_receptor4.py -r '+receptor_folder+''+lst[index]+'_H.pdb -A hydrogens -o '+receptor_folder+''+lst[index]+'_H.pdbqt')
		print "Done."
	return

def lig(lst,pyshell,ligand_folder,utilities):
	for index in range(len(lst)):
		print "Preparing ligands "+lst[index]+"."
		os.system(pyshell+'pythonsh '+utilities+'prepare_ligand4.py -l '+ligand_folder+''+lst[index]+'.pdb -A hydrogens -U lps -B amide -o '+ligand_folder+''+lst[index]+'.pdbqt')		
		print "Done."
	return

def zinc(lst,pyshell,utilities):
	for index in range(len(lst)):
		print "Preparing pseudo zinc of receptor "+lst[index]+"."
		os.system(pyshell+'pythonsh '+utilities+'zinc_pseudo.py -r '+lst[index]+'.pdbqt -o '+lst[index]+'_tz.pdbqt')
		print "Done."
	return

def gpf(lst1,lst2,x,y,z,xc,yc,zc,ligand_folder,receptor_folder,pyshell,utilities,spacing):
	indexr=0
	for indexr in range(len(lst1)):
		indexl=0	
		for indexl in range(len(lst2)):
			print "Creating gpf file for "+lst1[indexr]+"_"+lst2[indexl]+"."
			os.system('mkdir '+lst1[indexr]+'_'+lst2[indexl]+'/')
			os.system('cp '+receptor_folder+'AD4Zn.dat '+receptor_folder+''+lst1[indexr]+'_'+lst2[indexl]+'/')
			os.system('cp '+ligand_folder+''+lst2[indexl]+'.pdbqt '+receptor_folder+''+lst1[indexr]+'_'+lst2[indexl]+'/')
			os.system('cp '+receptor_folder+''+lst1[indexr]+'_tz.pdbqt '+receptor_folder+''+lst1[indexr]+'_'+lst2[indexl]+'/')
			os.system(pyshell+'pythonsh '+utilities+'prepare_gpf4zn.py -l '+ligand_folder+''+lst2[indexl]+'.pdbqt -r '+receptor_folder+''+lst1[indexr]+'_'+lst2[indexl]+'/'+lst1[indexr]+'_tz.pdbqt -o '+receptor_folder+''+lst1[indexr]+'_'+lst2[indexl]+'/'+lst1[indexr]+'_'+lst2[indexl]+'.gpf -p npts=%s,%s,%s -p gridcenter=%s,%s,%s -p parameter_file=./AD4Zn.dat -p spacing=%s' % (x,y,z,xc,yc,zc,spacing))
			print "Done."	
	return

def autogrid(lst,utilities,receptor_folder):
	for indexr in range(len(lst)):
		print "Preparing autogrid of "+lst[indexr]+"."
		os.chdir(receptor_folder+''+lst[indexr])
		os.system(utilities+'/autogrid4 -p '+lst[indexr]+'.gpf -l '+receptor_folder+''+lst[indexr]+'/'+lst[indexr]+'.glg')
		print "Done."
	return

def dpf(lst1,lst2,pyshell,utilities,ligand_folder,receptor_folder,torsdof, rmstol, extnrg, ga_pop_size, ga_num_evals, ga_num_generations, ga_elitism, ga_mutation_rate, ga_crossover_rate, ga_window_size, ga_cauchy_alpha, ga_cauchy_beta, sw_max_its, sw_max_succ, sw_max_fail, sw_rho, sw_lb_rho, ls_search_freq, ga_run,):
	indexr=0
	for indexr in range(len(lst1)):
		indexl=0	
		for indexl in range(len(lst2)):
			print "Creating dpf files of "+lst1[indexr]+"_"+lst2[indexl]+"."
			os.system(pyshell+"pythonsh "+utilities+"prepare_dpf42.py -l "+ligand_folder+''+lst2[indexl]+".pdbqt -r "+receptor_folder+lst1[indexr]+'_'+lst2[indexl]+"/"+lst1[indexr]+"_tz.pdbqt -o "+receptor_folder+lst1[indexr]+"_"+lst2[indexl]+"/"+lst1[indexr]+"_"+lst2[indexl]+".dpf -p torsdof="+torsdof+" -p rmstol="+rmstol+" -p extnrg="+extnrg+"-p ga_pop_size="+ga_pop_size+" -p ga_num_evals="+ga_num_evals+" -p ga_num_generations="+ga_num_generations+" -p ga_elitism="+ga_elitism+" -p ga_mutation_rate="+ga_mutation_rate+" -p ga_crossover_rate="+ga_crossover_rate+" -p ga_window_size="+ga_window_size+" -p ga_cauchy_alpha="+ga_cauchy_alpha+" -p ga_cauchy_beta="+ga_cauchy_beta+" -p sw_max_its="+sw_max_its+" -p sw_max_succ="+sw_max_succ+" -p sw_max_fail="+sw_max_fail+" -p sw_rho="+sw_rho+" -p sw_lb_rho="+sw_lb_rho+" -p ls_search_freq="+ls_search_freq+" -p ga_run="+ga_run+"")
			print "Done."	
	return

def dock(lst,utilities,receptor_folder,pyshell):
	for i in range(len(lst)):
		print "Docking "+lst[i]+"."
		os.chdir(receptor_folder+''+lst[i])
		os.system(utilities+'autodock4 -p '+lst[i]+'.dpf -l '+lst[i]+'.log')
		os.system(pyshell+'pythonsh '+utilities+'summarize_docking.py -a -l '+receptor_folder+''+lst[i]+'/'+lst[i]+'.log -o '+receptor_folder+'summary_docking.txt')		
		os.system('rm -rf AD4Zn.dat ')
		os.system('rm -rf *.pdbqt')
		print "Done."
	os.system('rm -rf '+receptor_folder+'AD4Zn.dat')
	return

def split_seq(line,ncpus):
	sublst = []
	splitsize = 1.0/ncpus*len(line)
	for i in range(ncpus):
		sublst.append(line[int(round(i*splitsize)):int(round((i+1)*splitsize))])
	return sublst

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
		print "prepare_gpf4zn.py has been found.\n"
	else:
		print "Sorry, prepare_gpf4zn.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'.\n"
		quit()
	if os.path.exists('./prepare_dpf42.py')==True:
		print "prepare_dpf42.py has been found.\n"
	else:
		print "Sorry, prepare_dpf42.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'.\n"
		quit()
	if os.path.exists('./prepare_receptor4.py')==True:
		print "prepare_receptor4.py has been found.\n"
	else:
		print "Sorry, prepare_receptor4.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'.\n"
		quit()
	if os.path.exists('./prepare_ligand4.py')==True:
			print "prepare_gpf4zn.py has been found.\n"
	else:
		print "Sorry, prepare_ligand4.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'.\n"
		quit()
	if os.path.exists('./zinc_pseudo.py')==True:
		print "zinc_pseudo.py has been found.\n"
	else:
		print "Sorry, zinc_pseudo.py isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'.\n"
		quit()
	if os.path.exists('./autodock4')==True:
		print "autodock4.py has been found.\n"
	else:
		print "Sorry, autodock4 isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'.\n"
		quit()
	if os.path.exists('./autogrid4')==True:
		print "autogrid4 has been found.\n"
	else:
		print "Sorry, autogrid4 isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'.\n"
		quit()
	if os.path.exists('./AD4Zn.dat')==True:
		print "AD4Zn.dat has been found.\n"
		print "Copying AD4Zn.dat to working directory: "+receptor_folder+"\n"	
		os.system('cp ./AD4Zn.dat '+receptor_folder)
	else:
		print "Sorry, AD4Zn.dat isn't in 'MGLTools-1.5.7rc1/MGLToolsPackgs/AutoDockTools/Utilities24/'.\n"
		quit()

	os.chdir(pyshell)
	#now we're inside ~/MGLTools-1.5.7rc1/bin/
	if os.path.exists('./pythonsh')==True:
		print "pythonsh has been found.\n"
	else:
		print "Sorry, pythonsh isn't in 'MGLTools-1.5.7rc1/bin/'.\n"
		quit()
	return True

#################################################################
#			   RECEPTOR				#
#################################################################

def prep_receptor(receptor_folder,ncpus,pyshell,utilities):
	# Creating a list
	lines=[]
	sublines=[]

	# Listing all pdb targets

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

	print "\nGreat, pdbqt of receptor files are ready! :D\n"
	return True
#################################################################
#			   LIGAND				#
#################################################################

def prep_ligand(ligand_folder,ncpus,pyshell,utilities):
	lines=[]
	sublines=[]

	# Listing all pdb bullets

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

	print "\nGreat, pdbqt of ligand files are ready! :D\n"
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

	print "\nGreat, pseudo zinc atoms were added to receptor pdbqt files! :D\n"
	return True

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

	print "\nGreat, grid point files (*.gpf) are ready! :D\n"
	return True
#################################################################
#		GENERATING THE MAPS PARAMETER FILE		#
#################################################################

def map_par(receptor_folder,ncpus,utilities):
	lines=[]
	sublines=[]
	linesr=[]
	os.chdir(receptor_folder)

	# Getting lists of gpf files to work with
	
	linesr = listgpf(receptor_folder)

	sublines = split_seq(linesr,int(ncpus)) 
	
	# Creating map files
####################################################################################################################################################################	
#		MUST TO KEEP DEBBUGING FROM HERE, IN POSITIVE CASE, YOU SHALL CHANGE ALL REMAINING CODE SIMILAR TO THE OLD ONE IN HERE 				   #
####################################################################################################################################################################	
	
	# It has been changed by Manoel in 03/04/2017 for the sake of memory usage. Before, it was building up until it explodes to a higher amount of memory than
	# that acceptable by the host computer. Thus, it was chocking to notice an absence of a memory cleaning routine. It shall be tested for futher debugging.

	if __name__ == '__main__':
		jobs = []
		for i in range(int(ncpus)):
			p = multiprocessing.Process(target=autogrid, args=(sublines[i],utilities,receptor_folder,))
			jobs.append(p)
			p.start()
		for job in jobs:
			job.join()

####################################################################################################################################################################	

	# Ends everything

	print "\nGreat, maps files are ready! :D\n"
	return True

#################################################################
#		GENERATING DOCKING PARAMETER FILE		#
#################################################################

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

	print "\nGreat, docking parameter files are ready! :D\n"
	return True


#################################################################
#	GETTING BACK FROM THE LAST ENDING POINT & DOCKING	#
#################################################################

def docking(ligand_folder,receptor_folder,ncpus,utilities,pyshell):
	lines=[]
	sublines=[]
	linesr=[]
	linesl=[]
	os.chdir(receptor_folder)

	# Checking for old runs. Listlog is going to get any old runs and update a backup list named listlog.txt. Thus, listlog.txt is going to be read once again
	# to ensure that it has something in it, in positive case (there is more than 1 line in it) it will restart from it has stopped
	# if not positive (listlog.txt has less than 0 line) it will begin a new wholy new docking campaign

	linesr = listlog(receptor_folder)

	if len(linesr) < 1:
		del linesr[:]
		del sublines[:]
		linesr = listdpf(receptor_folder)
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

	os.system('rm -rf *_tz.pdbqt')
	os.system('rm *pdbqt')
	os.system('rm *pdb')
	os.system('rm -rf '+ligand_folder+'*.pdbqt')
	return True

#################################################################
#		       WHERE EVERYGTHING BEGINS		        #
#################################################################

if __name__ == '__main__':
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'ihF:',["file="])
	if not opt_list:
		print "\nNo options supplied.\n"
		usage()
		sys.exit()
    except getopt.GetoptError, msg:
        print '\nError: %s\n' % msg
	usage()
        sys.exit(2)

    for o, a in opt_list:
	if o in ('-h', '--Help','-H','--help'):
            usage()
            sys.exit()
	if o in ('-F', '--File','--F','--file','-f', '--f'):
	    print "\nTaking parameter file.\n"
	    par_file = a
	    value = False
	    while value == False:
		value = automatic(par_file)
	else:
	    usage()
	    sys.exit()

# Ends everything

	print "\nGreat, everything is done, now it's time to statiscally evaluate your docking results! Good luck and ~take care~! ;D\n"
	quit()
