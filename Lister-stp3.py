#!/usr/bin/python
#################################################################
# Lister-stp3.py - created by Manoel Barrionuevo - 2016         #
#################################################################

import sys, time, os, shutil

# Defining yes or no query verifier

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
            print "You have hit [ENTER] thus I'm going to proceed as default."
            print " "
	    return False
        elif choice in valid:
            if choice == 'no' or choice == 'n':
		print " "
		print "So you have chosen '%s', thus you may like to change 'list.txt' by your own." % choice
		print " "
		exiting = open('exit.txt', 'wr+')
		exiting.close()
		quit()
	    elif choice == 'yes' or choice == 'y':
		print " "
		print "You have chosen '%s'." %choice		
		print " "
		return False 
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

# Getting arguments to work with

path = sys.argv[1]
dirs = os.listdir(path)

# Creating txt file to write down a new list

if not os.path.exists('./listgpf.txt'):
	print " "
	print "Hey, I'm going do write down a new file name 'listgpf.txt' that will be our reference file to all gpf files herein."
	print " "
	lst = open("./listgpf.txt","wr+")
else:
	print " "
	query_yes_no("Sorry, but you already have a 'listgpf.txt' file herein. Would you like me to delete it? ")
	print "I'm going to delete 'listgpf.txt' and create a new file."
	print " "
	os.remove('./listgpf.txt')
	lst = open("./listgpf.txt","wr+")

# Creating lista and writing it down to a txt file

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
quit()
