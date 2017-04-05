# GremmlenzDOCK - V2.0

Docking algorithms for VHTS of zinc enzymes.

This script is a wholy new version of GremmlenzDOCK - V1.0.0, in which new functions were added in order to make docking pipeline faster with lower memory need. Interactive methods were all removed and now the wholy script is going to work only with a protocol file (also provided herein as protocol_example_file), which must be carefully written before submmiting to GremmlenzDOCK. As stated before, GremmlenzDOCK were made as an attempt to create an easy to use routine of VHTS for zinc enzymes. GremmlenzDOCK should be able to tackle the herculean job of docking thousands of compounds (bullets) to a huge group of proteins (targets), even better, proteins whose active site has zinc coordinated to it, such the B group of beta-lactamase hydrolases.

Actually, herein these scripts heavily rely on two major works:

1) Michel F. Sanner. Python: A Programming Language for Software Integration and Development. J. Mol. Graphics Mod., 1999, Vol 17, February. pp57-61

2) AutoDock4Zn: An Improved AutoDock Force Field for Small-Molecule Docking to Zinc Metalloproteins Santos-Martins, D., Forli, S., João Ramos, M., Olson, A., J. J.Chem.Info.Mod. 2014

GremmlenzDOCK still under development and may behave unproperly under different OS. We just can guarantee that GremmlenzDOCK will work at Linux based systems with python and glibc libraries installed.

There are some observed disadvantages:

i) As GremmlenzDOCK relies on AutoDock4ZN, which generate a lot of files to dock any bullet into any target (such as .xyz, .fld, .map, .gpf, .dpf etc), the amount of disk space required can be awfully huge if you are going to perform VHTS using huge pool of compounds (bullets);

ii) As GremmlenzDOCK is based on Python language, multithreading is done by python multiprocess internal routine and it can be trick since the user can set a higher amount of CPUs to use, though this number can be greater than the available number of physical cores. By doing so, one can observe that there is no advantage to set a number of CPUs higher than the avilalbe number of physical cores.

As future work, we have proposed to compile a new version of AutoDock4.2 that work with multilevel parallelization as shown by:

3) Multilevel Parallelization of AutoDock 4.2 Norgan, A.P., Coffman, P.K., Kocher,J-P.A, Katzmann, D.J. and Sosa, C.P. J.Chem.Info 2011

However, since there is a lack of information and doubts about computer architecture where such compiled version could work properly, this task was delayed and will be considered to be abondoned (I'm alone with this and would like to get some help). Any sort of corrections that you may would like to point out will be welcomed. I hope to make this work useful for whoever may like to use it.
