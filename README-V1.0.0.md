# GremmlenzDOCK - V.1.0.0

Docking algorithms to do VHTS of zinc enzymes

These scripts were made as an attempt to create an easy to use routine of VHTS for zinc enzymes. By following the stepts provided by GremmlenzDOCK.py one should be able to tackle the herculean job of docking thousands of compounds to a huge group of target proteins, even better, proteins whose active site has zinc coordinated to it, such the B group of beta-lactamase hydrolases.

Actually, herein these scripts heavily rely on two major works:

1) Michel F. Sanner. Python: A Programming Language for Software Integration and Development. J. Mol. Graphics Mod., 1999, Vol 17, February. pp57-61

2) AutoDock4Zn: An Improved AutoDock Force Field for Small-Molecule Docking to Zinc
Metalloproteins Santos-Martins, D., Forli, S., João Ramos, M., Olson, A., J. J.Chem.Info.Mod. 2014

These scripts are still under development and may behave unproperly under different OS.

As a disadvantage to be tackled GremmlenzDOCK can not handle multithreading and this makes things awfully slow. Thus, for the next version to be realeased in a few weeks, I'm looking forward to test a new AutoDock tool integrated with OpenMPI by using the following workt:

3) Multilevel Parallelization of AutoDock 4.2 Norgan, A.P., Coffman, P.K., Kocher,J-P.A, Katzmann, D.J. and Sosa, C.P. J.Chem.Info 2011

Any sort of corrections that you may would like to point out will be welcomed. I hope to make this work useful for whoever may like to use it.
