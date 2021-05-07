# BoltzmannSolver

Solver of Boltzmann equations for general Beyond the Standard Model (BSM) physics.

## Handling git

### Basic manipulations
Before doing any modification, the remote branch must be pulled:

    > git pull
Once done, modifications can be done locally. Once finished, commit and push on remote using:

    > git add .
    > git commit -m "Message describing the commit"
    > git push

### Merge conflicts
Merge conflicts can happen when remote and local branches are not properly synchronized. In that case, all conflicts must be resolved one by one. Once done, replace the `git add .` command by the same but specifying this time not the entire directory but only files with conflicts:

    > git add file_conflict

Once all conflicts have been resolved that way, changes can finally be pushed on remote branch using again:

    > git commit -m "Message describing the commit"
    > git push

## Compiling
Compiling a program in script/:

    > make program.x

Compiling all programs in script/:

    > make

## Executing programs
Executable files are stores in the bin/ directory and should be launched from
the root directory using:

    > bin/program.x

## The generated libary
The generated library solving the Boltzmann equation is then generated in 
BRparity/. 

**Warning:**
Each time the MARTY script is executed, the previous library is removed. This 
means that additional files (in particular BEv1.h/cpp) must be stored in the
directories include_lib/ src_lib/, and exported when wanted to the library 
typing:

    > ./export_to_lib.sh

**Note:**
The MARTY program actually exports itself these files, thus the only subtlety
left is to modify files in include_lib/ and src_lib/, not directly in the 
library BRparity/.


## Clang format
To apply clang format simply type the following command:

    > ./apply_clang_format.sh
