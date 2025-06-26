A variant of MSEIS (QSEIS), for calculating synthetic f-k spectrum (dispersion curves) based on a layered elastic halfspace model with water as the surface layer.

For Windows user, the executable file is provided under folder "WindowsEXE". Linux user may compile the source codes with "gfortran" via a single command like, e.g.,

~>cd .../SourceCode

~>gfortran -o mseisfk *.f -O3

to get the excutable code mseisfk.

After start the executable code, the program ask for an input file in the ASCII format. An example input file is provided under folder "InputFile". You may change the input data included in this file for your own applications.

A similar application example is given in QSEISFK_2011.

References

Wang, R., (1999), A simple orthonormalization method for stable and efficient computation of Green's functions, Bulletin of the Seismological Society of America, 89(3), 733-741.
