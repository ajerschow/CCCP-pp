# CCCP++
Calculation of Coherence Pathways for NMR Spectroscopy

https://wp.nyu.edu/jerschow/resources/cccp-complete-calculation-of-coherence-pathways/

A. Jerschow and R. Kumar, Calculation of Coherence Pathway Selection and Cogwheel Cycles,
J. Magn. Reson. 160, 59-64, (2003).
http://dx.doi.org/10.1016/S1090-7807(02)00031-9

Complete Calculation of Coherence Pathways – CCCP++

Some background and calculations are described in

A. Jerschow and R. Kumar, Calculation of Coherence Pathway Selection and Cogwheel Cycles,
J. Magn. Reson. 160, 59-64, (2003).
http://dx.doi.org/10.1016/S1090-7807(02)00031-9

Please cite this article when you use the program.

Quick Start

Download CCCP++v1.2.tgz, current version is 1.2. Works with Linux and OSX (10.4 Tiger). The older versions are available below.
Unpack by “tar -xvzf CCCP++v1.2.tgz” in an appropriate directory.
Type “make”
run the tests by “make test”.
read comments in the supplied “*.par” files for information on how to write such files for new experiments.
Compilation Issues

The compilation was tested on SuSe Linux 8.2 and MacOS 10.3 and 10.4. If you experience compilation problems you could try out the following binaries which were compiled on different Unix systems:
v. 1.0 compiled on Mandrake 8.0
v. 1.1 compiled on SuSe Linux 8.2
v. 1.2 compiled on OSX 10.4 (Tiger)
older src distributions:
CCCP++v1.1.tgz
CCCP++v1.0.tgz

If you get an error related to “parse.*” try the following: “make parse.cc”, then “make”.

The error ” ‘yy_current_buffer’ was not declared in this scope ” can be fixed by adding a line in the Makefile saying “USE_GMAKE=Yes” (fix reported for Fedora 8 OS and Ubuntu by Rangeet Bhattacharyya).

If you get an error related to “flex++” not found then you need to install the proper version of flex++ (available under the GNU license).

The program was written in ANSI C++ using the Standard Template Library. It is expected to compile without difficulties (using the supplied Makefile) on any GNU g++ compatible compiler of version 2.95 or higher. Some lower versions probably would work also, as well as c++ compilers as long as they are ANSI compliant. It has been compiled on Linux (Mandrake 8.2, i586), SGI (Irix 6.5, mips R12000). The parsing module was created using flex++ version 2.5.1 and 2.5.4, it is possible that these versions or later are required, although this could not be verified. (flex++ is available under the GNU license). The program itself including all source code is given away under the GNU General Public License, and is available for download at http://www.nyu.edu/projects/jerschow. If you get it to compile on other systems, please let me know, I can then update the information here.

Features

Calculate the coherence pathways selected by a given phase cycle
Find an appropriate cogwheel phase cycle
Calculate the coherence pathways selected by cogwheel cycles
Multiple nuclei species supported but not tested, basically for lack of need. If you have an interesting problem, let me know, I can work with you if the program does not work in this case.
Invocation

CCCP++ < parameter_file_name.par

Parameter File Description

It will probably be easier to take an existing file and modify it to your needs, but here is some description:

Normal operation:

The parameters may be input in any sequence
Comments start by # and run until the end of the current line and may be placed anywhere in the text.
A scalar is input as name=number.
The phase cycle for one pulse is input as “ph=number number …..” and can span several lines. A next ph= statement specifies the next phase program (for the next pulse). The phase programs need not have equal length, they will be automatically filled up to the largest number of phases. The last phase designates the receiver phase. Phases of pulses on different nuclei can be input consecutively.
base indicates the basis of the phases that are input. base=n means that 0, 1, 2, … refer to 0, 360/n, 360/2n… degrees. base=n1 n2…. specifies different bases for pulses 1, 2,….
which_nucl specifies which of the pulses act on which nuclei. This vector determines the number of pulses in the sequence against which all other variables will be compared. If they are not compatible, an error message will be generated.
allow_paths specifies the allowed and forbidden coherence pathways. This is a matrix containing zeros and ones for the respective coherence order. It has size (2p+1)x(2n+1), where p is the maximum coherence order, and n is the number of pulses (best to look at an example). A new row of the matrix is started after each newline character. A new allow_paths= statement allows one to input the allowed and forbidden coherences for the next nuclear species (this one may have its own maximum coherence order and a different number of pulses, according to the which_nucl vector).
cutoff=number is the signal intensity that should be considered zero.
Cogwheel operation, to find a cogwheel cycle:

This mode will be chosen when the parameters COGpw, COGmax_pw, COGmax_N are set.
All parameters are the same except that the ph= statements are not needed.
COGpw= is a list of vectors specifying the desired coherence pathway. A new row specifies the pathway for the next nucleus.
COGmax_pw= is the maximum number of tolerated coherence pathways (including the desired one).
COGmin_N, COGmax_N, COGinc_N specify the minimum and the maximum of the cogwheel cycle sizes that should be tested, as well as the stepsize. If the values are not given, COGmin_N defaults to 2 and COGinc_N to 1.
COGwdg is a list of vectors that specifies the winding numbers for each pulse and each nuclear species. When used for finding a proper cogwheel cycle the only information that’s used is whether the winding number is zero or not. If it is zero, then the winding number for this pulse will not be searched (this may speed up the calculation significantly, and is justified in many cases as mentioned in the paper).
Cogwheel cycle testing:

This mode will be chosen when COGN and COGwdg are set.
instead of inputting the phase cycle by ph= one may specify a cogwheel cycle by:
COGN specifies the cogwheel cycle size.
COGwdg specifies the winding numbers for every pulse and every nucleus.
Note on Programming Style

Well, this program is written in C++ but does not use all the fancy features of encapsulation, inheritance, etc, mainly because I didn’t quite know what classes to define for coherence pathways, etc, and whether that would actually be of any advantage. The main motivation for writing it in C++ were the availability of the vector and list classes in the STL which make life much easier. (well, there are also these nice little things like ‘<<‘, etc).
