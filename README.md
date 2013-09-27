README: HOLO_MULTIP version 1.0 C++

AUTHOR:	Marion Baumgartner

	This program is based on an erlier wesion writen in Fortran 77 by
		A.Stuck, Univ. Freiburg, Switzerland
		J. Osterwalder, Univ. Freiburg, 
		J. Wider (jw), Univ. Zurich, Switzerland
		A. Muntwiler (mm), Univ. Zurich, Switzerland
THANKS:

CHANGELOG:

	The program uese the C++ STD library as well as the Boost library (for calculation of special functions such asthe 		associated legendre Polynomials). For a detailed descriprion of Boost see http://www.boost.org/.

NEWS:


COMPILATION/INSTALL:

	for compilation call make
		-> compilation with g++
		-> using the C++98 standart

COPYING / LICENSE:

BUGS: Report Bugs to ...

MAIN IDEAS:

	1) read measured data from a data file formated as following
		#Intensity	theta	phi-angle
	2) from the measured intensity calculate a maximum amount of coefficients a_lm of the multipole epansion.

	3) from the coefficiets determine the intensity values at every theta and phi

	4) write claculated data to a file using a specified format

	
	special options:
		Calculate a Fermi functin like apodizing function for input data 
		
		calculate chi function from the measured hologram ...



