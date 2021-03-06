/***************************************************************************
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*! \file options.cpp 
  \brief Read options from command line.

    Command line options are:

    - **-h** : help (visualize the list of options)

    - **-i** : inputfile

    - **-t** : instance type:
            -# orLibrary
            -# Avella

    - **-v** : version:
            -# Single-source
            -# Multi-source
            -# Ellipsoidal
            -# Polyhedral

    - **-u** : uncertainty set
            -# Box uncertainty set
            -# Budget uncertanty set

    - **-r** : read from disk
            -# 0 No: A new Budget set $B_l$ is generated and stored
            -# 1 Yes: The Budget set is read from disk
*/

#include <iostream>
#include <cstdlib>
/**********************************************************/
#define   _TIMELIMITdef  18000   //!< default wall-clock time limit
#define   _VERSIONdef    1      //!< single source by default
#define   _FROMDISKdef    1      //!< single source by default
/**********************************************************/

using namespace std;

extern char* _FILENAME; 	//!< name of the instance file
extern int timeLimit;		//!< wall-clock time limit
extern int fType;           //!< instance type (1-2)
extern int version;         //!< 1-SS; 2-MS; 3-Ellipsoidal; 4-Polyhedral
extern int support;         //!< 1-Box; 2-Budget
extern int readFromDisk;    //!< 0-No; (Generate a new Budget set B_l); 1-Yes
extern string instanceType;
extern string versionType;
extern string supportType;


extern double _Omega_input; //!< value for omega directly given as input parameter (it overwrites the one read in the parameter file)
extern double _epsilon_input; //!< value for epsilon directly given as input parameter (it overwrites the one read in the parameter file)
extern double _delta_input; //!< value for delta directly given as input parameter (it overwrites the one read in the parameter file)
extern double _gamma_input; //!< value for gamma directly given as input parameter (it overwrites the one read in the parameter file)
extern int    L_input;  //!< value for L directly given as input parameter (it overwrites the one read in the parameter file)


/// Parse command line options.
/** Use -h to visualize the list of options.
 */
int parseOptions(int argc, char* argv[])
{
   bool setFile = false;
   bool setType = false;
   timeLimit    = _TIMELIMITdef;
   version      = _VERSIONdef;
   readFromDisk= _FROMDISKdef;
   cout <<endl << "R-CLSP v1.0 " << endl;
   if (argc == 1)
   {
      cout << "No options specified. Try -h " << endl;
      return -1;
   }  

   _Omega_input= -1;
   _epsilon_input= -1;
   _delta_input= -1;
   _gamma_input= -1;
   L_input= -1;
 
   int i = 0;
   while (++i < argc)
   {
      const char *option = argv[i];
      if (*option != '-')
	 return i;
      else if (*option == '\0')
	 return i;
      else if (*option == '-')
      {
	 switch (*++option)
	 {
	    case '\0':
	       return i + 1;
	    case 'i':
	       _FILENAME = argv[i+1];
	       setFile = true;
	       i++;
	       break;
	    case 't':
	       fType = atol(argv[i+1]);
           setType = true;
	       i++;
	       break;
	    case 'l':
	       timeLimit = atol(argv[i+1]);
	       i++;
	       break;
	    case 'v':
	       version = atol(argv[i+1]);
	       i++;
	       break;
        case 'u':
	       support = atol(argv[i+1]);
	       i++;
	       break;
        case 'r':
	       readFromDisk = atol(argv[i+1]);
	       i++;
	       break;



        case 'o':
	       _Omega_input = atof(argv[i+1]);
	       i++;
	       break;
        case 'e':
	       _epsilon_input = atof(argv[i+1]);
	       i++;
	       break;
        case 'd':
	       _delta_input = atof(argv[i+1]);
	       i++;
	       break;
        case 'g':
	       _gamma_input = atof(argv[i+1]);
	       i++;
	       break;
        case 'L':
	       L_input = atof(argv[i+1]);
	       i++;
	       break;


	    case 'h':
	       cout << "OPTIONS :: " << endl;
	       cout << "-i : problem instance file" << endl;
	       cout << "-l : time limit (real)" << endl;
	       cout << "-v : problem version (1-SS; 2-MS; 3-SOCP; 4- Poly)" << endl;
	       cout << "-t : instance type (1-OR Library; 2-Avella)" << endl;
	       cout << "-u : support type (1-Box; 2-Budget)" << endl;
	       cout << "-r : read Budget support set from disk (0-No; 1-Yes)" << endl;
	       cout << endl;
	       return -1;
	 }
      }
   }
 
   if (setFile && setType)
   {
        if (fType == 1)
            instanceType = "OR Library";
        else if (fType == 2)
            instanceType = "Avella";
        else
            instanceType = "***";

        if (version == 1)
            versionType = "Single-source-Nominal";
        else if (version == 2)
            versionType = "Multi-source-Nominal";
        else if (version == 3)
            versionType = "Multi-source-Ellipsoidal";
        else if (version == 4)
            versionType = "Polyhedral Uncertainty (Wd <= h)";

        if (support == 1)
            supportType = "Box Uncertainty Set";
        else if (support == 2)
            supportType = "Budget Uncertainty Set";

        return 0;
   }
   else
   {
      cout <<"Options -i and -t are mandatory. Try ./rcflp -h" << endl;
      return -1;
   }
}
