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

/*! \mainpage Robust Capacitated Facility Location Problem

  Description here.

  \authors
  \version v. 1.0.0
  \date Begins: 25.05.18
  \date Ends:

  The project is compiled using a makefile and run via command line:
  ~~~
  make
  ./bin/ScenarioProfile
  ~~~

  This function generates new scenario instances based on a nominal instance 


*/

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <limits> 

#include <chrono>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <string>
#include <cmath>
#include <climits>
#include <cfloat>
#include <cassert>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <sstream>

/* #include "timer.h" */
#include "options_se.h"

using namespace std;

double INFTY = std::numeric_limits<double>::infinity();
const long _MAXRANDOM  = 2147483647;
const double ZERO      = 0.0e0;
const double EPSI      = 0.00001;
const int seed = 27; //!< Initialization of Random Number Generator
mt19937_64 gen(seed); //!< 64-bit Mersenne Twister by Matsumoto and Nishimura, 2000
/* mt19937_64 gen(seed()); */

/****************** VARIABLES DECLARATION ***************************/
char * _FILENAME;		//!< Instance name file
char * _SOLNAME;		//!< Instance name file
char * _OUTNAME;
int fType;              //!< instance type (1-4)
string instanceType;

/// Structure used to define the instance data
// NOTE: Change the same structure in the file inout.cpp !!!
struct INSTANCE {
    int nF;        //!< Number of facilities
    int nC;        //!< Number of customers
    double  *f;    //!< Fixed costs
    double  *s;    //!< Capacity
    double  *d;    //!< Demand
    double **c;    //!< Allocation costs
    double   totS; //!< Total supply
    double   totD; //!< Total demand

    int     nR;    //!< Number of constraints polyhedron uncertainty set
    double  *h;    //!< Rhs of polyhedron definining support
    int     *W;    //!< Matrix W in column major format
    int *index;    //!< Index of column major format for w
    int *start;    //!< Starting position for elements of column j

};
INSTANCE inp; //!< Instance data

/// Optimal solution and obj function value
struct SOLUTION {
    int nOpen;
    int  *ySol;
    double **xSol;
    double zStar;
    IloAlgorithm::Status zStatus;
    IloNum startTime;
    IloNum cpuTime;
};
SOLUTION opt; //!< Solution data structure


/**** CPLEX DEFINITION ****/
typedef IloArray <IloNumVarArray> TwoD;
IloEnv env;
IloModel model(env, "cflp");
/* IloCplex cplex(model); */
TwoD x_ilo;
IloNumVarArray y_ilo;
IloNumVarArray q_ilo;
IloNumVar w_ilo;
TwoD psi_ilo;
IloNumVarArray u_ilo;
IloNumVarArray delta_ilo;
int solLimit     = 9999;
int displayLimit = 4;
int timeLimit;

/****************** FUNCTIONS DECLARATION ***************************/
int readProblemData(char * _FILENAME, int fType, INSTANCE & inp);
int readSolution(char * _SOLNAME, SOLUTION & opt, INSTANCE & inp);
void printOptions(char * _FILENAME, char * _SOLNAME, INSTANCE inp, int timeLimit);
double ComputeValue(SOLUTION & opt, INSTANCE & inp);
//double ComputeInfeasibility(SOLUTION & opt, INSTANCE & inp);
void printSolution(char * _FILENAME, INSTANCE inp, SOLUTION opt, bool toDisk,int fullOutput);
/****************** FUNCTIONS DECLARATION ***************************/

/************************ main program ******************************/
/// main program
/************************ main program ******************************/
int main(int argc, char *argv[])
{ 
 
	int err = parseOptions(argc, argv);
	if (err != 0) exit(1);

	readProblemData(_FILENAME, fType, inp);
	readSolution(_SOLNAME, opt, inp);

	printOptions(_FILENAME,  _SOLNAME, inp, timeLimit);
	printSolution(_FILENAME,inp, opt, 0,0);

	double OFvalue;
	double Infeasibility;
	OFvalue = ComputeValue(opt,inp);
//	Infeasibility = ComputeInfeasibility(opt,inp);
	

	return 0;
}
/************************ main program ******************************/
/// END main program
/************************ main program ******************************/

/****************** FUNCTIONS DEFINITION ***************************/
double ComputeValue(SOLUTION & opt, INSTANCE & inp){

	double const_part = 0.0;
	double variable_part = 0.0;

//	cout << "d = "<< endl;
//	for (int j = 0; j < inp.nC; j++) cout << inp.d[j]<< " ";
//	cout << endl;

//	cout << "c ="<< endl;
//	for (int i = 0; i < inp.nF; i++) {for (int j = 0; j < inp.nC; j++) cout << inp.c[i][j]<< " "; cout << endl;}

//	cout << "inp.nF" << inp.nF<< endl;

//	for (int i = 0 ; i < inp.nF ; i++) cout<<"y_"<<i <<" : " << opt.ySol[i] <<endl;
//	for (int i = 0 ; i < inp.nF ; i++) cout<<"f_"<<i <<" : " << inp.f[i] <<endl;

	for (int i = 0 ; i < inp.nF ; i++) const_part+= opt.ySol[i] * inp.f[i];


	for (int i = 0; i < inp.nF; i++) for (int j = 0; j < inp.nC; j++) variable_part += opt.xSol[i][j]* inp.d[j]*inp.c[i][j];

	cout << setprecision(15)<< const_part << " + " << variable_part << " = " << const_part + variable_part << endl;

	double * infeasibility_vector= new double[inp.nF];
	double infeasibility_max = 0; 
	for (int i = 0 ; i < inp.nF ; i++) infeasibility_vector[i] = opt.ySol[i] * inp.s[i];
	for (int i = 0; i < inp.nF; i++) for (int j = 0; j < inp.nC; j++) infeasibility_vector[i] -= opt.xSol[i][j]* inp.d[j];
	


	double infeasibility_tot = 0.0;

	for (int i = 0 ; i < inp.nF ; i++) {
		infeasibility_tot+= min(0.0,infeasibility_vector[i]); 
		if (infeasibility_max<-min(0.0,infeasibility_vector[i])) infeasibility_max=-min(0.0,infeasibility_vector[i]);
	}
	
	cout << "infeasibility_tot = "<< infeasibility_tot<< endl;
	cout << "infeasibility_max = "<< infeasibility_max<< endl;

//	for (int i = 0 ; i < inp.nF ; i++) cout <<opt.ySol[i]*inp.s[i]<<" -\t- " << infeasibility_vector[i]<< endl;

	cout << "inp.totD = "<< inp.totD<< endl;


	ofstream fWriter(_OUTNAME, ios::out);
	fWriter << _FILENAME << ";"<< _SOLNAME << ";" << const_part << ";" << variable_part << ";" << infeasibility_tot << ";" << infeasibility_max << endl;
	
	fWriter.close();

	return infeasibility_tot;

}



void printSolution(char * _FILENAME, INSTANCE inp, SOLUTION opt, bool toDisk, 
int fullOutput)
{
    cout << endl << "** SOLUTION **" << endl;
    cout << " ..z*     \t= " << setprecision(15) << opt.zStar << endl;
    cout << " ..time   \t= " << opt.cpuTime << endl;
    cout << " ..status \t= " << opt.zStatus << endl;

    ofstream fWriter("solution.txt", ios::out);
    fWriter << _FILENAME << "\t" << instanceType << "\t"
            << setprecision(15) << opt.zStar << "\t" 
            << opt.zStatus << "\t" << opt.cpuTime <<endl;
    fWriter.close();

    if (fullOutput >= 1)
    {
        cout << "Open Facilities: ";
        for (int i = 0; i < inp.nF; i++)
            if (opt.ySol[i] == 1)
                cout << setw(4) << i;
        cout << endl;
        if (fullOutput >= 2)
        {
            cout << "Allocation variables : " << endl;
            for (int i = 0; i < inp.nF; i++)
                for (int j = 0; j < inp.nC; j++)
                    if (opt.xSol[i][j] >= EPSI)
                        cout << "x(" << i << "," << j << ") = " << setprecision(3) 
                        << opt.xSol[i][j] << endl;
	    for (int i = 0; i < inp.nF; i++){
		for (int j = 0; j < inp.nC; j++) cout << setprecision(3) << opt.xSol[i][j] << " ";
		cout << endl;
	    }

        }
    }
}
