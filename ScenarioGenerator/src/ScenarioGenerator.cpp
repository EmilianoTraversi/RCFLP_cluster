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
#include <functional> //without .h


/* #include "timer.h" */
#include "options_sg.h"

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
int fType;              //!< instance type (1-4)
string instanceType;
double _epsilon;
int _seed = 0;
int _quantity = 1;

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
void printOptions(char * _FILENAME, INSTANCE inp, int timeLimit);
void GenerateDemand(int ind_seed);
/****************** FUNCTIONS DECLARATION ***************************/

/************************ main program ******************************/
/// main program
/************************ main program ******************************/
int main(int argc, char *argv[])
{

	int err = parseOptions(argc, argv);
	if (err != 0) exit(1);

	readProblemData(_FILENAME, fType, inp);
	printOptions(_FILENAME, inp, timeLimit);
	for (int ind_seed= _seed; ind_seed <_seed+_quantity; ind_seed++){
		GenerateDemand(ind_seed);
	}

	

	return 0;
}
/************************ main program ******************************/
/// END main program
/************************ main program ******************************/

/****************** FUNCTIONS DEFINITION ***************************/
void GenerateDemand(int ind_seed){

        string  s1      = string(_FILENAME);
        s1              = s1.substr(s1.find_last_of("\\/"), 100);
        string filename = "scenarios"+s1 + "_" + to_string((int)(_epsilon*1000)) +"_" + to_string(ind_seed) + ".box"; 

	ofstream fWriter(filename, ios::out);

	cout <<filename<< endl;
	
	

    // read OR Library instances
    if (fType == 1)
    {
	double somma_base = 0.0;
	double somma_nuovo = 0.0;
        //fWriter << filename<<endl;
        fWriter << inp.nF <<" "<< inp.nC<< endl;
        for (int i = 0; i < inp.nF; i++) fWriter << inp.s[i] <<" "<< inp.f[i]<< endl;
	double * domanda_temp =new double[inp.nC];
	mt19937 mt_rand(ind_seed);
        for (int j = 0; j < inp.nC; j++){
		double base_demand = (double)(inp.d[j]);

		int lb_dem = (int)ceil((1.0-_epsilon)* ( (double) base_demand));
		int ub_dem = (int)ceil((1.0+_epsilon)* ( (double) base_demand));

		double new_dem;
		new_dem = (double)(lb_dem+mt_rand()% (ub_dem-lb_dem));

//		cout << "*** " <<endl;
//		cout << "lb = " << lb_dem << " ub = " << ub_dem<< endl;
//		cout << "base_demand = " << base_demand<< endl;
//		cout << "new_dem = " << new_dem<< endl;
//		cout << "*** " << endl;
		

		//cout<< "inp.d[j] = " << inp.d[j] << " base_demand =" << base_demand << " coef = " << coef << " base_demand*coef = " << base_demand*coef << endl; 
		fWriter << new_dem<<" ";
		domanda_temp[j]=new_dem;

		somma_base+=base_demand;
		somma_nuovo+=ceil(new_dem);

	}

	cout << "somma_base = "<< somma_base << " somma_nuovo = " <<somma_nuovo << " ratio = " << somma_base/somma_nuovo << endl;
	
        fWriter<< endl;

        for (int i = 0; i < inp.nF; i++){
            for (int j = 0; j < inp.nC; j++) fWriter << inp.c[i][j]*domanda_temp[j]<<" ";
		fWriter<< endl;
        }

    }
    // read Avella instances
    else if (fType == 2)
    {
        cout <<"TODO AVELLA"<< endl;
    }
    else
        cout << "Problem type not defined (-t option). Use '-h' for help. " << endl;












	fWriter.close();
	


}
