
#include <iostream>
#include <cstdlib>
/**********************************************************/
#define   _TIMELIMITdef  3600   //!< default wall-clock time limit
#define   _VERSIONdef    1      //!< single source by default
#define   _FROMDISKdef    1      //!< single source by default
/**********************************************************/

using namespace std;

extern char* _FILENAME; 	//!< name of the instance file
extern char* _SOLNAME;
extern string instanceType;
extern char* _OUTNAME;
extern int fType;           //!< instance type (1-2)


int parseOptions(int argc, char* argv[]) 
{
   bool setFile = false;
   bool setType = false;

   cout <<endl << "R-CLSP v1.0 " << endl;
   if (argc == 1)
   {
      cout << "No options specified. Try -h " << endl;
      return -1;
   }  

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
	    case 'o':
	       _OUTNAME = argv[i+1];
	       setFile = true;
	       i++;
	       break;
	    case 's':
	       _SOLNAME = argv[i+1];
	       setFile = true;
	       i++;
	       break;
	    case 't':
	       fType = atol(argv[i+1]);
               setType = true;
	       i++;
	       break;
	    case 'h':
	       cout << "OPTIONS :: " << endl;
	       cout << "-i : problem instance file" << endl;
	       cout << "-o : output name" << endl;
	       cout << "-s : solution file" << endl;
	       cout << "-t : instance type (1-OR Library; 2-Avella)" << endl;
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

        return 0;
   }
   else
   {
      cout <<"Options -i and -t are mandatory. Try ./rcflp -h" << endl;
      return -1;
   }
}
