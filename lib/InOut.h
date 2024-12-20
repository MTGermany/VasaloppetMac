#ifndef IN_OUT_H
#define IN_OUT_H

#include <vector>


//#######################################################
// Check if input of arbitrary type T is a number (int,float, double etc)
//!!! does not work in the usual InOut.cpp, InOut.h style!! 

/*
template<typename T>
bool isNumber(T x){
   std::string s;
   std::stringstream ss; 
   ss << x;
   ss >>s;
   if(s.empty() || std::isspace(s[0]) || std::isalpha(s[0])) return false ;
   char * p ;
   double result=strtod(s.c_str(), &p);// use only &p of strtod=string2digit, 
                         // not return val of double => ignore warning
   if(false){cout<<"Result of strtod="<<result<<endl;}// just stop the whining
   return (*p == 0) ;
}
*/


//#######################################################


/**
   Provides various helper methods for input and output.  Examples are
   <ul>
   <li> reading columns from a file
   <li> reading a nXm matrix from a file
   <li> writing output to a file
   </ul>
   Furthermor, it defines the chars '%' and '#' as comments in any input file.
*/


//(jun 04) ACHTUNG: InOut gibt nun n_array als Zahl der Eintraege raus!
// max. Index also n_array-1  !!

class InOut{

 public:

  InOut();
  static const int NXMAX=1000;
  static const int NYMAX=201;
  static const int NTMAX=1441;// 2501;//1441;

  static const int LINEMAX=5000;
  static const int COLMAX=50;
  static const int MAXSTR=256;
  static const char COMMENTCHAR = '%';
  static const char COMMENTCHAR2 = '#';

  bool fileExists(const char* fname);
  int get_colnumber(const char* fname, string header);
  int getNumberOfDataLines(const char* fname);
  int getNumberOfDataLines(string fname);
  int getNumberOfLines(const char* fname);
  int getNumberOfCols(const char* fname);

  string getFirstLine(const char* fname);
  string getFirstDataLine(const char* fname);
  string getLastDataLine(const char* fname);
	
	
  void  get_array (char* fname, int & n_array, double array1[]);
  void  get_array (char* fname, int & n_array, int array1[]);
  void  get_array (char* fname, int & n_array, 
  		   double array1[], double array2[]); 
  void  get_array (char* fname, int & n_array, 
		   double array1[], double array2[], double array3[]);

  void  get_array (char* fname, int & n_array, 
		   int array1[], double array2[], double array3[]);

  //MT2018-12: string as filename parameter just does not work!! SUCK!!

  void  get_col_last(const char* fname, int col, int nlast,
		      int & n_array, int array1[]);
  void  get_col(const char* fname, int col,  int & n_array, int array1[]);
  void  get_col(const char* fname, int col,  int & n_array, int array1[],
		 bool log);
  void  get_col(const char* fname, int col,  int & n_array, double array1[]);
  void  get_col(const char* fname, int col,  int & n_array, double array1[],
		 bool log);
  void  get_col(const char* fname, int col,  int & n_array, string array1[],
		 bool log);
  void  get_col(const char* fname, int col,  int & n_array, string array1[]);

  void  getChar_col(const char* fname, int col,  int & n_array, char array1[]);

  //MT2023-02: get_col with vector argument

  void  get_col(const char* fname,int col,int & n_array,vector<double>& data);
  void  get_col(const char* fname,int col,int & n_array,vector<int>& data);
  void  get_col(const char* fname,int col,int & n_array,vector<string>& data);

  vector<double> get_col(string fname,int col);
  vector<int> get_intcol(string fname,int col);
  vector<string> get_stringcol(string fname,int col);
  

  
  void testHallo();
  
  void getGridnumbers_n1n2(char* fname, int & n1, int & n2);

  void get_matrix (char* fname, int n1, int n2, bool swap, 
		   int icol, double** matrix);

  
  void  get_array2d (char* fname, 
                     double unitTimeInDat_s, double unitSpaceInDat_m,
                     int & nt, int & nx, double & dt, double & dx, 
                     double array1[][NTMAX+1],
                     double array2[][NTMAX+1],
		     double array3[][NTMAX+1]);

  
  void  get_array2d_col (const char* fname, int col, 
		       double& xmin, double& dx, int& nx, 
		       double& ymin, double& dy, int& ny, 
		       double array[NXMAX][NYMAX]);



  void  write_array (char* fname, int nData, const double data_col1[],
                      char* name_col1);
  void  write_array (char* fname, int nData, const double data_col1[],
                     const double data_col2[], 
                     char* name_col1, char* name_col2);
  void  write_array (char* fname, int nData, const int data_col1[],
                     const int data_col2[], 
                     char* name_col1, char* name_col2);
  //<new apr19>
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     char* titleString);
  void  write_array (char* fname, int nData, 
                     const int data_col1[], const double data_col2[], 
                     char* titleString);
  //</new>
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],
                     char* name_col1, char* name_col2, char* name_col3);
  //<new aug16>
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],
                     char* titleString);
  //</new>

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4);
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const int data_col2[], 
                     const double data_col3[],const double data_col4[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4);
  //<new nov17>
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     char* titleString);
  //</new>


  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5);

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const int data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5);
  //<new nov17>
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],
                     char* titleString);
  //</new>

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6);
  //<new nov17>
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     char* titleString);
  //</new>

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7);
  //<new mar17>
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],
                     char* titleString);

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8);


  //<new mar17>
 void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     char* titleString);

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     const double data_col9[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8, 
                     char* name_col9);

 //<new mar17>  
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     const double data_col9[],
                     char* titleString);

 void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     const double data_col9[],const double data_col10[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8, 
                     char* name_col9, char* name_col10);

  //<new mar17>
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     const double data_col9[],const double data_col10[],
                     char* titleString);

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     const double data_col9[],const double data_col10[],
                     const double data_col11[],const double data_col12[],
                     char* titleString);


  
  void  write_array_gen(char* fname, int nLines, int nCols,  
			const double data[][LINEMAX], const char* headers);

  void  write_array_gen(char* fname, int nLines, int nCols,  
			const int data[][LINEMAX], const char* headers);


  //########################################################
  // general-purpose function  write_array2d 
  // writes gnu-plottable file with col1=x value, col2=y value, 
  // col3=array=z(z,y) 
  //########################################################


  //HINT: usage of arrays when calling:
  // double array1[NX][inout.NTMAX];
  // where inout=instance of InOut (InOut.NTMAX not working) and
  // NX any const value <=NXMAX
  // !! heineous BUG: static const int NTMAX=InOut.NTMAX and then define
  // array1[NX][NTMAX] does NOT (!!!) work !!

  // see also write_array2d_fromarray1d

  void  write_array2d (char* fname,
		       double xmin, double xmax, int nx,
		       double ymin, double ymax, int ny,
		       double array[][NYMAX],
		       char* titleString);


//##########################################################
// general-purpose function  write_array2d_fromarray1d 
// writes gnu-plottable file with col1=x value, col2=y value, 
// col3=array=z(z,y) where z[i][j]=array1d[nx*i+j]
//##########################################################

  void  write_array2d_fromarray1d (char* fname,
		       double xmin, double xmax, int nx,
		       double ymin, double ymax, int ny,
		       double array1d[],
		       char* titleString);
  
  

  // special-purpose; usage of arrays when calling: see HINT above
  
  void  write_array2d (char* fname, double dxout, double dtout, 
                       double xmin, double xmax, 
                       double tmin, double tmax, 
                       double array1[NXMAX+1][NYMAX]);

  // special-purpose; usage of arrays when calling: see HINT above

  void  write_array2d (char* fname, int nt, int nx, double dt, double dx,
  double array1[NXMAX+1][NYMAX],
  double array2[NXMAX+1][NYMAX],
  double array3[NXMAX+1][NYMAX]);


  // special-purpose; usage of arrays when calling: see HINT above

  void  write_array2d (char* fname, double dxout, double dtout, 
    double xmin, double xmax, double tmin, double tmax,
                       double array1[][NYMAX],
                       double array2[][NYMAX],
                       double array3[][NYMAX],
		       double array4[][NYMAX]);

  void  write_array2d (char* fname, double dxout, double dtout, //!! bisher int
		       double xmin, double xmax, double tmin, double tmax,
                       double array1[][NYMAX+1],
                       double array2[][NYMAX+1],
                       double array3[][NYMAX+1],
		       double array4[][NYMAX+1], bool data3d);


  static void getvar(FILE *fp, int* pint);
  static void getvar(FILE *fp, long* plong);
  static void getvar(FILE *fp, double* pdouble);

  int is_not_white(char line[], int nchar);
  bool is_data_line(char line[]);
  bool is_data_line(string line);
  void prepare_line(char line[]);
  void prepare_line(string & line);


};

#endif // IN_OUT_H
