/* Definition of Input Function */
double radxx(double, double, double);
double radyy(double, double, double);
double radzz(double, double, double);
double radxy(double, double, double);
double radyz(double, double, double);
double radxz(double, double, double);
int     double2int(double);
int     arrondi(double);
int     imp2i(int, int, int, int);
int     i2imp(int, int, int);
int     i2icpu(int, int, int);

double staggardv4(double, double, double,
	double, double, double, double,
	double, double, double, double,
	double, double, double, double);
double staggards4(double, double, double, double,
        double, double, double, double,
        double, double, double, double,
        double, double, double, double);
double staggardt4(double, double, double,
        double, double, double, double,
        double, double, double, double);

double staggardv2(double, double, double,
	double, double,  double, double, double, double);
double staggards2(double, double, double, double,
        double, double,  double, double, double, double);
double staggardt2(double, double, double,
        double, double, double, double);


double PMLdump(double, double, double, double,
	        double, double, double );
void dumpfactor(int, int, int, double, double *, double * );

double hardrock(double);
double softrock(double);
