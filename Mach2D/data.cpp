#include "stdafx.h"

// GRID PARAMETERS AND VARIABLES
int nx;      //< Number of volumes in the csi direction (real+fictitious)
int ny;      //< Number of volumes in the eta direction (real+fictitious)
int nxy;     //< nxy = nx * ny
int kg;      //< Kind of grid (1=uniform, 2=geometric progression, 3=power law)
int kcm;     //< Kind of centroid mean (1=simple mean, 2=weighted mean)
int coord;   //< Kind of coord. system ( 1=cylindrical, 0 = cartesian)
double la;   //< Length of the elliptical x semi-axis (m)
double lb;   //< Length of the elliptical y semi-axis (m)
double a1;   //< Width of the volume closer to the wall (m)
double akn;  //< Exponent of the power law for the north boundary
double aks;  //< Exponent of the power law for the south boundary

double* x;      //< Coord. x at the northeast corner of the volume P (m)
double* y;      //< Coord. y at the northeast corner of the volume P (m)
double* xp;     //< Coord. x of the centroid of volume P (m)
double* yp;     //< Coord. y of the centroid of volume P (m)
double* xe;     //< x_eta at the center of east face of volume P (m)
double* ye;     //< y_eta at the center of east face of volume P (m)
double* xen;    //< x_eta at the center of north face of volume P (m)
double* yen;    //< y_eta at the center of north face of volume P (m)
double* xk;     //< x_csi at the center of north face  of volume P (m)
double* yk;     //< y_csi at the center of north face of volume P (m)
double* xke;    //< x_csi at the center of east face of volume P (m)
double* yke;    //< y_csi at the center of east face of volume P (m)
double* Jp;     //< Jacobian at the center of volume P (1/m2)
double* Je;     //< Jacobian at the center of east face of volume P (1/m2)
double* Jn;     //< Jacobian at the center of north face of volume P (1/m2)
double* alphae; //< (metric tensor component) alpha at the center of east face of volume P (m2)
double* gamman; //< (metric tensor component) gamma at the center of north face of volume P (m2)
double* betae;  //< (metric tensor component) beta  at the center of east face of volume P (m2)
double* betan;  //< (metric tensor component) beta  at the center of north face of volume P (m2)
double* re;     //< Radius at the center of east face of volume P (m)
double* rn;     //< Radius at the center of north face of volume P (m)
double* rp;     //< Radius at the center of volume P (m)

// NUMERIC PARAMETERS AND VARIABLES
int it;          //< Iteractions counter for time evolution cycle
int iti;         //< Iteractions counter for iterative correction cycle
int it1;         //< Number of iteractions up to which dt = dt1
int it2;         //< Number of iteractions from which dt = dt2
int itb1;        //< Number of iteractions up to which beta = beta1
int itb2;        //< Number of iteractions from which beta = beta2
int itm;         //< Iteractions counter for the mass cycle
int imax;        //< Maximum number of iteractions for mass cycle
int itmax;       //< Maximum number of iteractions for time evolution cycle
int itimax;      //< Maximum number of iteractions for the correction cycle
int nitm_u;      //< Maximum number of iteractions for solving the linear systems for u, v and T
int nitm_p;      //< Maximum number of iteractions for solving the linear system for p
int wlf;         //< Frequency of printing in the listing file
int sem_a;       //< 1 = do not open result files, 0 = open
int sem_g;       //< 0 = visualize the plot, 1 = do not visualize
int w_g;         //< Frequency of writing data for graphics
int w_cam;       //< 1 = write the fields, 0 = do not
int it_stop;     //< Number of iteractions before break time evolution cycle
int nit_res;     //< Number of iteractions to calculate the mean of the residuals
int res_check;   //< Checks if the it_stop was defined. 0 = false, 1 = true.
double tol_u;    //< Tolerance in the MSI for solving the linear systems for u, v and T
double tol_p;    //< Tolerance in the MSI for solving the linear system for p
double beta;     //< UDS/CDS mixing constant (0=UDS, 1=CDS)
double beta1;    //< Initial beta (UDS/CDS mixing constant (0=UDS, 1=CDS))
double beta2;    //< Final beta (UDS/CDS mixing constant (0=UDS, 1=CDS))
double dt;       //< Time step (s)
double dt1;      //< Initial time step (s)
double dt2;      //< Final time step (s)
double res_u;    //< Relative residual of the linear system for u
double res_v;    //< Relative residual of the linear system for v
double res_T;    //< Relative residual of the linear system for T
double res_p;    //< Residual of the linear system for p
double res;      //< Sum of all residuals of the linear systems for u, v, T and p.
double res_mean; //< Residual mean
double tol_res;  //< Tolerance for the sum of residuals
double norm;     //< Norm of the residuals of all linear systems
//double tcpu1;    //< First time measurement (s)
//double tcpu2;    //< Second time measurement (s)
//double tcpu;     //< Total cpu time (s)
double RAM;      //< RAM memory (MB)

double* bu;  //< Source of the linear system for u (N)
double* bv;  //< Source of the linear system for v (N)
double* bt;  //< Source of the linear system for T (J/s)
double* bp;  //< Source of the linear system for pl (kg/s)
double* pl;  //< Pressure deviation at center o volume P (Pa)
double* due; //< SIMPLEC coefficients for ue (m2.s/kg)
double* dve; //< SIMPLEC coefficients for ve (m2.s/kg)
double* dun; //< SIMPLEC coefficients for un (m2.s/kg)
double* dvn; //< SIMPLEC coefficients for vn (m2.s/kg)
double* de;  //< SIMPLEC coefficients for Uce (m3.s/kg)
double* dn;  //< SIMPLEC coefficients for Vcn (m3.s/kg)

double** au;  //< Coefficients of the linear system for u (kg/s)
double** av;  //< Coefficients of the linear system for v (kg/s)
double** at;  //< Coefficients of the linear system for T (J/(K.s))
double** ap;  //< Coefficients of the linear system for pl (m.s)
double** dl9; //< Lower matrix of the MSI method for 9 diagonals
double** du9; //< Upper matrix of the MSI method for 9 diagonals
double** dl5; //< Lower matrix of the MSI method for 5 diagonals
double** du5; //< Upper matrix of the MSI method for 5 diagonals

// GEOMETRIC PARAMETERS
double lr;  //< length of the body (m)
double rb;  //< base radius/semi-height of the body (m)

// FILE ID NUMBERS
int const tfid = 10;  //< Temporary file id
int const lid = 100;  //< Listing file id
int const rid = 101;  //< Residual file id

std::string date_s; //< System date
std::string time_s; //< System time
std::string sim_id; //< Simulation identification
std::string input_file_parameters; //< Input parameters data file

// GAS PROPERTIES AND VARIABLES
double GF; //< GF = gamma = Cp / Cv (for the free stream)  (dimensionless)
double Rg; //< Perfect gas constant (J/(kg.K))
double PF; //< Free stream pressure (Pa)
double TF; //< Free stream temperature (K)
double MF; //< Free stream Mach number (dimensionless)
double UF; //< Free stream speed (m/s)

double Cdfi; //< Pressure foredrag coefficient (dimensionless)

double* p;    //< Pressure at center o volume P (Pa)
double* pa;   //< Pressure of the previous time step (Pa)
double* T;    //< Temperature at center of volume P (K)
double* Ta;   //< Temperature of the previous time step (K)
double* ro;   //< Specific mass (absolute density) at center of volume P (kg/m3)
double* roa;  //< Specific mass of the previous time step (kg/m3)
double* roe;  //< Specific mass at east face of volume P (kg/m3)
double* ron;  //< Specific mass at north face of volume P (kg/m3)
double* g;    //< ro / p (kg/J)
double* u;    //< Cartesian velocity at the center of volume P (m/s)
double* ua;   //< u of the previous time step (m/s)
double* v;    //< Cartesian velocity at the center of volume P (m/s)
double* va;   //< v of the previous time step (m/s)
double* ue;   //< Cartesian velocity u at center of east face of volume P (m/s)
double* un;   //< Cartesian velocity u at center of north face of volume P (m/s)
double* ve;   //< Cartesian velocity v at center of east face of volume P (m/s)
double* vn;   //< Cartesian velocity v at center of north face of volume P (m/s)
double* Uce;  //< Contravariant velocity U at east face of volume P (m2/s)
double* Vcn;  //< Contravariant velocity V at north face of volume P (m2/s)
double* cp;   //< Specific heat at const pressure (J/(kg.K))
double* gcp;  //< gcp = gamma = Cp/Cv at center of CV P (dimensionless)
double* cup;  //< Term of deferred correction for u (N)
double* cvp;  //< Term of deferred correction for v (N)
double* uea;  //< ue of the previous time step (m/s)
double* vea;  //< ve of the previous time step (m/s)
double* una;  //< un of the previous time step (m/s)
double* vna;  //< vn of the previous time step (m/s)
double* Ucea; //< Uce of the previous time step (m2/s)
double* Vcna; //< Vcn of the previous time step (m2/s)
double* Vcw;  //< Contravariant velocity V on west face of fictitious volumes of east boundary (m2/s)

// SOLVER
double* r;	  //< residual vector
double* w;	  //< vector used to solve linear system Lw = -r
double* z;	  //< vector used to solve linear system Uz =  w

void get_parameters()
{
	//> Reads the parameters from an input data file, which name is in ./mach2d_input/mach2d_input.txt.
	//// Input and output variables are not explicited.

	//time_t now = time(0);
	//char* dt = ctime(&now);
	//std::cout << "The local date and time is: " << dt << std::endl;
	
	//std::time_t const now_c = std::time(NULL);
	//auto s = std::put_time(std::localtime(&now_c), "%F %T");
	//std::cout << s << std::endl;

	size_t strpos;
	std::ifstream inputFile;
	std::string line;
	std::vector<std::string> lines;
	inputFile.open("./mach2d_input/mach2d_input.txt");
	// Reading file name from which parameters will be read

	if (inputFile.is_open())
	{
		getline(inputFile, line);
		strpos = line.find(" ");
		line = line.substr(0, strpos);
		line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
		//std::cout << line << '\n';
	}
	else std::cout << "Unable to open file" << std::endl;
	inputFile.close();

	input_file_parameters = line;

	inputFile.open("./mach2d_input/" + input_file_parameters);
	if (inputFile.is_open())
	{
		while (getline(inputFile, line))
		{
			strpos = line.find(" ");
			line = line.substr(0, strpos);
			line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
			//std::cout << line << '\n';
			lines.push_back(line);
		}
	}
	else std::cout << "Unable to open file" << std::endl;
	inputFile.close();

	sim_id = lines[0];
	nx = std::stoi(lines[1]);
	ny = std::stoi(lines[2]);
	la = std::stod(lines[3]);
	lb = std::stod(lines[4]);
	lr = std::stod(lines[5]);
	rb = std::stod(lines[6]);
	kg = std::stoi(lines[7]);
	kcm = std::stoi(lines[8]);
	coord = std::stoi(lines[9]);
	a1 = std::stod(lines[10]);
	akn = std::stod(lines[11]);
	aks = std::stod(lines[12]);
	itmax = std::stoi(lines[13]);
	itimax = std::stoi(lines[14]);
	imax = std::stoi(lines[15]);
	it1 = std::stoi(lines[16]);
	it2 = std::stoi(lines[17]);
	dt1 = std::stod(lines[18]);
	dt2 = std::stod(lines[19]);
	nitm_u = std::stoi(lines[20]);
	nitm_p = std::stoi(lines[21]);
	nit_res = std::stoi(lines[22]);
	tol_u = std::stod(lines[23]);
	tol_p = std::stod(lines[24]);
	tol_res = std::stod(lines[25]);
	wlf = std::stoi(lines[26]);
	sem_a = std::stoi(lines[27]);
	sem_g = std::stoi(lines[28]);
	w_g = std::stoi(lines[29]);
	w_cam = std::stoi(lines[30]);
	itb1 = std::stoi(lines[31]);
	itb2 = std::stoi(lines[32]);
	beta1 = std::stod(lines[33]);
	beta2 = std::stod(lines[34]);
	GF = std::stod(lines[35]);
	Rg = std::stod(lines[36]);
	PF = std::stod(lines[37]);
	TF = std::stod(lines[38]);
	MF = std::stod(lines[39]);

	nx = nx + 2; // Real + Fictitious
	ny = ny + 2; // Real + Fictitious

	nxy = nx * ny;
}

void write_parameters(std::ofstream& outputFile)
{
	//> Writes the parameters, read from the input data file, to the file.
	//// Input and output variables are not explicited.
	if (outputFile.is_open())
	{
		outputFile << "" << std::endl;
		outputFile << "Date: " << std::endl;
		outputFile << "Time: " << std::endl;
		outputFile << "" << std::endl;
		outputFile << "PARAMETERS" << std::endl;
		outputFile << "" << std::endl;
		outputFile << sim_id << std::endl;
		outputFile << nx << std::endl;
		outputFile << ny << std::endl;
		outputFile << la << std::endl;
		outputFile << lb << std::endl;
		outputFile << lr << std::endl;
		outputFile << rb << std::endl;
		outputFile << kg << std::endl;
		outputFile << kcm << std::endl;
		outputFile << coord << std::endl;
		outputFile << a1 << std::endl;
		outputFile << akn << std::endl;
		outputFile << aks << std::endl;
		outputFile << itmax << std::endl;
		outputFile << itimax << std::endl;
		outputFile << imax << std::endl;
		outputFile << it1 << std::endl;
		outputFile << it2 << std::endl;
		outputFile << dt1 << std::endl;
		outputFile << dt2 << std::endl;
		outputFile << nitm_u << std::endl;
		outputFile << nitm_p << std::endl;
		outputFile << nit_res << std::endl;
		outputFile << tol_u << std::endl;
		outputFile << tol_p << std::endl;
		outputFile << tol_res << std::endl;
		outputFile << wlf << std::endl;
		outputFile << sem_a << std::endl;
		outputFile << sem_g << std::endl;
		outputFile << w_g << std::endl;
		outputFile << w_cam << std::endl;
		outputFile << itb1 << std::endl;
		outputFile << itb2 << std::endl;
		outputFile << beta1 << std::endl;
		outputFile << beta2 << std::endl;
		outputFile << GF << std::endl;
		outputFile << Rg << std::endl;
		outputFile << PF << std::endl;
		outputFile << TF << std::endl;
		outputFile << MF << std::endl;
	}
	else std::cout << "Unable to open file" << std::endl;
}

void allocate_initialize_variables()
{
	//> Reads the parameters from an input file, allocates the vectors and initialise them.
	//// Input and output variables are not explicited.

	// Reading data from input data file
	get_parameters();

	// Openning listing file
	std::ofstream outputFile;
	CreateFolder("./mach2d_output");
	outputFile.open("./mach2d_output/mach2d_" + sim_id + ".txt", std::fstream::in | std::fstream::out | std::fstream::trunc);

	if (outputFile.is_open())
	{
		outputFile << "" << std::endl;
		outputFile << "LISTING FILE OF MACH2D" << std::endl;
		outputFile << "" << std::endl;

		write_parameters(outputFile);

		outputFile.close();
	}
	else std::cout << "Unable to open file" << std::endl;

	//allocate_variables();
	
	Vcw = new double[ny]();

	x = new double[nxy];
	y = new double[nxy];
	xp = new double[nxy]();
	yp = new double[nxy]();
	xe = new double[nxy]();
	ye = new double[nxy]();
	xen = new double[nxy]();
	yen = new double[nxy]();
	xk = new double[nxy]();
	yk = new double[nxy]();
	xke = new double[nxy]();
	yke = new double[nxy]();
	Jp = new double[nxy]();
	Je = new double[nxy]();
	Jn = new double[nxy]();
	alphae = new double[nxy]();
	gamman = new double[nxy]();
	betae = new double[nxy]();
	betan = new double[nxy]();
	re = new double[nxy];
	rn = new double[nxy];
	rp = new double[nxy];
	p = new double[nxy];
	T = new double[nxy];
	ro = new double[nxy];
	roe = new double[nxy];
	ron = new double[nxy];
	u = new double[nxy];
	v = new double[nxy]();
	ue = new double[nxy];
	un = new double[nxy];
	ve = new double[nxy]();
	vn = new double[nxy]();
	Uce = new double[nxy];
	Vcn = new double[nxy];
	cp = new double[nxy];
	gcp = new double[nxy];
	ua = new double[nxy];
	va = new double[nxy];
	cup = new double[nxy];
	cvp = new double[nxy];
	g = new double[nxy];
	bu = new double[nxy];
	bv = new double[nxy];
	roa = new double[nxy];
	due = new double[nxy];
	dve = new double[nxy];
	dun = new double[nxy];
	dvn = new double[nxy];
	de = new double[nxy];
	dn = new double[nxy];
	uea = new double[nxy];
	vea = new double[nxy];
	una = new double[nxy];
	vna = new double[nxy];
	Ucea = new double[nxy];
	Vcna = new double[nxy];
	bp = new double[nxy];
	pl = new double[nxy]();
	Ta = new double[nxy];
	bt = new double[nxy];
	pa = new double[nxy];

	r = new double[nxy];
	w = new double[nxy];
	z = new double[nxy];

	// Vectors of vectors: reserve/resize just on inner dimension
	// VERIFY, check if replaces the global vars
	au  = new double*[nxy];
	av  = new double*[nxy];
	at  = new double*[nxy];
	ap  = new double*[nxy];
	dl9 = new double*[nxy];
	du9 = new double*[nxy];
	dl5 = new double*[nxy];
	du5 = new double*[nxy];

	for (int i = 0; i < nxy; i++)
	{
		au[i]  = new double[9];
		av[i]  = new double[9];
		at[i]  = new double[9];
		ap[i]  = new double[9];
		dl9[i] = new double[5];
		du9[i] = new double[4];
		dl5[i] = new double[4];
		du5[i] = new double[3];
	}
}

void CreateFolder(std::string path)
{
	if (!CreateDirectory(path.c_str(), NULL))
	{
		return;
	}
}