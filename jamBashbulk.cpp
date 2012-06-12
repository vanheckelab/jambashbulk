//g++ -O3 -o jam2d jamBashbulk2.cpp

// std libraries
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

using namespace std;

// definitions:
//ascii code for the escape key
#define ESCAPE 27
//permutation procedure
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
//loops
#define iloop(upperbound) for(int i = 0; i < upperbound; i++)
#define jloop(upperbound) for(int j = 0; j < upperbound; j++)

// global variables and constants:

// simulation parameters
int N = -1; // number of particles (per unit cell) in the packing
int Ncorrected;
const long double k = 1.0; // spring const. with respect to particle overlap

bool screenOutput = false;

char stop; // this variable waits for user input and halts the
// computation.

const bool report = false;
const bool test = true; // if true debugging output is shown
// if(test) cout << "test0000" << endl;

// graphics and output

// degrees of freedom of the periodic boundary unit cell in
long double alpha = 0.0; // lattice vector angle (simple shear)
long double delta = 0.0; // lattice vector aspect-ratio (pure shear)
long double L = 5.0; // lattice vector length
long double shear = 0.0;

bool alphaOnOff = false, deltaOnOff = false, pressOnOff = false;
bool dofOnOff = false;
int dofcriterium = 0;
// the previous parameters indicate in which way the simulated region
// can be deformed (these values change during execution of the
// program):
// alphaOnOff allows simple shear (alpha).
// deltaOnOff allows pure shear (delta).
// pressOnOff allows volume changes to reach given target pressure.

bool alphaOnOffInit = true, deltaOnOffInit = true, pressOnOffInit = true;
// these parameters in indicate the user's choice of the degrees of
// freedom.

string nameOfWorkingDirectory = "";
int particleNumberLength = 0;

// start values:
long double alphainit = 0.0; // the user's choice for the initial shear angle
long double deltainit = 0.0; // the initial aspect-ratio
long double P0init = 0.0;
long double P0 = 0.0; // the target pressure
long double P0_mem = 0.0;
long double phiinit = 0.86; // the initial fill fraction (0.869 ^= P~0.01)

int countAlphaFlip = 0;
int countDeltaFlip = 0;
int countPressFlip = 0;

// straining:
bool doSimpleShear = false; // switches on/off forced simple shear
bool doCompression = false;
bool fixedStepSize = false;
long double alphaBeforeDeformation; // the shear angle of the relaxed packing

int goalNumberOfContactChanges = 10; // the total number of contact changes we are interested in
long double goalStrain = 0.1; // the strain range we are interested in
int fixedStepNumber; // number of fixed size strain steps;

// OpenGL variables and parameters:
int window; //the number of our GLUT window

// iteration counters in various algorithms:
int iterationcountSimStep = 0;
int iterationcountfire = 0;
int maxIterationCountFire = 2e6;

int iterationcountmnbrak = 0;
int iterationcountbrent = 0;
int iterationcountfrprmn = 0;
int iterationcountfrprmnCUMULATIVE = 0;
int totaliterationcount = 0;
bool frprmnconverged = false;
bool fireconverged = false;
bool shearconverged = false;
bool converged = false;
int menumode = 0;
int programmode = 0;
int programmodeOld = 0;
int distributioncase = 0;
bool endprogram = false;
char *filename;
string filenameString;
int currentPackingNumber = 0;
int numPackingsToProcess = 0;
int firstPackingNumber = 0;
int packingNameNumber = 0;
bool redo = false;
long double distanceCalcs = 0.0;
long double Rneighbor = 100.0;
long double RneighborFrprmnLast;
bool onlydisplay = false;

long double energyDiffStep = 1e10;
long double energyDiffStepOld = 1e10;

long double enthalpieDiffStep = 1e10;
long double enthalpieDiffStepOld = 1e11;

long double UhelperLastFunctionCall = 1e20;
long double Uold = 1e10;
long double energyBeforeDeformation = 0.0;
long double sxxBeforeDeformation, sxyBeforeDeformation, syxBeforeDeformation,
		syyBeforeDeformation; // stress components

long double energywrite[50000];
long double enthalpiewrite[50000];
long double sxywrite[50000];

long double H = 1e9;
long double Hold = 1e10;
long double HLastFunctionCall = 1e11;

long double U; // potential energy packing (per unit cell)
long double phi; // fill fraction
long double Z; // average number of neighbors
long double P; // pressure
long double Pold;
long double sxx, sxy, syx, syy; // stress components
long double G;
long double sxyOld;

// overlap and distances
vector<long double> dij; // NxN matrix of particle overlap
vector<long double> rij; // NxN matrix of particle center distance
vector<long double> xij; // NxN matrix of particle x-position difference
vector<long double> yij; // NxN matrix of particle y-position difference

vector<bool> neighbors;
vector<bool> trueneighbors;
vector<bool> trueneighborsOld;
vector<bool> trueneighborsLast;
vector<int> trueneighborChanges;
vector<int> numberOfDirectNeighbors;
vector<bool> isRattler;
vector<bool> wasRattler;

int consideredNeighborNumber = 0;
int trueneighborNumber = 0;
int trueneighborNumberOld = 0;
int cumulativeNeighborchanges;
bool updateNeighbors;
long double maxgradient;

//generalized coordinates and gradients
vector<long double> p; // Positions of particles
vector<long double> pLast; // Positions at end of last simulationstep()
vector<long double> xi; // Potential gradient
vector<long double> g; // gradient helper
vector<long double> h; // gradient helper
vector<long double> R; // Radii of particles
long double Rmax; // Maximum particle radius

// fire algorithm variables
vector<long double> v; // Effective velocity
vector<long double> M; // Effective masses
long double beta;
long double power;
const int Nmin = 5;
const long double finc = 1.1;
const long double fdec = 0.5;
const long double betastart = 0.1;
const long double fbeta = 0.99;
long double dt = 1e-1;
long double dtmaxinit = 1e-1;
long double dtmax = dtmaxinit;
long double dtmin = 0.0;
long double damp = 1.0;
long double dampalpha = 0.9;
long double dampdelta = 0.9;
long double damppress = 0.999;
long double CG_step = 1;

//unit cell properties
long double lxx, lxy; // x-/y- component of L_x (1st unit cell vector)
long double lyx, lyy; // x-/y- component of L_y (2nd unit cell vector)
long double lxabs, lyabs; // lengths of L_x and L_y
vector<int> nx; // NxN matrix for periodicity calculation
vector<int> ny; // NxN matrix for periodicity calculation

//helper variables for energy and gradient calculation of 'hypothetical' configuration
long double Uhelper; // energy
vector<long double> phelper; // particle positions (x_i,y_i)
vector<long double> xihelper; // potential gradient (dU/dx_i,dU/dy_i)
long double alphahelper, deltahelper, Lhelper, Phelper;
long double lxxhelper, lxyhelper;
long double lyxhelper, lyyhelper;

//long double gg; // gradient squared
long double gg, vv;

// mathematical and program constants
const long double PI = 3.141592653589793;
const long double angle2rad = PI / 180.0;
const long double tiny = 1.0e-50; // small constant to avoid e.g. divion by zero
const long double gold = 0.5 * (1.0 + sqrt(5.0)); // golden ratio = 1.618033988749894885
const long double glimit = 100.0; // maximum magnification for parabolic-fit step in function mbrak

const long double CGOLD = 0.3819660; // golden ratio for brent
const long double ZEPS = 1e-25; // for brent
int ITMAXBRENT = 50; // maximum of iterations in brent
int ITMAX = 12; // maximum of iterations in frprmn
const long double TOL = 1e-7; // tolerance passed to brent by linmin
const long double AMIN = 1e-7; // starting step in linmin

int ncom; // dimension in linmin + f1dim
vector<long double> pcom, xicom; // positions and gradient
long double (*nrfunc)(); // function place-holder
int iterfrprmn; // number of iterations performed in frprmn
long double fret; // return minimum value of frprmn + linmin
long double ftol = 1e-17; // tolerance passed to frprmn()
long double ftolFrprmnBeforeFIRE = 1e-2;
long double ftolFIRE = 1e-17; // tolerance passed to fire()
int endcount = 0;

time_t starttime, endtime;
long double timediff1 = 0, timediff2 = 0;

long double dU, dH;

// function declarations:
void execute();
void initializeSimulation();
void initializeArrays();
void initializeOpenGL(int argc, char **argv);
void simulationstep(); // this is where the simulation is performed
long double testfunction(long double value);
long double energyfunc();
void particledistance(int i, int j);
void resethelpervars();
long double energy();
void gradientcalc();
void mnbrak(long double *ax, long double *bx, long double *cx, long double *fa,
		long double *fb, long double *fc, long double (*func)(long double));
long double brent(long double ax, long double bx, long double cx,
		long double (*f)(long double), long double tol, long double *xmin);
long double SIGN(long double a, long double b);
void linmin(int n, long double *fret, long double (*func)());
long double f1dim(long double x);
void frprmn(int n, long double *fret, long double (*func)());
void calcSysPara();
void menu();
void readPositionFile();
void writePositionFile();
void writeMultiplePackings(string name);
void packIntoBoundaries();
void createFileName();
void fire();
void calcShearModulus();
void calcBulkModulus();
void checkNeighborChanges(int& addedcontacts, int& removedcontacts,
		int& neighborChanges, int& neighborChangesLast, string& contactChanges);
void resolveRattler(int rattlerElement);
string ParticleAndPressureString();
void extractNandP(string foldername);
void checkFolderName(string foldername);

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	if ((argc == 2) and (strcmp(argv[1], "-screen") == 0)) {
		screenOutput = true;
	}

	starttime = time(NULL);

	menu();

	while (!endprogram)
		execute(); // either this

	return 0;
} // end main()

void extractNandP(string foldername) {

	int i = 0;
	int digit = 0;
	N = 0;
	P0 = 0;

	string comp2 = foldername.substr(i, 2);

	size_t pos = foldername.rfind("Packings"); // Name of the folder containing all the other folders with packings

	i = pos + 9;

	int i0 = i;

	if (foldername[i] != 'N') {
		cout << "Unexpected filename ERROR 1 : couldn't find N!" << endl;
		return;
	} else {
		i++;
		while (foldername[i] != '~') {
			N = N * 10;
			N = N + (foldername[i] - 48);
			i++;
			if (i - i0 > 10) {
				cout << "Unexpected filename ERROR 2 : N has too much digits!"
						<< endl;
				return;
			}
		}

		particleNumberLength = i - i0;

		if (foldername[i + 1] != 'P') {
			cout << "Unexpected filename ERROR 3 : no P in title" << endl;
			return;
		}

		int j = i + 2;
		int power = 1;
		while (foldername[j] != 'e') {
			digit = digit * 10;
			digit = digit + foldername[j] - 48;
			power = power - 1;
			j++;
		}

		if (foldername[j + 1] == '-') {
			P0 = digit * 1.0
					* pow(10, -1.0 * (-power + foldername[j + 2] - 48));
		} else if (foldername[j + 1] == '+') {
			P0 = digit * 1.0
					* pow(10, +1.0 * (-power + foldername[j + 2] - 48));
		}

	}
} // extractNandP

string ParticleAndPressureString() {

	char buffer;
	string bufferString;
	long double number, teni, teni1;
	int digit;

	long double pressureHelper = P0;
	int pressureExponent = 0;
	bool leadingZero = true;

	string Name = "";
	string ParticleNumberString = "N";
	string PressureString1 = "P";
	string PressureString2 = "e-";

	int pressuredigit;

	number = N;
	for (int i = 5; i > 0; i--) {
		teni = pow(10, i * 1.0);
		teni1 = pow(10, i - 1.0);
		digit = fmod(number, teni) / (teni1);
		if (digit > 0)
			leadingZero = false;
		buffer = digit + 48;
		if (digit >= 0 && !leadingZero)
			ParticleNumberString.push_back(buffer);
	}

	while (pressureHelper * 1.01 < 1.0) {
		pressureExponent++;
		pressureHelper *= 10.0;
	}

	pressuredigit = pressureHelper / 1;

	buffer = pressureExponent + 48;

	PressureString1.push_back(pressuredigit + 48);
	PressureString2.push_back(buffer);

	Name.append(ParticleNumberString);
	Name.append("~");
	Name.append(PressureString1);
	Name.append(PressureString2);

	return Name;
} // ParticleAndPressureString

////////////////////////////////////////////////////////////////////////////////
// execute()
void execute() {
	bool goodfile;

	if (!endprogram) {
		if (!converged) {
			simulationstep();
		} else {

			if ((iterationcountfire > maxIterationCountFire)
					|| fabs((Phelper - P0) / P0) > 0.1 || Z < 3.5 || Z > 10)
				goodfile = false;
			else
				goodfile = true;

			endtime = time(NULL); // clock function runtime
			timediff1 = endtime - starttime;
			starttime = endtime;
			if (screenOutput)
				cout << "Total runtime is " << timediff1 << " seconds." << endl;

			frprmnconverged = false;
			fireconverged = false;
			converged = false;

			if (programmode == 5) {
				if (!doSimpleShear && !doCompression) {
					writePositionFile();
				}
				iterationcountSimStep = 0;
				iterationcountfire = 0;
				totaliterationcount = 0;
				iterationcountfrprmnCUMULATIVE = 0;
				if (goodfile) {
					if (doSimpleShear)
						if (!redo)
							calcShearModulus();
					if (doCompression)
						if (!redo)
							calcBulkModulus();
				}

				if (currentPackingNumber
						< numPackingsToProcess - 1 + firstPackingNumber) {
					currentPackingNumber++;
					redo = false;
					if (distributioncase == 2)
						initializeSimulation();
					if (distributioncase == 3 || distributioncase == 4)
						readPositionFile();
					converged = false;
					fireconverged = false;

				}
				else
					endprogram = true;
			}
			else
				endprogram = true;

		}
	} // if (!endprogram)

	return;
} // execute

////////////////////////////////////////////////////////////////////////
// calcShearModulus()
void calcShearModulus() {
	int programmodeOld = programmode;
	int num = 0;

	long double shearfactor = 10; //= sqrt(10); // for fast calculation = 10 else = sqrt(10)
	long double dstrain = goalStrain / fixedStepNumber;
	int numberOfContactChanges = 0; // the number of contact changes at a given moment
	bool reachedGoal = false;

	int addedContacts = 0;
	int removedContacts = 0;

	int neighborChanges = 0; // sum of the number of added and removed contacts at a certain time
	int neighborChangesOld = neighborChanges;
	int neighborChangesLast = 0;
	int neighborChangesLastCumulative = 0;
	bool pastContactChange = false;
	bool sufficientAccuracy = false;
	int numberOfDataPoints = 0;

	vector<long double> extraPositionArray;
	extraPositionArray.reserve(2 * N + 3);

	long double shearLast = 0.0, sxyLast = sxy;

	time_t rawtime;
	struct tm *timeinfo;
	char timebuffer[80];

	string contactChanges = "";
	ofstream outG;
	ofstream outLog;
	ofstream outFirst;
	ofstream eraseFile;
	string dataFileName = filenameString;
	string logFileName = filenameString;
	string GpositionFile = filenameString;
	string Appendix = "";

	long double goalStrainHelper = goalStrain;
	int goalStrainExponent = 0, goalStraindigit = 0;
	string goalStrainString = "";

	while (goalStrainHelper * 1.01 < 1.0) {
		goalStrainExponent++;
		goalStrainHelper *= 10.0;
	}

	goalStraindigit = goalStrainHelper / 1;

	goalStrainString.push_back(goalStraindigit + 48);
	goalStrainString.append("e-");
	goalStrainString.push_back(goalStrainExponent + 48);

	if (!fixedStepSize) {
		Appendix = "~SR";
		Appendix.push_back(((goalNumberOfContactChanges / 100) % 10) + 48);
		Appendix.push_back(((goalNumberOfContactChanges / 10) % 10) + 48);
		Appendix.push_back(((goalNumberOfContactChanges / 1) % 10) + 48);
		Appendix.append("~step");
		Appendix.push_back(
				(((2 * goalNumberOfContactChanges + 1) / 100) % 10) + 48);
		Appendix.push_back(
				(((2 * goalNumberOfContactChanges + 1) / 10) % 10) + 48);
		Appendix.push_back(
				(((2 * goalNumberOfContactChanges + 1) / 1) % 10) + 48);
	} else {

		Appendix = "~SS";
		Appendix.append(goalStrainString);
		Appendix.append("~step");
		Appendix.push_back((((fixedStepNumber) / 1000) % 10) + 48);
		Appendix.push_back((((fixedStepNumber) / 100) % 10) + 48);
		Appendix.push_back((((fixedStepNumber) / 10) % 10) + 48);
		Appendix.push_back((((fixedStepNumber) / 1) % 10) + 48);
	}

	GpositionFile.insert(particleNumberLength + 6, Appendix);
	GpositionFile.insert(0, nameOfWorkingDirectory + "/");

	dataFileName.insert(particleNumberLength + 6, Appendix);
	dataFileName.insert(0, "data");
	dataFileName.insert(0, nameOfWorkingDirectory + "/");

	logFileName.insert(particleNumberLength + 6, Appendix);
	logFileName.insert(0, "log");
	logFileName.insert(0, nameOfWorkingDirectory + "/");

	if (!fixedStepSize)
		shear = 1e-9; // for fast calculation = 1e-12 else = 1e-16
	else
		shear = 0.0;

	outG.open((char*) dataFileName.c_str(), ios::trunc);
	outG.setf(ios::scientific, ios::floatfield);
	outG.precision(16);
	outG << "eta	s_xy	Ncontacts	Nchanges	N+	N-	P	Z" << endl;
	outG.close();

	outLog.open((char*) logFileName.c_str(), ios::trunc);
	outLog.setf(ios::scientific, ios::floatfield);
	outLog.precision(16);
	outLog << "step#" << "	N" << "	P0" << "	P" << "	alpha" << "	delta";
	outLog << "	L" << "	phi" << "	Z" << "	#rattler" << "	s_xx" << "	s_yy";
	outLog << "	s_xy" << "	U" << "	dU" << "	H" << "	dH" << "	t_run";
	outLog << "	#FIRE" << "	#CG" << "	gg" << " creation-date" << endl;
	outLog.close();

	alphaBeforeDeformation = p[2 * N];
	iloop(N) {
		jloop(N) {
			trueneighborsLast[j * N + i] = trueneighborsOld[j * N + i] =
					trueneighbors[j * N + i]; // the contacts before shearing
			trueneighborChanges[j * N + i] = 0;
		}
		wasRattler[i] = isRattler[i];
	}

	programmode = 3; // program-mode for shearing

	cumulativeNeighborchanges = 0;

	jloop(2*N+3) {
		pLast[j] = p[j];
	} // backup the particle position before the shear step

	energyBeforeDeformation = Uhelper;
	sxxBeforeDeformation = sxx;
	sxyBeforeDeformation = sxy;
	syxBeforeDeformation = syx;
	syyBeforeDeformation = syy;

	packIntoBoundaries();

	eraseFile.open((char*) GpositionFile.c_str(), ios::trunc);
	eraseFile.close();
	writeMultiplePackings(GpositionFile);

	while (!reachedGoal) { // shear < 0.15 && numberOfDataPoints < 5000 &&
		neighborChangesLastCumulative += neighborChangesLast;

		iloop(N) {
			jloop(N) {
				trueneighborsLast[j * N + i] = trueneighbors[j * N + i]; // the contacts before shearing step

			}
			wasRattler[i] = isRattler[i];
		}

		neighborChangesOld = neighborChanges;

		while (!sufficientAccuracy) {
			if (!fixedStepSize) {
				if (pastContactChange) {
					num++;

					jloop(2*N+3) {
						extraPositionArray[j] = p[j];
					} // save positions just after rearrangement
					jloop(2*N+3) {
						p[j] = pLast[j];
					} // go back to last particle positions
					shear = shear / shearfactor; // go back to previous shear
					shearfactor = sqrt(shearfactor);

				} else if (num > 0) {
					shearfactor = sqrt(shearfactor);
				}

				shear = shear * shearfactor; // shear is increased by shearfactor...

			} else {
				shear = shear + dstrain;
			}

			jloop(2*N+3) {
				pLast[j] = p[j];
			} // backup the particle position before the shear step
			p[2 * N] = alphaBeforeDeformation + shear; // ... and added to the shear of the relaxed packing

			iterationcountSimStep = 0;
			iterationcountfrprmn = 0;
			iterationcountfire = 0;
			iterationcountfrprmnCUMULATIVE = 0;

			G = 0.0;

			cumulativeNeighborchanges = 0;

			shearconverged = false;
			converged = false;
			frprmnconverged = false;
			fireconverged = false;

			while (!shearconverged) {
				simulationstep();
			}
			calcSysPara();
			if (iterationcountfire > maxIterationCountFire) {
				programmode = programmodeOld;
				return;
			}

			checkNeighborChanges(addedContacts, removedContacts,
					neighborChanges, neighborChangesLast, contactChanges);

			if (neighborChangesLast != 0)
				pastContactChange = true;
			else
				pastContactChange = false;

			if (!fixedStepSize) {
				if (((neighborChanges - neighborChangesOld)
						* (neighborChanges - neighborChangesOld) == 1
						&& (shearfactor - 1.0) < 0.0001)
						|| ((shearfactor - 1.0) < 1e-6))
					sufficientAccuracy = true;
			}

			endtime = time(NULL); // clock function runtime
			timediff1 = endtime - starttime;
			starttime = endtime;

			// fast finding of 1st contact change:
			if ((numberOfDataPoints < 1
					&& (neighborChangesLastCumulative + neighborChangesLast)
							!= 0) && !fixedStepSize) {
				numberOfContactChanges = 0;
				shear = 1e-16;
				sufficientAccuracy = false;
				shearfactor = 10;
				num = 0;
				pastContactChange = false;
			} //////////////////////////////////////////

			numberOfDataPoints++;

			energy();
			calcSysPara();

			long double maxGrad = 0;
			iloop(2*N) {
				if (fabs(xihelper[i]) > maxGrad)
					maxGrad = fabs(xihelper[i]);
			}

			G = (sxy - sxyLast) / (shear - shearLast);

			time(&rawtime);
			timeinfo = localtime(&rawtime);
			strftime(timebuffer, 80, "%Y-%m-%d_%H-%M-%S", timeinfo);

			outLog.open((char*) logFileName.c_str(), ios::app);

			outLog << numberOfDataPoints << "	" << N << "	" << P0 << "	" << P
					<< "	" << alpha << "	" << delta;
			outLog << "	" << L << "	" << phi << "	" << Z << "	"
					<< N - Ncorrected << "	" << sxx << "	" << syy;
			outLog << "	" << sxy << "	" << Uhelper << "	" << dU << "	" << H
					<< "	" << dH << "	" << timediff1;
			outLog << "	" << iterationcountfire << "	"
					<< iterationcountfrprmnCUMULATIVE << "	" << maxGrad << "	"
					<< timebuffer << endl;

			outLog.close();

			outG.open((char*) dataFileName.c_str(), ios::app);
			outG.setf(ios::scientific, ios::floatfield);
			outG.precision(16);
			outG << shear << "	" << sxy << "	" << trueneighborNumber << "	"
					<< (neighborChangesLastCumulative + neighborChangesLast)
					<< "	" << addedContacts << "	" << removedContacts << "	"
					<< Phelper << "	" << Z << endl;

			outG.close();

			shearLast = shear;
			sxyLast = sxy;

			if (fixedStepSize) {
				if (numberOfDataPoints < fixedStepNumber) {
					reachedGoal = false;
					writeMultiplePackings(GpositionFile);
				}

				else {
					reachedGoal = true;
					sufficientAccuracy = true;
					writeMultiplePackings(GpositionFile);
					programmode = programmodeOld;

					return;
				}
			}
		} // end while(!sufficientAccuracy)

		jloop(2*N+3) {
			p[j] = pLast[j];
		}
		writeMultiplePackings(GpositionFile);

		jloop(2*N+3) {
			p[j] = extraPositionArray[j];
		}
		writeMultiplePackings(GpositionFile);

		if (!fixedStepSize) {
			pastContactChange = false;
			sufficientAccuracy = false;
			num = 0;
			shearfactor = sqrt(sqrt(sqrt(sqrt(10)))); //sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(10)))))));
			numberOfContactChanges++;
		}

		if (!fixedStepSize) {
			if (numberOfContactChanges < goalNumberOfContactChanges)
				reachedGoal = false;
			else
				reachedGoal = true;
		}
	} // end while

	programmode = programmodeOld;
	return;
} // end calcShearModulus

////////////////////////////////////////////////////////////////////////
// calcBulkModulus()
void calcBulkModulus() {

	int programmodeOld = programmode;
	int num = 0;

	long double shearfactor = 10; //= sqrt(10); // for fast calculation = 10 else = sqrt(10)
	long double dstrain = goalStrain / fixedStepNumber;
	int numberOfContactChanges = 0; // the number of contact changes at a given moment
	bool reachedGoal = false;

	int addedContacts = 0;
	int removedContacts = 0;

	int neighborChanges = 0; // sum of the number of added and removed contacts at a certain time
	int neighborChangesOld = neighborChanges;
	int neighborChangesLast = 0;
	int neighborChangesLastCumulative = 0;
	bool pastContactChange = false;
	bool sufficientAccuracy = false;
	int numberOfDataPoints = 0;

	vector<long double> extraPositionArray;
	extraPositionArray.reserve(2 * N + 3);

	long double shearLast = 0.0, sxyLast = sxy;

	time_t rawtime;
	struct tm *timeinfo;
	char timebuffer[80];

	string contactChanges = "";
	ofstream outG;
	ofstream outLog;
	ofstream outFirst;
	ofstream eraseFile;
	string dataFileName = filenameString;
	string logFileName = filenameString;
	string GpositionFile = filenameString;
	bool compress;
	string Appendix = "";

	long double goalStrainHelper = goalStrain;
	int goalStrainExponent = 0, goalStraindigit = 0;
	string goalStrainString = "";

	if (goalStrain < 0 || goalNumberOfContactChanges < 0) {
		compress = true;
		if (goalStrain < 0)
			goalStrainHelper *= -1.0;
		if (goalNumberOfContactChanges < 0)
			goalNumberOfContactChanges *= -1;
	} else {
		compress = false;

	}

	while (goalStrainHelper * 1.01 < 1.0) {
		goalStrainExponent++;
		goalStrainHelper *= 10.0;
	}

	goalStraindigit = goalStrainHelper / 1;

	goalStrainString.push_back(goalStraindigit + 48);
	goalStrainString.append("e-");
	goalStrainString.push_back(goalStrainExponent + 48);

	if (compress) {
		if (!fixedStepSize) {
			Appendix = "~CR";
			Appendix.push_back(((goalNumberOfContactChanges / 100) % 10) + 48);
			Appendix.push_back(((goalNumberOfContactChanges / 10) % 10) + 48);
			Appendix.push_back(((goalNumberOfContactChanges / 1) % 10) + 48);
			Appendix.append("~step");
			Appendix.push_back(
					(((2 * (goalNumberOfContactChanges) + 1) / 100) % 10) + 48);
			Appendix.push_back(
					(((2 * (goalNumberOfContactChanges) + 1) / 10) % 10) + 48);
			Appendix.push_back(
					(((2 * (goalNumberOfContactChanges) + 1) / 1) % 10) + 48);
		} else {

			Appendix = "~CS";
			Appendix.append(goalStrainString);
			Appendix.append("~step");
			Appendix.push_back((((fixedStepNumber) / 1000) % 10) + 48);
			Appendix.push_back((((fixedStepNumber) / 100) % 10) + 48);
			Appendix.push_back((((fixedStepNumber) / 10) % 10) + 48);
			Appendix.push_back((((fixedStepNumber) / 1) % 10) + 48);
		}
	} else {
		if (!fixedStepSize) {
			Appendix = "~DR";
			Appendix.push_back(((goalNumberOfContactChanges / 100) % 10) + 48);
			Appendix.push_back(((goalNumberOfContactChanges / 10) % 10) + 48);
			Appendix.push_back(((goalNumberOfContactChanges / 1) % 10) + 48);
			Appendix.append("~step");
			Appendix.push_back(
					(((2 * goalNumberOfContactChanges + 1) / 100) % 10) + 48);
			Appendix.push_back(
					(((2 * goalNumberOfContactChanges + 1) / 10) % 10) + 48);
			Appendix.push_back(
					(((2 * goalNumberOfContactChanges + 1) / 1) % 10) + 48);
		} else {

			Appendix = "~DS";
			Appendix.append(goalStrainString);
			Appendix.append("~step");
			Appendix.push_back((((fixedStepNumber) / 1000) % 10) + 48);
			Appendix.push_back((((fixedStepNumber) / 100) % 10) + 48);
			Appendix.push_back((((fixedStepNumber) / 10) % 10) + 48);
			Appendix.push_back((((fixedStepNumber) / 1) % 10) + 48);
		}
	}

	GpositionFile.insert(particleNumberLength + 6, Appendix);
	GpositionFile.insert(0, nameOfWorkingDirectory + "/");

	dataFileName.insert(particleNumberLength + 6, Appendix);
	dataFileName.insert(0, "data");
	dataFileName.insert(0, nameOfWorkingDirectory + "/");

	logFileName.insert(particleNumberLength + 6, Appendix);
	logFileName.insert(0, "log");
	logFileName.insert(0, nameOfWorkingDirectory + "/");

	if (!fixedStepSize) {
		if (compress)
			shear = -1e-9; // for fast calculation = 1e-12 else = 1e-16
		else
			shear = 1e-9;
	} else
		shear = 0.0;

	outG.open((char*) dataFileName.c_str(), ios::trunc);
	outG.setf(ios::scientific, ios::floatfield);
	outG.precision(16);
	outG << "eta_V	s_xy	Ncontacts	Nchanges	N+	N-	P	Z" << endl;
	outG.close();

	outLog.open((char*) logFileName.c_str(), ios::trunc);
	outLog.setf(ios::scientific, ios::floatfield);
	outLog.precision(16);
	outLog << "step#" << "	N" << "	P0" << "	P" << "	alpha" << "	delta";
	outLog << "	L" << "	phi" << "	Z" << "	#rattler" << "	s_xx" << "	s_yy";
	outLog << "	s_xy" << "	U" << "	dU" << "	H" << "	dH" << "	t_run";
	outLog << "	#FIRE" << "	#CG" << "	gg" << " creation-date" << endl;
	outLog.close();

	alphaBeforeDeformation = p[2 * N + 2];
	iloop(N) {
		jloop(N) {
			trueneighborsLast[j * N + i] = trueneighborsOld[j * N + i] =
					trueneighbors[j * N + i]; // the contacts before shearing
			trueneighborChanges[j * N + i] = 0;
		}
		wasRattler[i] = isRattler[i];
	}

	programmode = 3; // program-mode for shearing

	cumulativeNeighborchanges = 0;

	jloop(2*N+3) {
		pLast[j] = p[j];
	} // backup the particle position before the shear step

	energyBeforeDeformation = Uhelper;
	sxxBeforeDeformation = sxx;
	sxyBeforeDeformation = sxy;
	syxBeforeDeformation = syx;
	syyBeforeDeformation = syy;

	packIntoBoundaries();

	eraseFile.open((char*) GpositionFile.c_str(), ios::trunc);
	eraseFile.close();
	writeMultiplePackings(GpositionFile);

	while (!reachedGoal) {

		neighborChangesLastCumulative += neighborChangesLast;

		iloop(N) {
			jloop(N) {
				trueneighborsLast[j * N + i] = trueneighbors[j * N + i]; // the contacts before shearing step

			}
			wasRattler[i] = isRattler[i];
		}

		neighborChangesOld = neighborChanges;

		while (!sufficientAccuracy) {

			if (!fixedStepSize) {
				if (pastContactChange) {
					num++;

					jloop(2*N+3) {
						extraPositionArray[j] = p[j];
					} // save positions just after rearrangement
					jloop(2*N+3) {
						p[j] = pLast[j];
					} // go back to last particle positions
					shear = shear / shearfactor; // go back to previous shear
					shearfactor = sqrt(shearfactor);

				} else if (num > 0) {
					shearfactor = sqrt(shearfactor);
				}

				shear = shear * shearfactor; // shear is increased by shearfactor...

			} else {
				shear = shear + dstrain;
			}

			jloop(2*N+3) {
				pLast[j] = p[j];
			} // backup the particle position before the shear step
			p[2 * N + 2] = alphaBeforeDeformation * (1 + shear); // ... and added to the shear of the relaxed packing

			iterationcountSimStep = 0;
			iterationcountfrprmn = 0;
			iterationcountfire = 0;
			iterationcountfrprmnCUMULATIVE = 0;

			G = 0.0;

			cumulativeNeighborchanges = 0;

			shearconverged = false;
			converged = false;
			frprmnconverged = false;
			fireconverged = false;

			while (!shearconverged) {
				simulationstep();
			}
			calcSysPara();
			if (iterationcountfire > maxIterationCountFire) {
				programmode = programmodeOld;
				return;
			}

			checkNeighborChanges(addedContacts, removedContacts,
					neighborChanges, neighborChangesLast, contactChanges);

			if (neighborChangesLast != 0)
				pastContactChange = true;
			else
				pastContactChange = false;

			if (!fixedStepSize) {
				if (((neighborChanges - neighborChangesOld)
						* (neighborChanges - neighborChangesOld) == 1
						&& fabs(shearfactor - 1.0) < 0.0001)
						|| (fabs(shearfactor - 1.0) < 1e-6))
					sufficientAccuracy = true;
			}

			endtime = time(NULL); // clock function runtime
			timediff1 = endtime - starttime;
			starttime = endtime;

			// fast finding of 1st contact change:
			if ((numberOfDataPoints < 1
					&& (neighborChangesLastCumulative + neighborChangesLast)
							!= 0) && !fixedStepSize) {
				numberOfContactChanges = 0;
				if (compress)
					shear = -1e-16;
				else
					shear = 1e-16;
				sufficientAccuracy = false;
				shearfactor = 10;
				num = 0;
				pastContactChange = false;
			} //////////////////////////////////////////

			numberOfDataPoints++;

			energy();
			calcSysPara();

			long double maxGrad = 0;
			iloop(2*N) {
				if (fabs(xihelper[i]) > maxGrad)
					maxGrad = fabs(xihelper[i]);
			}

			G = (sxy - sxyLast) / (shear - shearLast);

			time(&rawtime);
			timeinfo = localtime(&rawtime);
			strftime(timebuffer, 80, "%Y-%m-%d_%H-%M-%S", timeinfo);

			outLog.open((char*) logFileName.c_str(), ios::app);

			outLog << numberOfDataPoints << "	" << N << "	" << P0 << "	" << P
					<< "	" << alpha << "	" << delta;
			outLog << "	" << L << "	" << phi << "	" << Z << "	"
					<< N - Ncorrected << "	" << sxx << "	" << syy;
			outLog << "	" << sxy << "	" << Uhelper << "	" << dU << "	" << H
					<< "	" << dH << "	" << timediff1;
			outLog << "	" << iterationcountfire << "	"
					<< iterationcountfrprmnCUMULATIVE << "	" << maxGrad << "	"
					<< timebuffer << endl;

			outLog.close();

			outG.open((char*) dataFileName.c_str(), ios::app);
			outG.setf(ios::scientific, ios::floatfield);
			outG.precision(16);
			outG << ((1.0 + shear) * (1.0 + shear) - 1.0) << "	" << sxy << "	"
					<< trueneighborNumber << "	"
					<< (neighborChangesLastCumulative + neighborChangesLast)
					<< "	" << addedContacts << "	" << removedContacts << "	"
					<< Phelper << "	" << Z << endl;

			outG.close();

			shearLast = shear;
			sxyLast = sxy;

			if (fixedStepSize) {
				if (numberOfDataPoints < fixedStepNumber) {
					reachedGoal = false;
					writeMultiplePackings(GpositionFile);
				}

				else {
					reachedGoal = true;
					sufficientAccuracy = true;
					writeMultiplePackings(GpositionFile);
					programmode = programmodeOld;

					return;
				}
			}
		} // end while(!sufficientAccuracy)

		jloop(2*N+3) {
			p[j] = pLast[j];
		}
		writeMultiplePackings(GpositionFile);

		jloop(2*N+3) {
			p[j] = extraPositionArray[j];
		}
		writeMultiplePackings(GpositionFile);

		if (!fixedStepSize) {
			pastContactChange = false;
			sufficientAccuracy = false;
			num = 0;
			shearfactor = sqrt(sqrt(sqrt(sqrt(10)))); //sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(10)))))));

			if (numberOfContactChanges < goalNumberOfContactChanges)
				reachedGoal = false;
			else
				reachedGoal = true;

			numberOfContactChanges++;
		}
	} // end while

	programmode = programmodeOld;
	return;
} // end calcBulkModulus

////////////////////////////////////////////////////////////////////////
// checkNeighborChanges()
void checkNeighborChanges(int& addedContacts, int& removedContacts,
		int& neighborChanges, int& neighborChangesLast,
		string& contactChanges) {
	stringstream out(stringstream::out);

	int trueneighborChangesLast = 0;

	addedContacts = 0; // added contacts (with respect to contacts of the relaxed packing)
	removedContacts = 0; // removed contacts (-"-)
	neighborChanges = 0; // total number of contact changes (-"-)
	neighborChangesLast = 0; // contact changes since last step

	calcSysPara();

	iloop(N) {
		jloop(i) {
			trueneighborChanges[j * N + i] = trueneighbors[j * N + i]
					- trueneighborsOld[j * N + i];
			trueneighborChangesLast = trueneighbors[j * N + i]
					- trueneighborsLast[j * N + i];

			if (!isRattler[j] && !isRattler[i]) {
				if (trueneighborChangesLast > 0) {
					neighborChangesLast++;
					addedContacts++;
				}
			}

			if (!wasRattler[j] && !wasRattler[i]) {
				if (trueneighborChangesLast < 0) {
					neighborChangesLast++;
					removedContacts++;
				}
			}
		}
	}
} // end checkNeighborChanges

////////////////////////////////////////////////////////////////////////
// Simulation step
void simulationstep() {
	long double fret;
	long double frac;
	int neighborIterationCorrection;

	time_t rawtime1;

	timediff2 = time(NULL);
	distanceCalcs = 0.0;

	alpha = p[2 * N];
	delta = p[2 * N + 1];
	L = p[2 * N + 2];

	lxx = L / (1.0 + delta);
	lxy = L * 0.0;
	lyx = L * alpha;
	lyy = L * (1.0 + delta);

	packIntoBoundaries();

	if (programmode == 3 && iterationcountfire == 0)
		cumulativeNeighborchanges = 0;

	if ((programmode == 5) && frprmnconverged && Rneighbor / Rmax < 2.41) {
		if ((iterationcountfire % 10000) == 0) {

			trueneighborNumberOld = trueneighborNumber;
			trueneighborNumber = 0;
			consideredNeighborNumber = 0;

			iloop(N) {
				jloop(i) {
					neighbors[j * N + i] = true;

				}
				neighbors[i * N + i] = false;
			}

			iloop(2*N+3) {
				phelper[i] = p[i];
			}
			energy();
			gradientcalc();

			iloop(N) {
				jloop(i) {
					if (rij[j * N + i] < Rneighbor) {
						neighbors[j * N + i] = true;
						consideredNeighborNumber++;
					} else
						neighbors[j * N + i] = false;
					if (trueneighbors[j * N + i] && !isRattler[i]
							&& !isRattler[j])
						trueneighborNumber++;
				}
			}
		}
	} else {
		trueneighborNumberOld = trueneighborNumber;
		trueneighborNumber = 0;
		consideredNeighborNumber = 0;

		iloop(N) {
			jloop(i) {
				neighbors[j * N + i] = true;

			}
			neighbors[i * N + i] = false;
		}

		iloop(2*N+3) {
			phelper[i] = p[i];
		}
		energy();
		gradientcalc();

		iloop(N) {
			jloop(i) {
				if (rij[j * N + i] < Rneighbor) {
					neighbors[j * N + i] = true;
					consideredNeighborNumber++;
				} else
					neighbors[j * N + i] = false;
				if (trueneighbors[j * N + i] && !isRattler[i] && !isRattler[j])
					trueneighborNumber++;
			}
		}
	}

	iloop(2*N+3) {
		phelper[i] = p[i];
	}
	energy();
	gradientcalc();

	if (programmode == 1 && programmode == 2) {
		Rneighbor = 3.5 * Rmax
				+ 0.3 * L / pow(2.0, iterationcountSimStep * 2.0);
	}
	if (programmode == 5) {
		if (!frprmnconverged) {
			Rneighbor = 3.5 * Rmax
					+ 0.3 * L / pow(2.0, iterationcountSimStep * 0.1);
			RneighborFrprmnLast = Rneighbor;
		} else {
			Rneighbor = 3.5 * Rmax
					+ (RneighborFrprmnLast - 2.4 * Rmax)
							/ pow(2.0, (iterationcountfire) * 0.001);

		}
	}

	if (programmode == 7) {
		Rneighbor = 3.5 * Rmax;
	}

	if (distributioncase == 3 || distributioncase == 4) {
		Rneighbor = 3.5 * Rmax;
		dofOnOff = true;
	}

	if (programmode == 3) {
		Rneighbor = 3.5 * Rmax;
	}

	cumulativeNeighborchanges += fabs(
			trueneighborNumber - trueneighborNumberOld);

	if (programmode == 1 && programmode == 2) {
		frprmn(N, &fret, energy);
		if (frprmnconverged)
			converged = true;
	}
	if (programmode == 5) {
		if (!frprmnconverged) {
			alphaOnOff = false;
			deltaOnOff = false;
			pressOnOff = false;

			frprmn(N, &fret, energy);
		} else {
			dofOnOff = true;

			Pold = Phelper;

			if (dofOnOff) {

				if (alphaOnOffInit)
					alphaOnOff = true;
				if (deltaOnOffInit)
					deltaOnOff = true;
				if (pressOnOffInit)
					pressOnOff = true;
			} else {
				alphaOnOff = false;
				deltaOnOff = false;
				pressOnOff = false;
			} // end else

			fire();

			if (2.0 * fabs(H - HLastFunctionCall)
					< ftolFIRE * (fabs(H) + fabs(HLastFunctionCall) + ZEPS)) {
				if (fabs(sxy) < 1e-15) {

					if (screenOutput)
						cout << "FIRE algorithm converged!" << endl;
					fireconverged = true;
				}
			} else
				endcount = 0;
		} // end else
		if (fireconverged)
			converged = true;
	}

	if (programmode == 3) {
		dofOnOff = false;
		alphaOnOff = false;
		deltaOnOff = false;
		pressOnOff = false;

		shearconverged = false;
		converged = false;
		frprmnconverged = false;
		fireconverged = false;

		fire();

		if (2.0 * fabs(H - HLastFunctionCall)
				< 1e-13 * (fabs(H) + fabs(HLastFunctionCall) + ZEPS)) {

			if (screenOutput)
				cout << "FIRE algorithm converged!" << endl;
			fireconverged = true;
		}

		if (fireconverged) {
			converged = true;
			shearconverged = true;
		}
	}

	calcSysPara();

	if (fireconverged) {

		trueneighborNumberOld = trueneighborNumber;
		trueneighborNumber = 0;
		consideredNeighborNumber = 0;

		iloop(N) {
			jloop(i) {
				neighbors[j * N + i] = true;

			}
			neighbors[i * N + i] = false;
		}

		iloop(2*N+3) {
			phelper[i] = p[i];
		}
		energy();
		gradientcalc();

		iloop(N) {
			jloop(i) {
				if (rij[j * N + i] < Rneighbor) {
					neighbors[j * N + i] = true;
					consideredNeighborNumber++;
				} else
					neighbors[j * N + i] = false;
				if (trueneighbors[j * N + i] && !isRattler[i] && !isRattler[j])
					trueneighborNumber++;
			}
		}
	} // end if

	time(&rawtime1);

	if (programmode != 3) {
		energywrite[iterationcountSimStep] = Uhelper;
		enthalpiewrite[iterationcountSimStep] = Uhelper + P0 * L * L;
		sxywrite[iterationcountSimStep] = sxy;
	}

	if (dofOnOff) {
		int lowerswitch = 10;
		int upperswitch = 50;
		double dampDOF = 0.85;

		if (countAlphaFlip < lowerswitch) {
			dampalpha = 1.01 * (dampalpha);
		}
		if (countAlphaFlip > upperswitch)
			dampalpha = 0.99 * dampalpha;

		if (countDeltaFlip < lowerswitch) {
			dampdelta = 1.01 * (dampdelta);
		}
		if (countDeltaFlip > upperswitch)
			dampdelta = 0.99 * dampdelta;

		if (countPressFlip < lowerswitch) {
			damppress = 1.01 * (damppress);
		}
		if (countPressFlip > upperswitch)
			damppress = 0.99 * damppress;

		if (dampalpha > 1.0)
			dampalpha = 1.0;
		if (dampdelta > 1.0)
			dampdelta = 1.0;
		if (damppress > 1.0)
			damppress = 1.0;

		if (dampalpha < dampDOF)
			dampalpha = dampDOF;
		if (dampdelta < dampDOF)
			dampdelta = dampDOF;
		if (damppress < dampDOF)
			damppress = dampDOF;

		damp = 1.0;

	}

	M[N] = M[N + 1] = N * N;
	M[N + 2] = sqrt(N);

	if (2.0 * fabs(H - HLastFunctionCall)
			< 1e-14 * (fabs(H) + fabs(HLastFunctionCall) + ZEPS)) {
		dtmax = dt * 0.95;
	} else
		dtmax = 0.6;

	string filepath = nameOfWorkingDirectory + "/" + "errorLog.txt";
	if (iterationcountfire > maxIterationCountFire) {
		fireconverged = true;

		ofstream errorLog;
		errorLog.open((char*) filepath.c_str(), ios::app);

		createFileName();
		errorLog << filenameString
				<< " did not converge within maximum amount of iterations."
				<< endl;
	}

	countAlphaFlip = 0;
	countDeltaFlip = 0;
	countPressFlip = 0;

	iterationcountmnbrak = iterationcountbrent = iterationcountfrprmn = 0;
	iterationcountSimStep++;

	dU = Uhelper - UhelperLastFunctionCall;
	dH = H - HLastFunctionCall;

	UhelperLastFunctionCall = Uhelper;
	HLastFunctionCall = H;
	return;
}

// end simulation step

////////////////////////////////////////////////////////////////////////
// menu
void menu() {

	iterationcountSimStep = 0;
	iterationcountfrprmn = 0;
	iterationcountfire = 0;

	char answer;

	converged = false;
	frprmnconverged = false;
	fireconverged = false;
	alphaOnOff = deltaOnOff = false;

	doSimpleShear = false;

	string fileToOpen = "";

	/*
	 cout << endl;
	 cout << "Please choose a program mode." << endl;
	 cout << "1 = conjugate gradient reordering, FIXED unit cell."<< endl;
	 cout << "2 = conjugate gradient reordering, DEFORMABLE unit cell."<< endl;
	 cout << "3 = apply simple shear"<< endl;
	 cout << "4 = FIRE algorithm" << endl;

	 cout << "5 = create numerous (deformed) relaxed packings." << endl;
	 cout << "6 = create random packing." << endl;
	 cout << "7 = DISPLAY system state." << endl;
	 cout << "8 = WRITE system state to file." << endl;

	 cout << "9 = EXIT program"<< endl;
	 */

	if (screenOutput)
		cout << "Please specify WORKING DIRECTORY." << endl;
	cin >> nameOfWorkingDirectory;
	checkFolderName(nameOfWorkingDirectory);
	extractNandP(nameOfWorkingDirectory);

	if (screenOutput) {
		cout << endl;
		cout << "Please choose a program mode." << endl;
		cout << "1 = create DEFAULT jammed packing(s)." << endl;
		cout << "2 = create NON-default jammed packing(s)." << endl;
		cout << "3 = apply simple shear (rearrangements)." << endl;
		cout << "4 = apply simple shear (fixed step size)." << endl;
		cout << "5 = apply compression (rearrangements)." << endl;
		cout << "6 = apply compression (fixed step size)." << endl;
		cout << "7 = DISPLAY system state." << endl << endl;

		cout << "9 = EXIT program" << endl;
	}

	cin >> menumode;
	switch (menumode) {

	case 1:
		programmode = 5;
		ftol = ftolFrprmnBeforeFIRE;
		alphaOnOff = false;
		deltaOnOff = false;
		if (screenOutput)
			cout << "Please insert number of packings to be created!" << endl;
		cin >> numPackingsToProcess;
		if (screenOutput)
			cout << "Please insert starting number of packing names!" << endl;
		cin >> firstPackingNumber;
		currentPackingNumber = firstPackingNumber;
		break;

	case 2:
		programmode = 5;
		ftol = ftolFrprmnBeforeFIRE;
		if (screenOutput)
			cout << "Is simple shear (alpha) a degree of freedom? Y/N :";
		cin >> answer;
		if (answer == 'Y' || answer == 'y')
			alphaOnOffInit = true;
		else
			alphaOnOffInit = false;
		if (screenOutput)
			cout << "Is pure shear (delta) a degree of freedom? Y/N :";
		cin >> answer;
		if (answer == 'Y' || answer == 'y')
			deltaOnOffInit = true;
		else
			deltaOnOffInit = false;
		if (screenOutput)
			cout << "Equilibrate at target pressure? Y/N :";
		cin >> answer;
		if (answer == 'Y' || answer == 'y')
			pressOnOffInit = true;
		else {
			pressOnOffInit = false;
			if (screenOutput)
				cout << "Please insert desired fill fraction:";
			cin >> phiinit;
		}

		alphaOnOff = false;
		deltaOnOff = false;
		if (screenOutput)
			cout << "Please insert number of packings to be created!" << endl;
		cin >> numPackingsToProcess;
		currentPackingNumber = 0;
		break;

	case 3:
		programmode = 5;
		doSimpleShear = true;
		fixedStepSize = false;
		ftol = ftolFrprmnBeforeFIRE;
		alphaOnOff = false;
		deltaOnOff = false;
		if (screenOutput)
			cout << "Please insert number of packings to be created!" << endl;
		cin >> numPackingsToProcess;
		if (screenOutput)
			cout << "Please insert starting number of packing names!" << endl;
		cin >> firstPackingNumber;
		if (screenOutput)
			cout << "Please insert goal number of rearrangements!" << endl;
		cin >> goalNumberOfContactChanges;

		currentPackingNumber = firstPackingNumber;
		break;

	case 4:
		programmode = 5;
		doSimpleShear = true;
		fixedStepSize = true;
		ftol = ftolFrprmnBeforeFIRE;
		alphaOnOff = false;
		deltaOnOff = false;
		if (screenOutput)
			cout << "Please insert number of packings to be created!" << endl;
		cin >> numPackingsToProcess;
		if (screenOutput)
			cout << "Please insert starting number of packing names!" << endl;
		cin >> firstPackingNumber;
		if (screenOutput)
			cout << "Please insert strain goal!" << endl;
		cin >> goalStrain;
		if (screenOutput)
			cout << "Please insert number of equidistant strain steps!" << endl;
		cin >> fixedStepNumber;
		currentPackingNumber = firstPackingNumber;
		break;

	case 5:
		programmode = 5;
		doCompression = true;
		fixedStepSize = false;
		ftol = ftolFrprmnBeforeFIRE;
		alphaOnOff = false;
		deltaOnOff = false;
		if (screenOutput)
			cout << "Please insert number of packings to be created!" << endl;
		cin >> numPackingsToProcess;
		if (screenOutput)
			cout << "Please insert starting number of packing names!" << endl;
		cin >> firstPackingNumber;
		if (screenOutput)
			cout << "Please insert goal number of rearrangements!" << endl;
		if (screenOutput)
			cout << "(neg. values for compression, pos. decompression)" << endl;
		cin >> goalNumberOfContactChanges;
		currentPackingNumber = firstPackingNumber;
		break;

	case 6:
		programmode = 5;
		doCompression = true;
		fixedStepSize = true;
		ftol = ftolFrprmnBeforeFIRE;
		alphaOnOff = false;
		deltaOnOff = false;
		if (screenOutput)
			cout << "Please insert number of packings to be created!" << endl;
		cin >> numPackingsToProcess;
		if (screenOutput)
			cout << "Please insert starting number of packing names!" << endl;
		cin >> firstPackingNumber;
		if (screenOutput)
			cout << "Please insert strain goal!" << endl;
		if (screenOutput)
			cout << "(neg. values for compression, pos. decompression)" << endl;
		cin >> goalStrain;
		if (screenOutput)
			cout << "Please insert number of equidistant strain steps!" << endl;
		cin >> fixedStepNumber;
		currentPackingNumber = firstPackingNumber;
		break;

	case 7:
		onlydisplay = true;
		break;

	case 9:
		endprogram = true;
		//glutDestroyWindow(window);
		return;
		break;

	}
	if (screenOutput) {
		cout << "Please choose a system state." << endl;
		cout << "1 = use CURRENT particle distribution." << endl;
		cout << "2 = use RANDOM particle distribution." << endl;
		cout << "3 = read particle distribution from file." << endl;
		cout
				<< "4 = read particle distribution from file (NEW target pressure)."
				<< endl;
		cout << "5 = open file by name." << endl;

		cout << "9 = EXIT program" << endl;
	}
	cin >> distributioncase;

	switch (distributioncase) {
	case 1:
		if (screenOutput)
			cout << "Please insert target pressure P0: ";
		cin >> P0;
		break;
	case 2:

		initializeSimulation();
		resethelpervars();
		// iloop(N) cout << "[" << p[i] << ", " << p[N+i] << "]" << endl;
		break;
	case 3:
		readPositionFile();
		break;

	case 4:
		if (screenOutput)
			cout << "Please specify target pressure: ";
		cin >> P0init;
		readPositionFile();
		break;

	case 5:
		if (screenOutput)
			cout << "Filename: ";
		cin >> fileToOpen;
		filenameString = fileToOpen;
		if (screenOutput)
			cout << filenameString << endl;
		readPositionFile();
		simulationstep();
		break;

	case 9:
		endprogram = true; //glutDestroyWindow(window);
		break;

	}

	if (programmode == 3) {
		energy();
		energyBeforeDeformation = Uhelper;
		sxxBeforeDeformation = sxx;
		sxyBeforeDeformation = sxy;
		syxBeforeDeformation = syx;
		syyBeforeDeformation = syy;
	}

	return;
}
// end menu

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// fire
void fire() {
	// Velocity Verlet algorithm
	//v += 0.5*(aold+a)*dt
	//p += v*dt+0.5*aold*dt*dt;

	// if velocities are irrelevant and we only want the equilibrium position...
	//v += a*dt
	//p += v*dt

	long double iterPosPower = 0;
	int countAlpha = 0;
	int countDelta = 0;
	int countPress = 0;
	int itercount = 0;
	int endcount = 0;
	long double frac;
	long double fracAlpha;
	long double fracDelta;

	beta = betastart;

	//cout << iterationcountfire << endl;

	while (itercount < 1000) {

		totaliterationcount++;

		Uold = Uhelper;
		Hold = Uhelper + P0 * Lhelper * Lhelper;

		sxyOld = sxy;

		energy(); // calculate overlaps of new configuration
		gradientcalc(); // calulate the gradient for the new configuration

		H = Uhelper + P0 * Lhelper * Lhelper;
		//if (iterationcountfire == 0) cout << H << " fire begin   " <<H- HLastFunctionCall  << endl;

		if (programmode == 3)
			calcSysPara();

		power = 0;
		gg = 0;
		vv = 0;

		iloop(N) { // iloop(2*N+3) to include 3 DoF !
			power -= (xihelper[i] * v[i] + xihelper[N + i] * v[N + i]);
			gg += (xihelper[i] * xihelper[i])
					+ (xihelper[N + i] * xihelper[N + i]);
			vv += (v[i] * v[i]) + (v[N + i] * v[N + i]);
		} // end iloop

		gg = sqrt(gg) + 1e-30;
		vv = sqrt(vv);

		if (power > 0 && iterPosPower > Nmin) {
			dt = fmin(dt * finc, dtmax);
			//dt = fmax(dt, dtmin);
			beta *= fbeta;
		}
		if (power < 0) {
			//dt = fmin(dt*fdec, dtmax);
			if (fabs(sxy) < 1e-16)
				dt = fmax(dt * 0.2 * fdec, dtmin);
			else
				dt = fmax(dt * fdec, dtmin);
			iloop(2*N) {
				v[i] = 0.0;
			}
			beta = betastart;
			iterPosPower = 0;
			//damp = 0.1;
		}

		double radiusFraction = 1e-5;
		//cout << xihelper[2*N+2]<<endl;
		if (alphaOnOff && (fabs(Phelper - P0) / P0 < 5e-1)) {
			v[2 * N] = v[2 * N] * dampalpha - xihelper[2 * N] / M[N] * dt; // damping prevents too large changes
			if (v[2 * N] * Lhelper * dt > radiusFraction)
				v[2 * N] = radiusFraction / Lhelper / dt;
			else if (v[2 * N] * Lhelper * dt < -radiusFraction)
				v[2 * N] = -radiusFraction / Lhelper / dt;
			phelper[2 * N] = p[2 * N] = p[2 * N] + v[2 * N] * dt;
			fracAlpha = p[2 * N] - alphahelper;
		} else {
			fracAlpha = 0.0;
			v[2 * N] = 0.0;
		}

		if (deltaOnOff && (fabs(Phelper - P0) / P0 < 5e-1)) {
			v[2 * N + 1] = v[2 * N + 1] * dampdelta
					- xihelper[2 * N + 1] / M[N + 1] * dt;
			if (v[2 * N + 1] * Lhelper * dt > radiusFraction)
				v[2 * N + 1] = radiusFraction / Lhelper / dt;
			else if (v[2 * N + 1] * Lhelper * dt < -radiusFraction)
				v[2 * N + 1] = -radiusFraction / Lhelper / dt;
			phelper[2 * N + 1] = p[2 * N + 1] = p[2 * N + 1]
					+ v[2 * N + 1] * dt;
			fracDelta = (1 + p[2 * N + 1]) / (1 + deltahelper);
		} else {
			fracDelta = 1.0;
			v[2 * N + 1] = 0.0;
		}

		if (pressOnOff) {
			v[2 * N + 2] = v[2 * N + 2] * damppress
					- xihelper[2 * N + 2] / M[N + 2] * dt;
			if (v[2 * N + 2] * dt > radiusFraction)
				v[2 * N + 2] = radiusFraction / dt;
			else if (v[2 * N + 2] * dt < -radiusFraction)
				v[2 * N + 2] = -radiusFraction / dt;
			phelper[2 * N + 2] = p[2 * N + 2] = p[2 * N + 2]
					+ v[2 * N + 2] * dt;
			frac = p[2 * N + 2] / Lhelper;
		} else {
			frac = 1.0;
			v[2 * N + 2] = 0.0;
		}

		//cout << v[2*N+2] << ", " << xihelper[2*N+2] << ", " << damppress << ", " << M[N+2] << ", " << dt << endl;
		//if(stop != 'x') cin >> stop;

		iloop(N) {

			v[i] = v[i] * damp - xihelper[i] / M[i] * dt; // integrate accelerations to find velocities
			v[N + i] = v[N + i] * damp - xihelper[N + i] / M[i] * dt; // damping added

			v[i] = (1.0 - beta) * v[i] - beta * (xihelper[i] / gg * vv); // FIRE step
			v[N + i] = (1.0 - beta) * v[N + i]
					- beta * (xihelper[N + i] / gg * vv); // FIRE step

			phelper[i] = phelper[i] + phelper[N + i] * fracAlpha;
			phelper[i] = phelper[i] / fracDelta;
			phelper[N + i] = phelper[N + i] * fracDelta;

			phelper[i] = p[i] = (phelper[i] + v[i] * dt) * frac; // integrate velocities to find positions
			phelper[N + i] = p[N + i] = (phelper[N + i] + v[N + i] * dt) * frac;

		}
		long double vimaxloc;
		long double indexmaxloc;
		//		iloop(2*N+1){
		//			if(fabs(v[i+1])>fabs(v[i])) {
		//			vimaxloc=v[i+1];
		//		indexmaxloc=i+1;
		//				}
		//}
		//cout << "max gradient component " << indexmaxloc ;
		//cout << "  , value : " << vimaxloc <<endl;

		////////////

		if (dofOnOff) {

			if (xihelper[2 * N] * v[2 * N] > 0) {
				v[2 * N] = 0.0;
				countAlphaFlip++;
			} else {
				countAlpha++;
			}

			if (xihelper[2 * N + 1] * v[2 * N + 1] > 0) {
				v[2 * N + 1] = 0.0;
				countDeltaFlip++;
			} else {
				countDelta++;
			}

			if (xihelper[2 * N + 2] * v[2 * N + 2] > 0) {
				//v[2*N+2] = 0.0;
				countPressFlip++;
			} else {
				countPress++;
			}
		}
		////////////

		//enthalpieDiffStepOld = enthalpieDiffStep;
		enthalpieDiffStep = (H - Hold);
		energyDiffStep = (Uhelper - Uold);

		//if(energyDiffStep > 0.0) dt *= 0.05;

		iterPosPower++;
		iterationcountfire++;
		itercount++;
		//cout <<"loop end : H = "<< H << "  dH = "<< H-Hold << endl;

	} // end while

	//cout << "In FIRE: H-Hold = " << H - Hold << endl;

} // end FIRE

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// mnbrak
void mnbrak(long double *ax, long double *bx, long double *cx, long double *fa,
		long double *fb, long double *fc, long double (*func)(long double)) {

	long double ulim, u, r, q, fu, dum;
	int numberOfIterations = 0;

	*fa = (*func)(*ax);
	*fb = (*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum, *ax, *bx, dum)
		SHFT(dum, *fb, *fa, dum)
	} // end if
	*cx = (*bx) + gold * (*bx - *ax);
	*fc = (*func)(*cx);

	while (*fb > *fc) {
		numberOfIterations++;
		iterationcountmnbrak++;

		r = (*bx - *ax) * (*fb - *fc);
		q = (*bx - *cx) * (*fb - *fa);
		u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) / (2.0 * (q - r));
		// possible div by 0 !!!!!

		ulim = (*bx) + glimit * (*cx - *bx); // maximum step
		// "Where is u?"
		if ((*bx - u) * (u - *cx) > 0.0) { // IF_1
			fu = (*func)(u);
			if (fu < *fc) { // IF_2
				*ax = (*bx);
				*bx = u;
				*fa = (*fb);
				*fb = fu;
				//                      if(screenOutput) cout << "Number of iterations in mnbrak():  " << numberOfIterations << endl;
				return;
			} // end IF_2
			else if (fu > *fb) {
				*cx = u;
				*fc = fu;
				//                     if(screenOutput) cout << "Number of iterations in mnbrak(): " << numberOfIterations << endl;
				return;
			} // end ELSE IF
			u = (*cx) + gold * (*cx - *bx);
			fu = (*func)(u);
		} // end IF_1
		else if ((*cx - u) * (u - ulim) > 0.0) { // ELSE IF_1
			fu = (*func)(u);
			if (fu < *fc) {
				SHFT(*bx, *cx, u, *cx+gold*(*cx-*bx))
				SHFT(*fb, *fc, fu, (*func)(u))
			} // end IF
		} // end ELSE IF_1
		else if ((u - ulim) * (ulim - *cx) > 0.0) { // ELSE IF_2
			u = ulim;
			fu = (*func)(u);
		} // end ELSE IF_2
		else {
			u = (*cx) + gold * (*cx - *bx);
			fu = (*func)(u);
		} // end ELSE
		SHFT(*ax, *bx, *cx, u)
		SHFT(*fa, *fb, *fc, fu)
	} // end while
	  //      if(screenOutput) cout << "Number of iterations in mnbrak(): " << numberOfIterations << endl;

}
// end mnbrak

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// brent
long double brent(long double ax, long double bx, long double cx,
		long double (*f)(long double), long double tol, long double *xmin) {
	int iter;
	long double a, b, d, etemp, fu, fv, fw, fx, s, q, r, tol1, tol2, u, v, w, x,
			xm;
	long double e = 0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = (*f)(x);
	for (iter = 1; iter <= ITMAXBRENT; iter++) {
		iterationcountbrent++;
		//cout<< iterationcountbrent <<endl;
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			s = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0)
				s = -s;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(s) >= fabs(0.5 * q * etemp) || s <= q * (a - x)
					|| s >= q * (b - x))
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			else {
				d = s / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm - x);

			}
		} else {
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = (*f)(u);
		if (fu <= fx) {
			if (u >= x)
				a = x;
			else
				b = x;
			SHFT(v, w, x, u)
			SHFT(fv, fw, fx, fu)

		} else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;

			} else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}

	}
	//	if(screenOutput) cout << "Too many iterations in brent" << endl;
	*xmin = x;
	return fx;
} // end brent

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
long double SIGN(long double a, long double b) {
	if (b < 0.0)
		return -a;
	else
		return a;
} // end SIGN

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// linmin
void linmin(int n, long double *fret, long double (*func)()) {
	int j;
	long double xx, xmin, fx, fb, fa, bx, ax;

	//ncom = n;
	// pcom, xicom are global row-matrices
	nrfunc = func;

	//pcom = p;
	//xicom = xi;

	ax = 0.0;
	xx = AMIN;
	mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
	*fret = brent(ax, xx, bx, f1dim, TOL, &xmin);

	// cout << "xmin = " << xmin << "             (in linmin())" << endl;
	iloop(2*N+3) {
		p[i] += xi[i] * xmin;
	} // move by gradient*xmin to the 1D-minimum
	CG_step = xmin;
	return;
} // end linmin

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// f1dim
long double f1dim(long double x) {
	long double Uloc;

	iloop(2*N+3) {
		phelper[i] = p[i] + x * xi[i];
	} // move by gradient*x ...
	energy();
	Uloc = Uhelper;
	if (pressOnOff)
		Uloc = Uloc + P0 * Lhelper * Lhelper;

	return Uloc; // ... and return the energy of the new configuration
} // end f1dim

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// frprmn
void frprmn(int n, long double *fret, long double (*func)()) {
	// most input variable are global, dfunc is done in gradienU
	int j, its;
	int endcount = 0;
	long double gam, fp, dgg;

	iloop(2*N+3) {
		phelper[i] = p[i];
	}
	(*func)(); // func = energy does not need inputvector
	fp = Uhelper;

	if (pressOnOff)
		fp = fp + P0 * Lhelper * Lhelper;
	energy();
	gradientcalc();

	iloop(2*N+3) {
		xi[i] = h[i] = g[i] = -xihelper[i]; // the 'force' is in the negative gradient direction
	}

	for (its = 1; its <= ITMAX; its++) {
		packIntoBoundaries();
		iterationcountfrprmn++;
		iterationcountfrprmnCUMULATIVE++;
		totaliterationcount++;

		iterfrprmn = its;

		linmin(n, fret, func);
		//cout << "H = " << fp << endl;
		energyDiffStepOld = energyDiffStep;
		energyDiffStep = (*fret - fp);

		/*
		 if(energyDiffStep < -1e-4){
		 alphaOnOff = false;
		 deltaOnOff = false;
		 pressOnOff = false;
		 }
		 else if(energyDiffStep < -1e-8){
		 alphaOnOff = true;
		 deltaOnOff = true;
		 pressOnOff = false;
		 }
		 else{
		 alphaOnOff = true;
		 deltaOnOff = true;
		 pressOnOff = true;
		 }
		 */

		if (2.0 * fabs(*fret - fp) <= ftol * (fabs(*fret) + fabs(fp) + ZEPS)) {
			endcount++;

			if (endcount > 3) {
				if (screenOutput) {
					cout << "Conjugate gradient converged!" << endl;
					//				   cout << "EnergyDiffStep/U= "<< energyDiffStep/(*fret) << ", last: " << energyDiffStepOld/(*fret) << endl;
				}

				frprmnconverged = true;
				return;
			}
		} else
			endcount = 0;
		fp = *fret;
		energy();
		gradientcalc();
		iloop(2*N+3) {
			xi[i] = xihelper[i];
		}
		dgg = gg = 0.0;

		iloop(2*N+3) {
			gg += g[i] * g[i];
			dgg += xi[i] * xi[i]; // Fletcher-Reeves
			//dgg += (xi[i] + g[i])*xi[i]; // Polak-Ribiere
		}

		if (gg == 0.0) {
			if (screenOutput)
				cout << "Gradient is ZERO!" << endl;
			frprmnconverged = true;
			return;
		}
		gam = dgg / gg;
		iloop(2*N+3) {
			g[i] = -xi[i];
			xi[i] = h[i] = g[i] + gam * h[i];
		}
		H = Uhelper + P0 * Lhelper * Lhelper;

	}

	//	if(screenOutput) cout << "Too many iterations in frprmn!" << endl;
	return;
} // end frprmn

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// particle distance
void particledistance(int i, int j) {

	xij[j * N + i] = phelper[j] - phelper[i];
	yij[j * N + i] = phelper[N + j] - phelper[N + i];

	ny[j * N + i] = -floor((yij[j * N + i] + lyyhelper * 0.5) / lyyhelper);

	xij[j * N + i] = phelper[j] - phelper[i] + ny[j * N + i] * lyxhelper;
	yij[j * N + i] = phelper[N + j] - phelper[N + i]
			+ ny[j * N + i] * lyyhelper;
	/*
	 xij[j*N+i] = phelper[j]-phelper[i]+nx[j*N+i]*lxxhelper+ny[j*N+i]*lyxhelper;
	 yij[j*N+i] = phelper[N+j]-phelper[N+i]+nx[j*N+i]*lxyhelper+ny[j*N+i]*lyyhelper;
	 */
	nx[j * N + i] = -floor((xij[j * N + i] + lxxhelper * 0.5) / lxxhelper);

	xij[j * N + i] = phelper[j] - phelper[i] + nx[j * N + i] * lxxhelper
			+ ny[j * N + i] * lyxhelper;
	yij[j * N + i] = phelper[N + j] - phelper[N + i] + nx[j * N + i] * lxyhelper
			+ ny[j * N + i] * lyyhelper;

	rij[j * N + i] = sqrt(
			xij[j * N + i] * xij[j * N + i] + yij[j * N + i] * yij[j * N + i]);

	return;
} // end particledistance

/////////////////////////////////////////////////////////////////////////////////////////
// resetvars
void resethelpervars() {
	// reset variables

	alphahelper = phelper[2 * N];
	deltahelper = phelper[2 * N + 1];
	Lhelper = phelper[2 * N + 2];
	alpha = p[2 * N];
	delta = p[2 * N + 1];
	L = p[2 * N + 2];
	lxxhelper = Lhelper / (1.0 + deltahelper);
	lxyhelper = Lhelper * 0.0;
	lyxhelper = Lhelper * alphahelper;
	lyyhelper = Lhelper * (1.0 + deltahelper);

	Uhelper = 0.0;

	return;
} // end resetvars

///////////////////////////////////////////////////////////////////////////////////////////
// energy
long double energy() {

	resethelpervars();

	// calculate distances and overlaps
	iloop(N) {
		jloop(i) // only calculate j<i, i.e., upper right corner of matrices
		{

			if (neighbors[j * N + i]) { // avoid selfinteraction

				particledistance(i, j);

				// symmetrize/ antisymmetrize matrices
				/*
				 xij[i*N+j] = -xij[j*N+i];
				 yij[i*N+j] = -yij[j*N+i];
				 nx[i*N+j] = -nx[j*N+i];
				 ny[i*N+j] = -ny[j*N+i];
				 rij[i*N+j] = rij[j*N+i];
				 */

				dij[j * N + i] = R[i] + R[j] - rij[j * N + i];

				if (dij[j * N + i] > 0.0) { // this saves about 5-10 % in runtime
					Uhelper += 0.5 * k * (dij[j * N + i] * dij[j * N + i]);
					trueneighbors[j * N + i] = true;
				} else {
					trueneighbors[j * N + i] = false;
					dij[j * N + i] = 0.0;
				}

				distanceCalcs += 1.0;

			} //end if
		} // end jloop
		trueneighbors[i * N + i] = false;
	} // end iloop

	//Uhelper += 1e-30;
	return Uhelper;
} // end energy

//////////////////////////////////
// gradientcalc
void gradientcalc() {
	long double frac;
	long double component;
	long double term1;
	// calculate gradient and energy        cout << "xi[0]= " << xi[0] << endl;

	Phelper = 0.0;

	xihelper[2 * N] = 0.0;
	xihelper[2 * N + 1] = 0.0;
	xihelper[2 * N + 2] = 0.0;

	iloop(N) {

		xihelper[i] = 0.0;
		xihelper[N + i] = 0.0;

		jloop(i) {
			if (trueneighbors[j * N + i]) // this saves about 5-10 % in runtime
			{
				term1 = k * dij[j * N + i] / rij[j * N + i];
				// gradient expressions have been derived analytically
				component = term1 * xij[j * N + i]; //gradient component i with respect to particle j
				xihelper[i] += component;
				xihelper[j] -= component; // opposite sign in xji vs. xij

				component = term1 * yij[j * N + i];
				xihelper[N + i] += component;
				xihelper[N + j] -= component;

				if (alphaOnOff)
					xihelper[2 * N] += -2 * term1 * Lhelper * xij[j * N + i]
							* ny[j * N + i];
				if (deltaOnOff)
					xihelper[2 * N + 1] += -2 * term1 * Lhelper
							* (xij[j * N + i] * nx[j * N + i]
									* (-1.0 / ((1.0 + delta) * (1.0 + delta)))
									+ yij[j * N + i] * ny[j * N + i]);
				Phelper += 2 * term1 * rij[j * N + i] * rij[j * N + i];
			} // end if
		} // end jloop
	}

	// end iloop

	// normalize
	Phelper = Phelper / (Lhelper * Lhelper) / 4.0;

	if (fabs(p[2 * N + 1] - deltainit) > 0.2)
		deltaOnOff = false;

	if (pressOnOff) {
		//if(!frprmnconverged){
		//cout << "P0 = " << P0 << "  , and P = "<<Phelper << endl     ;

		//xihelper[2*N+2] += 2*Lhelper*(P0 - Phelper);
		xihelper[2 * N + 2] = 2 * Lhelper * (P0 - Phelper);
		//if(xihelper[2*N+2] > 0.2) xihelper[2*N+2] = 0.2;

		/*
		 frac = (Lhelper-xihelper[2*N+2])/Lhelper;

		 iloop(2*N){ // alpha=p[2*N] and delta=p[2*N+1] do not scale and L=p[2*N+2] is updated by its gradient
		 phelper[i] *= frac;
		 xihelper[i] *= frac;				// scale all the particle positions accordingly
		 }
		 */

	}

	return;
} // end gradientcalc

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// initializeSimulation

void initializeSimulation() {

	long double randLx, randLy;
	long double frac;

	dtmax = 0.1;
	dt = dtmax;

	if (currentPackingNumber == firstPackingNumber) {
		// size all the dynamic 'vector' arrays according to the particle number N
		dij.reserve(N * N);
		rij.reserve(N * N);
		xij.reserve(N * N);
		yij.reserve(N * N);
		nx.reserve(N * N);
		ny.reserve(N * N);
		neighbors.reserve(N * N);
		trueneighbors.reserve(N * N);
		trueneighborsOld.reserve(N * N);
		trueneighborsLast.reserve(N * N);
		trueneighborChanges.reserve(N * N);

		p.reserve(2 * N + 3);
		phelper.reserve(2 * N + 3);
		pcom.reserve(2 * N + 3);
		pLast.reserve(2 * N + 3);
		xi.reserve(2 * N + 3);
		xihelper.reserve(2 * N + 3);
		xicom.reserve(2 * N + 3);
		g.reserve(2 * N + 3);
		h.reserve(2 * N + 3);
		v.reserve(2 * N + 3);
		R.reserve(N);
		M.reserve(N + 3);

		numberOfDirectNeighbors.reserve(N);
		isRattler.reserve(N);
		wasRattler.reserve(N);
	} // end if

	// set initial unit cell properties
	alpha = alphainit;
	delta = deltainit;

	dampalpha = dampdelta = 0.9;
	damppress = 0.9;

	p[2 * N] = alpha;
	p[2 * N + 1] = delta;
	p[2 * N + 2] = L;

	lxx = L / (1.0 + delta);
	lxy = L * 0.0;
	lyx = L * alpha;
	lyy = L * (1.0 + delta);

	/*
	 //packingNameNumber++;
	 srand(numberOfPackings+1);
	 cout << numberOfPackings << ", " << rand() << endl;
	 cin >> stop;
	 */

	srand(currentPackingNumber + 1);

	// set initial particle positions and properties
	iloop(N) {
		randLx = (rand() * 1.0) / (RAND_MAX * 1.0); // random position along L_x
		randLy = (rand() * 1.0) / (RAND_MAX * 1.0); // random position along L_y
		p[i] = randLx * lxx + randLy * lyx; // set random x-position for particle
		p[N + i] = randLy * lxy + randLy * lyy; // set random y-position for particle
		if (i < 0.5 * N)
			R[i] = 1.0; // set particle radii
		else
			R[i] = 1.4;
		M[i] = 1.0; // set particle 'masses' for FIRE-algorithm
	}
	M[N] = M[N + 1] = M[N + 2] = N;

	// rescale to fit desired fill fraction
	iloop(2*N+2)
		phelper[i] = p[i];

	energy();
	calcSysPara();

	// calculate the current fill fraction phi
	frac = sqrt(phi / (phiinit)); // sqrt-ratio of desired fill fraction and actual
	// frac = sqrt(P0/P);
	L *= frac; // scale the boxlength accordingly
	iloop(2*N) {
		p[i] *= frac; // scale all the particle positions accordingly
		v[i] = 0.0;
	}
	p[2 * N + 2] = L;

	v[2 * N] = v[2 * N + 1] = v[2 * N + 2] = 0;

	energy();
	calcSysPara();

	Rmax = 0.0;
	iloop(N) {
		if (R[i] > Rmax)
			Rmax = R[i]; // determine the largest particle radius in the packing
		// (for neighbor determination)
	}

	iloop(N) {
		jloop(N) {
			neighbors[j * N + i] = true; // initially all particles are concidered
			if (j == i)
				neighbors[i * N + i] = false; // neighbors, except if j==i
		}
	}

	int ITMAXBRENT = max(100, 2 * N); // maximum of iterations in brent
	int ITMAX = 5; // maximum of iterations in frprmn

	dtmax = dtmaxinit;
	dofOnOff = false;

	return;
} // end initializeSimulation()

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// calcSysPara: calculate fillfracton & number of neighbors & pressure
void calcSysPara() {
	long double term1;

	//vector <int> numberOfDirectNeighbors;
	//vector <bool> isRattler;
	int numberOfRattlerChanges = 1;

	long double Z2 = 0;

	// reset sums
	phi = 0.0;
	Z = 0.0;
	P = 0.0;
	sxx = sxy = syx = syy = 0.0;

	// sum all contributions
	iloop(N) {
		phi += PI * R[i] * R[i]; // sum of particle surfaces
		jloop(i) {
			if (trueneighbors[j * N + i]) {

				Z++;

				term1 = k * dij[j * N + i] / rij[j * N + i];
				// if overlap > 0 --> one more contact
				P += term1 * rij[j * N + i] * rij[j * N + i]; // pressure
				sxx += -term1 * xij[j * N + i] * xij[j * N + i]; // xx-stress
				syy += -term1 * yij[j * N + i] * yij[j * N + i]; // yy-stress
				sxy += -term1 * xij[j * N + i] * yij[j * N + i]; // xy-/ yx-stress

			}
		}
	}

	if (fireconverged) {
		// Ncorrected
		iloop(N) {
			isRattler[i] = false;
			numberOfDirectNeighbors[i] = 0;
			jloop(i)
				trueneighbors[i * N + j] = trueneighbors[j * N + i];
		}

		while (numberOfRattlerChanges != 0) {
			numberOfRattlerChanges = 0;
			iloop(N) {
				numberOfDirectNeighbors[i] = 0;
				jloop(N) {

					if (trueneighbors[j * N + i]) {
						if (!isRattler[j]) {
							numberOfDirectNeighbors[i] += 1; //if(!isRattler[j] && !isRattler[i])
							//numberOfDirectNeighbors[j] += 1;
							Z2 += 1;
						}

					}
					//cout << trueneighbors[j*N+i];

				}
				//cout << endl;
				if (numberOfDirectNeighbors[i] < 3 && !isRattler[i]) {
					isRattler[i] = true;
					numberOfRattlerChanges++;
					//cout << "Rattler found!" << endl;
				}

			}
			//cout << endl << endl;
		} // end while
		Ncorrected = 0;
		iloop(N) {
			if (!isRattler[i])
				Ncorrected++;
			//cout << numberOfDirectNeighbors[i] << endl;
		}
		//    if(screenOutput) cout << "N = " << N << ", Ncorrected = " << Ncorrected << ", Z = " << 2*Z/(N*1.0) << ", Zcorrected = " << 2*Z/(Ncorrected*1.0 + 1e-16) << endl;
		//    if(screenOutput) cout << "Z2 = " << Z2 << endl;

		/*
		 iloop(N){
		 //cout << "Point 0" << endl;
		 if(isRattler[i]) resolveRattler(i);	// locally equilibrate all rattlers and get rid of contacs
		 }
		 */

	} // end if(fireconverged)
	else
		Ncorrected = N;

	// normalize
	phi = phi / (L * L);
	//Z = 2*Z/(N*1.0);
	Z = 2 * Z / (Ncorrected * 1.0 + 1e-16);
	P = P / (L * L) / 2.0;
	sxx = sxx / (L * L);
	syy = syy / (L * L);
	sxy = sxy / (L * L);
	syx = syx;

	//cout << "Point 2" << endl;

	return;
} // end calcSysPara

void resolveRattler(int rattlerElement) {

	vector<bool> bufferneighbors;
	bufferneighbors.reserve(N * N);

	iloop(N) {
		jloop(i) {
			bufferneighbors[j * N + i] = neighbors[j * N + i];

			if ((i == rattlerElement || j == rattlerElement)
					&& neighbors[j * N + i] == true) {
				neighbors[j * N + i] = true;
			} else
				neighbors[j * N + i] = false;
		}
	}

	//cout << "Point 1" << endl;

	while (numberOfDirectNeighbors[rattlerElement] > 0
			&& numberOfDirectNeighbors[rattlerElement] < 3) {

		//cout << xihelper[rattlerElement] << ", " << numberOfDirectNeighbors[rattlerElement]  << endl;

		phelper[rattlerElement] = p[rattlerElement];
		phelper[rattlerElement + N] = p[rattlerElement + N];

		energy();
		gradientcalc();

		p[rattlerElement] -= xihelper[rattlerElement];
		p[rattlerElement + N] -= xihelper[rattlerElement + N];

	}

	//cout << "Point 4" << endl;

	iloop(N) {
		jloop(i) {
			neighbors[j * N + i] = bufferneighbors[j * N + i];
		}
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// readPositionFile
void readPositionFile() {

	long double frac;
	char c = 'x', clast = 'x', clastlast = 'x';
	//long double pos[2*N+2], rad[N];
	int i = 0, j = 0;
	bool initializeNow = false;

	int decimal = 0;
	long double helperchar = 0.0;
	bool isnegative = false;
	int numofrightbrackets = 0;
	bool posRead = false;
	bool radRead = false;
	bool Nread = false;
	bool L1read = false;
	bool L2read = false;
	bool P0read = false;
	string filepath = nameOfWorkingDirectory + "/";

	long double L1x = 0.0, L1y = 0.0, L2x = 0.0, L2y = 0.0;

	ifstream infile;
	//infile.open("jamstate.txt");

	if (distributioncase != 5) {
		if (screenOutput)
			cout << "Creating file name." << endl;
		createFileName();
		filepath.append(filenameString);
		if (screenOutput)
			cout << "Filename: " << filepath << endl;
		infile.open((char*) filepath.c_str());
	} else
		infile.open((char*) filepath.c_str());

	if (!infile.is_open()) {
		currentPackingNumber = 0;
		if (screenOutput) {
			cout << "Input file did NOT OPEN!" << endl;
			cout << "Press any key + ENTER to proceed." << endl;
			cin >> stop;
		}
		return;
	} else if (screenOutput)
		cout << "Input file is opened succesfully!" << endl;

	while (infile.good()) {
		clastlast = clast;
		clast = c;
		c = infile.get();

		if (clastlast == 'N' && c == '=') {
			Nread = true;
			decimal = 0;
			helperchar = 0.0;

		}
		if (Nread) {
			if (c > 47 && c < 58) {
				helperchar = helperchar * 10 + (c - 48);
				decimal++;
			}
			if (c == ',') {
				N = helperchar;
				helperchar = 0.0;
				Nread = false;
				initializeNow = true;
				//cout << "N = " << N << endl;
			}

		}

		if (initializeNow) {
			initializeArrays();
			initializeNow = false;
		}

		if (clastlast == 'L' && clast == '1' && c == '=') {
			L1read = true;
			i = 0;
			decimal = 0;
			helperchar = 0.0;

		}
		if (L1read) {
			if (c == '-')
				isnegative = true;
			if (c > 47 && c < 58) {
				helperchar = helperchar * 10.0 + (c - 48) * 1.0;
				decimal++;
			}
			if (c == '.')
				decimal = 0;
			if (c == ',' || c == '}') {

				for (int j = 0; j < decimal; j++) {
					helperchar *= 0.1;
				}
				if (!isnegative && i == 0)
					L1x = helperchar;
				if (isnegative && i == 0)
					L1x = -helperchar;
				if (!isnegative && i == 1)
					L1y = helperchar;
				if (isnegative && i == 1)
					L1y = -helperchar;

				helperchar = 0.0;
				isnegative = false;
				i++;

			}
			if (i == 2) {
				L1read = false;
				//cout << "L1= (" << L1x << ", " << L1y << ") " << endl;
				//cin >> stop;
			}
		}

		if (clastlast == 'L' && clast == '2' && c == '=') {
			L2read = true;
			i = 0;
			decimal = 0;
			helperchar = 0.0;

		}
		if (L2read) {
			if (c == '-')
				isnegative = true;
			if (c > 47 && c < 58) {
				helperchar = helperchar * 10.0 + (c - 48) * 1.0;
				decimal++;
			}
			if (c == '.')
				decimal = 0;
			if (c == ',' || c == '}') {

				for (int j = 0; j < decimal; j++) {
					helperchar *= 0.1;
				}
				if (!isnegative && i == 0)
					L2x = helperchar;
				if (isnegative && i == 0)
					L2x = -helperchar;
				if (!isnegative && i == 1)
					L2y = helperchar;
				if (isnegative && i == 1)
					L2y = -helperchar;

				helperchar = 0.0;
				isnegative = false;
				i++;

				//cout << "L= " << L << endl;
			}
			if (i == 2) {
				L2read = false;
				//cout << "L2= (" << L2x << ", " << L2y << ") " << endl;
				//cin >> stop;
			}
		}

		if (clastlast == 'P' && clast == '0' && c == '=') {
			P0read = true;
			decimal = 0;
			helperchar = 0.0;

		}
		if (P0read) {
			if (c > 47 && c < 58) {
				helperchar = helperchar * 10.0 + (c - 48) * 1.0;
				decimal++;
			}
			if (c == '.')
				decimal = 0;
			if (c == ',') {

				for (int j = 0; j < decimal; j++) {
					helperchar *= 0.1;
				}
				P0 = helperchar;
				helperchar = 0.0;
				P0read = false;
				//cout << "P0= " << P0 << endl;

			}
		}

		if (!L1read && !L2read && c == '{') {
			posRead = true;
			decimal = 0;
			helperchar = 0.0;

		}
		if (posRead) {
			if (c == '-')
				isnegative = true;
			if (clastlast == '{') {
				i = 0;
				decimal = 0;
			}
			if (c > 47 && c < 58) {
				helperchar = helperchar * 10.0 + (c - 48) * 1.0;
				decimal++;
			}
			if (c == '.')
				decimal = 0;
			if (c == ',' || c == '}') {

				for (int j = 0; j < decimal; j++) {
					helperchar *= 0.1;
				}
				if (i < 3 * N && !isnegative && i % 3 == 0)
					p[i / 3] = helperchar;
				if (i < 3 * N && isnegative && i % 3 == 0)
					p[i / 3] = -helperchar;

				if (i < 3 * N && !isnegative && i % 3 == 1)
					p[N + i / 3] = helperchar;
				if (i < 3 * N && isnegative && i % 3 == 1)
					p[N + i / 3] = -helperchar;

				if (i < 3 * N && !isnegative && i % 3 == 2)
					R[i / 3] = helperchar;
				if (i < 3 * N && isnegative && i % 3 == 2)
					R[i / 3] = -helperchar;

				//cout << p[i/3] << " " << p[N+i/3] << " " << R[i/3] << endl;

				i++;

				isnegative = false;
				//cout << helperchar << endl;
				helperchar = 0.0;
			}

			if (c == '}') {
				posRead = false;
				//cin >> stop;
			}

		}

		/*
		 if(clastlast == 'R' && c == '='){
		 radRead = true;
		 decimal=0;
		 helperchar=0.0;

		 }
		 if(radRead){
		 if(c=='-') isnegative=true;
		 if(clast=='['){
		 i=0;
		 decimal=0;
		 helperchar=0.0;
		 }
		 if(c > 47 && c < 58) {
		 helperchar = helperchar*10.0+(c-48)*1.0;
		 decimal++;
		 }
		 if(c=='.') decimal=0;
		 if(c==',' || c==']'){

		 for(int j=0; j<decimal; j++){
		 helperchar *=0.1;
		 }
		 if(i<N && !isnegative) R[i] = helperchar;
		 if(i<N && isnegative) R[i] = -helperchar;

		 //cout << R[i] << endl;
		 i++;
		 isnegative=false;
		 //cout << helperchar << endl;
		 helperchar=0.0;
		 }
		 if(c==']') {
		 //posRead = false;
		 radRead = false;
		 //numofrightbrackets++;
		 }
		 }
		 */
	} // end while

	infile.close();

	Rmax = 0.0;
	iloop(N) {
		if (R[i] > Rmax)
			Rmax = R[i]; // determine the largest particle radius in the packing
		// (for neighbor determination)
	}

	L = p[2 * N + 2] = sqrt(L1x * L2y);
	alpha = p[2 * N] = L2x / L;
	delta = p[2 * N + 1] = sqrt(L2y / L1x) - 1.0;

	//cout << "L = " << L << ", alpha = " << alpha << ", delta = " << delta << endl;
	//cin >> stop;

	iloop(N) {
		jloop(N) {
			neighbors[j * N + i] = true; // initially all particles are concidered
			if (j == i)
				neighbors[i * N + i] = false; // neighbors, except if j==i
		}
	}

	if (test)
		cout << "test0000: readPositionFile-END: " << " N = " << N << ", L = "
				<< L << ", P0= " << P0 << endl;

	//    if(screenOutput) cout << "x = %, y = %, R = %:" << endl;
	//iloop(N){ cout << p[i] <<" " << p[N+i] <<" "<< R[i] << endl;}
	//    if(screenOutput) cout << "alpha = " << p[2*N] << ", delta = " << p[2*N+1] << ", P0 = " << P0 << endl;
	//    if(screenOutput) cout << " N = " << N << ", L = " << L << " = " << p[2*N+2] << endl;
	//cin >> stop;

	alpha = p[2 * N];
	delta = p[2 * N + 1];

	if (distributioncase == 4)
		P0 = P0init;

	return;
} // end readPositionFile

////////////////////////////////////////////////////////////////////////
// writePositionFile
inline void writePositionFile() {

	//char filename[256];
	time_t rawtime;
	struct tm *timeinfo;
	char timebuffer[80];
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(timebuffer, 80, "%Y-%m-%d_%H-%M-%S", timeinfo);

	string filepath = nameOfWorkingDirectory;

	filepath.append("/");

	createFileName();

	//	cout << "filenameString at step 1 : "<< endl;
	//	cout << filenameString << endl;

	filepath.append(filename);

	ofstream outfile;
	//	cout << filepath << endl;
	//outfile.open( "jamstate.txt", ios::trunc);
	if (programmode == 5)
		outfile.open((char*) filepath.c_str(), ios::trunc);

	//else outfile.open( "jamstate.txt" , ios::trunc);
	if (!outfile) {
		if (screenOutput)
			cout << "Cannot open output file!" << endl;
		outfile.close();
		ofstream outfile;
		createFileName();
		if (programmode == 5)
			outfile.open(filename, ios::trunc);
		//else outfile.open( "jamstate.txt" , ios::trunc);
		if (!outfile) {
			if (screenOutput)
				cout << "Cannot open output file!" << endl;
			return;
		} else if (screenOutput)
			cout << "Output file is opened succesfully (2nd attempt)!" << endl;
	} else if (screenOutput)
		cout << "Output file is opened succesfully!" << endl;
	outfile.setf(ios::fixed, ios::floatfield);
	outfile.precision(16);

	alpha = p[2 * N];
	delta = p[2 * N + 1];
	L = p[2 * N + 2];
	//    if(screenOutput) cout << "L = " << L << endl;
	lxx = L / (1.0 + delta);
	lxy = L * 0.0;
	lyx = L * alpha;
	lyy = L * (1.0 + delta);

	outfile << "N = " << N << " ,L = " << L << " ,L1= { " << lxx << " , " << lxy
			<< " } " << " ,L2= { " << lyx << " , " << lyy << " } " << " ,P = "
			<< Phelper << " ,P0= " << P0 << " ," << endl;

	//outfile << "L = " << Lhelper << ", phi = " << phi << endl;
	/*
	 outfile << "[";
	 for(int i=0; i<2*N+2; i++){
	 outfile << p[i] << ", ";
	 }
	 outfile << p[2*N+2] << "],";
	 outfile << endl << endl;
	 outfile << "R = [";
	 for(int i=0; i<N-1; i++){
	 outfile << R[i] << ", ";
	 }
	 outfile << R[N-1] << "],";
	 */
	outfile << "{" << endl;
	iloop(N) {
		outfile << p[i] << " ,	" << p[N + i] << " ,	" << R[i] << " ,	" << endl;
	}
	outfile << "}" << endl;
	outfile << endl << endl;

	outfile.close();
	outfile.close();

	//string filenameString;
	char buffer;
	string bufferString;
	long double number, teni, teni1;
	int digit;
	char *filenametime;

	/*
	 filenameString = "timefile";

	 number = N;
	 for(int i=5; i>0; i--){
	 teni = pow(10,i*1.0);
	 teni1 = pow(10,i-1.0);
	 digit = fmod(number,teni)/(teni1);
	 buffer = digit+48;
	 bufferString = buffer;
	 filenameString.append(bufferString);
	 }
	 filenameString.append(".txt");
	 filenametime = (char*)filenameString.c_str();
	 */

	/*
	 ofstream timefile;
	 timefile.open(filenametime, ios::app);
	 timefile << "t = " << timediff1 << ", fire = " << totaliterationcount << ", alpha = " << alphahelper << ", delta = " << deltahelper << ", phi = " << phi << ", Z = " << Z << ", dH/H = " << (H-HLastFunctionCall)/H << ", s_xy = " << sxy << ", P = " << P << endl;
	 timefile.close();
	 */

	string logFileNamePackings;
	logFileNamePackings = filenameString;
	//	size_t sz=logFileNamePackings.size();
	//	logFileNamePackings.erase(sz-9,5);

	logFileNamePackings.insert(0, "log");

	logFileNamePackings.insert(0, nameOfWorkingDirectory + "/");

	ofstream logfile;

	//    if(numberOfPackings == numberOfPackingsStartValue){
	//	  logfile.open((char*)logFileNamePackings.c_str(), ios::app);
	//		int position=logfile.tellp();
	//		if(position==0){
	//		cout << "position" << endl;
	//		cout << logfile.tellp() << endl;

	//	  logfile << "seed" << "	N" << "	P0" << "	P" << "	alpha" << "	delta";
	////	  logfile << "	L" << "	phi" << "	Z" << "	#rattler" << "	s_xx" << "	s_yy";
	//	  logfile << "	s_xy" << "	U" << "	dU" << "	H" << "	dH" << "	t_run";
	//	  logfile << "	#FIRE" << "	#CG" << "	gg" << "total runtime" <<" creation-date"  << endl;
	//		}
	//      logfile.close();
	//   }
	//
	long double maxGrad = 0;
	iloop(2*N) {
		if (fabs(xihelper[i]) > maxGrad)
			maxGrad = fabs(xihelper[i]);
	}

	logfile.open((char*) logFileNamePackings.c_str(), ios::app);
	//logfile.open((char*)filenameString2.c_str(), ios::app);

	logfile << currentPackingNumber << "	" << N << "	" << P0 << "	" << P << "	"
			<< alpha << "	" << delta;
	logfile << "	" << L << "	" << phi << "	" << Z << "	" << N - Ncorrected
			<< "	" << sxx << "	" << syy;
	logfile << "	" << sxy << "	" << Uhelper << "	" << dU << "	" << H << "	"
			<< dH << "	" << timediff1;
	logfile << "	" << iterationcountfire << "	"
			<< iterationcountfrprmnCUMULATIVE << "	" << maxGrad << "	"
			<< timediff1 << "       " << timebuffer << endl;

	logfile.close();

	return;
}

// end writePositionFile

////////////////////////////////////////////////////////////////////////
// writeMultiplePackings
inline void writeMultiplePackings(string name) {
	time_t rawtime;
	struct tm *timeinfo;
	char timebuffer[80];
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(timebuffer, 80, "%Y-%m-%d_%H-%M-%S", timeinfo);

	ofstream outfile;

	outfile.open((char*) name.c_str(), ios::app);

	if (!outfile) {
		if (screenOutput)
			cout << "Cannot open output file!" << endl;
		outfile.close();
		ofstream outfile;
		createFileName();
		if (programmode == 5)
			outfile.open((char*) name.c_str(), ios::app);
		if (!outfile) {
			if (screenOutput)
				cout << "Cannot open output file!" << endl;
			return;
		} else if (screenOutput)
			cout << "Output file is opened succesfully (2nd attempt)!" << endl;
	} else if (screenOutput)
		cout << "Output file is opened succesfully!" << endl;
	outfile.setf(ios::fixed, ios::floatfield);
	outfile.precision(16);

	alpha = p[2 * N];
	delta = p[2 * N + 1];
	L = p[2 * N + 2];

	lxx = L / (1.0 + delta);
	lxy = L * 0.0;
	lyx = L * alpha;
	lyy = L * (1.0 + delta);

	outfile << "N = " << N << " ,L = " << L << " ,L1= { " << lxx << " , " << lxy
			<< " } " << " ,L2= { " << lyx << " , " << lyy << " } " << " ,P = "
			<< Phelper << " ,P0= " << P0 << " ," << endl;

	outfile << "{" << endl;
	iloop(N) {
		outfile << p[i] << " ,	" << p[N + i] << " ,	" << R[i] << " ,	" << endl;
	}
	outfile << "}" << endl;
	outfile << endl << endl;

	outfile.close();

	return;
}
// end writeMultiplePackings

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/// void packIntoBoundaries()
void packIntoBoundaries() {
	iloop(N) {
		while (p[N + i] < (0)) {
			p[N + i] += lyy;
		}
		while (p[N + i] > (lyy)) {
			p[N + i] -= lyy;
		}

		while (p[i] < (0 + p[N + i] * lyx / lyy))
			p[i] += lxx;
		while (p[i] > (lxx + p[N + i] * lyx / lyy))
			p[i] -= lxx;

	}

	energy();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
// createFileName
inline void createFileName() {

	char buffer;
	string bufferString;
	long double number, teni, teni1;
	int digit;
	int phiInt;

	time_t rawtime;
	struct tm *timeinfo;
	char timebuffer[80];
	char namebuffer[4];

	filenameString = "";

	namebuffer[0] = 48 + (((currentPackingNumber) / 1000) % 10);
	namebuffer[1] = 48 + (((currentPackingNumber) / 100) % 10);
	namebuffer[2] = 48 + (((currentPackingNumber) / 10) % 10);
	namebuffer[3] = 48 + ((currentPackingNumber) % 10);

	if (screenOutput)
		cout << "FilenameString: " << filenameString << endl;

	int i = 0;

	string filebase = nameOfWorkingDirectory;

	while (nameOfWorkingDirectory[i] != 'N') {
		i++;
		filebase = nameOfWorkingDirectory.substr(i);

		if (i > nameOfWorkingDirectory.size()) {
			cout << "Unexpected filename ERROR 2!" << endl;
			return;
		}
	}

	filenameString.append(filebase);
	filenameString.append("~");


	filenameString.push_back(namebuffer[0]);
	filenameString.push_back(namebuffer[1]);
	filenameString.push_back(namebuffer[2]);
	filenameString.push_back(namebuffer[3]);


	filenameString.append(".txt");

	filename = (char*) filenameString.c_str();
	if (screenOutput)
		cout << filename << endl;

	return;
}
// end createFileName

void initializeArrays() {

	// size all the dynamic 'vector' arrays according to the particle number N
	if (test)
		cout << "test0002: N = " << N << endl;
	dij.reserve(N * N);
	rij.reserve(N * N);
	xij.reserve(N * N);
	yij.reserve(N * N);
	nx.reserve(N * N);
	ny.reserve(N * N);
	neighbors.reserve(N * N);
	trueneighbors.reserve(N * N);
	trueneighborsOld.reserve(N * N);
	trueneighborsLast.reserve(N * N);
	trueneighborChanges.reserve(N * N);

	p.reserve(2 * N + 3);
	phelper.reserve(2 * N + 3);
	pcom.reserve(2 * N + 3);
	pLast.reserve(2 * N + 3);
	xi.reserve(2 * N + 3);
	xihelper.reserve(2 * N + 3);
	xicom.reserve(2 * N + 3);
	g.reserve(2 * N + 3);
	h.reserve(2 * N + 3);
	v.reserve(2 * N + 3);
	R.reserve(N);
	M.reserve(N + 3);

	numberOfDirectNeighbors.reserve(N);
	isRattler.reserve(N);
	wasRattler.reserve(N);

	// set initial particle positions and properties
	iloop(N) {
		M[i] = 1.0; // set particle 'masses' for FIRE-algorithm
		M[N] = M[N + 1] = M[N + 2] = N;
	}

	iloop(2*N) {
		v[i] = 0.0;
	}
	v[2 * N] = v[2 * N + 1] = v[2 * N + 2] = 0;

	ITMAXBRENT = max(100, 2 * N); // maximum of iterations in brent
	ITMAX = 5; // maximum of iterations in frprmn

	dtmax = dtmaxinit;
	dofOnOff = false;

	if (test)
		cout << "test0003: initializeArrays-END" << endl;
	return;
}

void checkFolderName(string foldername) {

	string filename = "/runLog.txt";
	;
	bool folderexists = false;

	ofstream testFolderName;
	filename.insert(0, foldername);
	testFolderName.open((char*) filename.c_str(), ios::app);

	if (testFolderName.is_open())
		folderexists = true;
	else
		cout << "ERROR" << endl;
	if (folderexists)
		testFolderName.close();

	return;
}
