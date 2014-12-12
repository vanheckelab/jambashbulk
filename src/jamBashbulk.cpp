// jamBashbulk: 2d particle simulation using conjugate gradient & FIRE algorithms
//
// indentation / style:
//     astyle --mode=c --style=linux --indent=spaces=4 --break-blocks=all --pad-oper --unpad-paren --add-brackets --convert-tabs --align-pointer=middle --lineend=linux
//     vim: set expandtab shiftwidth=4 softtabstop=4 tabstop=8:
//

// std libraries
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>

#include "fheader.h"

using namespace std;

// definitions:

//loops
#define iloop(upperbound) for(int i = 0; i < upperbound; i++)
#define jloop(upperbound) for(int j = 0; j < upperbound; j++)

#define ALPHA (p[2 * N])
#define DELTA (p[2 * N + 1])
#define LENGTH (p[2 * N + 2])

#define PRECISION (19)
// global variables and constants:

// simulation parameters
static int N = -1; // number of particles (per unit cell) in the packing
static int Ncorrected;
static const LDBL k = 1.0; // spring const. with respect to particle overlap

static bool screenOutput = false;
static bool debug = false;

// degrees of freedom of the periodic boundary unit cell in
static LDBL shear = 0.0;

static bool alphaOnOff = false, deltaOnOff = false, pressOnOff = false;
static bool dofOnOff = false;

// the previous parameters indicate in which way the simulated region
// can be deformed (these values change during execution of the
// program):
// alphaOnOff allows simple shear (alpha).
// deltaOnOff allows pure shear (delta).
// pressOnOff allows volume changes to reach given target pressure.

static bool alphaOnOffInit = true, deltaOnOffInit = true, pressOnOffInit = true;
// these parameters in indicate the user's choice of the degrees of
// freedom.

static string nameOfWorkingDirectory = "";
static int particleNumberLength = 0;

// start values:
static const LDBL alphainit = 0.0; // the user's choice for the initial shear angle
static const LDBL deltainit = 0.0; // the initial aspect-ratio
static const LDBL Linit = 5.0;
static LDBL P0init = 0.0;
static LDBL P0 = 0.0; // the target pressure
static LDBL phiinit = 0.8; // the initial fill fraction (0.869 ^= P~0.01)

static int countAlphaFlip = 0;
static int countDeltaFlip = 0;
static int countPressFlip = 0;

// straining:
static bool doSimpleShear = false; // switches on/off forced simple shear
static bool doCompression = false;
static bool fixedStepSize = false;

static int goalNumberOfContactChanges = 10; // the total number of contact changes we are interested in
static LDBL goalStrain = 0.1; // the strain range we are interested in
static int fixedStepNumber; // number of fixed size strain steps;

// iteration counters in various algorithms:
static int iterationcountSimStep = 0;
static int iterationcountfire = 0;
static int maxIterationCountFire = 2e9;

static int iterationcountfrprmnCUMULATIVE = 0;
static bool frprmnconverged = false;
static bool fireconverged = false;
static bool shearconverged = false;
static bool converged = false;
static int menumode = 0;

static int distributioncase = 0;
static bool endprogram = false;
static char * filename;
static string filenameString;
static int currentPackingNumber = 0;
static int numPackingsToProcess = 0;
static int firstPackingNumber = 0;

static bool redo = false;
static LDBL Rneighbor = 100.0;
static LDBL RneighborFrprmnLast;

static LDBL UhelperLastFunctionCall = 1e20;

static LDBL H = 1e9;
static LDBL HLastFunctionCall = 1e11;

static LDBL phi; // fill fraction
static LDBL Z; // average number of neighbors
static LDBL P; // pressure
static LDBL sxx, sxy, syy; // stress components

// overlap and distances
static vector<LDBL> dij; // NxN matrix of particle overlap
static vector<LDBL> rij; // NxN matrix of particle center distance
static vector<LDBL> xij; // NxN matrix of particle x-position difference
static vector<LDBL> yij; // NxN matrix of particle y-position difference

static vector<bool> neighbors;
static vector<bool> trueneighbors;
static vector<bool> trueneighborsOld;
static vector<bool> trueneighborsLast;
static vector<int> trueneighborChanges;
static vector<int> numberOfDirectNeighbors;
static vector<bool> isRattler;
static vector<bool> wasRattler;

static int trueneighborNumber = 0;

//generalized coordinates and gradients
static vector<LDBL> p; // Positions of particles
static vector<LDBL> pLast; // Positions at end of last simulationstep()
static vector<LDBL> xi; // Potential gradient
static vector<LDBL> g; // gradient helper
static vector<LDBL> h; // gradient helper
static vector<LDBL> R; // Radii of particles
static LDBL Rmax; // Maximum particle radius

// fire algorithm variables
static vector<LDBL> v; // Effective velocity
static vector<LDBL> M; // Effective masses
static LDBL FIRE_alpha;
static LDBL power;
static const int Nmin = 5;
static const LDBL finc = 1.1;
static const LDBL fdec = 0.5;
static const LDBL FIRE_alpha_start = 0.1; // called alpha_{start} in Jo's thesis
static const LDBL f_FIRE_alpha = 0.99;    // called f_{alpha} in Jo's thesis
static LDBL dt = 1e-1;
static const LDBL dtmaxinit = 1e-1;
static LDBL dtmax = dtmaxinit;
static const LDBL dtmin = 0.0;
static const LDBL damp = 1.0;
static LDBL dampalpha = 0.9;
static LDBL dampdelta = 0.9;
static LDBL damppress = 0.999;

//unit cell properties
static LDBL lxx, lxy; // x-/y- component of L_x (1st unit cell vector)
static LDBL lyx, lyy; // x-/y- component of L_y (2nd unit cell vector)
static vector<int> nx; // NxN matrix for periodicity calculation
static vector<int> ny; // NxN matrix for periodicity calculation

//helper variables for energy and gradient calculation of 'hypothetical' configuration
static LDBL Uhelper; // energy
static vector<LDBL> phelper; // particle positions (x_i,y_i)
static vector<LDBL> xihelper; // potential gradient (dU/dx_i,dU/dy_i)
static LDBL alphahelper, deltahelper, Lhelper, Phelper;
static LDBL lxxhelper, lxyhelper;
static LDBL lyxhelper, lyyhelper;

//LDBL gg; // gradient squared
static LDBL gg, vv;

// mathematical and program constants
static const LDBL PI = 3.141592653589793;
static const LDBL gold = 0.5 * (1.0 + sqrt(5.0)); // golden ratio = 1.618033988749894885
static const LDBL glimit = 100.0; // maximum magnification for parabolic-fit step in function mbrak

static const LDBL CGOLD = 0.3819660; // golden ratio for brent
static const LDBL ZEPS = 1e-25; // for brent
static int ITMAXBRENT = 50; // maximum of iterations in brent

static int ITMAX = 12; // maximum of iterations in frprmn
static const LDBL TOL = 1e-7; // tolerance passed to brent by linmin

static const LDBL ftol = 1e-2; // tolerance passed to frprmn()
static const LDBL ftolFIRE = 1e-17; // tolerance passed to fire()
static int endcount = 0;

static time_t starttime, endtime;
static LDBL timediff1 = 0;

static LDBL dU, dH;

// function declarations:
static void execute();
static void initializeSimulation();
static void initializeArrays();
static void simulationstep(); // this is where the simulation is performed
static void particledistance(int i, int j);
static void resethelpervars();
static LDBL energy();
static void gradientcalc();
static void mnbrak(LDBL * ax, LDBL * bx, LDBL * cx, LDBL * fa,
                   LDBL * fb, LDBL * fc, LDBL(*func)(LDBL));
static LDBL brent(LDBL ax, LDBL bx, LDBL cx,
                         LDBL(*f)(LDBL), LDBL tol, LDBL * xmin);
static LDBL SIGN(LDBL a, LDBL b);
static void linmin(int n, LDBL * fret, LDBL(*func)());
static LDBL f1dim(LDBL x);
static void frprmn(int n, LDBL * fret, LDBL(*func)());
static void calcSysPara();
static void menu();
static void readPositionFile();
static void writePositionFile();
static void writeMultiplePackings(string name);
static void saveDebugState();
static void packIntoBoundaries();
static void createFileName();
static void fire();
static void calcShearModulus();
static void checkNeighborChanges(int & addedcontacts, int & removedcontacts,
                                 int & neighborChanges, int & neighborChangesLast);
static void extractNandP(string foldername);
static void checkFolderName(string foldername);
static bool pnpoly(int nvert, LDBL * vertx, LDBL * verty, LDBL testx, LDBL testy);

enum PROGRAMMODE {
    PROGRAMMODE_DEFORM_PACKING = 3,
    PROGRAMMODE_CREATE_PACKING = 5
};

static PROGRAMMODE programmode;
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
    cout.precision(PRECISION);
    screenOutput = false;
    debug = false;

    for (int i=1; i<argc; i++) {
        if (strcmp(argv[i], "-screen") == 0) {
            cout << "> screen output activated \n";
            screenOutput = true;
        }
        else if (strcmp(argv[i], "-v") == 0) {
            cout << FILE_HEADER;
            return 0;
        }
        else if (strcmp(argv[i], "-debug") == 0) {
            cout << "> debug output activated \n";
            debug = true;
        }
   }

    starttime = time(NULL);

    menu();

    while(!endprogram) {
        execute();    // either this
    }

    return 0;
} // end main()

void extractNandP(string foldername)
{
    // Correct format:
    // <anything>/Packings/N<N>~P1e<P>
    // where <N> is an integer
    // and <P> is +<integer> or -<integer>

    int i = 0;
    int digit = 0;
    N = 0;
    P0 = 0;

    string comp2 = foldername.substr(i, 2);

    size_t pos = foldername.rfind("Packings"); // Name of the folder containing all the other folders with packings

    i = pos + 9;

    int i0 = i;

    if(foldername[i] != 'N') {
        cout << "Unexpected filename ERROR 1 : couldn't find N!" << endl;
        return;

    } else {
        i++;

        while(foldername[i] != '~') {
            N = N * 10;
            N = N + (foldername[i] - 48);
            i++;

            if(i - i0 > 10) {
                cout << "Unexpected filename ERROR 2 : N has too much digits!"
                     << endl;
                return;
            }
        }

        particleNumberLength = i - i0;

        if(foldername[i + 1] != 'P') {
            cout << "Unexpected filename ERROR 3 : no P in title" << endl;
            return;
        }

        int j = i + 2;
        int power = 1;

        while(foldername[j] != 'e') {
            digit = digit * 10;
            digit = digit + foldername[j] - 48;
            power = power - 1;
            j++;
        }

        if(foldername[j + 1] == '-') {
            P0 = digit * 1.0
                 * pow(10, -1.0 * (-power + foldername[j + 2] - 48));

        } else if(foldername[j + 1] == '+') {
            P0 = digit * 1.0
                 * pow(10, +1.0 * (power + foldername[j + 2] - 48));
        }

        if (screenOutput) {
            cout << "N=" << N << "; P0=" << P0;
        }

    }
} // extractNandP

////////////////////////////////////////////////////////////////////////////////
// execute()
void execute()
{
    bool goodfile;

    if(!endprogram) {
        if (debug) {
            saveDebugState();
        }

        if(!converged) {
            simulationstep();
            if (debug) {
                saveDebugState();
            }
        } else {
            if((iterationcountfire > maxIterationCountFire)
               || fabs((Phelper - P0) / P0) > 0.1 || Z < 3.5 || Z > 10) {
                goodfile = false;

            } else {
                goodfile = true;
            }

            endtime = time(NULL); // clock function runtime
            timediff1 = (LDBL)(endtime - starttime);
            starttime = endtime;

            if(screenOutput) {
                cout << "Total runtime is " << timediff1 << " seconds." << endl;
            }

            frprmnconverged = false;
            fireconverged = false;
            converged = false;

            if(programmode == PROGRAMMODE_CREATE_PACKING) {
                if(!doSimpleShear && !doCompression) {
                    writePositionFile();
                }

                iterationcountSimStep = 0;
                iterationcountfire = 0;
                iterationcountfrprmnCUMULATIVE = 0;

                if(goodfile) {
                    if(doSimpleShear)
                        if(!redo) {
                            int olditercount = maxIterationCountFire;
                            maxIterationCountFire = 1e6;
                            calcShearModulus();
                            maxIterationCountFire = olditercount;
                        }

                    if(doCompression) {
                        cout << "NOT IMPLEMENTED";
                        exit(2);
                    }
                }

                if(currentPackingNumber
                   < numPackingsToProcess - 1 + firstPackingNumber) {
                    currentPackingNumber++;
                    redo = false;

                    if(distributioncase == 2) {
                        initializeSimulation();
                    }

                    if(distributioncase == 3 || distributioncase == 4) {
                        readPositionFile();
                    }

                    converged = false;
                    fireconverged = false;

                } else {
                    endprogram = true;
                }

            } else {
                endprogram = true;
            }

        }
    } // if (!endprogram)

    return;
} // execute

void saveShearSystemState(string logFileName, int numberOfDataPoints,
                          string dataFileName, int neighborChangesLastCumulative,
                          int neighborChangesLast, int addedContacts, int removedContacts,
                          string GpositionFile, int numberOfContactChanges)
{

    char timebuffer[80];
    time_t rawtime;
    struct tm * timeinfo;
    ofstream outG;
    ofstream outLog;
    LDBL maxGrad = 0;

    iloop(2 * N) {
        if(fabs(xihelper[i]) > maxGrad) {
            maxGrad = fabs(xihelper[i]);
        }
    }

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(timebuffer, 80, "%Y-%m-%d_%H-%M-%S", timeinfo);
    outLog.open((char *)(logFileName.c_str()), ios::app);
    outLog.setf(ios::scientific, ios::floatfield);
    outLog.precision(PRECISION);

    outLog << numberOfDataPoints << "	" << N << "	" << P0 << "	" << P << "	"
           << ALPHA << "	" << DELTA;
    outLog << "	" << LENGTH << "	" << phi << "	" << Z << "	" << N - Ncorrected << "	"
           << sxx << "	" << syy;
    outLog << "	" << sxy << "	" << Uhelper << "	" << dU << "	" << H << "	" << dH
           << "	" << timediff1;
    outLog << "	" << iterationcountfire << "	" << iterationcountfrprmnCUMULATIVE
           << "	" << maxGrad << "	" << timebuffer << " " << numberOfContactChanges << endl;
    outLog.close();

    outG.open((char *)(dataFileName.c_str()), ios::app);
    outG.setf(ios::scientific, ios::floatfield);
    outG.precision(PRECISION);
    outG << shear << "	" << sxy << "	" << trueneighborNumber << "	"
         << (neighborChangesLastCumulative + neighborChangesLast) << "	"
         << addedContacts << "	" << removedContacts << "	" << Phelper << "	"
         << Z << endl;
    outG.close();
    writeMultiplePackings(GpositionFile);
}

void saveDebugState() {
    // Save the system state to {debug.log, debug.data, debug.positions} in the current
    // directory. This uses the Shear output routines, so it may include some bogus data
    // (e.g. number of rearrangements)
    static int counter = 0;
    counter += 1;
    saveShearSystemState("debug.log", counter, "debug.data", 0, 0, 0, 0, "debug.positions", 0);
}

void gotoAlphaShear(LDBL targetAlpha) {
    ALPHA = targetAlpha;

    iterationcountSimStep = 0;
    iterationcountfire = 0;
    iterationcountfrprmnCUMULATIVE = 0;

    shearconverged = false;
    converged = false;
    frprmnconverged = false;
    fireconverged = false;

    while(!shearconverged) {
        if (debug) {
            saveDebugState();
        }
        simulationstep();
    }
    if (debug) {
        saveDebugState();
    }

    calcSysPara();
}

////////////////////////////////////////////////////////////////////////
// calcShearModulus()
void calcShearModulus()
{
    PROGRAMMODE programmodeOld = programmode;
    int num = 0;

    LDBL shearfactor = 10; //= sqrt(10); // for fast calculation = 10 else = sqrt(10)
    LDBL dstrain = goalStrain / fixedStepNumber;
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

    vector<LDBL> pAfterChange;
    pAfterChange.reserve(2 * N + 3);

    ofstream outG;
    ofstream outLog;
    ofstream outFirst;
    ofstream eraseFile;
    string dataFileName = filenameString;
    string logFileName = filenameString;
    string GpositionFile = filenameString;
    string Appendix = "";

    LDBL goalStrainHelper = goalStrain;
    int goalStrainExponent = 0, goalStraindigit = 0;
    string goalStrainString = "";

    // Reset shear value
    shear = 0.0;

    while(goalStrainHelper * 1.01 < 1.0) {
        goalStrainExponent++;
        goalStrainHelper *= 10.0;
    }

    goalStraindigit = (int)goalStrainHelper;

    goalStrainString.push_back((char)(goalStraindigit + 48));
    goalStrainString.append("e-");
    goalStrainString.push_back((char)(goalStrainExponent + 48));

    if(!fixedStepSize) {
        Appendix = "~SR";
        Appendix.push_back((char)(((goalNumberOfContactChanges / 100) % 10) + 48));
        Appendix.push_back((char)(((goalNumberOfContactChanges / 10) % 10) + 48));
        Appendix.push_back((char)(((goalNumberOfContactChanges / 1) % 10) + 48));
        Appendix.append("~step");
        Appendix.push_back((char)(
            (((2 * goalNumberOfContactChanges + 1) / 100) % 10) + 48));
        Appendix.push_back((char)(
            (((2 * goalNumberOfContactChanges + 1) / 10) % 10) + 48));
        Appendix.push_back((char)(
            (((2 * goalNumberOfContactChanges + 1) / 1) % 10) + 48));

    } else {

        Appendix = "~SS";
        Appendix.append(goalStrainString);
        Appendix.append("~step");
        Appendix.push_back((char)((((fixedStepNumber) / 1000) % 10) + 48));
        Appendix.push_back((char)((((fixedStepNumber) / 100) % 10) + 48));
        Appendix.push_back((char)((((fixedStepNumber) / 10) % 10) + 48));
        Appendix.push_back((char)((((fixedStepNumber) / 1) % 10) + 48));
    }

    GpositionFile.insert(particleNumberLength + 6, Appendix);
    GpositionFile.insert(0, "particles");
    GpositionFile.insert(0, nameOfWorkingDirectory + "/");

    dataFileName.insert(particleNumberLength + 6, Appendix);
    dataFileName.insert(0, "data");
    dataFileName.insert(0, nameOfWorkingDirectory + "/");

    logFileName.insert(particleNumberLength + 6, Appendix);
    logFileName.insert(0, "log");
    logFileName.insert(0, nameOfWorkingDirectory + "/");

    outG.open((char *) dataFileName.c_str(), ios::trunc);
    outG.setf(ios::scientific, ios::floatfield);
    outG.precision(PRECISION);
    outG << FILE_HEADER;
    outG << "gamma_alpha	s_xy	Ncontacts	Nchanges	N+	N-	P	Z" << endl;
    outG.close();

    outLog.open((char *) logFileName.c_str(), ios::trunc);
    outLog.setf(ios::scientific, ios::floatfield);
    outLog.precision(PRECISION);
    outLog << FILE_HEADER;
    outLog << "step#" << "	N" << "	P0" << "	P" << "	alpha" << "	delta";
    outLog << "	L" << "	phi" << "	Z" << "	#rattler" << "	s_xx" << "	s_yy";
    outLog << "	s_xy" << "	U" << "	dU" << "	H" << "	dH" << "	t_run";
    outLog << "	#FIRE" << "	#CG" << "	gg" << " creation-date" << " ccnum" << endl;
    outLog.close();

    LDBL alphaBeforeDeformation = ALPHA;
    iloop(N) {
        jloop(N) {
            trueneighborsLast[j * N + i] = trueneighborsOld[j * N + i] =
                                               trueneighbors[j * N + i]; // the contacts before shearing
            trueneighborChanges[j * N + i] = 0;
        }
        wasRattler[i] = isRattler[i];
    }

    programmode = PROGRAMMODE_DEFORM_PACKING; // program-mode for shearing

    jloop(2 * N + 3) {
        pLast[j] = p[j];
    } // backup the particle position before the shear step

    packIntoBoundaries();

    eraseFile.open((char *) GpositionFile.c_str(), ios::trunc);
    eraseFile <<  FILE_HEADER;
    eraseFile.close();

    saveShearSystemState(logFileName, 0,
                         dataFileName, 0,
                         0, 0, 0,
                         GpositionFile, numberOfContactChanges);

    if(!fixedStepSize) {
        // set shear to an initial (small) value > 0, or multiplication will yield 0.
        // for fast calculation = 1e-12 else = 1e-16
        shear = 1e-9;
    }

    while(!reachedGoal) {  // shear < 0.15 && numberOfDataPoints < 5000 &&
        neighborChangesLastCumulative += neighborChangesLast;

        iloop(N) {
            jloop(N) {
                trueneighborsLast[j * N + i] = trueneighbors[j * N + i]; // the contacts before shearing step

            }
            wasRattler[i] = isRattler[i];
        }

        neighborChangesOld = neighborChanges;

        while(!sufficientAccuracy) {
            if(!fixedStepSize) {
                if(pastContactChange) {
                    num++;

                    jloop(2 * N + 3) {
                        pAfterChange[j] = p[j];
                    } // save positions just after rearrangement
                    jloop(2 * N + 3) {
                        p[j] = pLast[j];
                    } // go back to last particle positions

                    if (numberOfDataPoints == 1) {
                        // if this is the first contact change, we have to manually go back to
                        // a very low shear value - or the smallest contact change shear we can
                        // detect is shear/shearfactor!
                        // instead, we reset the shear to a low value and reset the simulation
                        shear = 1e-16;
                        num = 0;
                        numberOfContactChanges = addedContacts = removedContacts = 0;
                        neighborChanges = neighborChangesOld = neighborChangesLast = neighborChangesLastCumulative = 0;
                        pastContactChange = false;
                        continue;
                    } else {
                        shear = shear / shearfactor; // go back to previous shear
                        shearfactor = sqrt(shearfactor);
                    }

                } else if(num > 0) {
                    shearfactor = sqrt(shearfactor);
                }

                shear = shear * shearfactor; // shear is increased by shearfactor...

            } else {
                shear = shear + dstrain;
            }

            jloop(2 * N + 3) {
                pLast[j] = p[j];
            } // backup the particle position before the shear step

            gotoAlphaShear(alphaBeforeDeformation + shear); // ... and added to the shear of the relaxed packing

            if(iterationcountfire > maxIterationCountFire) {
                programmode = programmodeOld;
                return;
            }

            checkNeighborChanges(addedContacts, removedContacts,
                                 neighborChanges, neighborChangesLast);

            if(neighborChangesLast != 0) {
                pastContactChange = true;

            } else {
                pastContactChange = false;
            }

            if(!fixedStepSize) {
                if(((neighborChanges - neighborChangesOld)
                    * (neighborChanges - neighborChangesOld) == 1
                    && (shearfactor - 1.0) < 0.0001)
                   || ((shearfactor - 1.0) < 1e-6)) {
                    sufficientAccuracy = true;
                }
            }

            endtime = time(NULL); // clock function runtime
            timediff1 = (LDBL)(endtime - starttime);
            starttime = endtime;

            numberOfDataPoints++;
            energy();
            calcSysPara();

            saveShearSystemState(logFileName, numberOfDataPoints,
                                 dataFileName, neighborChangesLastCumulative,
                                 neighborChangesLast, addedContacts, removedContacts,
                                 GpositionFile, numberOfContactChanges);



            if(fixedStepSize) {
                if(numberOfDataPoints < fixedStepNumber) {
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





        if(!fixedStepSize) {
            // IF we end BEFORE the CC, we need to reset to the POST-CC state
            // WE CANNOT DO THIS OTHERWISE, as pAfterChange is set in the BEGINNING of the loop
            // but I'm afraid moving will break other stuff...

            // also, we must NOT relax the packing further, as that might throw us into a different
            // state... (if the original state was not as far relaxed as it could have been)
            bool moved_to_after_cc = false;
            if (!pastContactChange) {
                moved_to_after_cc = true;
                // pAfterChange is saved... some 100 lines above here
                jloop(2 * N + 3) {
                    p[j] = pAfterChange[j];
                    phelper[j] = pAfterChange[j];
                }

                shear = ALPHA - alphaBeforeDeformation;
                fireconverged = true;
                energy();
                calcSysPara();
                checkNeighborChanges(addedContacts, removedContacts,
                                    neighborChanges, neighborChangesLast);
            }

            cout << "\nConverged to CC #" << numberOfContactChanges
                 << " (step #" << numberOfDataPoints << ")"
                 << " with " << neighborChangesLast << " changed contact" << (neighborChangesLast == 1 ? "" : "s")
                 << " at strain " << shear
                 << (moved_to_after_cc ? " [Jump]" : " [Direct]");

            // we should now be post contact change!
            if (neighborChangesLast == 0) {
                cout << "\n\nERROR: Not post contact change. Giving up.";
                reachedGoal = true;
                return;
            }

            pastContactChange = false;
            sufficientAccuracy = false;
            num = 0;
            shearfactor = sqrt(sqrt(sqrt(sqrt(10)))); //sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(10)))))));
            numberOfContactChanges++;
        }

        if(!fixedStepSize) {
            if(numberOfContactChanges < goalNumberOfContactChanges) {
                reachedGoal = false;

            } else {
                reachedGoal = true;
            }
        }
    } // end while

    programmode = programmodeOld;
    return;
} // end calcShearModulus

////////////////////////////////////////////////////////////////////////
// checkNeighborChanges()

// NB calcSysPara() needs to be called before this, and fireconverged must
// be true to correctly count rattlers!
void checkNeighborChanges(int & addedContacts, int & removedContacts,
                          int & neighborChanges, int & neighborChangesLast)
{
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

            if(!isRattler[j] && !isRattler[i]) {
                if(trueneighborChangesLast > 0) {
                    neighborChangesLast++;
                    addedContacts++;
                }
            }

            if(!wasRattler[j] && !wasRattler[i]) {
                if(trueneighborChangesLast < 0) {
                    neighborChangesLast++;
                    removedContacts++;
                }
            }
        }
    }
} // end checkNeighborChanges

void calculate_neighbors()
{
    trueneighborNumber = 0;
    iloop(N) {
        jloop(i) {
            neighbors[j * N + i] = true;

        }
        neighbors[i * N + i] = false;
    }
    iloop(2 * N + 3) {
        phelper[i] = p[i];
    }
    energy();
    gradientcalc();
    iloop(N) {
        jloop(i) {
            if(rij[j * N + i] < Rneighbor) {
                neighbors[j * N + i] = true;

            } else {
                neighbors[j * N + i] = false;
            }

            if(trueneighbors[j * N + i] && !isRattler[i] && !isRattler[j]) {
                trueneighborNumber++;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////
// Simulation step
void simulationstep()
{
    LDBL fret;

    time_t rawtime1;

    lxx = LENGTH / (1.0 + DELTA);
    lxy = LENGTH * 0.0;
    lyx = LENGTH * ALPHA;
    lyy = LENGTH * (1.0 + DELTA);

    packIntoBoundaries();

    if((programmode == PROGRAMMODE_CREATE_PACKING) && frprmnconverged && Rneighbor / Rmax < 2.41) {
        if((iterationcountfire % 10000) == 0) {
            calculate_neighbors();
        }

    } else {
        calculate_neighbors();
    }

    iloop(2 * N + 3) {
        phelper[i] = p[i];
    }
    energy();
    gradientcalc();

    if(programmode == PROGRAMMODE_CREATE_PACKING) {
        if(!frprmnconverged) {
            Rneighbor = 3.5 * Rmax
                        + 0.3 * LENGTH / pow(2.0, iterationcountSimStep * 0.1);
            RneighborFrprmnLast = Rneighbor;

        } else {
            Rneighbor = 3.5 * Rmax
                        + (RneighborFrprmnLast - 2.4 * Rmax)
                        / pow(2.0, (iterationcountfire) * 0.001);

        }
    }

    if(distributioncase == 3 || distributioncase == 4) {
        Rneighbor = 3.5 * Rmax;
        dofOnOff = true;
    }

    if(programmode == PROGRAMMODE_DEFORM_PACKING) {
        Rneighbor = 3.5 * Rmax;
    }

    if(programmode == PROGRAMMODE_CREATE_PACKING) {
        if(!frprmnconverged) {
            alphaOnOff = false;
            deltaOnOff = false;
            pressOnOff = false;

            frprmn(N, &fret, energy);

        } else {
            dofOnOff = true;

            if(dofOnOff) {

                if(alphaOnOffInit) {
                    alphaOnOff = true;
                }

                if(deltaOnOffInit) {
                    deltaOnOff = true;
                }

                if(pressOnOffInit) {
                    pressOnOff = true;
                }

            } else {
                alphaOnOff = false;
                deltaOnOff = false;
                pressOnOff = false;
            } // end else

            fire();

            if(2.0 * fabs(H - HLastFunctionCall)
               < ftolFIRE * (fabs(H) + fabs(HLastFunctionCall) + ZEPS)) {
                if(fabs(sxy) < 1e-15) {
                    fireconverged = true;
                }

            } else {
                endcount = 0;
            }
        } // end else

        if(fireconverged) {
            converged = true;
        }

        energy();
        gradientcalc();
    }

    if(programmode == PROGRAMMODE_DEFORM_PACKING) {
        dofOnOff = false;
        alphaOnOff = false;
        deltaOnOff = false;
        pressOnOff = false;

        shearconverged = false;
        converged = false;
        frprmnconverged = false;
        fireconverged = false;

        fire();

        if(2.0 * fabs(H - HLastFunctionCall)
           < 1e-13 * (fabs(H) + fabs(HLastFunctionCall) + ZEPS)) {
            fireconverged = true;
        }

        if(fireconverged) {
            converged = true;
            shearconverged = true; // NB fireconverged MUST be true before setting shearconverged=true
                                   // or contact counting will not work correctly
        }
    }

    calcSysPara();

    if(fireconverged) {
        calculate_neighbors();
    } // end if

    time(&rawtime1);

    if(dofOnOff) {
        int lowerswitch = 10;
        int upperswitch = 50;
        double dampDOF = 0.85;

        if(countAlphaFlip < lowerswitch) {
            dampalpha = 1.01 * (dampalpha);
        }

        if(countAlphaFlip > upperswitch) {
            dampalpha = 0.99 * dampalpha;
        }

        if(countDeltaFlip < lowerswitch) {
            dampdelta = 1.01 * (dampdelta);
        }

        if(countDeltaFlip > upperswitch) {
            dampdelta = 0.99 * dampdelta;
        }

        if(countPressFlip < lowerswitch) {
            damppress = 1.01 * (damppress);
        }

        if(countPressFlip > upperswitch) {
            damppress = 0.99 * damppress;
        }

        if(dampalpha > 1.0) {
            dampalpha = 1.0;
        }

        if(dampdelta > 1.0) {
            dampdelta = 1.0;
        }

        if(damppress > 1.0) {
            damppress = 1.0;
        }

        if(dampalpha < dampDOF) {
            dampalpha = dampDOF;
        }

        if(dampdelta < dampDOF) {
            dampdelta = dampDOF;
        }

        if(damppress < dampDOF) {
            damppress = dampDOF;
        }
    }

    M[N] = M[N + 1] = N * N;
    M[N + 2] = sqrt(N);

    if(2.0 * fabs(H - HLastFunctionCall)
       < 1e-14 * (fabs(H) + fabs(HLastFunctionCall) + ZEPS)) {
        dtmax = dt * 0.95;

    } else {
        dtmax = 0.6;
    }

    string filepath = nameOfWorkingDirectory + "/" + "errorLog.txt";

    if(iterationcountfire > maxIterationCountFire) {
        fireconverged = true;
        shearconverged = true;

        ofstream errorLog;
        errorLog.open((char *) filepath.c_str(), ios::app);

        createFileName();
        errorLog << filenameString
                 << " did not converge within maximum amount of iterations."
                 << endl;
    }

    countAlphaFlip = 0;
    countDeltaFlip = 0;
    countPressFlip = 0;

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
void menu()
{
    iterationcountSimStep = 0;
    iterationcountfire = 0;
    programmode = PROGRAMMODE_CREATE_PACKING;
    char answer;

    converged = false;
    frprmnconverged = false;
    fireconverged = false;
    alphaOnOff = deltaOnOff = false;

    doSimpleShear = false;

    string fileToOpen = "";

    if(screenOutput) {
        cout << "Please specify WORKING DIRECTORY." << endl;
    }

    cin >> nameOfWorkingDirectory;
    checkFolderName(nameOfWorkingDirectory);
    extractNandP(nameOfWorkingDirectory);

    if(screenOutput) {
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

    switch(menumode) {

    case 1:
        alphaOnOff = false;
        deltaOnOff = false;

        if(screenOutput) {
            cout << "Please insert number of packings to be created!" << endl;
        }

        cin >> numPackingsToProcess;

        if(screenOutput) {
            cout << "Please insert starting number of packing names!" << endl;
        }

        cin >> firstPackingNumber;
        currentPackingNumber = firstPackingNumber;
        break;

    case 2:
        if(screenOutput) {
            cout << "Is simple shear (alpha) a degree of freedom? Y/N :";
        }

        cin >> answer;

        if(answer == 'Y' || answer == 'y') {
            alphaOnOffInit = true;

        } else {
            alphaOnOffInit = false;
        }

        if(screenOutput) {
            cout << "Is pure shear (delta) a degree of freedom? Y/N :";
        }

        cin >> answer;

        if(answer == 'Y' || answer == 'y') {
            deltaOnOffInit = true;

        } else {
            deltaOnOffInit = false;
        }

        if(screenOutput) {
            cout << "Equilibrate at target pressure? Y/N :";
        }

        cin >> answer;

        if(answer == 'Y' || answer == 'y') {
            pressOnOffInit = true;

        } else {
            pressOnOffInit = false;

            if(screenOutput) {
                cout << "Please insert desired fill fraction:";
            }

            cin >> phiinit;
        }

        alphaOnOff = false;
        deltaOnOff = false;

        if(screenOutput) {
            cout << "Please insert number of packings to be created!" << endl;
        }

        cin >> numPackingsToProcess;
        currentPackingNumber = 0;
        break;

    case 3:
        doSimpleShear = true;
        fixedStepSize = false;
        alphaOnOff = false;
        deltaOnOff = false;

        if(screenOutput) {
            cout << "Please insert number of packings to be created!" << endl;
        }

        cin >> numPackingsToProcess;

        if(screenOutput) {
            cout << "Please insert starting number of packing names!" << endl;
        }

        cin >> firstPackingNumber;

        if(screenOutput) {
            cout << "Please insert goal number of rearrangements!" << endl;
        }

        cin >> goalNumberOfContactChanges;

        currentPackingNumber = firstPackingNumber;
        break;

    case 4:
        doSimpleShear = true;
        fixedStepSize = true;
        alphaOnOff = false;
        deltaOnOff = false;

        if(screenOutput) {
            cout << "Please insert number of packings to be created!" << endl;
        }

        cin >> numPackingsToProcess;

        if(screenOutput) {
            cout << "Please insert starting number of packing names!" << endl;
        }

        cin >> firstPackingNumber;

        if(screenOutput) {
            cout << "Please insert strain goal!" << endl;
        }

        cin >> goalStrain;

        if(screenOutput) {
            cout << "Please insert number of equidistant strain steps!" << endl;
        }

        cin >> fixedStepNumber;
        currentPackingNumber = firstPackingNumber;
        break;

    case 5:
        doCompression = true;
        fixedStepSize = false;
        alphaOnOff = false;
        deltaOnOff = false;

        if(screenOutput) {
            cout << "Please insert number of packings to be created!" << endl;
        }

        cin >> numPackingsToProcess;

        if(screenOutput) {
            cout << "Please insert starting number of packing names!" << endl;
        }

        cin >> firstPackingNumber;

        if(screenOutput) {
            cout << "Please insert goal number of rearrangements!" << endl;
        }

        if(screenOutput) {
            cout << "(neg. values for compression, pos. decompression)" << endl;
        }

        cin >> goalNumberOfContactChanges;
        currentPackingNumber = firstPackingNumber;
        break;

    case 6:
        doCompression = true;
        fixedStepSize = true;
        alphaOnOff = false;
        deltaOnOff = false;

        if(screenOutput) {
            cout << "Please insert number of packings to be created!" << endl;
        }

        cin >> numPackingsToProcess;

        if(screenOutput) {
            cout << "Please insert starting number of packing names!" << endl;
        }

        cin >> firstPackingNumber;

        if(screenOutput) {
            cout << "Please insert strain goal!" << endl;
        }

        if(screenOutput) {
            cout << "(neg. values for compression, pos. decompression)" << endl;
        }

        cin >> goalStrain;

        if(screenOutput) {
            cout << "Please insert number of equidistant strain steps!" << endl;
        }

        cin >> fixedStepNumber;
        currentPackingNumber = firstPackingNumber;
        break;

    case 7:
        break;

    case 9:
        endprogram = true;
        return;
        break;

    }

    if(screenOutput) {
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

    switch(distributioncase) {
    case 1:

        if(screenOutput) {
            cout << "Please insert target pressure P0: ";
        }

        cin >> P0;
        break;
    case 2:

        initializeSimulation();
        resethelpervars();
        break;
    case 3:
        readPositionFile();
        break;

    case 4:

        if(screenOutput) {
            cout << "Please specify target pressure: ";
        }

        cin >> P0init;
        readPositionFile();
        break;

    case 5:

        if(screenOutput) {
            cout << "Filename: ";
        }

        cin >> fileToOpen;
        filenameString = fileToOpen;

        if(screenOutput) {
            cout << filenameString << endl;
        }

        readPositionFile();
        simulationstep();
        break;

    case 9:
        endprogram = true;
        break;

    default:
        cout << "No value given; aborting";
        exit(3);

    }

    if(programmode == PROGRAMMODE_DEFORM_PACKING) {
        energy();
    }

    return;
}
// end menu

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// fire
void fire()
{
    // Velocity Verlet algorithm
    //v += 0.5*(aold+a)*dt
    //p += v*dt+0.5*aold*dt*dt;

    // if velocities are irrelevant and we only want the equilibrium position...
    //v += a*dt
    //p += v*dt

    LDBL iterPosPower = 0;
    int countAlpha = 0;
    int countDelta = 0;
    int countPress = 0;
    int itercount = 0;

    LDBL frac;
    LDBL fracAlpha;
    LDBL fracDelta;

    FIRE_alpha = FIRE_alpha_start;

    while(itercount < 1000) {
        energy(); // calculate overlaps of new configuration
        gradientcalc(); // calulate the gradient for the new configuration

        H = Uhelper + P0 * Lhelper * Lhelper;

        if(programmode == PROGRAMMODE_DEFORM_PACKING) {
            calcSysPara();
        }

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

        if(power > 0 && iterPosPower > Nmin) {
            dt = min(dt * finc, dtmax);
            FIRE_alpha *= f_FIRE_alpha;
        }

        if(power < 0) {
            if(fabs(sxy) < 1e-16) {
                dt = max(dt * 0.2 * fdec, dtmin);

            } else {
                dt = max(dt * fdec, dtmin);
            }

            iloop(2 * N) {
                v[i] = 0.0;
            }
            FIRE_alpha = FIRE_alpha_start;
            iterPosPower = 0;
        }

        double radiusFraction = 1e-5;

        // only vary alpha and delta IF
        bool varyAlphaDelta =
        // L has converged far enough OR
            (fabs(Phelper - P0) / P0 < 5e-1) ||
        // L cannot change
            (!pressOnOff);


        if(alphaOnOff && varyAlphaDelta) {
            v[2 * N] = v[2 * N] * dampalpha - xihelper[2 * N] / M[N] * dt; // damping prevents too large changes

            if(v[2 * N] * Lhelper * dt > radiusFraction) {
                v[2 * N] = radiusFraction / Lhelper / dt;

            } else if(v[2 * N] * Lhelper * dt < -radiusFraction) {
                v[2 * N] = -radiusFraction / Lhelper / dt;
            }

            phelper[2 * N] = ALPHA = ALPHA + v[2 * N] * dt;
            fracAlpha = ALPHA - alphahelper;

        } else {
            fracAlpha = 0.0;
            v[2 * N] = 0.0;
        }

        if(deltaOnOff && varyAlphaDelta) {
            v[2 * N + 1] = v[2 * N + 1] * dampdelta
                           - xihelper[2 * N + 1] / M[N + 1] * dt;

            if(v[2 * N + 1] * Lhelper * dt > radiusFraction) {
                v[2 * N + 1] = radiusFraction / Lhelper / dt;

            } else if(v[2 * N + 1] * Lhelper * dt < -radiusFraction) {
                v[2 * N + 1] = -radiusFraction / Lhelper / dt;
            }

            phelper[2 * N + 1] = DELTA = DELTA + v[2 * N + 1] * dt;
            fracDelta = (1 + DELTA) / (1 + deltahelper);

        } else {
            fracDelta = 1.0;
            v[2 * N + 1] = 0.0;
        }

        if(pressOnOff) {
            v[2 * N + 2] = v[2 * N + 2] * damppress
                           - xihelper[2 * N + 2] / M[N + 2] * dt;

            if(v[2 * N + 2] * dt > radiusFraction) {
                v[2 * N + 2] = radiusFraction / dt;

            } else if(v[2 * N + 2] * dt < -radiusFraction) {
                v[2 * N + 2] = -radiusFraction / dt;
            }

            phelper[2 * N + 2] = p[2 * N + 2] = p[2 * N + 2]
                                                + v[2 * N + 2] * dt;
            frac = p[2 * N + 2] / Lhelper;

        } else {
            frac = 1.0;
            v[2 * N + 2] = 0.0;
        }

        iloop(N) {

            v[i] = v[i] * damp - xihelper[i] / M[i] * dt; // integrate accelerations to find velocities
            v[N + i] = v[N + i] * damp - xihelper[N + i] / M[i] * dt; // damping added

            v[i] = (1.0 - FIRE_alpha) * v[i] - FIRE_alpha * (xihelper[i] / gg * vv); // FIRE step
            v[N + i] = (1.0 - FIRE_alpha) * v[N + i]
                       - FIRE_alpha * (xihelper[N + i] / gg * vv); // FIRE step

            phelper[i] = phelper[i] + phelper[N + i] * fracAlpha;
            phelper[i] = phelper[i] / fracDelta;
            phelper[N + i] = phelper[N + i] * fracDelta;

            phelper[i] = p[i] = (phelper[i] + v[i] * dt) * frac; // integrate velocities to find positions
            phelper[N + i] = p[N + i] = (phelper[N + i] + v[N + i] * dt) * frac;

        }

        ////////////

        if(dofOnOff) {

            if(xihelper[2 * N] * v[2 * N] > 0) {
                v[2 * N] = 0.0;
                countAlphaFlip++;

            } else {
                countAlpha++;
            }

            if(xihelper[2 * N + 1] * v[2 * N + 1] > 0) {
                v[2 * N + 1] = 0.0;
                countDeltaFlip++;

            } else {
                countDelta++;
            }

            if(xihelper[2 * N + 2] * v[2 * N + 2] > 0) {
                countPressFlip++;

            } else {
                countPress++;
            }
        }

        ////////////

        iterPosPower++;
        iterationcountfire++;
        itercount++;
    } // end while
} // end FIRE

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// mnbrak
// Given two initial points and a function, the algorithm finds another point such that
// a bracketing triple is formed.
//
// Returns (ax, bx, cx) such that ax < bx < cx   and  f(ax) > f(bx) < f(cx)
// Returns (fa, fb, fc) = f(ax), f(bx), f(cx)
//
// IN: LDBL *ax, *bx               two initial points
//     LDBL (*func)(LDBL)  (i.e. a pointer to a function that takes
//                                          one LDBL argument and returns
//                                          one LDBL value)
// OUT: LDBL *ax, *bx, *cx
//      LDBL *fa, *fb, *fc
//

#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);

void mnbrak(LDBL * ax, LDBL * bx, LDBL * cx, LDBL * fa,
            LDBL * fb, LDBL * fc, LDBL(*func)(LDBL))
{

    //permutation procedure
    LDBL ulim, u, r, q, fu, dum;
    int numberOfIterations = 0;

    *fa = (*func)(*ax);
    *fb = (*func)(*bx);

    if(*fb > *fa) {
        SHFT(dum, *ax, *bx, dum)
        SHFT(dum, *fb, *fa, dum)
    } // end if

    *cx = (*bx) + gold * (*bx - *ax);
    *fc = (*func)(*cx);

    while(*fb > *fc) {
        numberOfIterations++;

        r = (*bx - *ax) * (*fb - *fc);
        q = (*bx - *cx) * (*fb - *fa);
        u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) / (2.0 * (q - r));
        // possible div by 0 !!!!!

        ulim = (*bx) + glimit * (*cx - *bx); // maximum step

        // "Where is u?"
        if((*bx - u) * (u - *cx) > 0.0) {  // IF_1
            fu = (*func)(u);

            if(fu < *fc) {  // IF_2
                *ax = (*bx);
                *bx = u;
                *fa = (*fb);
                *fb = fu;
                //                      if(screenOutput) cout << "Number of iterations in mnbrak():  " << numberOfIterations << endl;
                return;
            } // end IF_2

            else if(fu > *fb) {
                *cx = u;
                *fc = fu;
                //                     if(screenOutput) cout << "Number of iterations in mnbrak(): " << numberOfIterations << endl;
                return;
            } // end ELSE IF

            u = (*cx) + gold * (*cx - *bx);
            fu = (*func)(u);
        } // end IF_1

        else if((*cx - u) * (u - ulim) > 0.0) {  // ELSE IF_1
            fu = (*func)(u);

            if(fu < *fc) {
                SHFT(*bx, *cx, u, *cx + gold * (*cx - *bx))
                SHFT(*fb, *fc, fu, (*func)(u))
            } // end IF
        } // end ELSE IF_1

        else if((u - ulim) * (ulim - *cx) > 0.0) {  // ELSE IF_2
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
#undef SHFT

// end mnbrak

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// brent
// Brent's method (see http://en.wikipedia.org/wiki/Brent's_method )

#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
LDBL brent(LDBL ax, LDBL bx, LDBL cx,
                  LDBL(*f)(LDBL), LDBL tol, LDBL * xmin)
{

    int iter;
    LDBL a, b, d, etemp, fu, fv, fw, fx, s, q, r, tol1, tol2, u, v, w, x,
         xm;
    LDBL e = 0.0;

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = w = v = bx;
    fw = fv = fx = (*f)(x);

    for(iter = 1; iter <= ITMAXBRENT; iter++) {
        xm = 0.5 * (a + b);
        tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);

        if(fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
            *xmin = x;
            return fx;
        }

        if(fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            s = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);

            if(q > 0.0) {
                s = -s;
            }

            q = fabs(q);
            etemp = e;
            e = d;

            if(fabs(s) >= fabs(0.5 * q * etemp) || s <= q * (a - x)
               || s >= q * (b - x)) {
                d = CGOLD * (e = (x >= xm ? a - x : b - x));

            } else {
                d = s / q;
                u = x + d;

                if(u - a < tol2 || b - u < tol2) {
                    d = SIGN(tol1, xm - x);
                }

            }

        } else {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }

        u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
        fu = (*f)(u);

        if(fu <= fx) {
            if(u >= x) {
                a = x;

            } else {
                b = x;
            }

            SHFT(v, w, x, u)
            SHFT(fv, fw, fx, fu)

        } else {
            if(u < x) {
                a = u;

            } else {
                b = u;
            }

            if(fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;

            } else if(fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }

    }

    *xmin = x;
    return fx;
}

#undef SHFT

// end brent

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
LDBL SIGN(LDBL a, LDBL b)
{
    if(b < 0.0) {
        return -a;

    } else {
        return a;
    }
} // end SIGN

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// linmin
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void linmin(int n, LDBL * fret, LDBL(*func)())
{
    static const LDBL AMIN = 1e-7; // starting step in linmin
    LDBL xx, xmin, fx, fb, fa, bx, ax;

    ax = 0.0;
    xx = AMIN;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
    *fret = brent(ax, xx, bx, f1dim, TOL, &xmin);

    iloop(2 * N + 3) {
        p[i] += xi[i] * xmin;
    } // move by gradient*xmin to the 1D-minimum
    return;
} // end linmin
#pragma GCC diagnostic pop
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// f1dim
LDBL f1dim(LDBL x)
{
    LDBL Uloc;

    iloop(2 * N + 3) {
        phelper[i] = p[i] + x * xi[i];
    } // move by gradient*x ...
    energy();
    Uloc = Uhelper;

    if(pressOnOff) {
        Uloc = Uloc + P0 * Lhelper * Lhelper;
    }

    return Uloc; // ... and return the energy of the new configuration
} // end f1dim

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// frprmn
// Fletcher-Reeves-Polak-Ribiere minimization
// see http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
void frprmn(int n, LDBL * fret, LDBL(*func)())
{
    // most input variable are global, dfunc is done in gradienU
    int its;
    int endcount = 0;
    LDBL gam, fp, dgg;

    iloop(2 * N + 3) {
        phelper[i] = p[i];
    }
    (*func)(); // func = energy does not need inputvector
    fp = Uhelper;

    if(pressOnOff) {
        fp = fp + P0 * Lhelper * Lhelper;
    }

    energy();
    gradientcalc();

    iloop(2 * N + 3) {
        xi[i] = h[i] = g[i] = -xihelper[i]; // the 'force' is in the negative gradient direction
    }

    for(its = 1; its <= ITMAX; its++) {
        packIntoBoundaries();
        iterationcountfrprmnCUMULATIVE++;

        linmin(n, fret, func);

        if(2.0 * fabs(*fret - fp) <= ftol * (fabs(*fret) + fabs(fp) + ZEPS)) {
            endcount++;

            if(endcount > 3) {
                if(screenOutput) {
                    cout << "Conjugate gradient converged!" << endl;
                }

                frprmnconverged = true;
                return;
            }

        } else {
            endcount = 0;
        }

        fp = *fret;
        energy();
        gradientcalc();
        iloop(2 * N + 3) {
            xi[i] = xihelper[i];
        }
        dgg = gg = 0.0;

        iloop(2 * N + 3) {
            gg += g[i] * g[i];
            dgg += xi[i] * xi[i]; // Fletcher-Reeves
            //dgg += (xi[i] + g[i])*xi[i]; // Polak-Ribiere
        }

        if(gg == 0.0) {
            if(screenOutput) {
                cout << "Gradient is ZERO!" << endl;
            }

            frprmnconverged = true;
            return;
        }

        gam = dgg / gg;
        iloop(2 * N + 3) {
            g[i] = -xi[i];
            xi[i] = h[i] = g[i] + gam * h[i];
        }
        H = Uhelper + P0 * Lhelper * Lhelper;

    }

    return;
} // end frprmn

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// particle distance
void particledistance(int i, int j)
{

    xij[j * N + i] = phelper[j] - phelper[i];
    yij[j * N + i] = phelper[N + j] - phelper[N + i];

    ny[j * N + i] = -(int)floor((yij[j * N + i] + lyyhelper * 0.5) / lyyhelper);
    ny[i * N + j] = -ny[j * N + i];

    xij[j * N + i] = phelper[j] - phelper[i] + ny[j * N + i] * lyxhelper;
    yij[j * N + i] = phelper[N + j] - phelper[N + i]
                     + ny[j * N + i] * lyyhelper;

    nx[j * N + i] = -(int)floor((xij[j * N + i] + lxxhelper * 0.5) / lxxhelper);
    nx[i * N + j] = -nx[j * N + i];

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
void resethelpervars()
{
    // reset variables

    alphahelper = phelper[2 * N];
    deltahelper = phelper[2 * N + 1];
    Lhelper = phelper[2 * N + 2];
    lxxhelper = Lhelper / (1.0 + deltahelper);
    lxyhelper = Lhelper * 0.0;
    lyxhelper = Lhelper * alphahelper;
    lyyhelper = Lhelper * (1.0 + deltahelper);

    Uhelper = 0.0;

    return;
} // end resetvars

///////////////////////////////////////////////////////////////////////////////////////////
// energy
LDBL energy()
{
    resethelpervars();

    // calculate distances and overlaps
    iloop(N) {
        jloop(i) { // only calculate j<i, i.e., upper right corner of matrices

            if(neighbors[j * N + i]) {  // avoid selfinteraction

                particledistance(i, j);

                dij[j * N + i] = R[i] + R[j] - rij[j * N + i];

                if(dij[j * N + i] > 0.0) {  // this saves about 5-10 % in runtime
                    Uhelper += 0.5 * k * (dij[j * N + i] * dij[j * N + i]);
                    trueneighbors[j * N + i] = true;

                } else {
                    trueneighbors[j * N + i] = false;
                    dij[j * N + i] = 0.0;
                }
            } //end if
        } // end jloop
        trueneighbors[i * N + i] = false;
    } // end iloop

    return Uhelper;
} // end energy

//////////////////////////////////
// gradientcalc
void gradientcalc()
{
    LDBL component;
    LDBL term1;
    // calculate gradient and energy

    Phelper = 0.0;

    xihelper[2 * N] = 0.0;
    xihelper[2 * N + 1] = 0.0;
    xihelper[2 * N + 2] = 0.0;

    iloop(N) {

        xihelper[i] = 0.0;
        xihelper[N + i] = 0.0;

        jloop(i) {
            if(trueneighbors[j * N + i]) { // this saves about 5-10 % in runtime
                term1 = k * dij[j * N + i] / rij[j * N + i];
                // gradient expressions have been derived analytically
                component = term1 * xij[j * N + i]; //gradient component i with respect to particle j
                xihelper[i] += component;
                xihelper[j] -= component; // opposite sign in xji vs. xij

                component = term1 * yij[j * N + i];
                xihelper[N + i] += component;
                xihelper[N + j] -= component;

                if(alphaOnOff)
                    xihelper[2 * N] += -2 * term1 * Lhelper * xij[j * N + i]
                                       * ny[j * N + i];

                if(deltaOnOff)
                    xihelper[2 * N + 1] += -2 * term1 * Lhelper
                                           * (xij[j * N + i] * nx[j * N + i]
                                              * (-1.0 / ((1.0 + DELTA) * (1.0 + DELTA)))
                                              + yij[j * N + i] * ny[j * N + i]);

                Phelper += 2 * term1 * rij[j * N + i] * rij[j * N + i];
            } // end if
        } // end jloop
    }

    // end iloop

    // normalize
    Phelper = Phelper / (Lhelper * Lhelper) / 4.0;

    if(fabs(DELTA - deltainit) > 0.2) {
        deltaOnOff = false;
    }

    if(pressOnOff) {
        xihelper[2 * N + 2] = 2 * Lhelper * (P0 - Phelper);
    }

    return;
} // end gradientcalc

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// initializeSimulation

void initializeSimulation()
{

    LDBL randLx, randLy;
    LDBL frac;

    dtmax = 0.1;
    dt = dtmax;

    if(currentPackingNumber == firstPackingNumber) {
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
        pLast.reserve(2 * N + 3);
        xi.reserve(2 * N + 3);
        xihelper.reserve(2 * N + 3);
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
    ALPHA = alphainit;
    DELTA = deltainit;

    dampalpha = dampdelta = 0.9;
    damppress = 0.9;

    LENGTH = Linit;

    lxx = LENGTH / (1.0 + DELTA);
    lxy = LENGTH * 0.0;
    lyx = LENGTH * ALPHA;
    lyy = LENGTH * (1.0 + DELTA);

    srand(currentPackingNumber + 1);

    // set initial particle positions and properties
    iloop(N) {
        randLx = (rand() * 1.0) / (RAND_MAX * 1.0); // random position along L_x
        randLy = (rand() * 1.0) / (RAND_MAX * 1.0); // random position along L_y
        p[i] = randLx * lxx + randLy * lyx; // set random x-position for particle
        p[N + i] = randLy * lxy + randLy * lyy; // set random y-position for particle

        if(i < 0.5 * N) {
            R[i] = 1.0;    // set particle radii

        } else {
            R[i] = 1.4;
        }

        M[i] = 1.0; // set particle 'masses' for FIRE-algorithm

        // polydispersity :

//      R[i]=1+(1.4-1)*i/(N-1) ;
//      cout << "X[i] : "<< p[i] << "   Y[i]: "<< p[N+i] <<endl;
//      cout << "index : "<< i << " and radius : "<< R[i] <<endl;



    }
    M[N] = M[N + 1] = M[N + 2] = N;

    // rescale to fit desired fill fraction
    iloop(2 * N + 2)
    phelper[i] = p[i];

    energy();
    calcSysPara();

    // calculate the current fill fraction phi
    frac = sqrt(phi / (phiinit)); // sqrt-ratio of desired fill fraction and actual

    LENGTH *= frac; // scale the boxlength accordingly
    iloop(2 * N) {
        p[i] *= frac; // scale all the particle positions accordingly
        v[i] = 0.0;
    }

    v[2 * N] = v[2 * N + 1] = v[2 * N + 2] = 0;

    energy();
    calcSysPara();

    Rmax = 0.0;
    iloop(N) {
        if(R[i] > Rmax) {
            Rmax = R[i];    // determine the largest particle radius in the packing
        }

        // (for neighbor determination)
    }

    iloop(N) {
        jloop(N) {
            neighbors[j * N + i] = true; // initially all particles are concidered

            if(j == i) {
                neighbors[i * N + i] = false;    // neighbors, except if j==i
            }
        }
    }

    dtmax = dtmaxinit;
    dofOnOff = false;

    return;
} // end initializeSimulation()

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// calcSysPara: calculate fillfracton & number of neighbors & pressure
void calcSysPara()
{
    LDBL term1;

    int numberOfRattlerChanges = 1;

    LDBL Z2 = 0;

    // Rattlers check
    LDBL X_3n[3], Y_3n[3]; // in case there are 3 neighbors, we have to check their positions.
    int stable_particle;

    iloop(3) {
        X_3n[i] = 0;
        Y_3n[i] = 0;
    }
    int ip;

    lxx = LENGTH / (1.0 + DELTA);
    lxy = LENGTH * 0.0;
    lyx = LENGTH * ALPHA;
    lyy = LENGTH * (1.0 + DELTA);


    // reset sums
    phi = 0.0;
    Z = 0.0;
    P = 0.0;
    sxx = sxy = syy = 0.0;

    // sum all contributions
    iloop(N) {
        phi += PI * R[i] * R[i]; // sum of particle surfaces
        jloop(i) {
            if(trueneighbors[j * N + i]) {

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



    if(fireconverged) {



        // Ncorrected :
        // first loop defines the matrix "trueneighbors" and symmetrizes it.

        iloop(N) {
            isRattler[i] = false;
            numberOfDirectNeighbors[i] = 0;
            jloop(i)
            trueneighbors[i * N + j] = trueneighbors[j * N + i];
        }


        // now the algorithm of rattler detections :
        // First no particle is a rattler, then if this particle has less than 3 contacts, it is a rattler.
        // If at least one particle is a rattler within the loop, we have to start it again to check if this rattler did not make another rattler.
        // In the mean time, we count the number of contacts, excluding the rattlers.

        while(numberOfRattlerChanges != 0) {
            Z2 = 0;
            numberOfRattlerChanges = 0;
            iloop(N) {
//              cout << i << endl;
                numberOfDirectNeighbors[i] = 0;

                if(!isRattler[i]) {
                    jloop(N) {
                        if(trueneighbors[j * N + i]) {
                            if(!isRattler[j]) {
                                numberOfDirectNeighbors[i] += 1;
                                Z2 += 0.5;
                            }

                        }
                    }
                }

                if(numberOfDirectNeighbors[i] < 3 && !isRattler[i]) {
                    isRattler[i] = true;
                    numberOfRattlerChanges++;

                }

                if(numberOfDirectNeighbors[i] == 3 && !isRattler[i]) {

                    ip = 0;
                    jloop(N) {
                        if(trueneighbors[j * N + i] && !isRattler[j]) {
                            X_3n[ip] = p[j] + nx[j * N + i] * lxx  + ny[j * N + i] * lyx;
                            Y_3n[ip] = p[j + N] + nx[j * N + i] * lxy  + ny[j * N + i] * lyy;
                            ip += 1;

                        }
                    }


                    stable_particle = pnpoly(3, X_3n, Y_3n, p[i], p[i + N]);


                    if(!stable_particle) {
                        cout << "Particle " << i << " is not bounded by forces on all axes !!" << endl;
                        isRattler[i] = true;
                    }

                }




            }
        } // end while


        Ncorrected = 0;
        iloop(N) {
            if(!isRattler[i]) {
                Ncorrected++;
            }
        }
    } // end if(fireconverged)

    else {
        Ncorrected = N;
    }



    // normalize
    phi = phi / (LENGTH * LENGTH);


    Z2 = 2 * Z2 / (Ncorrected * 1.0 + 1e-16);
    Z = Z2;
    P = P / (LENGTH * LENGTH) / 2.0;
    sxx = sxx / (LENGTH * LENGTH);
    syy = syy / (LENGTH * LENGTH);
    sxy = sxy / (LENGTH * LENGTH);

    return;
} // end calcSysPara

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// readPositionFile
void readPositionFile()
{
    char c = 'x', clast = 'x', clastlast = 'x';

    int i = 0;
    bool initializeNow = false;

    int decimal = 0;
    LDBL helperchar = 0.0;
    bool isnegative = false;

    bool posRead = false;
    bool Nread = false;
    bool L1read = false;
    bool L2read = false;
    bool P0read = false;
    string filepath = nameOfWorkingDirectory + "/";


    LDBL L1x = 0.0, L1y = 0.0, L2x = 0.0, L2y = 0.0;
    // L1y is not used currently; suppress warning about this
    L1y = L1y;

    ifstream infile;

    if(distributioncase != 5) {
        if(screenOutput) {
            cout << "Creating file name." << endl;
        }

        createFileName();
        filepath.append(filenameString);

        if(screenOutput) {
            cout << "Filename: " << filepath << endl;
        }

        infile.open((char *) filepath.c_str());

    } else {
        infile.open((char *) filepath.c_str());
    }

    if(!infile.is_open()) {
        currentPackingNumber = 0;

        cout << "Input file did NOT OPEN!" << endl;
        exit(4);
    } else if(screenOutput) {
        cout << "Input file is opened succesfully!" << endl;
    }

    while(infile.good()) {
        clastlast = clast;
        clast = c;
        c = (char)infile.get();

        if(clastlast == 'N' && c == '=') {
            Nread = true;
            decimal = 0;
            helperchar = 0.0;

        }

        if(Nread) {
            if(c > 47 && c < 58) {
                helperchar = helperchar * 10 + (c - 48);
                decimal++;
            }

            if(c == ',') {
                N = (int)helperchar;
                helperchar = 0.0;
                Nread = false;
                initializeNow = true;
            }
        }

        if(initializeNow) {
            initializeArrays();
            initializeNow = false;
        }

        if(clastlast == 'L' && clast == '1' && c == '=') {
            L1read = true;
            i = 0;
            decimal = 0;
            helperchar = 0.0;

        }

        if(L1read) {
            if(c == '-') {
                isnegative = true;
            }

            if(c > 47 && c < 58) {
                helperchar = helperchar * 10.0 + (c - 48) * 1.0;
                decimal++;
            }

            if(c == '.') {
                decimal = 0;
            }

            if(c == ',' || c == '}') {

                for(int j = 0; j < decimal; j++) {
                    helperchar *= 0.1;
                }

                if(!isnegative && i == 0) {
                    L1x = helperchar;
                }

                if(isnegative && i == 0) {
                    L1x = -helperchar;
                }

                if(!isnegative && i == 1) {
                    L1y = helperchar;
                }

                if(isnegative && i == 1) {
                    L1y = -helperchar;
                }

                helperchar = 0.0;
                isnegative = false;
                i++;

            }

            if(i == 2) {
                L1read = false;
            }
        }

        if(clastlast == 'L' && clast == '2' && c == '=') {
            L2read = true;
            i = 0;
            decimal = 0;
            helperchar = 0.0;

        }

        if(L2read) {
            if(c == '-') {
                isnegative = true;
            }

            if(c > 47 && c < 58) {
                helperchar = helperchar * 10.0 + (c - 48) * 1.0;
                decimal++;
            }

            if(c == '.') {
                decimal = 0;
            }

            if(c == ',' || c == '}') {

                for(int j = 0; j < decimal; j++) {
                    helperchar *= 0.1;
                }

                if(!isnegative && i == 0) {
                    L2x = helperchar;
                }

                if(isnegative && i == 0) {
                    L2x = -helperchar;
                }

                if(!isnegative && i == 1) {
                    L2y = helperchar;
                }

                if(isnegative && i == 1) {
                    L2y = -helperchar;
                }

                helperchar = 0.0;
                isnegative = false;
                i++;
            }

            if(i == 2) {
                L2read = false;
            }
        }

        if(clastlast == 'P' && clast == '0' && c == '=') {
            P0read = true;
            decimal = 0;
            helperchar = 0.0;

        }

        if(P0read) {
            if(c > 47 && c < 58) {
                helperchar = helperchar * 10.0 + (c - 48) * 1.0;
                decimal++;
            }

            if(c == '.') {
                decimal = 0;
            }

            if(c == ',') {

                for(int j = 0; j < decimal; j++) {
                    helperchar *= 0.1;
                }

                P0 = helperchar;
                helperchar = 0.0;
                P0read = false;
            }
        }

        if(!L1read && !L2read && c == '{') {
            posRead = true;
            decimal = 0;
            helperchar = 0.0;

        }

        if(posRead) {
            if(c == '-') {
                isnegative = true;
            }

            if(clastlast == '{') {
                i = 0;
                decimal = 0;
            }

            if(c > 47 && c < 58) {
                helperchar = helperchar * 10.0 + (c - 48) * 1.0;
                decimal++;
            }

            if(c == '.') {
                decimal = 0;
            }

            if(c == ',' || c == '}') {

                for(int j = 0; j < decimal; j++) {
                    helperchar *= 0.1;
                }

                if(i < 3 * N && !isnegative && i % 3 == 0) {
                    p[i / 3] = helperchar;
                }

                if(i < 3 * N && isnegative && i % 3 == 0) {
                    p[i / 3] = -helperchar;
                }

                if(i < 3 * N && !isnegative && i % 3 == 1) {
                    p[N + i / 3] = helperchar;
                }

                if(i < 3 * N && isnegative && i % 3 == 1) {
                    p[N + i / 3] = -helperchar;
                }

                if(i < 3 * N && !isnegative && i % 3 == 2) {
                    R[i / 3] = helperchar;
                }

                if(i < 3 * N && isnegative && i % 3 == 2) {
                    R[i / 3] = -helperchar;
                }

                i++;

                isnegative = false;
                helperchar = 0.0;
            }

            if(c == '}') {
                posRead = false;
            }

        }
    } // end while

    infile.close();

    Rmax = 0.0;
    iloop(N) {
        if(R[i] > Rmax) {
            Rmax = R[i];    // determine the largest particle radius in the packing
        }

        // (for neighbor determination)
    }

    LENGTH = sqrt(L1x * L2y);
    ALPHA = L2x / LENGTH;
    DELTA = sqrt(L2y / L1x) - 1.0;


    iloop(N) {
        jloop(N) {
            neighbors[j * N + i] = true; // initially all particles are concidered

            if(j == i) {
                neighbors[i * N + i] = false;    // neighbors, except if j==i
            }
        }
    }

    if(screenOutput)
        cout << "test0000: readPositionFile-END: " << " N = " << N << ", L = "
             << LENGTH << ", P0= " << P0 << endl;

    if(distributioncase == 4) {
        P0 = P0init;
    }

    return;
} // end readPositionFile

////////////////////////////////////////////////////////////////////////
// writePositionFile
inline void writePositionFile()
{

    time_t rawtime;
    struct tm * timeinfo;
    char timebuffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(timebuffer, 80, "%Y-%m-%d_%H-%M-%S", timeinfo);

    string filepath = nameOfWorkingDirectory;

    filepath.append("/");

    createFileName();

    filepath.append(filename);

    ofstream outfile;

    if(programmode == PROGRAMMODE_CREATE_PACKING) {
        outfile.open((char *) filepath.c_str(), ios::trunc);
    }

    if(!outfile) {
        if(screenOutput) {
            cout << "Cannot open output file!" << endl;
        }

        outfile.close();
        ofstream outfile;
        createFileName();

        if(programmode == PROGRAMMODE_CREATE_PACKING) {
            outfile.open(filename, ios::trunc);
        }

        if(!outfile) {
            if(screenOutput) {
                cout << "Cannot open output file!" << endl;
            }

            return;

        } else if(screenOutput) {
            cout << "Output file is opened succesfully (2nd attempt)!" << endl;
        }
    }

    outfile.setf(ios::fixed, ios::floatfield);
    outfile.precision(PRECISION);

    outfile  << FILE_HEADER;

    lxx = LENGTH / (1.0 + DELTA);
    lxy = LENGTH * 0.0;
    lyx = LENGTH * ALPHA;
    lyy = LENGTH * (1.0 + DELTA);

    outfile << "N = " << N << " ,L = " << LENGTH << " ,L1= { " << lxx << " , " << lxy
            << " } " << " ,L2= { " << lyx << " , " << lyy << " } " << " ,P = "
            << Phelper << " ,P0= " << P0 << " ," << endl;

    outfile << "{" << endl;
    iloop(N) {
        outfile << p[i] << " ,	" << p[N + i] << " ,	" << R[i] << " ,	" << endl;
    }
    outfile << "}" << endl;
    outfile << endl << endl;

    outfile.close();
    outfile.close();

    string logFileNamePackings;

    logFileNamePackings = filenameString;

    logFileNamePackings.insert(0, "log");

    logFileNamePackings.insert(0, nameOfWorkingDirectory + "/");

    ofstream logfile;

    LDBL maxGrad = 0;
    iloop(2 * N) {
        if(fabs(xihelper[i]) > maxGrad) {
            maxGrad = fabs(xihelper[i]);
        }
    }

    logfile.open((char *) logFileNamePackings.c_str(), ios::app);
    logfile.setf(ios::scientific, ios::floatfield);
    logfile.precision(PRECISION);

    logfile << currentPackingNumber << "	" << N << "	" << P0 << "	" << P << "	"
            << ALPHA << "	" << DELTA;
    logfile << "	" << LENGTH << "	" << phi << "	" << Z << "	" << N - Ncorrected
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
inline void writeMultiplePackings(string name)
{
    time_t rawtime;
    struct tm * timeinfo;
    char timebuffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(timebuffer, 80, "%Y-%m-%d_%H-%M-%S", timeinfo);

    ofstream outfile;

    outfile.open((char *) name.c_str(), ios::app);

    if(!outfile) {
        if(screenOutput) {
            cout << "Cannot open output file!" << endl;
        }

        outfile.close();
        ofstream outfile;
        createFileName();

        if(programmode == PROGRAMMODE_CREATE_PACKING) {
            outfile.open((char *) name.c_str(), ios::app);
        }

        if(!outfile) {
            if(screenOutput) {
                cout << "Cannot open output file!" << endl;
            }

            return;

        } else if(screenOutput) {
            cout << "Output file is opened succesfully (2nd attempt)!" << endl;
        }
    }

    outfile.setf(ios::fixed, ios::floatfield);
    outfile.precision(PRECISION);

    lxx = LENGTH / (1.0 + DELTA);
    lxy = LENGTH * 0.0;
    lyx = LENGTH * ALPHA;
    lyy = LENGTH * (1.0 + DELTA);

    outfile << "N = " << N << " ,L = " << LENGTH << " ,L1= { " << lxx << " , " << lxy
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
void packIntoBoundaries()
{
    iloop(N) {
        // up/down boundaries
        while(p[N + i] < (0)) {
            p[i] += lyx;
            p[N + i] += lyy;
        }

        while(p[N + i] > (lyy)) {
            p[i] -= lyx;
            p[N + i] -= lyy;
        }

        // left/right boundaries
        // why is this a slanted boundary instead of a straight??
        while(p[i] < (0 + p[N + i] * lyx / lyy)) {
            p[i] += lxx;
            p[N+i] += lxy; // 0
        }

        while(p[i] > (lxx + p[N + i] * lyx / lyy)) {
            p[i] -= lxx;
            p[N+i] -= lxy; // 0
        }

    }

    energy();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
// createFileName
inline void createFileName()
{
    char namebuffer[4];
    filenameString = "";

    namebuffer[0] = (char)(48 + (((currentPackingNumber) / 1000) % 10));
    namebuffer[1] = (char)(48 + (((currentPackingNumber) / 100) % 10));
    namebuffer[2] = (char)(48 + (((currentPackingNumber) / 10) % 10));
    namebuffer[3] = (char)(48 + ((currentPackingNumber) % 10));

    if(screenOutput) {
        cout << "FilenameString: " << filenameString << endl;
    }

    unsigned int i = 0;

    string filebase = nameOfWorkingDirectory;

    while(nameOfWorkingDirectory[i] != 'N') {
        i++;
        filebase = nameOfWorkingDirectory.substr(i);

        if(i > nameOfWorkingDirectory.size()) {
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

    filename = (char *) filenameString.c_str();

    if(screenOutput) {
        cout << filename << endl;
    }

    return;
}
// end createFileName

void initializeArrays()
{

    // size all the dynamic 'vector' arrays according to the particle number N
    if(screenOutput) {
        cout << "test0002: N = " << N << endl;
    }

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
    pLast.reserve(2 * N + 3);
    xi.reserve(2 * N + 3);
    xihelper.reserve(2 * N + 3);
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

    iloop(2 * N) {
        v[i] = 0.0;
    }
    v[2 * N] = v[2 * N + 1] = v[2 * N + 2] = 0;

    ITMAXBRENT = max(100, 2 * N); // maximum of iterations in brent
    ITMAX = 5; // maximum of iterations in frprmn

    dtmax = dtmaxinit;
    dofOnOff = false;

    if(screenOutput) {
        cout << "test0003: initializeArrays-END" << endl;
    }

    return;
}

void checkFolderName(string foldername)
{

    string filename = "/runLog.txt";
    ;
    bool folderexists = false;

    ofstream testFolderName;
    filename.insert(0, foldername);
    testFolderName.open((char *) filename.c_str(), ios::app);

    if(testFolderName.is_open()) {
        folderexists = true;

    } else {
        cout << "ERROR" << endl;
    }

    if(folderexists) {
        testFolderName.close();
    }

    return;
}



//! Point Inclusion in Polygon Test
//! http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
//!
//! @param nvert         Number of vertices in polygon.
//! @param vertx, verty  Arrays containing nvert (x,y) coordinates of the
//!                      vertices of the polygons.
//! @param testx, testy  (x,y) coordinates of the test point
//!
//! @return true for strictly interior points
//!         false for strictly exterior points
//!         true or false for points on the vertices
bool pnpoly(int nvert, LDBL * vertx, LDBL * verty, LDBL testx, LDBL testy)
{
    int i, j;
    bool c = false;

    for(i = 0, j = nvert - 1; i < nvert; j = i++) {

        if(((verty[i] > testy) != (verty[j] > testy)) &&
           (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i])) {
            c = !c;
        }
    }

    return c;
}

extern "C" {
    struct packingparams {
        LDBL P;
        LDBL phi;
        LDBL Z;
        int Ncorrected;
        LDBL sxx;
        LDBL sxy;
        LDBL syy;
        LDBL U;
        LDBL H;
        LDBL gg;
    };

    void relax_packing(bool alphaFree, bool deltaFree, bool LFree) {
        screenOutput = true;

        alphaOnOffInit = alphaFree;
        deltaOnOffInit = deltaFree;
        pressOnOffInit = LFree;

        converged = frprmnconverged = fireconverged = false;

        calculate_neighbors(); // already calls energy() and gradientcalc()
        calcSysPara();


		programmode = PROGRAMMODE_DEFORM_PACKING;
        while(!converged) { // && !frprmnconverged) {
            cout << ".";
        
            simulationstep();

            /*
            cout << converged << ", " << frprmnconverged << ", " << fireconverged << ", " <<
            2.0 * fabs(H - HLastFunctionCall) << (ftolFIRE * (fabs(H) + fabs(HLastFunctionCall) + ZEPS)) << ", " <<
            ALPHA << ", " << fabs(Phelper - P0) / P0 << ", " << sxy << ", " << fabs(sxy) / 1e-15 << endl;
            */
        }
    }

    // import packing
    //
    // @param[in]  _N      Number of particles = len(x) = len(y) = len(r)
    // @param[in]  _P0     Base pressure (required to calculate entropy H)
    // @param[in]  _x[_N]  x coordinates of particles
    // @param[in]  _y[_N]  y coordinates of particles
    // @param[in]  _r[_N]  radii of particles
    // @param[in]  _alpha  alpha shear parameter (simple shear)
    // @param[in]  _delta  delta shear parameter (pure shear)
    // @param[in]  _L      length of box
    void import_packing(int _N, LDBL _P0,
                          LDBL *_x, LDBL *_y, LDBL *_r,
                          LDBL _alpha, LDBL _delta, LDBL _L) {
        fireconverged = false;
        // load packing
        N = _N;
        P = P0 = _P0;
        Rmax = 0;

        initializeArrays();
        initializeSimulation();

        iloop(N) {
            p[i] = _x[i];
            p[i+N] = _y[i];
            R[i] = _r[i];
            Rmax = max(Rmax, _r[i]);
        }

        Rneighbor = 3.5 * Rmax;

        ALPHA = _alpha;
        DELTA = _delta;
        LENGTH = _L;

        calculate_neighbors();
        fireconverged = true; // to trigger the correct part of calcSysPara()...
        calcSysPara();
    }

    // import packing
    //
    // @param[in]  _N      Number of particles = len(x) = len(y) = len(r)
    // @param[out]  _P0     Base pressure (required to calculate entropy H)
    // @param[out]  _x[_N]  x coordinates of particles
    // @param[out]  _y[_N]  y coordinates of particles
    // @param[out]  _r[_N]  radii of particles
    // @param[out]  _alpha  alpha shear parameter (simple shear)
    // @param[out]  _delta  delta shear parameter (pure shear)
    // @param[out]  _L      length of box

    void export_packing(int _N, LDBL * _P0,
                          LDBL *_x, LDBL *_y, LDBL *_r,
                          LDBL *_alpha, LDBL *_delta, LDBL *_L) {
        // write packing
        if(N != _N) {
            cout << "ERROR: N is not equal to the packing N. ABORTING.";
            return;
        }
        *_P0 = P0;

        iloop(N) {
            _x[i] = p[i];
            _y[i] = p[i+N];
            _r[i] = R[i];
        }

        *_alpha = ALPHA;
        *_delta = DELTA;
        *_L = LENGTH;
    }


    // @returns number of particles N
    int get_packing_size() {
        return N;
    }

    // calculate packing params for the packing that has been imported
    // @param[out] out     packingparams struct containing the calculated values
    void get_packing_data(packingparams *out) {
        // manually calculate H and gg (otherwise done in fire())
        out->H = Uhelper + P0 * Lhelper * Lhelper;

        LDBL gg_max=0;
        gg_max = 0;
        iloop(N) {
            gg_max = max(gg_max, xihelper[i]);
            gg_max = max(gg_max, xihelper[i+N]);
        }
        out->gg = gg_max;

        // fill in the other parameters

        out->P = P;
        out->phi = phi;
        out->Z = Z;
        out->Ncorrected = Ncorrected;
        out->sxx = sxx;
        out->sxy = sxy;
        out->syy = syy;
        out->U = Uhelper;
    }

}
