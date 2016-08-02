jamBashBulk
=============

The Van Hecke Lab shear-stabilized packing simulation code.

JamBashBulk is a C/C++ simulation code that provides three basic operations:
 * *Creating* shear stabilized packings,
 * *Shearing* packings with non-square boundary conditions, and
 * Arbitrary deformations via a C/Python interface.
 
These three operations will be described below.

For background and the physics behind the simulation, please see:
 * Simon Dagois-Bohy, Brian P. Tighe, Johannes Simon, Silke Henkes, and Martin van Hecke.
   _Soft-Sphere Packings at Finite Pressure but Unstable to Shear_.
   [Phys. Rev. Lett. **109**, 095703](http://dx.doi.org/10.1103/PhysRevLett.109.095703),
   [arXiv:1203.3364](http://arxiv.org/abs/1203.3364).

* Merlijn S. van Deen, Johannes Simon, Zorana Zeravcic, Simon Dagois-Bohy, Brian P. Tighe, and Martin van Hecke.
  _Contact changes near jamming_.
  [Phys. Rev. E **90** 020202(R)](http://dx.doi.org/10.1103/PhysRevE.90.020202), [arXiv:1404.3156](http://arxiv.org/abs/1404.3156).

* Merlijn S. van Deen, Brian P. Tighe, and Martin van Hecke.
  _Contact Changes of Sheared Systems: Scaling, Correlations, and Mechanisms_.
  [arXiv:1606.04799](https://arxiv.org/abs/1606.04799)

* Merlijn S. van Deen.
  _Mechanical Response of Foams: Elasticity, Plasticity, and Rearrangements_.
  PhD Thesis, Leiden University, 2016. [hdl:1887/40902](https://openaccess.leidenuniv.nl/handle/1887/40902)

Several data sets created using the simulation code are available via Zenodo:
 * http://dx.doi.org/10.5281/zenodo.59216 (stable packings)
 * http://dx.doi.org/10.5281/zenodo.59217 (sheared packings)

Compilation and installation
----------------------------
```bash
git clone https://github.com/vanheckelab/jambashbulk
cd jambashbulk
make
```

This will build:
 * `bin/jam2D`, the main simulation code,
 * `bin/jam2D_d`, the same, but with debugging information,
 * `bin/jamBashbulk.so`, the shared library.


Usage from the command line
---------------------------
Interactive usage is supported by passing the `-screen` flag:
```bash
bin/jam2D -screen
```

This will show a number of prompts. The following sections show the typical in/output for creating and shearing packings.

### Creating packings
First, create the directory where the packings will be stored. This has to be in the format `Packings/N16~P1e-3`,
where `Packings` is a required prefix, `N16` determines the number of particles (16), and `P1e-3` determines
the external pressure (10<sup>-3</sup>). For non-integer powers, the format `N16~P1234e-3` can be used, which specifies
an external pressure of 1.234·10<sup>-3</sup>.

To create packings, use program mode _1_ ('create DEFAULT jammed packing(s)') with
system state _2_ ('use RANDOM particle distribution'). For example, pass in the following commands:

```
$ mkdir -p Packings/N16~P1e-3   # directory needs to be created in advance
$ bin/jam2D -screen
(...)
Packings/N16~P1e-3
(... program mode ...)
1
(... number of packings ...)
1
(... starting number ...)
1
(... system state ...)
2
```

This will create a shear stabilized packing in `Packings/N16~P1e-3/N16~P1e-3~0001.txt`. The file format is described below.

### Shearing packings
To _shear_ packings (searching for contact changes), use program mode _3_ ('apply simple shear (rearrangements)') with system state _3_ ('read particle distribution from file').
For example, the following commands:
```
$ bin/jam2D -screen
Packings/N16~P1e-3
3
1
1
3
3
```

finds the first three contact changes for the packing created above.

### argvtolines.sh
To simplify automation, `argvtolines.sh` is provided. Simply pass each line of input as parameter. Example usage:
```
./argvtolines.sh Packings/N16~P1e-3 1 1 1 2 | bin/jam2D -screen
```
to create a packing as in the first example above.


### C/DLL interface
The simulation code can be used from other code by using the C/DLL interface. (jamBashBulk.cpp, in the
extern "C" block). Provided are:

 * `void import_packing(...)`, to load a packing in the simulation,
 * `void relax_packing(...)`, to relax the packing with given boundary conditions,
 * `void export_packing(...)`, to read out the packing
 
and the convenience methods
 * `int get_packing_size` to get the number of particles `N`, and
 * `void get_packing_data(packingparams *out)` to read out information on the current state.
are provided

#### import_packing
```C
void import_packing(int _N, LDBL _P0,
                    LDBL *_x, LDBL *_y, LDBL *_r,
                    LDBL _alpha, LDBL _delta, LDBL _L)
```
`import_packing` loads a packing (particle positions and boundary conditions) into the simulation. Parameters are:
* system size `_N`,
* external pressure `_P0`,
* particle positions and radii `_x`, `_y` and `_r` (`LDBL` is `double` on win32/64, `long double` on linux)
* simple shear and pure shear parameters `_alpha` and `_delta`,
* square root of area `_L`.

#### relax_packing
```C
void relax_packing(bool alphaFree, bool deltaFree, bool LFree)
```
Relaxes the currently loaded packing until converged (|∆H| < 10<sup>−13</sup> · H)). Parameters are:
* alphaFree (is the simple shear boundary condition free to change?),
* deltaFree (same, for pure shear), and
* LFree (same, for the area)

To apply a given deformation, call `import_packing` with the _new_ boundary conditions, then call `relax_packing`. For example:
```C
import_packing(_N, _P0, _x, _y, _r, _alpha + 0.01, _delta, _L);
relax_packing(false, false, false);
export_packing(_N, &_P0, _x, _y, _r, &_alpha, &_delta, &_L);
```
Will apply a simple shear γ=0.01, and read out the new system state.

#### export_packing
```C
void export_packing(int _N, LDBL * _P0,
                    LDBL *_x, LDBL *_y, LDBL *_r,
                    LDBL *_alpha, LDBL *_delta, LDBL *_L)
```
Exports the current packing into _user-provided_ arrays. This function _writes_ to all parameter except for `_N`; the parameters are the same as for `import_packing`. To get the correct value of `_N` (and to pre-allocate arrays), call `int get_packing_size()`.

#### get_packing_size
```C
int get_packing_size()
```
Returns the current number of particles `N`.

#### get_packing_data
```C
struct packingparams {
    LDBL P;          // internal pressure
    LDBL phi;        // packing fraction
    LDBL Z;          // contact number
    int Ncorrected;  // effective number of particles (=#particles-#rattlers)
    LDBL sxx;        // stress on xx-boundary
    LDBL sxy;        // simple shear stress
    LDBL syy;        // stress on yy-boundary
    LDBL U;          // energy
    LDBL H;          // enthalpy
    LDBL gg;         // maximum gradient
};

void get_packing_data(packingparams *out)
```
Fills a _user-provded_ `packingparams` struct with the current simulation state. 

### Python interface
A Python interface which wraps the C interface is provided. This interface takes care of allocating numpy arrays and calls the relevant C interface functions.

The Python interface provides:
  * direct `ctypes` based wrapping functions (jamBashBulk.DLL.*)
  * Python wrapping functions: `def import_packing`, `def export_packing`, `def relax_packing`,
                               `def get_packing_size`, `def get_packing_data`

#### import_packing
`import_packing` takes a single parameter; a dictionary containing:
```
{'L': 26.21776889062624,                       # square root of area,
 'L1': array([ 25.77606417,   0.        ]),    # Lxx, Lxy
 'L2': array([ -0.12999204,  26.66704277]),    # Lyx, Lyy
 'N': 128,                                     # number of particles
 'P': 0.01000000000000001,                     # internal pressure (not required for import)
 'P0': 0.01000000000000001,                    # external pressure
 'particles': array([(9.349917925369324, 0.0, 12.224353792603848, 0.0, 1.000000000000001),
       ...], 
      dtype=[('x', '<f8'), ('x_err', '<f8'), ('y', '<f8'), ('y_err', '<f8'), ('r', '<f8')])}
```
`particles` can also be a dictionary `{'x': [x0, x1, ...], 'y': [y0, y1, ...], 'r': [r0, r1, ...]}`.

#### relax_packing
See `relax_packing` of the C API.

#### export_packing
Exports the current packing in the same format as `import_packing`.

#### get_packing_size
Returns the number of particles `N`.

#### get_packing_data
See `get_packing_data`; converts the struct into a dict with keys P, phi, etc.



File format (creating a packing)
--------------------------------
_Creating_ a packing creates two files (see `examples/`):
 * `N16~P1e-3~0001.txt`, the relaxed, shear-stabilized packing,
 * `logN16~P1e-3~0001.txt`, with a number of packing properties.
 
The packing format is as follows:
```
# src/jamBashbulk.cpp (unix)
# commit 85670c21403a4f54092b9da783e3a906d5e2f869 (2016-08-01 17:09:16 +0200)
# RAND_MAX = 2147483647
# srand(0); rand() = 1804289383
N = 16 ,L = 9.5543369225244207275 ,L1= { 9.5393607291793076335 , 0.0000000000000000000 }  ,L2= { 0.0630575654921525210 , 9.5693366275464142399 }  ,P = 0.0010000000000000000 ,P0= 0.0010000000000000000 ,
{
8.8767069662298361795 , 7.2547583839685033783 , 1.0000000000000000000 ,
0.0594811234618037227 , 2.1385979775910763725 , 1.0000000000000000000 ,
2.4694173606288082063 , 4.7755001031424884048 , 1.0000000000000000000 ,
6.0292182826694455775 , 0.2223613971052526768 , 1.0000000000000000000 ,
5.5615954890760337039 , 5.2120391800267744696 , 1.0000000000000000000 ,
0.6032967432765678209 , 4.0617315606256442218 , 1.0000000000000000000 ,
2.0551790531853149613 , 2.2633237872532981267 , 1.0000000000000000000 ,
0.9221716315664489174 , 6.0388375891240512958 , 1.0000000000000000000 ,
3.6441763659566363022 , 0.4648779694119330627 , 1.3999999999999999112 ,
6.9595452979294494607 , 2.4322762560624367426 , 1.3999999999999999112 ,
4.2678726517533163713 , 3.1919539651122893701 , 1.3999999999999999112 ,
8.4215430695368404561 , 0.0490738802499536716 , 1.3999999999999999112 ,
3.7601085052506494589 , 6.7966197174681128482 , 1.3999999999999999112 ,
7.9547166101854404198 , 5.0440810567318968504 , 1.3999999999999999112 ,
1.4535012083627461571 , 8.3791982947573971164 , 1.3999999999999999112 ,
6.4855578360586323317 , 7.4264381068616523774 , 1.3999999999999999112 ,
}
```
The first lines are comments noting information on the platform, the version of the simulation code and the `rand()`
used in the standard C library.
The first line of data contains system properties: the number of particles `N`, the square root of the area, `L`,
the periodic boundary condition vectors `L1={Lxx,Lxy}` and `L2={Lyx,Lyy}`, the _internal_ pressure `P`
and the _external_ pressure `P0`. The following lines describe each particle, one per line. Each particle has an `x` and `y`
position and a radius `r`.


The `logN16~P1e-3~0001.txt` data format contains the following columns (`void writePositionFile()`):
* current packing number,
* number of particles N,
* external pressure P0,
* internal pressure P,
* simple shear α,
* pure shear δ,
* square root of area, L,
* packing fraction φ,
* contact number Z,
* number of rattlers #rattler,
* boundary stress σ_{xx},
* boundary stress σ_{yy},
* simple shear stress σ_{xy},
* energy U,
* energy change in last minimization step dU,
* enthalpy H,
* enthalpy change in last minimization step dH,
* simulation time,
* FIRE iteration count,
* conjugate gradient iteration count,
* maximum gradient,
* simulation time,
* date/time of simulation

File format (shearing a packing)
--------------------------------
_Shearing_ a packing creates three files:
 * `particlesN16~P1e-3~SR003~step007~0001.txt`, containing packing for each bisection step,
 * `logN16~P1e-3~SR003~step007~0001.txt`, with packing properties for each bisection step,
 * `dataN16~P1e-3~SR003~step007~0001.txt`, with contact change information for each step.
 
The `particles` file contains packing states for each bisection step. The format is the same as for `N16~P1e-3~0001.txt`; each
individual state is separated by two newlines.

The `log` file contains the following columns:
* step number of the shear simulation step#
* N, P0, P, alpha, delta, L, phi, Z, #rattler, s_xx, s,_yy, s_xy, U, dU, H, dH, t_run, #FIRE, #CG, maximum gradient, creating date as before
* contact change number (0 while converging to the first cc, 1 while converging to the second, etc)

The `data` file contains the following columns:
* applied shear strain (alpha-alpha(0)),
* shear stress s_xy,
* number of contacts Ncontacts
* number of changed contacts when compared to base state (initial state while converging to the first cc,
  just after the first cc when converging to the second cc, etc.)
* number of _created_ contacts when compared to base state,
* number of _broken_ contacts when compared to base state,
* internal pressure P,
* contact number Z.
