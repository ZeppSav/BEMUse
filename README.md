# BEMUse - A lightweight boundary element method library for hydrodynamics and aerodynamics.

The purpose of the BEMUse library is threefold:
* Provide users with an intuitive, lightweight boundary element method (BEM) library;
* As an alternative to proprietary software;
* As a teaching platform for numerical methods for aerodynamics and hydrodynamics.

BEMUse is written in C++ and can be compiled in two configurations depending on your application case:
* **Command-line configuration:** For code development or for the inclusion of BEMUse into a design loop, a command-line interface is available;
* **GUI configuration:** For visual inspection of meshes and visualisation of analysis results, a minimal graphical user interface (GUI) is available.

The core architecture of BEMUse is identical in both cases. 

## Why should you use BEMUse?
* You want a solver to carry out frequency-domain potential flow hydrodynamic calculations (1st order- radiation & diffraction);
* You want an open-source alternative to proprietary hydrodynamic solvers;
* Execution is fast: BEMUse leverages OpenMP in order to use multi-threading to speed up calculations;
* BEMUse can import, `WAMIT`®, `NEMOH` and `.stl` geometry formats; 
* BEMUse can export meshes to `WAMIT`® and `NEMOH` geometry formats; 
* You want a simple GUI to visualise and perform sanity checks on your mesh;
* You want a simple GUI to inspect the radiation/excitation potential solution on the boundary of your mesh;
* You want a simple GUI to visualise free surface displacement due to radiation or diffraction potentials;
* You want to generate your own geometries. It is easy to create your own meshes based on the template `Boundary` class.
* You want a solver to carry out potential flow aerodynamic calculations;
* Minimal dependencies: The only external dependency is the [eigen](https://eigen.tuxfamily.org) linear algebra library;
* You want a lightweight boundary element method library to form the basis of your own software;

<img src="/img/OC4.PNG" width=75% height=75%>

## Licensing and authorship
BEMUse is developed by Joseph Saverin and is distributed under the GNU General Public License Version 2.0, copyright © Joseph Saverin 2022.

## How to use the BEMUse command-line configuration
In order to carry out an analysis with BEMUse, you need to specify the geometry which you wish to analyse and the simulation parameters. 
The geometry is specified over the command line. There are two options for specifying a geometry. 

### Importing a Geometry
You can import an external geometry simply by specifying a `WAMIT`® (`.gdf`), `NEMOH` (`.mar`) or (`.stl`) geometry file. 
The specified path should be relative to the directory in which the executable is found. An example might be:
```
BEMUse <ExternalGeos/MyGeo.mar> <OutputPrefix>
```
BEMUse will automatically recognise the geometry type based on the file suffix. As the surface elements are externally defined, 
you do not need to specify the input files `Discretisation.bemin` or `Dimensions.bemin` (these are described below).

### Using a Template Geometry
A range of template geometries are available which are already parametrised in BEMuse, 
so that you only need to specify the discretisation or dimensions.
These geometries include analytical geometries as well as standard geometries for floating offshore wind turbine platforms. 
This list is regularly being extended and if you would like to have your geometry on this list, please get in touch. 
In order to see a list of the available geometries, type the following into the command line:
```
BEMUse Templates
```
In order to check what the necessary inputs (described below) are for your chosen template geometry, `<GeoType>`, you need simply to type:
```
BEMUse Inputs <GeoType>
```
You are required to specify a directory named `Outputs` within the directory of the executable where the output files will be written. 
In order to discern between multiple runs, you are also required to specify the output file prefix.
Provided you have provided the input parameters (described below), a complete call to execute BEMUse and write the outputs will look like this:
```
BEMUse <GeoType> <OutputPrefix>  
```
That simple! If there are insufficient inputs defined for the chosen geometry, the behaviour will be undefined.

### Specifying input files
A range of simulation parameters also need to be defined for a BEMUse simulation. These are defined in text-based input files. 
You are required to specify a directory named `Inputs` within the directory of the executable from where the input files will be read.  
These files are defined below along with the corresponding variables:
* `Frequencies.bemin`: Specify here the desired frequencies of analysis
	* One frequency per line, ideally in increasing order
* `Discretisation.bemin`: Specify here the discretisation for the template geometries. These have three formats:
	* [BDRY] Specify the discretisation of elements on the surface of the geometry being investigated;
	* [IFS] Specify the discretisation of elements on the interior free-surface of the geometry being investigated. This is relevant if the `IFR` feature has been activated;
	* [EFS] Specify the discretisation of elements on the exterior free-surface of the geometry being investigated. This is only revelant if you wish to inspect the free surface displacement around the geometry.
* `Dimensions.bemin`- Specify here the dimensions for the template geometries. These have three formats:
	* [BDRY] Specify the dimensions of elements on the surface of the geometry being investigated;
	* [EFS] Specify the dimensions of elements on the exterior free-surface of the geometry being investigated. This is only revelant if you wish to inspect the free surface displacement around the geometry.
* `Inertia.bemin`- Specify here the inertial parameters of the geometry. These are important for the calculation of the RAOs of the geometry. These include:
	* Centre of gravity of the geometry;
	* Mass matrix option: 1) Calculates the mass matrix (M_ii) assuming a uniform density equal to water density. 2) Imported;
	* Radii of Gyration: If option 1) above has been selected, these terms are used to calculate M_ii(4,4)-M_ii(6,6);
	* Imported Mass Matrix: If option 2) above has been selected, the elements of M_ii;
	* Hydrostatic Stiffness Matrix Option: 1) Calculates the hydrostatic stiffness matrix (H_ii) assuming a uniform density equal to water density. 2) Imported;
	* Imported Hydrostatic Stiffness Matrix: If option 2) above has been selected, the elements of H_ii.
* `SolverParams.bemin`- Specify solver flags
	* Irregular Frequency Removal- If this is specified as true, you need to specify the [IFS] values in the `Discretisation.bemin` file 
	* Kochin Angles [-]
	* Density of the fluid [kg/m^3]
	* Acceleration due to gravity [m/s^2]

## How to use the BEMUse GUI configuration
In this case you really don't need to do very much at all, just follow the steps in the GUI. 
Once you have specified a geometry, you can toggle any of the visualisation options at the top right of the window.
BEMUse solves a radiation potential for each of the six degrees of freedom (surge, sway, heave, roll, pitch, yaw) at each frequency. 
BEMUse solves a diffraction potential for each of the six degree of freedom for each incoming wave angle selected in the specification of the analysis.
The radiation potential and diffraction potential can only be visualised once an analysis has completed.
The free surface visualisation is only possible if you have specified an exterior free surface grid in the creation of the geometry.
### Modifying the view angle or visualisation
A number of options exist for modifying the visualisation:
* To rotate the view hold the left mouse button and shift the mouse;
* To translate the view hold the right mouse button and shift the mouse;
* To zoom in and out use the mouse scroll wheel;
* To snap to an isometric viewing angle double click the left mouse button;
* To toggle the degree of freedom which is visualised, hold the `Shift` button and use the scroll on your mouse;
* To toggle the frequency which is being visualised, hold the `Ctrl` button and use the scroll on your mouse;
* To toggle the degree of freedom which is visualised, hold the `Shift` button and use the scroll on your mouse; 
* To toggle the incoming wave angle which is being visualised, hold the `Alt` button and use the scroll on your mouse. 

<img src="/img/Wigley.PNG" width=75% height=75%>

## Citation information		
I am currently in the process of preparing a simple reference paper for BEMUse on ArXiv. 
In the meantime, if you use BEMUse please cite it in your publication as follows:
- Saverin, J., Grueter, L. **Quadrature Schemes for the Integration of the Free-Surface Greens Function within the Open-Source Boundary Element Method Library BEMUse**, 
Proceedings of the ASME 2022 41st International Conference on Ocean, Offshore and Arctic Engineering. June 2022, Hamburg, Germany.

## Compilation
The only dependency of BEMUse is the [eigen](https://eigen.tuxfamily.org) lineary algebra library. 
The compilation of BEMUse has been tested with GCC (v7.3).

Two options are available for compiling:
- qmake: BEMUse was prepared with the cross-platform development environment [Qt Creator](https://www.qt.io/product/development-tools). 
The .pro file required for compiling with qmake has been provided. 
You simply need to include the path to the eigen directory on your device.
For example:
```
INCLUDEPATH += -L$$PWD/../eigen/
```
As the GUI config relies on the Qt widget libraries, this configuration is probably best compiled with [Qt Creator]. 
In order to avoid compiling the GUI config the entire section in the `.pro` file below the line `# GUI configuration` should be commented out.

- CMake: The `CMakeLists.txt` file has been provided. 
You simply need to point the `EIGEN_DIR` variable to the eigen directory on your device.
For example:
```
set(EIGEN_DIR C:\\MyCodeisHere\\eigen\\)
```
Nb. The CMakeLists file is currently only configured to compile the command-line version of BEMUse.

## Pre-compiled libraries
Pre-built binaries are available on the [releases](https://github.com/ZeppSav/BEMUse/releases) page.

## Documentation

A [readthedocs] page is currently being prepared. 
