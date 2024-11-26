# vasaMac

a small macroscopic traffic flow
simulator for mass-sports events, in the present version specialized to the Swedish Vasaloppet event with a mass start of approx 16.000 skiers. At the beginning, the number of tracks (at the Vasaloppet main race, the classical style is mandatory) decreases from 49 down to 6 and there is a steep uphill section. This will lead to a traffic jam with a delay time exceeding one hour for the last runners.  After calibrating the model to the 2024 race (sim1), a wave start can be simulated (sim2) or the initial uphill section of the course can be widenened to a minimum of 20 tracks (sim3).
The results are about to be published in Physica A. Everything is text
based, there is no GUI.

## Compiling

Assume that you are on a linux system and gnu g++ is
implemented, go to ```trunk/src/``` 
and just enter 

```
make -f makefile_withLibs mic
```
When using another compiler, change the makefile
accordingly.

## Running the Simulation

The simulations are organized into `projects`. Each project has some input and output files to be explained below. Three example simulations (projects `sim1` to `sim3`) can be found in the sim/ subdirectory.

Just run a simulation by calling

```
vasaMac <project>
```
(without any extension)

## The vasaMac project files


### Input files

* `.proj`  The top-level input file controlling and specifying the
  simulation project.  The order of the entries is **fixed** (I did not want to bother with xml
  parsing). Its entries are the following 

  * `tmax`: Simulation time

 
* `.trackParam` Specification of the coursev

### Model specification

Two-phase model with dispersion-transport model for free traffic and a triangular LWR model for congested situations


### Initial and boundary conditions

This section specifies the traffic demand and its spatiotemporal
change. 


## Graphics

This simulator is purely text-based and has neither GUI nor
graphics. However, in the sample projects running under Linux, there
is also a gnuplot command file generating some plots. 
## References 

[1]  M. Treiber and A. Kesting. [_Traffic Flow Dynamics, Data, Models and Simulation_](http://www.traffic-flow-dynamics.org). Springer 2013. [Link](http://www.springer.com/physics/complexity/book/978-3-642-32459-8)

[2] Physica A submission
