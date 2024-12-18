# vasaMac

a small macroscopic traffic flow
simulator for mass-sports events, in the present version specialized to the Swedish Vasaloppet event with a mass start of approx 16.000 skiers. At the beginning, the number of tracks (at the Vasaloppet main race, the classical style is mandatory) decreases from 49 down to 6 and there is a steep uphill section. This will lead to a traffic jam with a delay time exceeding one hour for the last runners.  After calibrating the model to the 2024 race (sim1), a wave start can be simulated (sim2) or the initial uphill section of the course can be widenened to a minimum of 20 tracks (sim3).
The results are about to be published in Physica A. Everything is text
based, there is no GUI.

## Compiling

Assume that you are on a linux system and gnu g++ is
implemented, just enter 

```
make -f makefileLocal vasaMac
```
When using another compiler, change the makefile `makefileLocal`
accordingly.

## Running the Simulation

The simulations are organized into projects. Each project has some input and output files to be explained below. The three simulations of the published paper (projects `sim1` to `sim3`) can be found in the `sim/` subdirectory.

Just go to `sim/` and run a simulation by calling

```
vasaMac <project>
```
where the projects `sim1` to `sim3` are supplied as a starter

## The vasaMac project files


### Input files

* `.proj`  The top-level input file controlling and specifying the
  simulation project.  The order of the entries is **fixed** (I did not want to bother with xml
  parsing) and are explained in the comments (second row) of this file. We have the categories

  * grid specification (total simulation time, update time, length of the simulated section, cell size)
  * the four model parameters maximum density `rhomax` per track, wave velocity `w`, gliding friction coefficient `mu` and a dispersion variation coeficient `sigmarel`. The first two parameters pertain to the LWR model with a triangular FD. The third parameter, the free-flow speed V0, is calculated from the abilities of the athletes and the local uphill gradient using the assumption of constant power and the friction coefficient. `sigmarel` is for the free-flow dispersion phase of the proposed model
  * output specification. Every `dntout` 'th time step and every `dnxout` th cell is written to the output files
 
* `.trackParam` Specification of the course

* `.skiers` Specification of the athlete crowd and the starting configuration
  * `Group`: Starting groups (group 0=elite)
  * `size`: number of athletes in each group passing the first checkpoint "HoegstaPunkten" (I did not use the number in the starting lists because there were some no-shows and dropouts already to this point);
  * `paceLevel`: Median pace (inverse of the speed) of each group on level terrain (as estimated by the split speeds between two checkpoints with no significant hills in between)
  * `delayToStart` Median time interval for the athletes of the respective group between the starting shot and the passing time of the starting line (the simulation starts at the starting line, not at the starting blocks)
  * `delayWavestart` Only nonzero for the wavestart scenario: Time difference between the release time of the respective group and the release of the elite (starting shot). If `delayToStart` is included, the later waves cannot proceed to the starting line once the earlier ones are released (if set to zero, they can proceed)
  * common data file `HoegstaPunkten_splitTimes.csv` for all projects

### Output files

* `.macro`: The main macroscopic output

* `.HoegstaPunktenCalib`: Macroscopic output for direct comparison with the "single-athlete stationary detector data" (file `HoegstaPunkten_splitTimes.csv`) obtained from this checkpoint

* `.traj`: Synthetic (virtual) trajectories generated from the macroscopic output for comparison with the GNSS-trajectories ("floating-athlete data") of some athletes (not included)

### Miscellaneous files

* `.run` runs and plots the simulations using the gnuplot script `.gnu`


### Model specification

Two-phase model with dispersion-transport model for free traffic and a triangular LWR model for congested situations


### Initial and boundary conditions

Block mass start or wave start specified by the `.skiers` files; homogeneous Von-Neumann boundary conditions for the outflowing boundary (open system)


## Graphics

This simulator is purely text-based and has neither GUI nor
graphics. However, in the sample projects running under Linux, there
is a `.gnu` gnuplot command file for each project generating some plots. 


## References 

[1] K. Bogenberger, M. Treiber, and P. Malcolm, "Traffic Flow Phenomena in Large Sporting Events -- Empirical Analysis and Macroscopic Simulation of the Vasaloppet", Physica A, submitted.

[2] M. Treiber,
"Crowd Flow Modeling of Athletes in Mass Sports Events - a Macroscopic Approach" in Traffic and Granular Flow '13 (Springer), pp. 21-29
DOI: 10.1007/978-3-319-10629-8_3 (2015).

[3]  M. Treiber and A. Kesting. [_Traffic Flow Dynamics, Data, Models and Simulation_](http://www.traffic-flow-dynamics.org). Springer 2013. [Link](http://www.springer.com/physics/complexity/book/978-3-642-32459-8).

