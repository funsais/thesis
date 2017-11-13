# Monte Carlo Simlation of H2 onto Graphene Sheet to develop an Isotherm.

### At a Glance Status:
This is my Master's thesis project, and is currently on hiatus due to my studies also being on hiatus for financial reasons. Currently the project is a stand-alone script that simulates the unit cell, generates, and stores the data as a csv, and visualizes the unit cell with ggplot.  
The goal is to develop an isotherm, compare the isotherm with experimental data to evaluate the accuracy of the simulation, and determine the effect of early stopping of the simulation on accuracy.

### Top Level Reasons for Research:
If you need separate or move a particular type of gas, often you'll use some sort of lattice or other material that gas will want to "stick" to. There are many materials like this, and there are also many more in development. It would be nice to know approximately how well these materials adsorb gases before spending all the time and money towards development only to find out it didnt work as well as intended. Simulations allow for approximations of effectiveness without all of that development and production, and the closer they are to reality the better. The main problem is computation time, which can be fairly expensive in other ways. YOu could theoretically run your simulation forever depending on the amount of "accuracy" you want. Our goal is to experiment on the effect that computation time has on accuracy and provide some theoretical stopping points other than "enough time has passed and seems equilibrated."

## The Simulation Process:

### 1: Define your System

#### Unit Cell
The core of a simulation like this is the notion of a "Unit Cell" that is a representation of the rest of your surface. Lattices repeat, and so building a unit cell that when stacked and repeated, represents the whole lattice framework without computing the whole thing. This assumes that any particle leaving the cell, re-enters the cell from the opposite direction if it is repeating in that direction. For our system, our unit cell repeats in the X and Y directions, so any H2 particle that leaves one side of our unit cell on the X direction re-enters the same cell as if the particle were coming from the cell next to it.
Our unit cell is a roughly 40Ax40A graphene sheet at Z = 0, the space spanning that sheet in the X and Y direction, and extending 70A above the graphene sheet in the Z direction.

#### Properties
In order to develop an isotherm, you need the amount of particles adsorbed at a particular Temperature and Pressure. You also need to know how the particles are interacting at a local level. For this project we are using a Lennard-Jones model.

### 2: The Trial Move:
At the start of the simulation proper, one of three trial moves could happen: Move, Add, or Remove. Once one of these is chosen, the simulation considers what the system would look like after the step takes place, calculates the change in potential energy do to this change, and determines a probability of this move happening would be. In theory, an equilibrium would occur when the graphene sheet is saturated with H2 and will add and remove particles at an equal rate.

#### Move:
A random H2 particle is chosen from the unit cell, and is moved some distance from its location to a new one (the average amount of movement is referred to as delx).
#### Add:
At a random location in the unit cell, an H2 particle is placed.
#### Remove:
A random H2 particle is chosen to be removed from the system.

The change in potential energy is determined by the particle's Lennard-Jones potential with each H2 and C particle in the unit cell. So not only should H2 particles congregate somewhere around the minimum potential energy between H2 and C (~3.61 A), but they should also disperse among themselves somewhere around the minimum potential energy between H2 and H2. 

These trial moves continue until an "equilibrium" is reached. Typically this is determined by assuming enough trial moves have elapsed in general, though there are many other potential equilibrium definitions that can be reached.

#### Probability:
In order to determine whether or not a trial move has taken place, a probability is assigned based on the change in potential energy caused by the trial move. For our purposes we are using probabilities derived from the Grand Canonical Thermodynamic model assuming an Ideal gas state.

Q: Why perform this simulation numerically? Aren’t there plenty of resources for pre-developed simulation environments to use?
A: Doing this simulation numerically gives me more power to target the variables I want and to bug-fix appropriately. We are also testing the numerical method after all, our goal is not production. The primary reason this simulation is performed in R is because that was the language I knew best moving into the drafting and initial work of this thesis. It also works fairly well with stochastic models and simulations. The plotting functions in R are in my opinion the easiest to use and interact with.

Q: What is your unit cell?
A: Currently the unit cell is a graphene sheet at Z = 0 that is roughly 40 A by 40 A.
Q: How long does it take to get a dataset? Isn't this method very computationally intensive?
A: It depends on the device obviously. With a unit cell of what we're working with now, it takes quite a few hours to run the system to equilibrium. Part of the work here is that some of that time can be cut off by knowing when to cut off the simulation. This is also why the project work is currently on hold. My device is not equipped to run scripts like this for longer than 2 or 3 hours, and since I'm not attending Rose at the moment, I don't have a password to the remote desktop connection to the cluster.
