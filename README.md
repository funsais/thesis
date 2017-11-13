# Monte Carlo Simlation of H2 onto Graphene Sheet to develop an Isotherm.

### At a Glance:
This is my Master's thesis project, and is currently on hiatus due to my studies also being on hiatus for financial reasons. Currently the project is a stand-alone script that simulates the unit cell, generates, and stores the data as a csv, and visualizes the unit cell with ggplot.  
The goal is to develop an isotherm, compare the isotherm with experimental data to evaluate the accuracy of the simulation, and determine the effect of early stopping of the simulation on accuracy.

Q: Why perform this simulation numerically? Aren’t there plenty of resources for pre-developed simulation environments to use?
A: Doing this simulation numerically gives me more power to target the variables I want and to bug-fix appropriately. We are also testing the numerical method after all, our goal is not production. The primary reason this simulation is performed in R is because that was the language I knew best moving into the drafting and initial work of this thesis. It also works fairly well with stochastic models and simulations. The plotting functions in R are in my opinion the easiest to use and interact with.

Q: What is your unit cell?
A: Currently the unit cell is a graphene sheet at Z = 0 that is roughly 40 A by 40 A.

Q: How long does it take to get a dataset? Isn't this method very computationally intensive?
A: It depends on the device obviously. With a unit cell of what we're working with now, it takes quite a few hours to run the system to equilibrium. Part of the work here is that some of that time can be cut off by knowing when to cut off the simulation. This is also why the project work is currently on hold. My device is not equipped to run scripts like this for longer than 2 or 3 hours, and since I'm not attending Rose at the moment, I don't have a password to the remote desktop connection to the cluster.
