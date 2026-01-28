## Optimal Substation Reconfiguration ## 

See OptimalGridConfiguration.pdf for details.

**an AMPL implementations using a MILP to solve the problem.** 

- Grid_MaxFlow_V(.run,.dat,.mod) without randomization

and randomized version

- Grid_MaxFlow_Vrandom(.run,.dat,.mod) with randomization

**A Julia implemnentation (using Gurobi)**

- Data : {bus,branch,breaker)_data.csv

- Julia code : MC_Max_Flow_Breakers_v2.jl
  
Maximize inter-area power exchange by optimally configuring substation
breakers, while respecting DC power flow physics, thermal limits,
and basic connectivity constraints.

More advanced **Objective function**:
- Lambda: maximize inter-area exchange
- Connectivity: mild preference for closed breakers
- Penalty: discourages unnecessary breaker power circulation

We are checking for bridges in the optimal topology. Bridges are single points of failure that we want to avoid.





