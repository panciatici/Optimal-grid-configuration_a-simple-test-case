################################################################################
# Optimal Grid Configuration with Breakers
# Monte Carlo Simulation – Full Script
#
# This code implements the simplified DC power flow test case described in:
#   "Optimal grid configuration – a simple test case"
#
# Goal:
#   Maximize inter-area power exchange by optimally configuring substation
#   breakers, while respecting DC power flow physics, thermal limits,
#   and basic connectivity constraints.
################################################################################

using CSV
using DataFrames
using Random
using Statistics
using JuMP
using Gurobi
using Graphs
using Printf
using MathOptInterface
const MOI = MathOptInterface

# ==============================================================================
# FIXED PARAMETERS (SCALING AND OBJECTIVE WEIGHTS)
# ==============================================================================

# Base power for per-unit normalization (DC model assumption)
const Pn = 1000.0

# Total available generation in per-unit
# Used to normalize Monte Carlo random injections
const T_Pg = 10000.0 / Pn

# Objective function weights
# - Lambda: maximize inter-area exchange
# - Connectivity: mild preference for closed breakers
# - Penalty: discourages unnecessary breaker power circulation
const W_LAMBDA       = 1.0
const W_CONNECTIVITY = 1e-4
const W_PENALTY      = 1e-7

# Big-M constants
# - M_flow: breaker flow when closed
# - M_angle: maximum phase angle difference
const M_flow  = T_Pg * 1.5
const M_angle = pi

# Base thermal limit for internal branches
# Interconnection branches are scaled up later
const BF_max_base = 400 / Pn

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

"""
    bridge_statistics(branch_df, breaker_df, Breaker_S)

Builds a graph consisting of:
- all permanent transmission branches
- only the breakers that are closed in the optimal solution

Returns:
- number of bridges (cut edges)
- list of bridges

Bridges correspond to single points of failure, which are discussed
as undesirable topologies in the paper.
"""
function bridge_statistics(branch_df, breaker_df, Breaker_S)
    G = SimpleGraph()

    # Permanent transmission branches
    for r in eachrow(branch_df)
        add_edge!(G, Int(r.From_Bus), Int(r.To_Bus))
    end

    # Add breakers only if closed
    for i in 1:nrow(breaker_df)
        if value(Breaker_S[i]) > 0.5
            r = breaker_df[i, :]
            add_edge!(G, Int(r.From_Bus), Int(r.To_Bus))
        end
    end

    br = bridges(G)
    br_list = [(src(e), dst(e)) for e in br]
    return length(br_list), br_list
end


"""
    max_interconnection_flow(branch_df, Angle, B_param)

Computes the total power exchanged between the two areas
using only interconnection branches.

This is the main performance indicator in the test case.
"""
function max_interconnection_flow(branch_df, Angle, B_param)
    sum(
        value(
            B_param[i] *
            (Angle[branch_df.From_Bus[i]] -
             Angle[branch_df.To_Bus[i]])
        ) * Pn
        for i in findall(branch_df.InterC .== 1)
    )
end

# ==============================================================================
# LOAD INPUT DATA
# ==============================================================================

# Bus data:
# - Bus index
# - Zone (1 or 2)
bus_df     = CSV.read("bus_data.csv", DataFrame)

# Branch data:
# - Permanent transmission lines
# - InterC = 1 indicates inter-area branches
branch_df  = CSV.read("branch_data.csv", DataFrame)

# Breaker data:
# - Switchable connections inside substations
breaker_df = CSV.read("breaker_data.csv", DataFrame)

nb = nrow(bus_df)     # number of buses
nl = nrow(branch_df)  # number of branches
nk = nrow(breaker_df) # number of breakers

# ==============================================================================
# BUILD OPTIMIZATION MODEL
# ==============================================================================

# Mixed-integer linear program
model = Model(Gurobi.Optimizer)
set_silent(model)

# ==============================================================================
# DECISION VARIABLES
# ==============================================================================

# Lambda scales the inter-area power exchange
@variable(model, Lambda >= 0)

# Voltage angles (DC power flow)
@variable(model, -M_angle <= Angle[bus_df.Bus] <= M_angle)

# Power flows on permanent branches
@variable(model, Branch_F[1:nl])

# Breaker status (1 = closed, 0 = open)
@variable(model, Breaker_S[1:nk], Bin)

# Power flows through breakers
@variable(model, -M_flow <= Breaker_F[1:nk] <= M_flow)

# Absolute value of breaker flows (for linear penalty)
@variable(model, Abs_Brk_F[1:nk] >= 0)

# ==============================================================================
# MODEL PARAMETERS
# ==============================================================================

@variable(model, B_param[1:nl]     in Parameter(0.0))  # susceptances
@variable(model, Fmax_param[1:nl]  in Parameter(0.0))  # thermal limits
@variable(model, Pgen_param[1:nb]  in Parameter(0.0))  # generation
@variable(model, Pload_param[1:nb] in Parameter(0.0))  # load

# Coupling parameters between zones (derived analytically)
@variable(model, alpha_param in Parameter(0.0))
@variable(model, beta_param  in Parameter(0.0))

# Reference (slack) bus
fix(Angle[1], 0.0; force = true)

# ==============================================================================
# CONSTRAINTS
# ==============================================================================

# Absolute value linearization
@constraint(model, [i=1:nk], Abs_Brk_F[i] >=  Breaker_F[i])
@constraint(model, [i=1:nk], Abs_Brk_F[i] >= -Breaker_F[i])

# DC power flow on permanent branches
@constraint(model, [i=1:nl],
    Branch_F[i] ==
        B_param[i] *
        (Angle[branch_df.From_Bus[i]] -
         Angle[branch_df.To_Bus[i]])
)

# Thermal limits
@constraint(model, [i=1:nl],  Branch_F[i] <=  Fmax_param[i])
@constraint(model, [i=1:nl], -Branch_F[i] <=  Fmax_param[i])

# Breaker flow only if breaker is closed
@constraint(model, [i=1:nk], Breaker_F[i] <=  M_flow * Breaker_S[i])
@constraint(model, [i=1:nk], Breaker_F[i] >= -M_flow * Breaker_S[i])

# Angle decoupling when breaker is open
@constraint(model, [i=1:nk],
    Angle[breaker_df.From_Bus[i]] -
    Angle[breaker_df.To_Bus[i]] <= M_angle * (1 - Breaker_S[i])
)
@constraint(model, [i=1:nk],
    Angle[breaker_df.From_Bus[i]] -
    Angle[breaker_df.To_Bus[i]] >= -M_angle * (1 - Breaker_S[i])
)

# ==============================================================================
# POWER BALANCE EQUATIONS (TWO-ZONE FORMULATION)
# ==============================================================================

for (i, r) in enumerate(eachrow(bus_df))
    b = r.Bus
    out_f = findall(branch_df.From_Bus .== b)
    in_f  = findall(branch_df.To_Bus   .== b)
    out_b = findall(breaker_df.From_Bus .== b)
    in_b  = findall(breaker_df.To_Bus   .== b)

    if r.Zone == 1
        # Exporting area
        @constraint(model,
            sum(Branch_F[j] for j in out_f) +
            sum(Breaker_F[j] for j in out_b) -
            sum(Branch_F[j] for j in in_f) -
            sum(Breaker_F[j] for j in in_b)
            ==
            Pgen_param[i] * Lambda - Pload_param[i]
        )
    else
        # Importing area
        @constraint(model,
            sum(Branch_F[j] for j in out_f) +
            sum(Breaker_F[j] for j in out_b) -
            sum(Branch_F[j] for j in in_f) -
            sum(Breaker_F[j] for j in in_b)
            ==
            Pgen_param[i] -
            Pload_param[i] * (alpha_param * Lambda + beta_param)
        )
    end
end

# Minimum connectivity: avoid isolated buses and radial cofiguration
for b in bus_df.Bus
    bk = findall((breaker_df.From_Bus .== b) .| (breaker_df.To_Bus .== b))
    br = findall((branch_df.From_Bus  .== b) .| (branch_df.To_Bus  .== b))
    @constraint(model, sum(Breaker_S[i] for i in bk) + length(br) >= 2)
end

# ==============================================================================
# OBJECTIVE FUNCTION
# ==============================================================================

@objective(model, Max,
    W_LAMBDA * Lambda +
    W_CONNECTIVITY * sum(Breaker_S) -
    W_PENALTY * sum(Abs_Brk_F)
)

# ==============================================================================
# RESULTS STORAGE
# ==============================================================================

results = DataFrame(
    seed = Int[],
    lambda = Float64[],
    n_bridges = Int[],
    bridges = Vector{Tuple{Int,Int}}[],
    n_closed_breakers = Int[],
    initial_exchange = Float64[],
    max_exchange = Float64[],
    solve_time = Float64[],
    node_count = Int[],
    simplex_iterations = Int[]
)

# ==============================================================================
# MONTE CARLO SIMULATION
# ==============================================================================

for seed in 100:10:200
    Random.seed!(seed)
    @printf("\nSeed: %d\n", seed)

    # Random generation and load profiles
    xgg = rand(nb) .* 500
    xll = rand(nb) .* 500
    # Imposing a total generation = Total load = T_Pg
    Pgen  = xgg ./ sum(xgg) .* T_Pg
    Pload = xll ./ sum(xll) .* T_Pg

    for i in 1:nb
        set_parameter_value(Pgen_param[i],  Pgen[i])
        set_parameter_value(Pload_param[i], Pload[i])
    end

    # Aggregate zone quantities
    TG1 = sum(Pgen[bus_df.Zone .== 1])
    TG2 = sum(Pgen[bus_df.Zone .== 2])
    TL1 = sum(Pload[bus_df.Zone .== 1])
    TL2 = sum(Pload[bus_df.Zone .== 2])

    # Analytical coupling coefficients ensure balancing.
	# When generated powers are scaled proportionally to their initial values. 
    set_parameter_value(alpha_param, TG1 / TL2)
    set_parameter_value(beta_param,  (TG2 - TL1) / TL2)

    initial_exchange = (TG1 - TL1) * Pn

    # Randomize branch parameters
    for i in 1:nl
        if branch_df.InterC[i] == 0
            Fmax = rand(0.99:0.0001:1.01) * BF_max_base
            B    = 1 / rand(0.99:0.0001:1.01)
        else
            Fmax = rand(0.99:0.0001:1.01) * 10 * BF_max_base
            B    = 1 / (rand(0.99:0.0001:1.01) * 2)
        end
        set_parameter_value(Fmax_param[i], Fmax)
        set_parameter_value(B_param[i],    B)
    end

    # Solve
    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL
	
		L = value(Lambda)
		@printf("\nOptimal Lambda: %.4f\n", L)
        
		n_br, br_list = bridge_statistics(branch_df, breaker_df, Breaker_S)
        n_closed = count(value(Breaker_S[i]) > 0.5 for i in 1:nk)

        push!(results, (
            seed,
            value(Lambda),
            n_br,
            br_list,
            n_closed,
            initial_exchange,
            max_interconnection_flow(branch_df, Angle, B_param),
            MOI.get(model, MOI.SolveTimeSec()),
            MOI.get(model, MOI.NodeCount()),
            MOI.get(model, MOI.SimplexIterations())
        ))
    end
end

# ==============================================================================
# EXPORT RESULTS
# ==============================================================================

CSV.write("grid_montecarlo_results.csv", results)
println("\nResults written to grid_montecarlo_results.csv")