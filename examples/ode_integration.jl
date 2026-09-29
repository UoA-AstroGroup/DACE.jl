# # ODE integration example
#
# This tutorial shows how to perform the numerical integration of an ordinary differential 
# equation (ODE) using DA objects as the state components' type.
# It also demonstrates how to extract the state transition matrix (STM) from the polynomial 
# expansion of the dynamics flow.
#
# In this case, the state vector represents an object in orbit around a central body and 
# subject only to the gravitational pull of the latter.
# Its motion is described by the differential equations for the Kepler problem expressed in 
# Cartesian coordinates.

# ## Install dependencies
#
# Make sure the required packages are installed

# ```julia
# using Pkg
# Pkg.add("DACE")
# ```

# ## Using DACE
#
# Load the required modules

using OrdinaryDiffEq
using DACE

# Define the parameters and intial conditions for a circular orbit with a normalized radius 
# equal to one
μ = 1.0
x0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]

# Set the integration time span equal to one revolution
t0 = 0.0
tf = 2π

# Define the equations of motion for the resticted two-body problem
function kepler_ode!(du, u, μ, _)
    du[1:3] .= u[4:6]
    du[4:6] .= -u[1:3] * μ ./ sum(u[1:3].^2)^(3/2)
end

# Compute the nominal solution
prob = ODEProblem(kepler_ode!, x0, [t0, tf], μ)
sol = solve(prob, Vern9(), abstol=1e-12, reltol=1e-12)

xf = sol.u[end]

# Initialize DACE to compute second-order expansions
DACE.init(2,6)

# Define the initial conditions as a DA vector
dx_dace = [DA(i,1.0) for i in 1:6]
x0_dace = x0 .+ dx_dace

println(x0_dace)

# Compute the polynomial expansion of the dynamics flow
sol_dace = solve(remake(prob, u0=x0_dace), Vern9(), abstol=1e-12, reltol=1e-12)

xf_dace = sol_dace.u[end]
xf_cons = DACE.cons.(xf_dace)

# Verify the (nominal) final state

println("Max abs error on final state: " * string(maximum(abs.(xf_cons - xf))))
@assert maximum(abs.(xf_cons - xf)) < 1e-9

# Extract the state transition matrix (STM)

stm_dace = DACE.jacobian(xf_dace)
stm_cons = DACE.cons.(stm_dace)
