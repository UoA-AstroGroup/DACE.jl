"""
Differential Algebra Monte Carlo (DAMC) method.

This example is taken from M. Valli, R. Armellin, P. Di Lizia, and M. R. Lavagna, 
“Nonlinear Mapping of Uncertainties in Celestial Mechanics,” 
Journal of Guidance, Control, and Dynamics, vol. 36, no. 1, pp. 48-63, Jan. 2013, doi: 10.2514/1.58068.

"""

using BenchmarkTools
using LinearAlgebra
using DACE

using CairoMakie

function lagrange_propagator(rv0::AbstractArray{T}, Δt::Float64, μ::Float64) where T<:Union{DA,Float64}

    # extract initial position and velocity vectors
    rr0 = rv0[1:3]
    vv0 = rv0[4:6]

    # compute angular momentum, initial radius and velocity
    hh = cross(rr0,vv0)
    h  = norm(hh)
    r0 = norm(rr0)
    v0 = norm(vv0)

    # compute semi-major axis and semilatus rectum
    a = μ/(2.0*μ/r0 - v0*v0)
    p = h*h/μ
    sigma0 = dot(rr0,vv0)/sqrt(μ)

    ord = DACE.getMaxOrder()
    iter = 0
    err = 1.0

    if (h isa DA)
        # println("Lagrange propagator for DA state")
        tol = 0.0
        maxiter = ord + 1
    else
        # println("Lagrange propagator for double state")
        tol = 1e-14
        maxiter = 20
    end

    if (DACE.cons(a) > 0.0)
        MmM0 = Δt * sqrt(μ/a^3)
        EmE0 = DACE.cons(MmM0)

        # println("Solving KeplerSolver equation...")
        while (err > tol && iter < maxiter)
            iter += 1
            fx0 = -MmM0 + EmE0 + sigma0/sqrt(a)*(1.0 - cos(EmE0)) - (1.0 - r0/a) * sin(EmE0)
            fxp = 1.0 + sigma0/sqrt(a) * sin(EmE0) - (1.0 - r0/a) * cos(EmE0)
            err = abs(DACE.cons(fx0 / fxp))
            EmE0 = EmE0 - fx0/fxp
        end

        theta = 2.0*atan(sqrt(a*p)*tan(EmE0/2), r0 + sigma0*sqrt(a)*tan(EmE0/2))
        r = p*r0/(r0 + (p - r0)*cos(theta) - sqrt(p)*sigma0*sin(theta))

        # compute the Lagrangian coefficients
        F = 1.0 - a/r0 * (1.0 - cos(EmE0))
        G = a*sigma0/sqrt(μ)*(1.0 - cos(EmE0)) + r0 * sqrt(a/μ) * sin(EmE0)
        Ft = -sqrt(μ*a)/(r*r0) * sin(EmE0)
        Gt = 1.0 - a/r * (1.0 - cos(EmE0))
    else
        println("Negative semimajor axis")
        NmN0 = Δt*sqrt(-μ/a^3)
        HmH0 = 0.0

        while (err > tol && iter < maxiter)
            iter += 1
            fx0 = -NmN0 - HmH0 + sigma0/sqrt(-a) * (-1.0 + cosh(HmH0)) + (1.0 - r0/a) * sinh(HmH0)
            fxp = -1.0 + sigma0/sqrt(-a)*sinh(HmH0) + (1.0 - r0/a)*cosh(HmH0)
            err = abs(DACE.cons(fx0 / fxp))
            HmH0 = HmH0 - fx0/fxp
        end

        F = 1.0 - a/r0 * (1.0 - cosh(HmH0))
        G = a*sigma0/sqrt(μ)*(1.0 - cosh(HmH0)) + r0 * sqrt(-a/μ) * sinh(HmH0)

        rv_tmp = F.*rr0 .+ G.*vv0
        r = norm(rv_tmp);
        Ft = - sqrt(μ*(-a))/(r*r0) * sinh(HmH0)
        Gt = 1.0 - a/r*(1.0 - cosh(HmH0))
    end

    rv_fin = zeros(typeof(h),6)
    rv_fin[1:3] = F.*rr0 .+ G.*vv0
    rv_fin[4:6] = Ft.*rr0 .+ Gt.*vv0

    return rv_fin
end;


# initial conditions
x0 = [-0.68787, -0.39713, 0.28448, -0.51331, 0.98266, 0.37611]
var0 = [1e-7, 1e-7, 1e-7, 1e-9, 1e-9, 1e-9]

# physical and propagation parameters
μ  = 1.0
tp = 2.0π
nr = 30
Δt = nr * tp

zs = 3.0
ns = 10000

# nominal solution
_ = lagrange_propagator(x0, Δt, μ)
@time x1 = lagrange_propagator(x0, Δt, μ)

# Monte Carlo simulation
u0_mc = randn(ns,length(x0))
δ0_mc = u0_mc.*(sqrt.(var0)')
x0_mc = x0' .+ δ0_mc

f_mc(xs) = stack([lagrange_propagator(x, Δt, μ) for x in eachrow(xs)], dims=1)

_ = f_mc(x0_mc)
@time x1_mc = f_mc(x0_mc)

DACE.init(8,length(var0))

# DAMC simulation (order 2)
DACE.pushTO(2)
x0_da_2nd = x0 .+ zs.*sqrt.(var0).*[DA(i,1) for i in eachindex(x0)]
_ = lagrange_propagator(x0_da_2nd, Δt, μ)
@time x1_da_2nd = lagrange_propagator(x0_da_2nd, Δt, μ)

x1_cp_2nd = DACE.compile(x1_da_2nd)
x1_ev_2nd = stack([DACE.eval(x1_cp_2nd, u0_mc[i,:]/zs) for i in 1:ns], dims=1)

# DAMC simulation (order 4)
DACE.pushTO(4)
x0_da_4th = x0 .+ zs.*sqrt.(var0).*[DA(i,1) for i in eachindex(x0)]
@time x1_da_4th = lagrange_propagator(x0_da_4th, Δt, μ)

x1_cp_4th = DACE.compile(x1_da_4th)
x1_ev_4th = stack([DACE.eval(x1_cp_4th, u0_mc[i,:]/zs) for i in 1:ns], dims=1)

# DAMC simulation (order 8)
DACE.pushTO(8)
x0_da_8th = x0 .+ zs.*sqrt.(var0).*[DA(i,1) for i in eachindex(x0)]
@time x1_da_8th = lagrange_propagator(x0_da_8th, Δt, μ)

x1_cp_8th = DACE.compile(x1_da_8th)
x1_ev_8th = stack([DACE.eval(x1_cp_8th, u0_mc[i,:]/zs) for i in 1:ns], dims=1)

# plot results
fig = Figure()
ax = Axis(fig[1,1], xlabel="x [-]", ylabel="y [-]")
scatter!(ax, x1_mc[:,1], x1_mc[:,2], marker=:diamond, color=:black, label="MC")
scatter!(ax, x1_ev_2nd[:,1], x1_ev_2nd[:,2], marker=:dtriangle, color="#666666", label="DAMC-2")
scatter!(ax, x1_ev_4th[:,1], x1_ev_4th[:,2], marker=:circle, color="#999999", label="DAMC-4")
scatter!(ax, x1_ev_8th[:,1], x1_ev_8th[:,2], marker=:utriangle, color="#CCCCCC", label="DAMC-8")
axislegend(ax, position=:rt)
save("damc_kepler.pdf", fig)
