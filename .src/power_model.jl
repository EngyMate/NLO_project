using JuMP
import Ipopt

# Import data from the data file
include("power_data.jl")

the_model = Model(Ipopt.Optimizer)

# ---- Variables -------------------------------------------------------------
# r^[a] + s^[a] = 1, s >= 0, r >= 0   <=>   0 <= r^[a] <= 1 (same for r^[r])
@variable(the_model, 0 <= ra[1:N] <= 1)   # r^[a]: fraction of active capacity used
@variable(the_model, 0 <= rr[1:N] <= 1)   # r^[r]: fraction of reactive capacity used

@variable(the_model, -Float64(π) <= θ[1:N] <=Float64(π))

# V + s^[v+] = 1.02, V - s^[v-] = 0.98, slacks >= 0   <=>   0.98 <= V <= 1.02
@variable(the_model, Vmin <= V[1:N] <= Vmax, start = 1.0)

# ---- Line flows (both directions of every edge) ----------------------------
# p_kl = V_k^2 g_kl - V_k V_l ( g_kl cos(θ_k-θ_l) + b_kl sin(θ_k-θ_l) )
# q_kl = -V_k^2 b_kl + V_k V_l ( b_kl cos(θ_k-θ_l) - g_kl sin(θ_k-θ_l) )
@expression(
    the_model,
    p[(k, l) in arcs],
    V[k]^2 * g[(k, l)] -
    V[k] * V[l] * (g[(k, l)] * cos(θ[k] - θ[l]) + b[(k, l)] * sin(θ[k] - θ[l])),
)
@expression(
    the_model,
    q[(k, l) in arcs],
    -V[k]^2 * b[(k, l)] +
    V[k] * V[l] * (b[(k, l)] * cos(θ[k] - θ[l]) - g[(k, l)] * sin(θ[k] - θ[l])),
)

# ---- Generation quantities (sums over neighbours l of node k) --------------
@expression(
    the_model,
    ga[k = 1:N],
    sum(p[(l, k)] for (l, kk) in arcs if kk == k) + ra[k] * gmax[k],
)
@expression(
    the_model,
    gr_plus[k = 1:N],
    sum(q[(l, k)] for (l, kk) in arcs if kk == k) + q_ratio * ra[k] * gmax[k],
)
@expression(the_model, gr_minus[k = 1:N], q_ratio * rr[k] * gmax[k])

# ---- Objective: f = Σ_i ΔP_i^2 + Δq_i^2, with
#      ΔP_i = g^[a]_i - c_i  and  Δq_i = g^[r+]_i - g^[r-]_i ------------------
@expression(the_model, ΔP[i = 1:N], ga[i] - c[i])
@expression(the_model, Δq[i = 1:N], gr_plus[i] - gr_minus[i])
@objective(the_model, Min, sum(ΔP[i]^2 + Δq[i]^2 for i in 1:N))

println(the_model)

optimize!(the_model)

println("")
println("Termination status: ", termination_status(the_model))
println("Optimal objective function value: ", objective_value(the_model))
println("r^[a]   = ", round.(value.(ra); digits = 4))
println("r^[r]   = ", round.(value.(rr); digits = 4))
println("V       = ", round.(value.(V); digits = 4))
println("θ [rad] = ", round.(value.(θ); digits = 4))
println("ΔP      = ", round.(value.(ΔP); digits = 5))
println("Δq      = ", round.(value.(Δq); digits = 5))
