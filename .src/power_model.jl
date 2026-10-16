using JuMP
import Ipopt

include("power_data.jl")

the_model = Model(Ipopt.Optimizer)

# ---- Variables -------------------------------------------------------------
@variable(the_model, 0 <= ra[1:NG] <= 1)   # r^[a]_j: fraction of generator j's active capacity used
@variable(the_model, 0 <= rr[1:NG] <= 1)   # r^[r]_j: fraction of generator j's reactive capacity used
@variable(the_model, -Float64(π) <= θ[1:N] <= Float64(π))
@variable(the_model, Vmin <= V[1:N] <= Vmax, start = 1.0)

# ---- Line flows (both directions of every edge) ----------------------------
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

# ---- Node quantities: sums over neighbours, generators and consumers at node k
@expression(
    the_model,
    ga[k = 1:N],
    sum(p[(l, k)] for (l, kk) in arcs if kk == k) +
    sum(ra[j] * gen_cap[j] for j in gens_at[k]; init = 0.0),
)
@expression(
    the_model,
    gr_plus[k = 1:N],
    sum(q[(l, k)] for (l, kk) in arcs if kk == k) +
    q_ratio * sum(ra[j] * gen_cap[j] for j in gens_at[k]; init = 0.0),
)
@expression(
    the_model,
    gr_minus[k = 1:N],
    q_ratio * sum(rr[j] * gen_cap[j] for j in gens_at[k]; init = 0.0),
)

# ---- System-wide balance: Σ generation - Σ consumption = 0 ------------------
@constraint(
    the_model,
    system_balance,
    sum(ra[j] * gen_cap[j] for j in 1:NG) - sum(con_dem[j] for j in 1:NC) == 0,
)

# ---- Objective: f = Σ_k [ ΔP_k² + Δq_k² + λ · (Σ_{j∈G(k)} r_j cap_j price_j) / cost_ref ]
@expression(the_model, ΔP[k = 1:N], ga[k] - sum(con_dem[j] for j in cons_at[k]; init = 0.0))
@expression(the_model, Δq[k = 1:N], gr_plus[k] - gr_minus[k])
@expression(
    the_model,
    gen_cost[k = 1:N],   # normalised cost of the generators at node k
    sum(ra[j] * gen_cap[j] * gen_price[j] for j in gens_at[k]; init = 0.0) / cost_ref,
)
@objective(the_model, Min, sum(ΔP[k]^2 + Δq[k]^2 + gen_cost[k] for k in 1:N))

println(the_model)
optimize!(the_model)

println("")
println("Termination status: ", termination_status(the_model))
println("Optimal objective function value: ", objective_value(the_model))
println("r^[a] per generator = ", round.(value.(ra); digits = 4))
println("r^[r] per generator = ", round.(value.(rr); digits = 4))
println("V       = ", round.(value.(V); digits = 4))
println("θ [rad] = ", round.(value.(θ); digits = 4))
println("ΔP      = ", round.(value.(ΔP); digits = 5))
println("Δq      = ", round.(value.(Δq); digits = 5))
println("Generation cost [SEK] = ", sum(value.(gen_cost*cost_ref)))
println("Weighted mismatch     = ", w_mismatch * sum(value.(ΔP).^2 .+ value.(Δq).^2))
println("Total generation = ", sum(value.(ra) .* gen_cap), "  total demand = ", sum(con_dem))
println("Dual of system balance = ", dual(system_balance))

using Plots

# Node positions (pixels read from your illustration; y points down)
px = Dict(1 => (215, 50),  2 => (157, 152), 11 => (288, 187), 3 => (94, 262),
          10 => (252, 287), 9 => (247, 344), 4 => (77, 355),  8 => (177, 392),
          5 => (117, 421),  7 => (228, 435), 6 => (160, 472))
xy = [(px[k][1] / 100, -px[k][2] / 100) for k in 1:N]

# Where the text block of each node is placed: (dx, dy, alignment)
label_pos = Dict(k => (-0.30, 0.0, :right) for k in 1:N)          # default: left of node
for k in (1, 11, 10, 9, 7)                                         # these go to the right
    label_pos[k] = (0.30, 0.0, :left)
end
label_pos[8] = (-0.30, 0.30, :right)                               # keep clear of edge 5–8

const MILLI = 1000   # pu -> mPU

half_disk(cx, cy, r, a0, a1; n = 40) =
    Shape([cx; cx .+ r .* cos.(range(a0, a1; length = n))],
          [cy; cy .+ r .* sin.(range(a0, a1; length = n))])

function plot_network(flow, quantity::String, reactive::Bool)
    # Net flow over each edge in mPU: average of what leaves one end and what leaves the other
    net  = Dict(e => MILLI * (value(flow[(e[1], e[2])]) - value(flow[(e[2], e[1])])) / 2 for e in edges)
    fmax = maximum(abs.(values(net)))
    cg   = cgrad(:YlOrRd)

    # Generator output: active = r^[a]·cap, reactive = q_ratio·r^[r]·cap
    rfrac = reactive ? value.(rr) : value.(ra)
    scale = reactive ? q_ratio : 1.0

    plt = plot(; aspect_ratio = :equal, legend = false, axis = false, grid = false,
                 framestyle = :none, xlims = (-1.0, 4.4), ylims = (-5.2, -0.2),
                 size = (1000, 850), title = "$quantity power at the optimum [mPU]")

    # Arrows
    R = 0.2   # shorten arrows so they start/end at the node edge
    for (k, l) in edges
        P = net[(k, l)]
        s, t = P >= 0 ? (k, l) : (l, k)
        dx, dy = xy[t][1] - xy[s][1], xy[t][2] - xy[s][2]
        len = hypot(dx, dy)
        ux, uy = dx / len, dy / len
        x0, y0 = xy[s][1] + R * ux, xy[s][2] + R * uy
        x1, y1 = xy[t][1] - R * ux, xy[t][2] - R * uy
        plot!(plt, [x0, x1], [y0, y1]; arrow = :closed,
              color = get(cg, abs(P) / fmax), lw = 3 + 4 * abs(P) / fmax)
        annotate!(plt, (x0 + x1) / 2 - 0.2 * uy, (y0 + y1) / 2 + 0.2 * ux,
                  text(string(round(abs(P); digits = 1)), 8, :black))
    end

    # Colour bar: 0 to the largest arrow value
    scatter!(plt, [NaN], [NaN]; marker_z = [0.0], clims = (0, fmax), color = cg,
             colorbar = true, colorbar_title = "|flow| [mPU]", label = false)

    # Nodes: green = generation, red = consumption, split circle = both
    for k in 1:N
        x, y, r = xy[k][1], xy[k][2], 0.17
        gen, load = gmax[k] > 0, c[k] > 0
        if gen && load
            plot!(plt, half_disk(x, y, r, 0, π);  fillcolor = :seagreen, linecolor = :black, label = false)
            plot!(plt, half_disk(x, y, r, π, 2π); fillcolor = :indianred, linecolor = :black, label = false)
        else
            col = gen ? :seagreen : load ? :indianred : :gray60
            plot!(plt, half_disk(x, y, r, 0, 2π); fillcolor = col, linecolor = :black, label = false)
        end
        annotate!(plt, x, y, text(string(k), 9, :white, :center))
    end

    # Values next to each node: generators as output / capacity [mPU] (% of capacity), consumers as demand
    for k in 1:N
        lines = String[]
        for (i, (node, cap, _)) in enumerate(generators)
            node == k || continue
            val  = MILLI * scale * rfrac[k] * cap
            vmax = MILLI * scale * cap
            push!(lines, "G$i: $(round(val; digits = 1)) / $(round(vmax; digits = 1)) mPU ($(round(100 * rfrac[k]; digits = 1))%)")
        end
        if !reactive
            for (j, (node, d)) in enumerate(consumers)
                node == k && push!(lines, "C$j: $(round(MILLI * d; digits = 1)) mPU")
            end
        end
        isempty(lines) && continue
        dx, dy, al = label_pos[k]
        annotate!(plt, xy[k][1] + dx, xy[k][2] + dy, text(join(lines, "\n"), 7, al))
    end

    # Small colour label (top left)
    scatter!(plt, [-0.8], [-0.5]; markersize = 9, markershape = :square, color = :seagreen, markerstrokecolor = :black)
    annotate!(plt, -0.65, -0.5, text("generation", 9, :left))
    scatter!(plt, [-0.8], [-0.8]; markersize = 9, markershape = :square, color = :indianred, markerstrokecolor = :black)
    annotate!(plt, -0.65, -0.8, text("consumption", 9, :left))

    return plt
end

plt_p = plot_network(p, "Active",   false)
plt_q = plot_network(q, "Reactive", true)
display(plt_p);  savefig(plt_p, "network_active.png")
display(plt_q);  savefig(plt_q, "network_reactive.png")

