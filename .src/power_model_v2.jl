using JuMP
import Ipopt
using Plots

include("power_data_v2.jl")

model = Model(Ipopt.Optimizer)

# ---- Variables -------------------------------------------------------------
@variable(model, 0 <= active_usage[1:num_generators] <= 1)
@variable(model, -1 <= reactive_usage[1:num_generators] <= 1, start = 0.0)
@variable(model, -Float64(π) <= angle[1:num_nodes] <= Float64(π))
@variable(model, voltage_min <= voltage[1:num_nodes] <= voltage_max, start = 1.0)

# ---- Line flows (both directions of every edge) ----------------------------
@expression(
    model,
    active_flow[(sender, receiver) in arcs],
    voltage[sender]^2 * conductance[(sender, receiver)] -
    voltage[sender] * voltage[receiver] * (
        conductance[(sender, receiver)] * cos(angle[sender] - angle[receiver]) +
        susceptance[(sender, receiver)] * sin(angle[sender] - angle[receiver])
    ),
)
@expression(
    model,
    reactive_flow[(sender, receiver) in arcs],
    -voltage[sender]^2 * susceptance[(sender, receiver)] +
    voltage[sender] * voltage[receiver] * (
        susceptance[(sender, receiver)] * cos(angle[sender] - angle[receiver]) -
        conductance[(sender, receiver)] * sin(angle[sender] - angle[receiver])
    ),
)

# ---- Node quantities (flows measured at the node, leaving into the lines) ---
@expression(
    model,
    active_mismatch[node = 1:num_nodes],
    sum(active_usage[g] * generator_capacity[g] for g in generators_at_node[node]; init = 0.0) -
    sum(consumer_demand[c] for c in consumers_at_node[node]; init = 0.0) -
    sum(active_flow[(node, neighbour)] for (sender, neighbour) in arcs if sender == node),
)
@expression(
    model,
    reactive_mismatch[node = 1:num_nodes],
    reactive_ratio * sum(reactive_usage[g] * generator_capacity[g] for g in generators_at_node[node]; init = 0.0) -
    sum(reactive_flow[(node, neighbour)] for (sender, neighbour) in arcs if sender == node),
)
@expression(
    model,
    normalized_cost[node = 1:num_nodes],
    sum(active_usage[g] * generator_capacity[g] * generator_price[g] for g in generators_at_node[node]; init = 0.0) / cost_reference,
)

@constraint(
    model,
    system_balance,
    sum(active_usage[g] * generator_capacity[g] for g in 1:num_generators) >= sum(consumer_demand),
)

tuning_param = 0.00001
# ---- Objective -------------------------------------------------------------
@objective(
    model,
    Min,
    sum(active_mismatch[node]^2 + reactive_mismatch[node]^2 + tuning_param*cost_weight * normalized_cost[node] for node in 1:num_nodes),
)

println(model)
optimize!(model)

println("")
println("Termination status: ", termination_status(model))
println("Optimal objective function value: ", objective_value(model))
println("Active usage per generator   = ", round.(value.(active_usage); digits = 4))
println("Reactive usage per generator = ", round.(value.(reactive_usage); digits = 4))
println("Voltage            = ", round.(value.(voltage); digits = 4))
println("Angle [rad]        = ", round.(value.(angle); digits = 4))
println("Active mismatch    = ", round.(value.(active_mismatch); digits = 5))
println("Reactive mismatch  = ", round.(value.(reactive_mismatch); digits = 5))
println("Generation cost [SEK] = ", cost_reference * sum(value.(normalized_cost)))
println("Total generation = ", sum(value.(active_usage) .* generator_capacity), "  total demand = ", sum(consumer_demand))
println("Total active loss = ", sum(value(active_flow[arc]) for arc in arcs))
println("Dual of system balance = ", dual(system_balance))

# ---- Plots -----------------------------------------------------------------
pixel_position = Dict(1 => (215, 50),  2 => (157, 152), 11 => (288, 187), 3 => (94, 262),
                      10 => (252, 287), 9 => (247, 344), 4 => (77, 355),  8 => (177, 392),
                      5 => (117, 421),  7 => (228, 435), 6 => (160, 472))
node_position = [(pixel_position[node][1] / 100, -pixel_position[node][2] / 100) for node in 1:num_nodes]
network_center = (sum(first.(node_position)) / num_nodes, sum(last.(node_position)) / num_nodes)

# Which side column each node's text goes to
left_nodes  = [2, 3, 4, 5, 6, 8]
right_nodes = [1, 11, 10, 9, 7]
left_anchor_x  = 0.35
right_anchor_x = 3.30

# Manual nudges of arrow labels, keyed by (smaller node, larger node), e.g. (5, 6) => (0.0, -0.1)
arrow_label_shift = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()

const MILLI = 1000   # pu -> mPU
show_mismatch = false   # set to true to also print ΔP_k / g^[r]_k in the node text

subscript_digits = Dict('0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
                        '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉')
subscript(number) = join(get(subscript_digits, character, character) for character in string(number))

function format_number(number, digits)
    rounded = round(number; digits = digits)
    return string(rounded == 0 ? abs(rounded) : rounded)
end

function half_disk(center_x, center_y, radius, angle_start, angle_end; points = 40)
    angles = range(angle_start, angle_end; length = points)
    return Shape([center_x; center_x .+ radius .* cos.(angles)],
                 [center_y; center_y .+ radius .* sin.(angles)])
end

# Vertical positions of the text blocks in one side column: close to their node, but never overlapping
function column_positions(nodes, text_blocks; line_height = 0.10, gap = 0.12)
    positions = Dict{Int,Float64}()
    previous_bottom = Inf
    for node in sort(nodes; by = node -> -node_position[node][2] * -1, rev = true)
        height = line_height * length(text_blocks[node])
        center = node_position[node][2]
        if center + height / 2 > previous_bottom - gap
            center = previous_bottom - gap - height / 2
        end
        positions[node] = center
        previous_bottom = center - height / 2
    end
    return positions
end

function plot_network(flow, is_reactive::Bool)
    quantity_name  = is_reactive ? "Reactive" : "Active"
    flow_symbol    = is_reactive ? "q" : "p"
    usage_symbol   = is_reactive ? "rʳ" : "rᵃ"
    usage          = is_reactive ? value.(reactive_usage) : value.(active_usage)
    capacity_scale = is_reactive ? reactive_ratio : 1.0
    digits         = is_reactive ? 2 : 1
    colors         = cgrad(:YlOrRd)

    # For every edge: sender, receiver (from the net flow direction) and the flow leaving the sender
    arrows = []
    for (node_a, node_b) in edges
        net_flow = value(flow[(node_a, node_b)]) - value(flow[(node_b, node_a)])
        sender, receiver = net_flow >= 0 ? (node_a, node_b) : (node_b, node_a)
        push!(arrows, (sender, receiver, MILLI * value(flow[(sender, receiver)])))
    end
    max_flow = maximum(abs(flow_value) for (_, _, flow_value) in arrows)

    # Text block of every node: V and θ, then its generators and consumers
    text_blocks = Dict{Int,Vector{String}}()
    for node in 1:num_nodes
        lines = String[
            "V$(subscript(node)) = $(format_number(value(voltage[node]), 4)) pu,  " *
            "θ$(subscript(node)) = $(format_number(value(angle[node]), 4)) rad",
        ]
        for generator in generators_at_node[node]
            output   = MILLI * capacity_scale * usage[generator] * generator_capacity[generator]
            capacity = MILLI * capacity_scale * generator_capacity[generator]
            push!(lines,
                  "G$(subscript(generator)): $usage_symbol$(subscript(generator)) = $(format_number(usage[generator], 3))" *
                  "  ($(format_number(output, digits)) / $(format_number(capacity, digits)) mPU)")
        end
        if !is_reactive
            for consumer in consumers_at_node[node]
                push!(lines, "C$(subscript(consumer)): d$(subscript(consumer)) = $(format_number(MILLI * consumer_demand[consumer], 1)) mPU")
            end
        end
        if show_mismatch
            mismatch = is_reactive ? value(reactive_mismatch[node]) : value(active_mismatch[node])
            symbol   = is_reactive ? "gʳ" : "ΔP"
            push!(lines, "$symbol$(subscript(node)) = $(format_number(MILLI * mismatch, 2)) mPU")
        end
        text_blocks[node] = lines
    end
    left_y  = column_positions(left_nodes,  text_blocks)
    right_y = column_positions(right_nodes, text_blocks)

    figure = plot(; aspect_ratio = :equal, legend = false, axis = false, grid = false,
                    framestyle = :none, xlims = (-3.2, 5.6), ylims = (-5.6, -0.1),
                    size = (1500, 950), title = "$quantity_name power flow at the optimum")

    # Arrows, labelled with the flow variable leaving the sending node
    arrow_gap = 0.2
    for (sender, receiver, flow_value) in arrows
        delta_x = node_position[receiver][1] - node_position[sender][1]
        delta_y = node_position[receiver][2] - node_position[sender][2]
        length_xy = hypot(delta_x, delta_y)
        unit_x, unit_y = delta_x / length_xy, delta_y / length_xy
        start_x = node_position[sender][1] + arrow_gap * unit_x
        start_y = node_position[sender][2] + arrow_gap * unit_y
        end_x   = node_position[receiver][1] - arrow_gap * unit_x
        end_y   = node_position[receiver][2] - arrow_gap * unit_y
        plot!(figure, [start_x, end_x], [start_y, end_y]; arrow = :closed,
              color = get(colors, abs(flow_value) / max_flow),
              linewidth = 3 + 4 * abs(flow_value) / max_flow)

        # Label on the side of the edge that points away from the middle of the network
        middle_x, middle_y = (start_x + end_x) / 2, (start_y + end_y) / 2
        normal_x, normal_y = -unit_y, unit_x
        side = sign(normal_x * (middle_x - network_center[1]) + normal_y * (middle_y - network_center[2]))
        side = side == 0 ? 1.0 : side
        shift_x, shift_y = get(arrow_label_shift, (min(sender, receiver), max(sender, receiver)), (0.0, 0.0))
        label = "$flow_symbol$(subscript(sender)),$(subscript(receiver)) = $(format_number(flow_value, digits))"
        annotate!(figure, middle_x + side * 0.25 * normal_x + shift_x, middle_y + side * 0.25 * normal_y + shift_y,
                  text(label, 7, :black))
    end

    scatter!(figure, [NaN], [NaN]; marker_z = [0.0], clims = (0, max_flow), color = colors,
             colorbar = true, colorbar_title = "|$flow_symbol| [mPU]", label = false)

    # Leader lines and text blocks in the side columns
    for node in left_nodes
        plot!(figure, [node_position[node][1] - 0.17, left_anchor_x + 0.05], [node_position[node][2], left_y[node]];
              color = :gray70, linewidth = 0.8, label = false)
        annotate!(figure, left_anchor_x, left_y[node], text(join(text_blocks[node], "\n"), 7, :right))
    end
    for node in right_nodes
        plot!(figure, [node_position[node][1] + 0.17, right_anchor_x - 0.05], [node_position[node][2], right_y[node]];
              color = :gray70, linewidth = 0.8, label = false)
        annotate!(figure, right_anchor_x, right_y[node], text(join(text_blocks[node], "\n"), 7, :left))
    end

    # Nodes: green = generation, red = consumption, split circle = both
    for node in 1:num_nodes
        x, y, radius = node_position[node][1], node_position[node][2], 0.17
        has_generation  = node_capacity[node] > 0
        has_consumption = node_demand[node] > 0
        if has_generation && has_consumption
            plot!(figure, half_disk(x, y, radius, 0, π);  fillcolor = :seagreen,  linecolor = :black, label = false)
            plot!(figure, half_disk(x, y, radius, π, 2π); fillcolor = :indianred, linecolor = :black, label = false)
        else
            fill_color = has_generation ? :seagreen : has_consumption ? :indianred : :gray60
            plot!(figure, half_disk(x, y, radius, 0, 2π); fillcolor = fill_color, linecolor = :black, label = false)
        end
        annotate!(figure, x, y, text(string(node), 9, :white, :center))
    end

    # Key (top left): colours and the meaning of the symbols
    scatter!(figure, [-3.1], [-0.4]; markersize = 9, markershape = :square, color = :seagreen, markerstrokecolor = :black)
    annotate!(figure, -2.95, -0.4, text("generation", 9, :left))
    scatter!(figure, [-3.1], [-0.7]; markersize = 9, markershape = :square, color = :indianred, markerstrokecolor = :black)
    annotate!(figure, -2.95, -0.7, text("consumption", 9, :left))
    key_lines = [
        "$(flow_symbol)ₖₗ: $(lowercase(quantity_name)) flow from node k to node l [mPU]",
        "Vₖ, θₖ: voltage [pu] and angle [rad] at node k",
        "$(usage_symbol)ⱼ: fraction of the $(lowercase(quantity_name)) capacity of generator j in use" *
            (is_reactive ? " (negative = absorbing)" : ""),
        is_reactive ? "output = 0.03 · rʳⱼ · capⱼ" : "dᵢ: demand of consumer i",
    ]
    annotate!(figure, -3.1, -1.2, text(join(key_lines, "\n"), 7, :left))

    return figure
end

active_figure   = plot_network(active_flow,   false)
reactive_figure = plot_network(reactive_flow, true)
display(active_figure);   savefig(active_figure,   "network_active.png")
display(reactive_figure); savefig(reactive_figure, "network_reactive.png")