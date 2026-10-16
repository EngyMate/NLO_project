num_nodes = 11

edges = [
    (1, 2), (1, 11), (2, 3), (2, 11), (3, 4), (3, 9), (4, 5), (5, 6),
    (5, 8), (6, 7), (7, 8), (7, 9), (8, 9), (9, 10), (10, 11),
]
edge_susceptance = [
    -20.1, -22.3, -16.8, -17.2, -11.7, -19.4, -10.8, -12.3,
    -9.2, -13.9, -8.7, -11.3, -7.7, -13.5, -26.7,
]
edge_conductance = [
    4.12, 5.67, 2.41, 2.78, 1.98, 3.23, 1.59, 1.26,
    1.71, 1.11, 1.32, 2.01, 4.41, 2.14, 5.06,
]

# Both directions of every edge, with symmetric parameters
arcs = Tuple{Int,Int}[]
susceptance = Dict{Tuple{Int,Int},Float64}()
conductance = Dict{Tuple{Int,Int},Float64}()
for (index, (node_a, node_b)) in enumerate(edges)
    for arc in ((node_a, node_b), (node_b, node_a))
        push!(arcs, arc)
        susceptance[arc] = edge_susceptance[index]
        conductance[arc] = edge_conductance[index]
    end
end

# (node, maximum capacity [pu], cost [SEK/pu])
generators = [
    (2, 0.02, 175),
    (2, 0.15, 100),
    (2, 0.08, 150),
    (3, 0.07, 150),
    (4, 0.04, 300),
    (5, 0.17, 350),
    (7, 0.17, 400),
    (9, 0.26, 300),
    (9, 0.05, 200),
]

# (node, active demand [pu])
consumers = [
    (1, 0.10),
    (4, 0.19),
    (6, 0.11),
    (8, 0.09),
    (9, 0.21),
    (10, 0.05),
    (11, 0.04),
]

num_generators = length(generators)
num_consumers  = length(consumers)

generator_node     = [generator[1] for generator in generators]
generator_capacity = [generator[2] for generator in generators]
generator_price    = [generator[3] for generator in generators]
consumer_node      = [consumer[1] for consumer in consumers]
consumer_demand    = [consumer[2] for consumer in consumers]

generators_at_node = [findall(==(node), generator_node) for node in 1:num_nodes]
consumers_at_node  = [findall(==(node), consumer_node) for node in 1:num_nodes]

# Node totals, used for colouring nodes in the plots
node_capacity = [sum(generator_capacity[g] for g in generators_at_node[node]; init = 0.0) for node in 1:num_nodes]
node_demand   = [sum(consumer_demand[c] for c in consumers_at_node[node]; init = 0.0) for node in 1:num_nodes]

voltage_min    = 0.98
voltage_max    = 1.02
reactive_ratio = 0.03     # reactive capacity = reactive_ratio * active capacity

# Cost normalisation: mean full-output cost of a generator [SEK]
cost_reference = sum(generator_capacity .* generator_price) / num_generators
cost_weight    = 1.0      # tuning factor for the normalised cost term