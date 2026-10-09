# Data for the power-system model

N = 11 # Number of nodes

# ---- Edges (Table 2) -------------------------------------------------------
edges = [
    (1, 2), (1, 11), (2, 3), (2, 11), (3, 4), (3, 9), (4, 5), (5, 6),
    (5, 8), (6, 7), (7, 8), (7, 9), (8, 9), (9, 10), (10, 11),
]
b_vals = [
    -20.1, -22.3, -16.8, -17.2, -11.7, -19.4, -10.8, -12.3,
    -9.2, -13.9, -8.7, -11.3, -7.7, -13.5, -26.7,
]
g_vals = [
    4.12, 5.67, 2.41, 2.78, 1.98, 3.23, 1.59, 1.26,
    1.71, 1.11, 1.32, 2.01, 4.41, 2.14, 5.06,
]

# Symmetric parameters: b_kl = b_lk and g_kl = g_lk.
# `arcs` contains both directions of every edge.
arcs = Tuple{Int,Int}[]
b = Dict{Tuple{Int,Int},Float64}()
g = Dict{Tuple{Int,Int},Float64}()
for (idx, (k, l)) in enumerate(edges)
    for arc in ((k, l), (l, k))
        push!(arcs, arc)
        b[arc] = b_vals[idx]
        g[arc] = g_vals[idx]
    end
end

# ---- Generators (Table 3): (node, max capacity [pu], cost [SEK/pu]) --------
generators = [
    (2, 0.02, 175),  # G1
    (2, 0.15, 100),  # G2
    (2, 0.08, 150),  # G3
    (3, 0.07, 150),  # G4
    (4, 0.04, 300),  # G5
    (5, 0.17, 350),  # G6
    (7, 0.17, 400),  # G7
    (9, 0.26, 300),  # G8
    (9, 0.05, 200),  # G9
]

# ---- Consumers (Table 4): (node, demand [pu]) ------------------------------
consumers = [
    (1, 0.10),  # C1
    (4, 0.19),  # C2
    (6, 0.11),  # C3
    (8, 0.09),  # C4
    (9, 0.21),  # C5
    (10, 0.05), # C6
    (11, 0.04), # C7
]

# ---- Aggregate per node (nodes without generator/consumer get 0) -----------
gmax = zeros(N)   # g_k^max : total generation capacity at node k
c = zeros(N)      # c_k     : total active demand at node k
for (node, cap, _) in generators
    gmax[node] += cap
end
for (node, d) in consumers
    c[node] += d
end

# ---- Other parameters ------------------------------------------------------
Vmin = 0.98
Vmax = 1.02
q_ratio = 0.03   # reactive capacity ratio used in g^[r+] and g^[r-]
