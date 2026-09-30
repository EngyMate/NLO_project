using JuMP
import Ipopt

# Import data from the data file
include("intro_data.jl")

# Create the model object
the_model = Model(Ipopt.Optimizer)

# Create the variables and set their common lower bound
@variable(
    the_model,
    x[1:n_vars_first_dimension, 1:n_vars_second_dimension] >= lb,
)

# Set the upper bounds for the variables as named constraints.
# Keeping these as constraints (rather than variable bounds) is to illustrate
# below how this can be used to query dual variables.
@constraint(
    the_model,
    ub_constr[i = 1:n_vars_first_dimension, j = 1:n_vars_second_dimension],
    x[i, j] <= ub[i, j],
)

# Create the nonlinear (quadratic) objective sum((x[i,j] - 2)^2), which we minimize.
@objective(
    the_model,
    Min,
    sum(
        (x[i, j] - 2)^2
        for i in 1:n_vars_first_dimension, j in 1:n_vars_second_dimension
    ),
)

# Create the nonlinear constraint that the sum of all variables squared should
# be less than or equal to a given constant. The constrain is given the name
# "SOS_constr", which can later be referenced in the code (see below).
@constraint(
    the_model,
    SOS_constr,
    sum(
        x[i, j]^2
        for i in 1:n_vars_first_dimension, j in 1:n_vars_second_dimension
    ) <= sum_bound,
)

# Print the optimization problem in the terminal
println(the_model)

# Solve the optimization problem
optimize!(the_model)

# Print selected results for further analysis
# NOTE: This is the type of output you need to analyise in your project.
#       You can use the output to check if the solver has found a solution, what the solution is, and what the dual variables are.
println("") # Printing white line after solver output, before printing
println("Termination status: ", termination_status(the_model))
println("Optimal objective function value: ", objective_value(the_model))
println("Optimal point: ", value.(x))
println("Dual variables/Lagrange multipliers corresponding to some constraints:")
println(dual(SOS_constr))
println(dual.(ub_constr))
println(dual.(LowerBoundRef.(x)))
