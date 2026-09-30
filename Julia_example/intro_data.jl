# Small data file to show how one can import julia files into other julia files.
# For example, this gives a convenient way to create files containing all the data for a problem.

n_vars_first_dimension = 2 # Number of variables in first dimension
n_vars_second_dimension = 2 # Number of variables in second dimension
lb = 0 # Lower bounds for the varaibles; same for all of them
ub = [[3, 2] [4, 1]] # Upper bounds for the variables
sum_bound = 7.75 # Constraint on sum of variables
