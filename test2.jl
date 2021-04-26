using JuMP
using Gurobi
using CSV
using DataFrames
using Bijections
using LinearAlgebra

# Read the data
DATA_dir = joinpath("hsa") # Path to the data
raw_df = CSV.read(joinpath(DATA_dir, "GSE10063.txt"), DataFrames.DataFrame, header=0)

# Create the set of genes
genes_set = Set(vcat(raw_df[:,1], raw_df[:,2]))
n_genes = length(genes_set) # Number of genes

# Create a bijection from gene number to gene index
b = Bijection()
i = 1
for g in genes_set
	b[i] = g
	i+=1
end

# Create the co-expression network matrix
threshold = 1
net_mat = zeros(Int8, n_genes, n_genes) # Initialize matrix to zeroes

for i in 1:DataFrames.nrow(raw_df)
	if raw_df[i, 3] >= threshold
		net_mat[b(raw_df[i, 1]), b(raw_df[i, 2])] = 1
		net_mat[b(raw_df[i, 2]), b(raw_df[i, 1])] = 1
	end
end

# Write a model and an algorithm that finds the largest module 
# where the module is defined as a set of genes in which 
# all pairwise co-expressions are present, except a predefined number of them

model = Model(Gurobi.Optimizer)

## Variables

# X_i = 1 if gene i is in the module
#     = 0 otherwise

@variable(model, x[i=1:n_genes], Bin)

# Y_ij = 1 when gene i and j don't have to be related
#      = 0 otherwise 

@variable(model, y[i=1:n_genes, j=1:n_genes; i > j], Bin)

## Objective

# Max sum x_i

@objective(model, Max, sum(x))

## Constraints

# X_i + X_j <= 1 + net_mat[i,j] + Y_ij (with i=/=j if diag is 0) for all i and j

@constraint(model, [i = 1:n_genes, j = 1:n_genes; i > j], x[i] + x[j] <= 1 + net_mat[i,j] + y[i,j])

# sum_i sum_j Y_ij = 2*p where p is the prefined number

p = 5
@constraint(model, sum(y) == p)

# Y_ij = Y_ji for all i and j

#@constraint(model, [i = 1:n_genes, j = 1:n_genes; i != j], y[i,j] == y[j,i])

## Solve
optimize!(model)

# Get solution
x_sol = value.(x)
objective_value(model)