using JuMP
using Gurobi
using CSV
using DataFrames
using Bijections
using LinearAlgebra

# Read the data
DATA_dir = joinpath("hsa") # Path to the data
file_name = "GSE10063"
raw_df = CSV.read(joinpath(DATA_dir, string(file_name, ".txt")), DataFrames.DataFrame, header=0)

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

# x_i = 1 if gene i is in the module
#     = 0 otherwise

@variable(model, x[i=1:n_genes], Bin)

# s_i = the number of additional connection allowed for gene i

@variable(model, s[i=1:n_genes], lower_bound=0, Int)

## Objective

# Max sum x_i

@objective(model, Max, sum(x))

## Constraints

# Number of edge linked to other gene in module for gene i must be >= size(module)-1
# Constraint not effective if gene i is not part of the module

@constraint(model, [i = 1:n_genes], x' * net_mat[i,:] + s[i] >= sum(x)-1 - (1-x[i])*n_genes)

# sum_s = 2*p where p is the prefined number

p = 0.2
@constraint(model, sum(s) <= 2*p*sum(x))

## Solve
optimize!(model)

# Get solution
x_sol = value.(x)

# Set of gene in the module (clique)
clique = Set()
for i in 1:n_genes
	if x_sol[i] == 1
		push!(clique, i)
	end
end

# Translate the genes back to their original number
module_gene = Set()
for g_index in clique
	push!(module_gene, b[g_index])
end