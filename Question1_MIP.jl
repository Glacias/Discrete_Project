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
# all pairwise co-expressions are present

model = Model(Gurobi.Optimizer)

## Variables

# X_i = 1 if gene i is in the module
#     = 0 otherwise

@variable(model, x[i=1:n_genes], Bin)

## Objective

# Max sum x_i

@objective(model, Max, sum(x))

## Constraints

# X_i + X_j <= 1 + net_mat[i,j] (with i=/=j if diag is 0)

@constraint(model, [i = 1:n_genes, j = 1:n_genes; i > j], x[i] + x[j] <= 1 + net_mat[i,j])

## Solve
optimize!(model)

# Get solution
x_sol = value.(x)

# Number of gene in the module (clique)
objective_value(model)

# Set of gene in the module (clique)
clique = Set()
for i in 1:n_genes
	if x_sol[i] == 1
		push!(clique, i)
	end
end

# Translate the genes back to their original number
# The results are stored in module_gene
module_gene = Set()
for g_index in clique
	push!(module_gene, b[g_index])
end



### Display the graph and the clique found

#using LightGraphs
#using GraphPlot
#using Colors
#using Cairo
#using Compose

# Create graph
#g = Graph(net_mat)

#membership = ones(Int64, n_genes)
#nodecolor = [colorant"lightseagreen", colorant"orange"]

#for k in clique
#	membership[k] = 2
#end

#nodefillc = nodecolor[membership]

# Save graph in pdf
#draw(PDF(string("graph-", file_name, ".pdf"), 16cm, 16cm), gplot(g, nodefillc=nodefillc))