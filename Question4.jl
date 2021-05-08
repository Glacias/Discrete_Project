using JuMP
using Gurobi
using CSV
using DataFrames
using Bijections
using LinearAlgebra
using LightGraphs

# Function that returns the maximum clique in a graph (represented by its adjacency matrix)
# Inputs :
# - net_mat -> the adjacency matrix of the graph
# Returns :
# - a set of vertices in the maximum clique
function find_maximum_clique(net_mat, set_ignore)
	n_genes = size(net_mat, 1)
	model = Model(Gurobi.Optimizer)
	@variable(model, x[i=1:n_genes], Bin)
	@objective(model, Max, sum(x))
	@constraint(model, [i = 1:n_genes, j = 1:n_genes; i > j], x[i] + x[j] <= 1 + net_mat[i,j])

	# Ignore genes that were already selected
	@constraint(model, [i in set_ignore], x[i] == 0)

	optimize!(model)
	x_sol = value.(x)
	clique = Set{Int64}()
	for i in 1:n_genes
		if x_sol[i] == 1
			push!(clique, i)
		end
	end
	return clique
end

# Function that returns the largest maximal clique in a graph (represented by its adjacency matrix)
# that covers a specified set
# Inputs :
# - net_mat -> the adjacency matrix of the graph
# - set_to_cover -> the set of vertices to cover
# Returns :
# - a set of vertices in the maximum clique
function find_maximal_clique_with_cover(net_mat, set_to_cover)
	n_genes = size(net_mat, 1)
	model = Model(Gurobi.Optimizer)
	@variable(model, x[i=1:n_genes], Bin)
	@objective(model, Max, sum(x))
	@constraint(model, [i = 1:n_genes, j = 1:n_genes; i > j], x[i] + x[j] <= 1 + net_mat[i,j])

	# Add constraint so that the set must be covered
	@constraint(model, [i in set_to_cover], x[i] == 1)

	optimize!(model)
	x_sol = value.(x)
	clique = Set{Int64}()
	for i in 1:n_genes
		if x_sol[i] == 1
			push!(clique, i)
		end
	end
	return clique
end

# Read the data
DATA_dir = joinpath("hsa") # Path to the data
file_name = "GSE10063"
raw_df = CSV.read(joinpath(DATA_dir, string(file_name, ".txt")), DataFrames.DataFrame, header=0)

# Create the set of genes
genes_set = Set{Int64}(vcat(raw_df[:,1], raw_df[:,2]))
n_genes = length(genes_set) # Number of genes

# Create a bijection from gene number to gene index
b = Bijection()
i = 1
for g in genes_set
	b[i] = g
	i+=1
end

# Create the co-expression network matrix
threshold_coexp = 1
net_mat = zeros(Int64, n_genes, n_genes) # Initialize matrix to zeroes

for i in 1:DataFrames.nrow(raw_df)
	if raw_df[i, 3] >= threshold_coexp
		net_mat[b(raw_df[i, 1]), b(raw_df[i, 2])] = 1
		net_mat[b(raw_df[i, 2]), b(raw_df[i, 1])] = 1
	end
end

# Set of the already selected genes
set_ignore = Set([])

# Store the solution
vect_clique = Vector{Set}(undef, 0)

# Find the first maximum clique
c = find_maximum_clique(net_mat, set_ignore)

# Add the clique to the set of genes to ignore
set_ignore = union(set_ignore, c)

# Add it to the solution
push!(vect_clique, c)

while length(set_ignore) < n_genes

	# Find the maximum clique in the graph where we ignore previously selected genes
	c_small = find_maximum_clique(net_mat, set_ignore)

	# Find the largest maximal clique in the initial graph that covers the clique found above
	c_sol = find_maximal_clique_with_cover(net_mat, c_small)

	# Add the previous clique to the set of genes to ignore
	set_ignore = union(set_ignore, c_small)

	# Save the solution
	push!(vect_clique, c_sol)
end

# Translate the genes back to their original number
# The results are stored in vect_mod that contrains ech module
vect_mod = Vector{Set}(undef, 0)

for clique in vect_clique
	module_gene = Set()
	for g_index in clique
		push!(module_gene, b[g_index])
	end
	push!(vect_mod, module_gene)
end