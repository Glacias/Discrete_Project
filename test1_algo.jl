using CSV
using DataFrames
using Bijections
using LinearAlgebra
using LightGraphs

# Function that returns :
# - true if a vertex is connected to all the members in the set provided
# - false otherwise
function is_related_to_entire_set(vertex, set, g)
	for v in set
		if !(has_edge(g, vertex, v))
			return false
		end
	end
	return true
end

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
net_mat = zeros(Int64, n_genes, n_genes) # Initialize matrix to zeroes

for i in 1:DataFrames.nrow(raw_df)
	if raw_df[i, 3] >= threshold
		net_mat[b(raw_df[i, 1]), b(raw_df[i, 2])] = 1
		net_mat[b(raw_df[i, 2]), b(raw_df[i, 1])] = 1
	end
end

# Create graph
g = Graph(net_mat)

### Greedy algo (sorted by degree)

# Sort gene by degree
degrees = vec(sum(net_mat, dims=1))
sorted_vertex = sortperm(degrees, rev=true)

# Start from the empty set
clique = Set()

# Try to add vertex if set remains a clique
for gene in sorted_vertex
	if is_related_to_entire_set(gene, clique, g)
		push!(clique, gene)
	end
end



### Display the graph and the clique found

using GraphPlot
using Colors
using Cairo
using Compose

membership = ones(Int64, n_genes)
nodecolor = [colorant"lightseagreen", colorant"orange"]

for k in clique
	membership[k] = 2
end

nodefillc = nodecolor[membership]

# Save graph in pdf
draw(PDF(string("graph-", file_name, "-greedy.pdf"), 16cm, 16cm), gplot(g, nodefillc=nodefillc))