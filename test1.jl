using JuMP
using Gurobi
using CSV
using DataFrames
using Bijections

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
threshold = 1.4
net_mat = zeros(Int8, n_genes, n_genes) # Initialize matrix to zeroes

for i in 1:DataFrames.nrow(raw_df)
	if raw_df[i, 3] >= threshold
		net_mat[b(raw_df[i, 1]), b(raw_df[i, 2])] = 1
		net_mat[b(raw_df[i, 2]), b(raw_df[i, 1])] = 1
	end
end