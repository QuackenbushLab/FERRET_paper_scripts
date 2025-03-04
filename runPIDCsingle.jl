using NetworkInference
using LightGraphs
using DelimitedFiles

dataset_name = string("/home/ubuntu/TCGA_single_cell_splits_pca/",ARGS[1])
algorithm = PIDCNetworkInference()
threshold = parse(Float64, ARGS[2])
@time genes = get_nodes(dataset_name);
@time network = InferredNetwork(algorithm, genes);
adjacency_matrix, labels_to_ids, ids_to_labels = get_adjacency_matrix(network, threshold)
graph = LightGraphs.SimpleGraphs.SimpleGraph(adjacency_matrix)
adjacency_list = LightGraphs.SimpleGraphs.SimpleGraphs.adj(graph)
writedlm(string("/home/ubuntu/PIDC_CPTAC/",ARGS[2], "_", ARGS[1]),adjacency_matrix, ',')
