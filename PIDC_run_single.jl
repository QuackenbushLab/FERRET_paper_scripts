using NetworkInference
using LightGraphs
using DelimitedFiles

dataset_name = ARGS[1]
algorithm = PIDCNetworkInference()
threshold = parse(Float64, ARGS[2])
println(dataset_name)
@time genes = get_nodes(dataset_name);
@time network = InferredNetwork(algorithm, genes);
adjacency_matrix, labels_to_ids, ids_to_labels = get_adjacency_matrix(network, threshold)
graph = LightGraphs.SimpleGraphs.SimpleGraph(adjacency_matrix)
adjacency_list = LightGraphs.SimpleGraphs.SimpleGraphs.adj(graph)
println(labels_to_ids)
println(ids_to_labels)
writedlm(ARGS[3],adjacency_matrix, ',')
