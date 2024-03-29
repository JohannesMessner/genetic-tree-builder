#########################################################
#######        COMPUTATIONAL BIOLOGY         ############
#######             HOMEWORK 3               ############
#########################################################
#                                                       #
# Reconstruct the phylogenetic tree of given sequences  #
# using UPGMA with hamming, JC69, K80 distances.        #
#                                                       #
#########################################################
#########################################################

#########################################################
######    Code below this should not be changed   #######
#########################################################

library(ape)

transform_to_phylo = function(sequences, named_edges, edge_lengths, node_description) {
    # Produce a tree of the phylo class from the matrix of edges, vector of lengths and
    # the dataframe  containing node descriptions.
    #    sequences: a list of the original sequences;
    #    named_edges: an Mx2 matrix of pairs of nodes connected by an edge,
    #                 where the M rows are the different edges and the 2 columns
    #                 are the parent node and the child node of an edge.
    #    edge_lengths: a vector of length M of the corresponding edge lengths.
    #    node_description: a data frame of the node descriptions as defined in
    #                      initialize_node_description and extended by the upgma code.
    
    edges = named_edges
    for (name in rownames(node_description)) {
        index = which(rownames(node_description) == name)
        edges[which(edges == name)] = as.numeric(index)
    }
    edges = matrix(as.numeric(edges), ncol = 2)
    
    edges[which(edges > length(sequences))] = - edges[which(edges > length(sequences))]
    root = setdiff(edges[,1], edges[,2])
    edges[which(edges==root)] = length(sequences) + 1
    
    k = length(sequences) + 2
    for(x in unique(edges[which(edges < 0)])) {
        edges[which(edges==x)] = k
        k = k + 1
    }
    
    tree = list()
    class(tree) = "phylo"
    tree$edge = edges
    tree$edge.length = edge_lengths
    tree$Nnode = as.integer(length(sequences) - 1)
    tree$tip.label = names(sequences)
    
    # Return the tree in the form of the phylo class from ape
    return(tree)
}

plot_tree = function(tree) {
    # Plot the phylogenetic tree with node labels and edge lengths.
    #    tree: an object of class phylo from the ape package

    plot(tree)
    edgelabels(format(tree$edge.length, digits = 2))
}

initialize_node_description = function(sequences) {
    # Initialize the structure that will hold node descriptions.
    # The created structure is a data frame where the rows are node names, and the columns are
    # node_height -- distance from the node to any tip in the ultrametric tree, and
    # node_size -- number of tips that this node is ancestral to.
    #    sequences: a list of the original sequences;

    N = length(sequences)
    node_names = names(sequences)
    node_sizes = rep(1, times = N)
    node_heights = rep(0, times = N)
    node_description = data.frame(node_sizes, node_heights)
    rownames(node_description) = node_names

    # Return a data frame that contains information on currently existing tip nodes.
    # node_description: a dataframe containing information on the currently existing nodes.
    #                   The row names are the names of the currently existing tips, i.e.
    #                   are the same as the names in the sequence list, node_height is
    #                   0 and node_size is 1 as the all the currently existing nodes are tips.
    return(node_description)
}

add_new_node = function(node_description, merging_nodes) {
    # Add new merged node to the node description data frame.
    # The new node is a combination of the nodes supplied in the merging_nodes,
    # e.g. if one needs to merge nodes "bird" and "fish", the new node in the
    # dataframe will be called "bird.fish".
    #    node_description: the dataframe created by initialize_node_description, containing
    #                      current node sizes and node heights
    #    merging_nodes: a vector of two names of the nodes being merged

    new_node_name = paste(merging_nodes, collapse = ".")
    new_node_row = data.frame(node_sizes = 0, node_heights = 0)
    rownames(new_node_row) = new_node_name
    new_node_description = rbind(node_description, new_node_row)
    
    # Return the node_description dataframe with a row for the new node added, and
    # the new node name.
    #    node_description: the dataframe where the rows are labelled by current nodes and columns
    #                      contain the node heights and sizes.
    #    new_node_name: the name of the newly added node, created from names in merging_nodes.
    return(list(node_description = new_node_description, new_node_name = new_node_name))
}

#########################################################
######    Code above this should not be changed   #######
#########################################################

get_hamming_distance = function(sequence1, sequence2) {
    # Compute the Hamming distance between two sequences.
    #    sequence1: first sequence 
    #    sequence2: second sequence

    seq1_split <- strsplit(sequence1, "")[[1]]
    seq2_split <- strsplit(sequence2, "")[[1]]
    distance <- 0
    # Both sequences are alignments -> they are assumed to have the same length
    for (i in (1:length(seq1_split))) {
        if (seq1_split[i] != seq2_split[i]) {
          distance <- distance+1
        }
    }

    # Return the numerical value of the distance
    return(distance)
}

get_JC69_distance = function(sequence1, sequence2) {
    # Compute the JC69 distance between two sequences.
    #    sequence1: first sequence
    #    sequence2: second sequence

    # Both sequences are alignments -> they are assumed to have the same length
    seq1_split <- strsplit(sequence1, "")[[1]]
    
    prop <- (get_hamming_distance(sequence1 = sequence1, sequence2 = sequence2))/length(seq1_split)
    distance <- (-3/4)*log(1-(4/3)*prop)

    # Return the numerical value of the distance
    return(distance)
}

get_K80_distance = function(sequence1, sequence2) {
    # Compute the K80 distance between two sequences.
    #    sequence1: first sequence
    #    sequence2: second sequence
    seq1_split <- strsplit(sequence1, "")[[1]]
    seq2_split <- strsplit(sequence2, "")[[1]]

    transit_diffs <- 0
    transvers_diffs <- 0
    for (i in (1:length(seq1_split))) {
        base1 <- seq1_split[i]
        base2 <- seq2_split[i]
        transit_pair1 <- c("C","T")
        transit_pair2 <- c("A", "G")
        if (base1 != base2) {
            if ((base1 %in% transit_pair1 && base2 %in% transit_pair1)
                || (base1 %in% transit_pair2 && base2 %in% transit_pair2)) {
                transit_diffs <- transit_diffs + 1
            } else {
                transvers_diffs <- transvers_diffs +1
            }
        }
    }
    
    S <- transit_diffs/length(seq1_split)
    V <- transvers_diffs/length(seq2_split)

    distance <- -(1/2)*log(1-(2*S)-V)-(1/4)*log(1-2*V)

    # Return the numerical value of the distance
    return(distance)
}

compute_initial_distance_matrix = function(sequences, distance_measure) {
    # Compute the initial distance matrix using one of the distance measures.
    # The matrix is of dimension NxN, where N is the number of sequences.
    # The matrix columns and rows should be labelled with tip names, each row and column
    # corresponding to the appropriate sequence.
    # The matrix can be filled completely (i.e symmetric matrix) or only the upper half
    # (as shown in the lecture).
    # The diagonal elements of the matrix should be Inf.
    #    sequences: the sequence alignment in the form of a list of species names and 
    #               the associated genetic sequences
    #    distance_measure: a string indicating whether the 'hamming', 'JC69' or 'K80' 
    #                      distance measure should be used
  
    N <- length(sequences)
    distance_matrix <- matrix(nrow = N, ncol = N)
    diag(distance_matrix) <- Inf
    rownames(distance_matrix) <- names(sequences)
    colnames(distance_matrix) <- names(sequences)

    # determine distance measure
    if (distance_measure == "hamming") {
        distance_func <- get_hamming_distance
    } else if (distance_measure == "JC69") {
        distance_func <- get_JC69_distance
    } else {
        distance_func <- get_K80_distance
    }
    
    # fill matrix
    for (i in 1:N) {
        for (j in 1:N) {
            if ( i != j) {
                seq1 <- sequences[[rownames(distance_matrix)[i]]]
                seq2 <- sequences[[colnames(distance_matrix)[j]]]
                distance_matrix[i,j] <- distance_func(seq1, seq2)
            }
        }
    }

    # Return the NxN matrix of inter-sequence distances with Inf on the diagonal
    return(distance_matrix)
}

get_merge_node_distance = function(node_description, distance_matrix, merging_nodes, existing_node) {
    # Compute the new distance between the newly created merge node and an existing old node in the tree
    #    node_description: a dataframe containing information on the currently existing nodes
    #    distance_matrix: the matrix of current distances between nodes
    #    merging_nodes: a vector of two node names that are being merged in this step
    #    existing_node: one of the previously existing nodes, not included in the new node
    
    n_i <- node_description[[merging_nodes[1],"node_sizes"]]
    n_j <- node_description[[merging_nodes[2],"node_sizes"]]
    
    d_i <- distance_matrix[merging_nodes[1], existing_node]
    d_j <- distance_matrix[merging_nodes[2], existing_node]
    
    new_distance <- (n_i*d_i + n_j*d_j)/(n_i + n_j)
  
    # Returns the distance between the newly created merge node and the existing node
    return(new_distance)
}

update_distance_matrix = function(node_description, distance_matrix, merging_nodes, new_node_name) {
    # Update the distance matrix given that two nodes are being merged.
    #    node_description: a dataframe containing information on the currently existing nodes
    #    distance_matrix: the current distance matrix that needs to be updated
    #    merging_nodes: a vector of two node names that need to be merged in this step
    #    new_node_name: the name with which the merged node should be labelled
    # The resulting matrix should be one column and one row smaller, i.e. if the given distance matrix
    # was MxM, then the updated matrix will be M-1xM-1, where the 2 rows and cols represent the separate
    # nodes undergoing the merge are taken out and a new row and col added that represents the new node.

    # create updated matrix
    # remove entries containing merging (old) nodes
    updated_distance_matrix <- distance_matrix[! rownames(distance_matrix) %in% merging_nodes,! colnames(distance_matrix) %in% merging_nodes]
    if (! is.matrix(updated_distance_matrix)) {
        updated_distance_matrix <- matrix(updated_distance_matrix)
    }
    
    # add row and column for merged (new) node
    previous_names <- rownames(updated_distance_matrix)
    if (is.null(previous_names)) {
        previous_names <- rownames(distance_matrix)[! rownames(distance_matrix) %in% merging_nodes]
    }
    updated_distance_matrix <- rbind(updated_distance_matrix, rep(Inf, ncol(updated_distance_matrix)))
    updated_distance_matrix <- cbind(updated_distance_matrix, rep(Inf, nrow(updated_distance_matrix)))
    names <- c(previous_names, new_node_name)
    rownames(updated_distance_matrix) <- names
    colnames(updated_distance_matrix) <- names

    # compute new distances to merged (new) node
    for (i in 1:ncol(distance_matrix)) {
        if (rownames(updated_distance_matrix)[i] == new_node_name) {
            break
        }
        new_dist <- get_merge_node_distance(node_description = node_description,
                                            distance_matrix = distance_matrix,
                                            merging_nodes = merging_nodes,
                                            existing_node = colnames(updated_distance_matrix)[i])
        
        updated_distance_matrix[nrow(updated_distance_matrix), i] <- new_dist
        updated_distance_matrix[i, ncol(updated_distance_matrix)] <- new_dist
    }

    # Returns the updated matrix of cluster distances
    return(updated_distance_matrix)
}

upgma_one_step = function(node_description, distance_matrix, edges, edge_lengths) {
    # Performs one step of the upgma algorithm, i.e. the nodes with the smallest distance are merged, 
    # the node height of the new node is calculated and the distance matrix is newly calculated.
    # Values that are expected to be returned are listed below.
    #    node_description: a dataframe containing information on the currently existing nodes
    #    distance_matrix: the current distance matrix that needs to be updated (OxO)
    #    edges: an Mx2 matrix of pairs of nodes connected by an edge, where the M rows are the different
    #           edges and the 2 columns are the parent node and the child node of an edge.
    #    edge_lengths: a vector of length M of the corresponding edge lengths.
  
    # determine which nodes will be merged (the ones with minimal distance)
    merging_nodes_indices <- which(distance_matrix == min(distance_matrix), arr.ind=T)
    merging_nodes_indices <- c(merging_nodes_indices[1,1], merging_nodes_indices[1,2])
    merging_nodes <- c(rownames(distance_matrix)[merging_nodes_indices[1]], rownames(distance_matrix)[merging_nodes_indices[2]])
    
    # obtain information regarding these merging nodes
    d_to_tip_1 <- node_description[merging_nodes[1], "node_heights"]
    d_to_tip_2 <- node_description[merging_nodes[2], "node_heights"]
    size_1 <- node_description[merging_nodes[1], "node_sizes"]
    size_2 <- node_description[merging_nodes[2], "node_sizes"]
    
    # calculate lenghts for newly merged node
    d_between_merging_nodes <- distance_matrix[merging_nodes[1], merging_nodes[2]]
    new_length_1 <- d_between_merging_nodes/2 - d_to_tip_1
    new_length_2 <- d_between_merging_nodes/2 - d_to_tip_2
    
    # update information about edge lengths with newly calculated lengths
    edge_lengths <- c(edge_lengths, new_length_1, new_length_2)
    
    # add new merged node to existing nodes
    desc_update <- add_new_node(node_description, merging_nodes)
    node_description <- desc_update[["node_description"]]
    new_node_name <- desc_update[["new_node_name"]]
    # calculate properties of new entry for the newly merged node
    node_description[new_node_name, "node_sizes"] <- size_1 + size_2
    node_description[new_node_name, "node_heights"] <- new_length_1 + d_to_tip_1
    # add new edges
    edges <- rbind(edges, c(new_node_name, merging_nodes[1]))
    edges <- rbind(edges, c(new_node_name, merging_nodes[2]))
    
    #update distance metrix with the new merged node
    distance_matrix <- update_distance_matrix(node_description = node_description,
                                             distance_matrix = distance_matrix,
                                             merging_nodes = merging_nodes,
                                             new_node_name = new_node_name)
                                    
    # Return the updated distance matrix, edge description matrix, edge length vector and 
    # node_description data frame
    #    node_description: data frame containing sizes and heights of all nodes 
    #        (to be updated using add_new_node())
    #    distance_matrix: the updated matrix of distances between nodes (O-1xO-1)
    #    edges: an Mx2 matrix of pairs of nodes connected by an edge, where the M rows are
    #    the different edges and the 2 columns are the parent node and the child node of an edge.
    #    edge_lengths: a vector of length M of the corresponding edge lengths.
    return(list(node_description = node_description, distance_matrix = distance_matrix,
                edges = edges, edge_lengths = edge_lengths))
}

build_upgma_tree = function(sequences, distance_measure) {
    # Build the tree from given sequences using the UPGMA algorithm.
    #    sequences: the sequences in the format of a list of species names and the associated genetic 
    #               sequences
    #    distance_measure: a string indicating whether the 'hamming', 'JC69' or 'K80' distance measure
    #                      should be used
    
    # initialize empty data structures
    N <- length(sequences)
    node_description <- initialize_node_description(sequences)
    edges <- matrix(nrow = 0, ncol = 2)
    edge_lengths <- vector(mode = "numeric", length = 0)
    # intial distance matrix
    distance_matrix <- compute_initial_distance_matrix(sequences, distance_measure)
    
    # perform UPGMA-steps (node-merging and information update) until root is formed
    while(nrow(distance_matrix) > 1 && ncol(distance_matrix) > 1) {
        updated_values <- upgma_one_step(node_description, distance_matrix, edges, edge_lengths)
        node_description <- updated_values[["node_description"]]
        distance_matrix <- updated_values[["distance_matrix"]]
        edges <- updated_values[["edges"]]
        edge_lengths <- updated_values[["edge_lengths"]]
    }
  
    # Return the UPGMA tree of sequences
    tree <- transform_to_phylo(sequences, edges, edge_lengths, node_description)
    return(tree)
}

test_tree_building = function() {
    sequences <- list(orangutan = "TCACACCTAAGTTATATCTATATATAATCCGGGCCGGG",
                     chimp = "ACAGACTTAAAATATACATATATATAATCCGGGCCGGG",
                     human = "AAAGACTTAAAATATATATATATATAATCCGGGCCGGG",
                     gorilla = "ACACACCCAAAATATATATATATATAATCCGGGCCGGG",
                     unicorn = "ACACACCCAAAATATATACGCGTATAATCCGGGCCGAA")
    distance_measure <- 'hamming'
    distance_measure <- 'JC69'
    distance_measure <- 'K80'

    tree <- build_upgma_tree(sequences, distance_measure)
    plot_tree(tree)
}

test_tree_building()