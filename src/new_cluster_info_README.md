# NLCE-Rust
Numerical Linked Cluster Expansion code written in Rust.

# How to interpret the data
The output of this code is formatted into three json files per order. Each file is setup to be imported into a HashMap similar to Python's Dictionary object. The keys to the dictionary are the unique identification number for each graph. This identification number is unique to every graph that is isomorphic to the graph given. Meaning, two graphs that are topologically identical will have the same identification number. 

The first file is the graph bond information. As with the other files, the keys here are the unique identification numbers for each graph. In this file, each identification number is linked to the graph bond information that corresponds to the cluster.

Take, for example, the following output from the NLCE data for a square lattice with clusters of order 4. In this example, there are three isomorphically distinct clusters (given by the three unique ids) and they each have their own set of bond information. The first entry in the first graph (graphID: 2664808976974672910, entry: [0, 1, 1]) has bond information that tells us that site 0 is connected to site 1 through bond type 1 (hence [0,1,1]). The other two entries tell us that site 1 is connected to site 2 and site 2 to site 3. Overall, this tells us that the first cluster is a straight line of four connected sites. Similarly, we can look at the other two to see what possible isomorphically distinct clusters exist in this type of geometry. 
{
    "2664808976974672910" : [[0,1,1],[1,2,1],[2,3,1]],
    "6110569970769299842" : [[0,2,1],[0,1,1],[1,3,1],[2,3,1]],
    "850428628311615989" : [[0,2,1],[1,2,1],[2,3,1]]
}

The second file is the graph multiplicity information. In this file, each identification number is linked to the graph's total multiplicity. This is the multiplicity of each of the underlying symmetrically distinct clusters summed up. 

For example, the below data is taken from the NLCE data for a square lattice with clusters of order 4. From the previous example about bond information, we know that the cluster with graphID 6110569970769299842 is a square with four sites. (Take a second to show this to yourself based on the bond information in the previous example) Based on this, we know that the total graph multiplicity of this cluster, the number of ways it and symmetrically distinct versions of itself can be oriented on the underlying lattice, should be one. There is only one way that you can orient this square that is unique. This is reflected in the example below where the graphID 6110569970769299842 has a graph multiplicity of one. 
{
    "6110569970769299842" : 1,
    "2664808976974672910" : 14,
    "850428628311615989" : 4
}

The last file is the subgraph multiplicity information. In this file, each identification number is linked to another hashmap that tells us the clusters that make up the original cluster, and their corresponding multiplicites.

For example, the below data is again taken from the NLCE data for a square lattice with clusters of order 4. Let us look at the square cluster again. The graphID for the square cluster is 6110569970769299842. Immediately, we can see that the square cluster is made up of three different subclusters. The first subcluster (4288071131319860613) can be seen in the square lattice clusters of order 2. Specifically, its the two site connection with bond type 1. (There should only be one graph in the graph_bond_info for order 2.)
Thus, this example file tells us that within the square cluster of order 4, there are 4 subclusters of order 2 with the shape given in the graph bond information from order 2.  
{
    "850428628311615989" : { 
                                "15130871412783076140" : 4,
                                "4288071131319860613" : 3,
                                "15478155942260055368" : 3
                           },
    "2664808976974672910" : {
                                "15478155942260055368" : 2,
                                "15130871412783076140" : 4,
                                "4288071131319860613" : 3
                            },
    "6110569970769299842" : {   
                                "4288071131319860613" : 4,
                                "15478155942260055368" : 4,
                                "15130871412783076140" : 4
                            }
}

