package jebl.evolution.trees;

import jebl.evolution.distances.DistanceMatrix;
import jebl.evolution.taxa.Taxon;

/**
 * A meeting point for tree building from sequence data. A very initial form which will develope to encompass more
 * methods and distances. Currently only pairwise distance methods are implemented.
 *
 * @author Joseph Heled
 * @version $Id: TreeBuilderFactory.java 853 2007-12-07 07:32:42Z twobeers $
 *
 */

public class TreeBuilderFactory {

    /**
     * Supported methods for tree building
     */
    public static enum Method { NEIGHBOR_JOINING("Neighbor-Joining"), UPGMA("UPGMA");
        Method(String name) { this.name = name; }
        public String toString() { return getName(); }
        public String getName() { return name; }
        private String name;
    }

    /**
     * Supported pairwise distance methods
     */
    public static enum DistanceModel { JukesCantor, F84, HKY, TamuraNei }

    /**
     * Supported consensus methods.
     */
    public static enum ConsensusMethod { GREEDY, MRCAC }

    /**
     *
     * @param method to check
     * @return Wheather method generates a rooted or unrooted tree.
     */
    public static boolean isRootedMethod(Method method) {
        switch( method ) {
            case UPGMA:
            {
                return true;
            }
            case NEIGHBOR_JOINING:
            default:
            {
                return false;
            }
        }
    }

    // TT: This method should probably have been called createBuilder() because it
    // creates a new builder every time instead of reusing an existing one. This
    // makes a difference because it means that progress listeners added to a
    // builder don't need to be removed afterwards if the builder is discarded.

    /**
     *
     * @param method build method to use.
     * @param distances Pre computed pairwise distances.
     * @return A tree builder using method and distance matrix
     */
    static public ClusteringTreeBuilder getBuilder(Method method, DistanceMatrix distances) {
        ClusteringTreeBuilder builder;
        switch( method ) {
            case UPGMA:
            {
                builder = new UPGMATreeBuilder(distances);
                break;
            }
            case NEIGHBOR_JOINING:
            default:
            {
                builder = new NeighborJoiningTreeBuilder(distances);
                break;
            }
        }
        return builder;
    }

    // trees must have the same taxa
    static public ConsensusTreeBuilder buildUnRooted(Tree[] trees, Taxon outGroup, double supportThreshold, ConsensusMethod method) {
        if( ! (supportThreshold >= 0 && supportThreshold <= 1) ) {
             throw new IllegalArgumentException("support not in [0..1]: " + supportThreshold);
        }
        switch( method ) {
            case GREEDY: {
                return new GreedyUnrootedConsensusTreeBuilder(trees, outGroup, supportThreshold);
            }
        }
        // bug
        throw new IllegalArgumentException(method.toString());
    }

    // trees must have the same taxa
    static public ConsensusTreeBuilder buildRooted(RootedTree[] trees, double supportThreshold, ConsensusMethod method) {
        if( ! (supportThreshold >= 0 && supportThreshold <= 1) ) {
             throw new IllegalArgumentException("Expected support value in [0..1], got " + String.format("%.3f", supportThreshold));
        }

        switch (method) {
            case GREEDY:
                return new GreedyRootedConsensusTreeBuilder(trees, supportThreshold);
            case MRCAC:
                return new MRCACConsensusTreeBuilder(trees, supportThreshold);
            default:
                throw new IllegalArgumentException("Unknown consensus method: " + method);
        }
    }

    /**
     * convenience method. Convert arrays of trees, guaranteed to be rooted to the array of the appropriate
     * type.
     * @param trees trees - all must be rooted
     * @param supportThreshold minimum required consensus support (in [01])
     * @param method which consensus method to use
     * @return consensus tree builder
     */
    static public ConsensusTreeBuilder buildRooted(Tree[] trees, double supportThreshold, ConsensusMethod method) {
        RootedTree[] rtrees = new RootedTree[trees.length];
        for(int i = 0; i < trees.length; ++i) {
            rtrees[i] = (RootedTree)trees[i];
        }
        return buildRooted(rtrees, supportThreshold, method);
    }
}
