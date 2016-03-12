/*
 * Tree.java
 *
 * (c) 2005 JEBL Development Team
 *
 * This package is distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package jebl.evolution.trees;

import java.util.Set;

/**
 * A rooted or unrooted tree. This interface is the common base class for all trees,
 * and contains only operations for unrooted trees. The subinterface RootedTree
 * contains additional methods that make sense only on rooted trees.
 *
 * Both interfaces contain no mutator methods. As of 2006-12-08, the only way
 * to mutate a tree after it has been built is to use its concrete class
 * instead of the Tree or RootedTree interface.
 *
 * @author rambaut
 * @author Alexei Drummond
 *
 * @version $Id: Tree.java 627 2007-01-15 03:50:40Z pepster $
 */
public interface Tree extends jebl.evolution.graphs.Graph {

    /**
     * @return a set of all nodes that have degree 1.
     * These nodes are often refered to as 'tips'.
     */
    Set<jebl.evolution.graphs.Node> getExternalNodes();

    /**
     * @return a set of all nodes that have degree 2 or more.
     * These nodes are often refered to as internal nodes.
     */
    Set<jebl.evolution.graphs.Node> getInternalNodes();

	/**
	 * @return a set of all edges that have a degree 1 node.
	 */
	Set<jebl.evolution.graphs.Edge> getExternalEdges();

	/**
	 * @return a set of all edges for which both nodes have degree 2 or more.
	 */
	Set<jebl.evolution.graphs.Edge> getInternalEdges();

    /**
     * @return the set of taxa associated with the external
     * nodes of this tree. The size of this set should be the
     * same as the size of the external nodes set.
     */
    Set<jebl.evolution.taxa.Taxon> getTaxa();

    /**
     * @param node the node whose associated taxon is being requested.
     * @return the taxon object associated with the given node, or null
     * if the node is an internal node.
     */
    jebl.evolution.taxa.Taxon getTaxon(jebl.evolution.graphs.Node node);

    /**
     * @param node the node
     * @return true if the node is of degree 1.
     */
    boolean isExternal(jebl.evolution.graphs.Node node);

    /**
     * @param taxon the taxon
     * @return the external node associated with the given taxon, or null
     * if the taxon is not a member of the taxa set associated with this tree.
     */
    jebl.evolution.graphs.Node getNode(jebl.evolution.taxa.Taxon taxon);

    void renameTaxa(jebl.evolution.taxa.Taxon from, jebl.evolution.taxa.Taxon to);
}