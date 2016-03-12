/*
 * Node.java
 *
 * (c) 2005 JEBL Development Team
 *
 * This package is distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package jebl.evolution.graphs;

import jebl.util.Attributable;

/**
 * Represents a node in a graph or tree. In general it is
 * used only as a handle to traverse a graph or tree structure and
 * it has no methods or instance variables.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: Node.java 904 2008-05-04 23:16:09Z matt_kearse $
 */
public interface Node extends Attributable {

    /**
     * Get the number of edges connected to this node. (this is also the number of nodes connected to this node).
     * <p/>
     * Note: In previous implementations, this method would return a different value on the root node of trees.
     * It would previously return 1 more than the correct value in that case.
     * @return the number of edges connected to this node
     */
    int getDegree();
}
