package jebl.evolution.trees;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jebl.evolution.graphs.Graph;
import jebl.evolution.graphs.Node;

/**
 * A basic implementation on an unrooted tree.
 *
 * @author Joseph Heled
 * @version $Id: SimpleTree.java 935 2008-07-22 16:52:04Z rambaut $
 *
 */

public final class SimpleTree implements Tree {

    /**
     * Tree (to be constructed by subsequent calls).
     */
    public SimpleTree() {}

    /**
     *  Duplicate a tree.
     *
     * @param tree
     */
    public SimpleTree(Tree tree) {
        try {
            createTree(tree, tree.getExternalNodes().iterator().next(), null);
        } catch (Graph.NoEdgeException e) {
            throw new IllegalArgumentException("BUG: invalid tree");
        }
    }

    /**
     * Creates a new external node with the given taxon. See createInternalNode
     * for a description of how to use these methods.
     * @param taxon the taxon associated with this node
     * @return the created node reference
     */
    public jebl.evolution.graphs.Node createExternalNode(jebl.evolution.taxa.Taxon taxon) {
        SimpleNode node = new SimpleNode(taxon);
        externalNodes.put(taxon, node);
        return node;
    }

    /**
     * Once a SimpleTree has been created, the node stucture can be created by
     * calling createExternalNode and createInternalNode. First of all createExternalNode
     * is called giving Taxon objects for the external nodes. Then these are put into
     * sets and passed to createInternalNode to create new internal nodes.
     *
     * It is the caller responsibility to insure no cycles are created.
     *
     * @param adjacencies the child nodes of this node
     * @return the created node.
     */
    public jebl.evolution.graphs.Node createInternalNode(List<jebl.evolution.graphs.Node> adjacencies) {
        SimpleNode node = new SimpleNode(adjacencies);

        internalNodes.add(node);

        for( jebl.evolution.graphs.Node c : adjacencies ) {
            ((SimpleNode)c).addAdjacency(node);
        }
        return node;
    }

    /**
     * Set edge distance between two adjacent nodes.
     * @param node1
     * @param node2
     * @param length
     */
    public void setEdgeLength(final jebl.evolution.graphs.Node node1, final jebl.evolution.graphs.Node node2, final double length) {
        assert getAdjacencies(node1).contains(node2) && getAdjacencies(node2).contains(node1) && length >= 0;

        final jebl.evolution.graphs.Edge edge = new SimpleEdge(node1, node2, length);

        edges.put(new HashPair<Node>(node1, node2), edge);
        edges.put(new HashPair<Node>(node2, node1), edge);
    }

    /**
     * Change length of an existing edge.
     * @param edge
     * @param length
     */
    public void setEdgeLength(final jebl.evolution.graphs.Edge edge, final double length) {
       ((SimpleEdge)edge).length = length;
    }

    /**
     * Add a new edge between two existing (non adjacent yet)  nodes.
     * @param node1
     * @param node2
     * @param length
     */
    public void addEdge(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2, double length) {
        assert !getAdjacencies(node1).contains(node2);

        ((SimpleNode)node1).addAdjacency(node2);
        ((SimpleNode)node2).addAdjacency(node1);
        setEdgeLength(node1, node2, length);
    }

    /**
     * Copy partition of source starting from node but skipping root.
     *
     * @param source to copy from
     * @param node  in source to copy
     * @param root  adjacent node already copied
     * @return copy of node inside tree
     * @throws Graph.NoEdgeException
     */
    private jebl.evolution.graphs.Node createTree(Tree source, jebl.evolution.graphs.Node node, jebl.evolution.graphs.Node root) throws Graph.NoEdgeException {
        jebl.evolution.graphs.Node h;
        if( source.isExternal(node) ) {
            h = createExternalNode(source.getTaxon(node));
        } else {
            h = createInternalNode(new ArrayList<jebl.evolution.graphs.Node>() );
        }

        final List<jebl.evolution.graphs.Node> adjacencies = source.getAdjacencies(node);

        for( jebl.evolution.graphs.Node c : adjacencies ) {
            if( c == root ) continue;
            final jebl.evolution.graphs.Node n = createTree(source, c, node);
            addEdge(h, n, source.getEdgeLength(node, c) );
        }
        return h;
    }

    /* Graph IMPLEMENTATION */

    /**
     * Returns a list of edges connected to this node
     *
     * @param node
     * @return the set of nodes that are attached by edges to the given node.
     */
    public List<jebl.evolution.graphs.Edge> getEdges(jebl.evolution.graphs.Node node) {
        //return null;
        List<jebl.evolution.graphs.Node> adjacencies = getAdjacencies(node);
        List<jebl.evolution.graphs.Edge> edges = new ArrayList<jebl.evolution.graphs.Edge>();
        for(jebl.evolution.graphs.Node adjNode : adjacencies){
            try{
                edges.add(getEdge(node,adjNode));
            }
            catch(Graph.NoEdgeException ex){/*do nothing*/}
        }
        return edges;
    }

    /**
     * @param node
     * @return the set of nodes that are attached by edges to the given node.
     */
    public List<jebl.evolution.graphs.Node> getAdjacencies(jebl.evolution.graphs.Node node) {
        return ((SimpleNode)node).getAdjacencies();
    }

    /**
     * Returns the Edge that connects these two nodes
     *
     * @param node1
     * @param node2
     * @return the edge object.
     * @throws Graph.NoEdgeException
     *          if the nodes are not directly connected by an edge.
     */
    public jebl.evolution.graphs.Edge getEdge(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) throws Graph.NoEdgeException {
        jebl.evolution.graphs.Edge edge = edges.get(new HashPair<Node>(node1, node2));
        if( edge == null ) {
            // not connected
            throw new Graph.NoEdgeException();
        }
        return edge;
    }

    /**
     * @return a set of all nodes that have degree 1.
     *         These nodes are often refered to as 'tips'.
     */
    public Set<jebl.evolution.graphs.Node> getExternalNodes() {
        return new LinkedHashSet<jebl.evolution.graphs.Node>(externalNodes.values());
    }

    /**
     * @return a set of all nodes that have degree 2 or more.
     *         These nodes are often refered to as internal nodes.
     */
    public Set<jebl.evolution.graphs.Node> getInternalNodes() {
        return new LinkedHashSet<jebl.evolution.graphs.Node>(internalNodes);
    }

	/**
	 * @return the set of taxa associated with the external
	 *         nodes of this tree. The size of this set should be the
	 *         same as the size of the external nodes set.
	 */
	public Set<jebl.evolution.taxa.Taxon> getTaxa() {
	    return new LinkedHashSet<jebl.evolution.taxa.Taxon>(externalNodes.keySet());
	}
    /**
     * @param node the node whose associated taxon is being requested.
     * @return the taxon object associated with the given node, or null
     *         if the node is an internal node.
     */
    public jebl.evolution.taxa.Taxon getTaxon(jebl.evolution.graphs.Node node) {
        return ((SimpleNode)node).getTaxon();
    }

    /**
     * @param node the node
     * @return true if the node is of degree 1.
     */
    public boolean isExternal(jebl.evolution.graphs.Node node) {
        return node.getDegree() == 1;
    }

	/**
	 * @param edge the edge
	 * @return true if the edge has a node of degree 1.
	 */
	public boolean isExternal(jebl.evolution.graphs.Edge edge) {
	    return ((SimpleEdge)edge).isExternal();
	}

    /**
     * @param taxon the taxon
     * @return the external node associated with the given taxon, or null
     *         if the taxon is not a member of the taxa set associated with this tree.
     */
    public jebl.evolution.graphs.Node getNode(jebl.evolution.taxa.Taxon taxon) {
        return externalNodes.get(taxon);
    }

    public void renameTaxa(jebl.evolution.taxa.Taxon from, jebl.evolution.taxa.Taxon to) {
        SimpleNode node = (SimpleNode)externalNodes.get(from);
        node.setTaxa(to);

        externalNodes.remove(from);
        externalNodes.put(to, node);
    }

    /**
     * @param node1
     * @param node2
     * @return the length of the edge connecting node1 and node2.
     * @throws Graph.NoEdgeException if the nodes are not directly connected by an edge.
     */
    public double getEdgeLength(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) throws Graph.NoEdgeException {
        return getEdge(node1, node2).getLength();
    }

	/**
	 * Returns an array of 2 nodes which are the nodes at either end of the edge.
	 *
	 * @param edge
	 * @return an array of 2 edges
	 */
	public jebl.evolution.graphs.Node[] getNodes(jebl.evolution.graphs.Edge edge) {
		return new jebl.evolution.graphs.Node[] { ((SimpleEdge)edge).getNode1(), ((SimpleEdge)edge).getNode2() };
	}

	/**
	 * @return the set of all nodes in this graph.
	 */
	public Set<jebl.evolution.graphs.Node> getNodes() {
	    Set<jebl.evolution.graphs.Node> nodes = new LinkedHashSet<jebl.evolution.graphs.Node>(internalNodes);
	    nodes.addAll(externalNodes.values());
	    return nodes;
	}

    /**
     * @return the set of all edges in this graph.
     */
    public Set<jebl.evolution.graphs.Edge> getEdges() {
        return new LinkedHashSet<jebl.evolution.graphs.Edge>(edges.values());
    }

    /**
     * @param degree the number of edges connected to a node
     * @return a set containing all nodes in this graph of the given degree.
     */
    public Set<jebl.evolution.graphs.Node> getNodes(int degree) {
        Set<jebl.evolution.graphs.Node> nodes = new LinkedHashSet<jebl.evolution.graphs.Node>();
        for (jebl.evolution.graphs.Node node : getNodes()) {
            if (((SimpleNode)node).getDegree() == degree) nodes.add(node);
        }
        return nodes;
    }

	/**
	 * The set of external edges. This is a pretty inefficient implementation because
	 * a new set is constructed each time this is called.
	 * @return the set of external edges.
	 */
	public Set<jebl.evolution.graphs.Edge> getExternalEdges() {
		Set<jebl.evolution.graphs.Edge> externalEdges = new LinkedHashSet<jebl.evolution.graphs.Edge>();
		for (jebl.evolution.graphs.Edge edge : getEdges()) {
			if (((SimpleEdge)edge).isExternal()) {
 				externalEdges.add(edge);
			}
		}
		return externalEdges;
	}

	/**
	 * The set of internal edges. This is a pretty inefficient implementation because
	 * a new set is constructed each time this is called.
	 * @return the set of internal edges.
	 */
	public Set<jebl.evolution.graphs.Edge> getInternalEdges() {
		Set<jebl.evolution.graphs.Edge> internalEdges = new LinkedHashSet<jebl.evolution.graphs.Edge>();
		for (jebl.evolution.graphs.Edge edge : getEdges()) {
			if (!((SimpleEdge)edge).isExternal()) {
 				internalEdges.add(edge);
			}
		}
		return internalEdges;
	}

    // Attributable IMPLEMENTATION

    public void setAttribute(String name, Object value) {
        if (helper == null) {
            helper = new jebl.util.AttributableHelper();
        }
        helper.setAttribute(name, value);
    }

    public Object getAttribute(String name) {
        if (helper == null) {
            return null;
        }
        return helper.getAttribute(name);
    }

    public void removeAttribute(String name) {
        if( helper != null ) {
            helper.removeAttribute(name);
        }
    }

    public Set<String> getAttributeNames() {
        if (helper == null) {
            return Collections.emptySet();
        }
        return helper.getAttributeNames();
    }

    public Map<String, Object> getAttributeMap() {
        if (helper == null) {
            return Collections.emptyMap();
        }
        return helper.getAttributeMap();
    }

    // PRIVATE members

    private jebl.util.AttributableHelper helper = null;
    private final Set<jebl.evolution.graphs.Node> internalNodes = new LinkedHashSet<jebl.evolution.graphs.Node>();
    private final Map<jebl.evolution.taxa.Taxon, jebl.evolution.graphs.Node> externalNodes = new LinkedHashMap<jebl.evolution.taxa.Taxon, jebl.evolution.graphs.Node>();
    /**
     * A mapping between edges and edge length.
     */
    Map<HashPair, jebl.evolution.graphs.Edge> edges = new LinkedHashMap<HashPair, jebl.evolution.graphs.Edge>();

    final class SimpleNode extends BaseNode {

        /**
         * A tip having a taxon
         * @param taxon
         */
        private SimpleNode(jebl.evolution.taxa.Taxon taxon) {
            this.adjacencies = Collections.unmodifiableList(new ArrayList<jebl.evolution.graphs.Node>());
            this.taxon = taxon;
        }

        /**
         * An internal node.
         * @param adjacencies set of adjacent noeds
         */
        private SimpleNode(List<jebl.evolution.graphs.Node> adjacencies) {
            this.adjacencies = Collections.unmodifiableList(adjacencies);
            this.taxon = null;
        }

        /**
         * Add an adjacency.
         * @param node
         */
        public void addAdjacency(jebl.evolution.graphs.Node node) {
            List<jebl.evolution.graphs.Node> a = new ArrayList<jebl.evolution.graphs.Node>(adjacencies);
            a.add(node);
            adjacencies = Collections.unmodifiableList(a);
        }

        public jebl.evolution.taxa.Taxon getTaxon() {
            return taxon;
        }

        public int getDegree() {
            return (adjacencies == null ? 0 : adjacencies.size());
        }

        public List<jebl.evolution.graphs.Node> getAdjacencies() {
            return adjacencies;
        }

        public void setTaxa(jebl.evolution.taxa.Taxon to) {
            taxon = to;
        }

        // PRIVATE members
        private List<jebl.evolution.graphs.Node> adjacencies;
        private jebl.evolution.taxa.Taxon taxon;

    }

	final class SimpleEdge extends BaseEdge {

		private SimpleEdge(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2, double length) {
			this.node1 = node1;
			this.node2 = node2;
			this.length = length;
		}

		public jebl.evolution.graphs.Node getNode1() {
			return node1;
		}

		public jebl.evolution.graphs.Node getNode2() {
			return node2;
		}

		public double getLength() {
			return length;
		}

		private boolean isExternal() {
			return (node1.getDegree() == 1 || node2.getDegree() == 1);
		}

		private double length;
		jebl.evolution.graphs.Node node1, node2;
	}
}