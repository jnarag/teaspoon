package jebl.evolution.trees;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jebl.evolution.graphs.Graph;

/**
 * Root an unrooted tree. This class works as a wrapper over any tree to root it. There are two
 * constructors, one which roots the tree at any internal node, the other roots the tree between any two
 * internal nodes. Be aware that rooting between nodes where one of them has less than 3 adjacencies may
 * be problematic when converting back from the Newick format.
 *
 * @author Joseph Heled
 * @version $Id: RootedFromUnrooted.java 936 2008-08-06 14:12:07Z rambaut $
 *
 */

public class RootedFromUnrooted implements RootedTree {
	/**
	 * The unrooted tree
	 */
	private Tree source;

	/**
	 * Root of rooted tree. Either an existing internal node or a new "synthetic" node.
	 */
	private jebl.evolution.graphs.Node root;
	/**
	 * Maps each nodes to its parent.
	 */
	private Map<jebl.evolution.graphs.Node, jebl.evolution.graphs.Node> parents;

	/**
	 *  Children of the synthetic root (when rooted between nodes)
	 */
	private jebl.evolution.graphs.Node topLeft, topRight;
	/**
	 * branch lengths from synthetic root to its children (when rooted between nodes)
	 */
	private double rootToLeft, rootToRight;
	private boolean intentUnrooted;

	/**
	 * Set <arg>parent</arg> as parent of <arg>node</arg>, and recursivly set parents for node subtree
	 * (whose root is parent)
	 * @param node
	 * @param parent
	 */
    private void setParent(jebl.evolution.graphs.Node node, jebl.evolution.graphs.Node parent) {
		parents.put(node, parent);
		for( jebl.evolution.graphs.Node adj : source.getAdjacencies(node) ) {
			if( adj != parent && ! (node == topLeft && adj == topRight) && !(node == topRight && adj == topLeft) ) {
                setParent(adj, node);
			}
		}
	}

	/**
	 * Root tree at some internal node.
	 *
	 * @param source tree to root
	 * @param root  internal node to root at
	 * @param intentUnrooted
	 */
	public RootedFromUnrooted(Tree source, jebl.evolution.graphs.Node root, boolean intentUnrooted) {
		this.source = source;
		this.root = root;
		this.intentUnrooted = intentUnrooted;
		topLeft = topRight  = null;
		rootToLeft = rootToRight = 0.0;
		parents = new LinkedHashMap<jebl.evolution.graphs.Node, jebl.evolution.graphs.Node>();
		for( jebl.evolution.graphs.Node adj : source.getAdjacencies(root) ) {
            setParent(adj, root);
			}
		}

	/**
	 * Root source by creating a new internal node whose children are (the adjacent) left and right.
	 * @param source
	 * @param left
	 * @param right
	 * @param fromLeft branch from new root to left node.
	 */
	public RootedFromUnrooted(Tree source, jebl.evolution.graphs.Node left, jebl.evolution.graphs.Node right, double fromLeft) {
		this.source = source;
		intentUnrooted = false;
		topLeft = left;
		topRight = right;
		rootToLeft = fromLeft;
		try {
			rootToRight = source.getEdgeLength(left, right) - rootToLeft;
		} catch (Graph.NoEdgeException e) {
			// bug
		}
		parents = new LinkedHashMap<jebl.evolution.graphs.Node, jebl.evolution.graphs.Node>();

		// This is just a handle used to refer to the root so create the simplest possible implementation...
        root = new BaseNode() { public int getDegree() { return 2; } };

		parents.put(root, null);
        setParent(left, root);
        setParent(right, root);
    }

	public List<jebl.evolution.graphs.Node> getChildren(jebl.evolution.graphs.Node node) {
		ArrayList<jebl.evolution.graphs.Node> s = new ArrayList<jebl.evolution.graphs.Node>(getAdjacencies(node));
		if( node != root ) {
			s.remove(getParent(node));
		}
		return s;
	}

	public boolean hasHeights() {
		return false;
	}

	private double findNodeHeightFromTips(jebl.evolution.graphs.Node node) {
		if( isExternal(node) ) return 0.0;

		double h = 0.0;
		for( jebl.evolution.graphs.Node n : getChildren(node) ) {
			h = Math.max(h, getLength(n) + findNodeHeightFromTips(n));
		}
		return h;
	}

	public double getHeight(jebl.evolution.graphs.Node node) {
		double hr = findNodeHeightFromTips(root);
		if( node == root ) {
			return hr;
		}

		double toRoot = 0.0;
		while( node != root ) {
			toRoot += getLength(node);
			node = getParent(node);
		}
		return hr - toRoot;
	}

	public boolean hasLengths() {
		return true;
	}

	public double getLength(jebl.evolution.graphs.Node node) {
		if( node == root ) return 0.0;
		if( node == topLeft ) return rootToLeft;
		if( node == topRight ) return rootToRight;
		double l = 0.0;
		try {
			l = source.getEdgeLength(node, getParent(node));
		} catch (Graph.NoEdgeException e) {
			// bug, should not happen
		}
		return l;
	}

	public jebl.evolution.graphs.Node getParent(jebl.evolution.graphs.Node node) {
		return parents.get(node);
	}

    public jebl.evolution.graphs.Node getRootNode() {
        return root;
    }

	public boolean conceptuallyUnrooted() {
		return intentUnrooted;
	}

	public Set<jebl.evolution.graphs.Node> getExternalNodes() {
		return source.getExternalNodes();
	}

	public Set<jebl.evolution.graphs.Node> getInternalNodes() {
		HashSet<jebl.evolution.graphs.Node> s = new LinkedHashSet<jebl.evolution.graphs.Node>(source.getInternalNodes());
		s.add(root);
		return s;
	}

	public Set<jebl.evolution.taxa.Taxon> getTaxa() {
		return source.getTaxa();
	}

	public jebl.evolution.taxa.Taxon getTaxon(jebl.evolution.graphs.Node node) {
		if( node == root ) return null;
		return source.getTaxon(node);
	}

	public boolean isExternal(jebl.evolution.graphs.Node node) {
		return node != root && source.isExternal(node);
	}

	public jebl.evolution.graphs.Node getNode(jebl.evolution.taxa.Taxon taxon) {
		return source.getNode(taxon);
	}

	public void renameTaxa(jebl.evolution.taxa.Taxon from, jebl.evolution.taxa.Taxon to) {
		source.renameTaxa(from, to);
	}

	/**
	 * Returns a list of edges connected to this node
	 *
	 * @param node
	 * @return the set of nodes that are attached by edges to the given node.
	 */
	public List<jebl.evolution.graphs.Edge> getEdges(jebl.evolution.graphs.Node node) {
		return source.getEdges(node);
	}

	public jebl.evolution.graphs.Node[] getNodes(jebl.evolution.graphs.Edge edge) {
		return source.getNodes(edge);
	}

	public List<jebl.evolution.graphs.Node> getAdjacencies(jebl.evolution.graphs.Node node) {
		// special case when syntetic root
		if( topLeft != null ) {
			if( node == root ) {
				jebl.evolution.graphs.Node[] d = {topLeft, topRight};
				return Arrays.asList(d);
			}
			if( node == topLeft || node == topRight ) {
				List<jebl.evolution.graphs.Node> s = new ArrayList<jebl.evolution.graphs.Node>(source.getAdjacencies(node));
				s.remove(node == topLeft ? topRight : topLeft);
				s.add(root);
				return s;
			}
		}
		return source.getAdjacencies(node);
	}

	public double getEdgeLength(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) throws Graph.NoEdgeException {
		// special case when syntetic root
		if( topLeft != null ) {
			if( node2 == root ) {
				jebl.evolution.graphs.Node tmp = node1;
				node1 = node2;
				node2 = tmp;
			}
			if( node1 == root ) {
				if( ! (node2 == topLeft || node2 == topRight) ) {
					throw new Graph.NoEdgeException();
				}
				return node2 == topLeft ? rootToLeft : rootToRight;
			}
		}
		return source.getEdgeLength(node1, node2);
	}

	public jebl.evolution.graphs.Edge getEdge(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) throws Graph.NoEdgeException {
		return source.getEdge(node1, node2);
	}

	public Set<jebl.evolution.graphs.Node> getNodes() {
		Set<jebl.evolution.graphs.Node> nodes = new LinkedHashSet<jebl.evolution.graphs.Node>(getInternalNodes());
		nodes.addAll(getExternalNodes());
		if( topLeft != null ) {
			nodes.add(root);
		}
		return nodes;
	}

	/**
	 * @return the set of all edges in this graph.
	 */
	public Set<jebl.evolution.graphs.Edge> getEdges() {
		return source.getEdges();
	}

	/**
	 * The set of external edges.
	 * @return the set of external edges.
	 */
	public Set<jebl.evolution.graphs.Edge> getExternalEdges() {
		return source.getExternalEdges();
	}

	/**
	 * The set of internal edges.
	 * @return the set of internal edges.
	 */
	public Set<jebl.evolution.graphs.Edge> getInternalEdges() {
		return source.getInternalEdges();
	}

	public Set<jebl.evolution.graphs.Node> getNodes(int degree) {
		Set<jebl.evolution.graphs.Node> nodes = source.getNodes(degree);
		if( degree == 2 ) {
			nodes.add(root);
		}
		return nodes;
	}

	public boolean isRoot(jebl.evolution.graphs.Node node) {
		return node == root;
	}

	// Attributable IMPLEMENTATION

	public void setAttribute(String name, Object value) {
		source.setAttribute(name, value);
	}

	public Object getAttribute(String name) {
		return source.getAttribute(name);
	}

	public void removeAttribute(String name) {
		source.removeAttribute(name);
	}

	public Set<String> getAttributeNames() {
		return source.getAttributeNames();
	}

	public Map<String, Object> getAttributeMap() {
		return source.getAttributeMap();
	}
}