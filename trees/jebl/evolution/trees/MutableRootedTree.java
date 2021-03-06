package jebl.evolution.trees;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jebl.evolution.graphs.Graph;

/**
 * A simple rooted tree providing some ability to manipulate the tree.
 *
 *   - Root an unrooted tree using an outgroup.
 *   - Remove internal node: all children of node are adopted by it's parent.
 *   - Split/Refine node by creating two new children and distributing the children to new nodes.
 *   - Re-root a rooted tree given an outgroup.

 * @author Joseph Heled
 * @version $Id: MutableRootedTree.java 935 2008-07-22 16:52:04Z rambaut $
 *
 */

public class MutableRootedTree implements RootedTree {
    MutableRootedTree() {  super(); }

    /**
     * Construct a rooted tree from unrooted.
     *
     * @param tree      Unrooted tree to root
     * @param outGroup  Node in tree assumed to be the outgroup
     */
    public MutableRootedTree(Tree tree, jebl.evolution.graphs.Node outGroup) {
        if( ! tree.isExternal(outGroup) ) throw new IllegalArgumentException("Outgroup must be a tip");

        // Adjacency of node to become new root.
        jebl.evolution.graphs.Node root = tree.getAdjacencies(outGroup).get(0);

        try {
            MutableRootedNode newSubtreeRoot = rootAdjaceincesWith(tree, root, outGroup);


            // Add the outgroup in
            MutableRootedNode out = (MutableRootedNode)createExternalNode( tree.getTaxon(outGroup) );
            setLength(out, tree.getEdgeLength(outGroup, root));
            // Create new root
            ArrayList<MutableRootedNode> rootChildren = new ArrayList<MutableRootedNode>();
            rootChildren.add(out);
            rootChildren.add(newSubtreeRoot);
            //MutableRootedNode newRoot =
            	this.createInternalNode( rootChildren );
            setLength(newSubtreeRoot,0);
        } catch (Graph.NoEdgeException e) {
            // bug
        }
    }


    /**
     *  Remove internal node. Move all children to their grandparent.
     *  @param node  to be removed
     */
    public void removeInternalNode(jebl.evolution.graphs.Node node) {
        assert ! isExternal(node) && getRootNode() != node;

        MutableRootedNode parent = (MutableRootedNode)getParent(node);
        for( jebl.evolution.graphs.Node n : getChildren(node) ) {
            parent.addChild((MutableRootedNode)n);
        }
        parent.removeChild(node);
        internalNodes.remove(node);
    }

    /**
     *
     * @param node     Node to refine
     * @param leftSet  indices of children in the left new subtree.
     */
    public void refineNode(jebl.evolution.graphs.Node node, int[] leftSet) {
        List<jebl.evolution.graphs.Node> allChildren = getChildren(node);

        List<jebl.evolution.graphs.Node> left = new ArrayList<jebl.evolution.graphs.Node>();
        List<jebl.evolution.graphs.Node> right = new ArrayList<jebl.evolution.graphs.Node>();

        for( int n : leftSet ) {
            left.add(allChildren.get(n));
        }
        for( jebl.evolution.graphs.Node n : allChildren ) {
            if( !left.contains(n) ) {
                right.add(n);
            }
        }
        internalNodes.remove(node);
        MutableRootedNode saveRoot = rootNode;

        MutableRootedNode lnode = (left.size() > 1) ? createInternalNode(left) : (MutableRootedNode)left.get(0);
        MutableRootedNode rnode = (right.size() > 1) ? createInternalNode(right) : (MutableRootedNode)right.get(0);

        List<MutableRootedNode> nodes = new ArrayList<MutableRootedNode>(2);
        nodes.add(lnode);
        nodes.add(rnode);
        ((MutableRootedNode)node).replaceChildren(nodes);

        rootNode = saveRoot;
    }

    /**
     *  Re-root tree using an outgroup.
     * @param outGroup
     * @param attributeNames Move those attributes (if they exist in node) to their previous parent. The idea is to
     * preserve "branch" attributes which we now store in the child since only "node" properties are supported.
     */
    public void reRootWithOutgroup(jebl.evolution.graphs.Node outGroup, Set<String> attributeNames) {
        assert isExternal(outGroup);
        reRoot((MutableRootedNode)getAdjacencies(outGroup).get(0), attributeNames);
    }

    /**
     * Construct a rooted sub-tree from unrooted. Done recursivly: Given an internal node N and one adjacency A to become
     * the new parent, recursivly create subtrees for all adjacencies of N (ommiting A) using N as parent, and return
     * an internal node with all subtrees as children. A tip simply creates an external node and returns it.
     *
     * @param tree   Unrooted source tree
     * @param node   span sub-tree from this node
     * @param parent adjacency of node which serves as the parent.
     * @return  rooted subtree.
     * @throws Graph.NoEdgeException
     */
    private MutableRootedNode rootAdjaceincesWith(Tree tree, jebl.evolution.graphs.Node node, jebl.evolution.graphs.Node parent) throws Graph.NoEdgeException {
        if( tree.isExternal(node) ) {
            return (MutableRootedNode)createExternalNode( tree.getTaxon(node) );
        }

        List<jebl.evolution.graphs.Node> children = new ArrayList<jebl.evolution.graphs.Node>();
        for( jebl.evolution.graphs.Node adj : tree.getAdjacencies(node) ) {
            if( adj == parent ) continue;
            MutableRootedNode rootedAdj = rootAdjaceincesWith(tree, adj, node);
            setLength(rootedAdj, tree.getEdgeLength(adj, node));
            children.add(rootedAdj);
        }
        return createInternalNode(children);
    }

    /**
     * Similar to  rootAdjaceincesWith.
     * @param node
     * @param attributeNames
     */
    private void reRoot(MutableRootedNode node, Set<String> attributeNames) {
        MutableRootedNode parent = (MutableRootedNode)getParent(node);
        if( parent == null) {
            return;
        }
        double len = getLength(node);
        parent.removeChild(node);
        reRoot(parent, attributeNames);
        if( parent == getRootNode() ) {
            rootNode = node;
        }

        if( parent.getChildren().size() == 1 ) {
            parent = (MutableRootedNode)parent.getChildren().get(0);
            len += parent.getLength();
        }

        node.addChild(parent);
        parent.setLength(len);
        node.setParent(null);

        if( attributeNames != null ) {
            for( String name : attributeNames ) {
                Object s = node.getAttribute(name);
                if( s != null ) {
                    parent.setAttribute(name, s);
                    node.removeAttribute(name);
                }
            }
        }
    }

    public jebl.evolution.graphs.Node detachChildren(jebl.evolution.graphs.Node node, List<Integer> split) {
        assert( split.size() > 1 );

        List<jebl.evolution.graphs.Node> allChildren = getChildren(node);

        List<jebl.evolution.graphs.Node> detached = new ArrayList<jebl.evolution.graphs.Node>();

        for( int n : split ) {
            detached.add(allChildren.get(n));
        }

        MutableRootedNode saveRoot = rootNode;

        for( jebl.evolution.graphs.Node n : allChildren ) {
            if( detached.contains(n) ) {
                ((MutableRootedNode)node).removeChild(n);
            }
        }

        MutableRootedNode dnode = createInternalNode(detached);
        ((MutableRootedNode)node).addChild(dnode);

        rootNode = saveRoot;

        return dnode;
    }

    /**
     * Creates a new external node with the given taxon. See createInternalNode
     * for a description of how to use these methods.
     * @param taxon the taxon associated with this node
     * @return the created node reference
     */
    public jebl.evolution.graphs.Node createExternalNode(jebl.evolution.taxa.Taxon taxon) {
        MutableRootedNode node = new MutableRootedNode(taxon);
        externalNodes.put(taxon, node);
        return node;
    }

    /**
     * Once a SimpleRootedTree has been created, the node stucture can be created by
     * calling createExternalNode and createInternalNode. First of all createExternalNode
     * is called giving Taxon objects for the external nodes. Then these are put into
     * sets and passed to createInternalNode to create a parent of these nodes. The
     * last node created using createInternalNode is automatically the root so when
     * all the nodes are created, the tree is complete.
     *
     * @param children the child nodes of this nodes
     * @return the created node reference
     */
    public MutableRootedNode createInternalNode(List<? extends jebl.evolution.graphs.Node> children) {
        MutableRootedNode node = new MutableRootedNode(children);

        for (jebl.evolution.graphs.Node child : children) {
            ((MutableRootedNode)child).setParent(node);
        }

        internalNodes.add(node);

        rootNode = node;
        return node;
    }


    /**
     * @param node the node whose height is being set
     * @param height the height
     */
    public void setHeight(jebl.evolution.graphs.Node node, double height) {
        lengthsKnown = false;
        heightsKnown = true;

        // If a single height of a single node is set then
        // assume that all nodes have heights and by extension,
        // branch lengths as well as these will be calculated
        // from the heights
        hasLengths = true;
        hasHeights = true;

        ((MutableRootedNode)node).setHeight(height);
    }

    /**
     * @param node the node whose branch length (to its parent) is being set
     * @param length the length
     */
    public void setLength(jebl.evolution.graphs.Node node, double length) {
        heightsKnown = false;
        lengthsKnown = true;

        // If a single length of a single branch is set then
        // assume that all branch have lengths and by extension,
        // node heights as well as these will be calculated
        // from the lengths
        hasLengths = true;
        hasHeights = true;

        ((MutableRootedNode)node).setLength(length);
    }

    /**
     * @param node the node whose children are being requested.
     * @return the list of nodes that are the children of the given node.
     *         The list may be empty for a terminal node (a tip).
     */
    public List<jebl.evolution.graphs.Node> getChildren(jebl.evolution.graphs.Node node) {
        return new ArrayList<jebl.evolution.graphs.Node>(((MutableRootedNode)node).getChildren());
    }

    /**
     * @return Whether this tree has node heights available
     */
    public boolean hasHeights() {
        return hasHeights;
    }

    /**
     * @param node the node whose height is being requested.
     * @return the height of the given node. The height will be
     *         less than the parent's height and greater than it children's heights.
     */
    public double getHeight(jebl.evolution.graphs.Node node) {
        if (!hasHeights) throw new IllegalArgumentException("This tree has no node heights");
        if (!heightsKnown) calculateNodeHeights();
        return ((MutableRootedNode)node).getHeight();
    }

    /**
     * @return Whether this tree has branch lengths available
     */
    public boolean hasLengths() {
        return hasLengths;
    }

    /**
     * @param node the node whose branch length (to its parent) is being requested.
     * @return the length of the branch to the parent node (0.0 if the node is the root).
     */
    public double getLength(jebl.evolution.graphs.Node node) {
        if (!hasLengths) throw new IllegalArgumentException("This tree has no branch lengths");
        if (!lengthsKnown) calculateBranchLengths();
        return ((MutableRootedNode)node).getLength();
    }

    /**
     * @param node the node whose parent is requested
     * @return the parent node of the given node, or null
     *         if the node is the root node.
     */
    public jebl.evolution.graphs.Node getParent(jebl.evolution.graphs.Node node) {
        return ((MutableRootedNode)node).getParent();
    }

	public jebl.evolution.graphs.Edge getParentEdge(jebl.evolution.graphs.Node node) {
	    return ((MutableRootedNode)node).getEdge();
	}
    /**
     * The root of the tree has the largest node height of
     * all nodes in the tree.
     *
     * @return the root of the tree.
     */
    public jebl.evolution.graphs.Node getRootNode() {
        return rootNode;
    }

	public boolean isRoot(jebl.evolution.graphs.Node node) {
		return node == rootNode;
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
        return ((MutableRootedNode)node).getTaxon();
    }

    /**
     * @param node the node
     * @return true if the node is of degree 1.
     */
    public boolean isExternal(jebl.evolution.graphs.Node node) {
        return ((MutableRootedNode)node).getChildren().size() == 0;
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
        MutableRootedNode node = (MutableRootedNode)externalNodes.get(from);
        node.setTaxa(to);

        externalNodes.remove(from);
        externalNodes.put(to, node);
    }

    /**
     * Returns a list of edges connected to this node
     *
     * @param node
     * @return the set of nodes that are attached by edges to the given node.
     */
    public List<jebl.evolution.graphs.Edge> getEdges(jebl.evolution.graphs.Node node) {
        List<jebl.evolution.graphs.Edge> edges = new ArrayList<jebl.evolution.graphs.Edge>();
        for (jebl.evolution.graphs.Node adjNode : getAdjacencies(node)) {
            edges.add(((MutableRootedNode)adjNode).getEdge());

        }
        return edges;
    }

	/**
	 * Returns an array of 2 nodes which are the nodes at either end of the edge.
	 *
	 * @param edge
	 * @return an array of 2 edges
	 */
	public jebl.evolution.graphs.Node[] getNodes(jebl.evolution.graphs.Edge edge) {
		for (jebl.evolution.graphs.Node node : getNodes()) {
			if (((MutableRootedNode)node).getEdge() == edge) {
				return new jebl.evolution.graphs.Node[] { node, ((MutableRootedNode)node).getParent() };
			}
		}
		return null;
	}


    /**
     * @param node
     * @return the set of nodes that are attached by edges to the given node.
     */
    public List<jebl.evolution.graphs.Node> getAdjacencies(jebl.evolution.graphs.Node node) {
        return ((MutableRootedNode)node).getAdjacencies();
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
        if (((MutableRootedNode)node1).getParent() == node2) {
            return ((MutableRootedNode)node1).getEdge();
        } else if (((MutableRootedNode)node2).getParent() == node1) {
            return ((MutableRootedNode)node2).getEdge();
        } else {
            throw new Graph.NoEdgeException();
        }
    }

    /**
     * @param node1
     * @param node2
     * @return the length of the edge connecting node1 and node2.
     * @throws Graph.NoEdgeException
     *          if the nodes are not directly connected by an edge.
     */
    public double getEdgeLength(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) throws Graph.NoEdgeException {
        if (((MutableRootedNode)node1).getParent() == node2) {
            if (heightsKnown) {
                return ((MutableRootedNode)node2).getHeight() - ((MutableRootedNode)node1).getHeight();
            } else {
                return ((MutableRootedNode)node1).getLength();
            }
        } else if (((MutableRootedNode)node2).getParent() == node1) {
            if (heightsKnown) {
                return ((MutableRootedNode)node1).getHeight() - ((MutableRootedNode)node2).getHeight();
            } else {
                return ((MutableRootedNode)node2).getLength();
            }
        } else {
            throw new Graph.NoEdgeException();
        }
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
        Set<jebl.evolution.graphs.Edge> edges = new LinkedHashSet<jebl.evolution.graphs.Edge>();
        for (jebl.evolution.graphs.Node node : getNodes()) {
            if (node != getRootNode()) {
                edges.add(((MutableRootedNode)node).getEdge());
            }

        }
        return edges;
    }

	/**
	 * The set of external edges. This is a pretty inefficient implementation because
	 * a new set is constructed each time this is called.
	 * @return the set of external edges.
	 */
	public Set<jebl.evolution.graphs.Edge> getExternalEdges() {
		Set<jebl.evolution.graphs.Edge> edges = new LinkedHashSet<jebl.evolution.graphs.Edge>();
		for (jebl.evolution.graphs.Node node : getExternalNodes()) {
			edges.add(((MutableRootedNode)node).getEdge());
		}
		return edges;
	}

	/**
	 * The set of internal edges. This is a pretty inefficient implementation because
	 * a new set is constructed each time this is called.
	 * @return the set of internal edges.
	 */
	public Set<jebl.evolution.graphs.Edge> getInternalEdges() {
		Set<jebl.evolution.graphs.Edge> edges = new LinkedHashSet<jebl.evolution.graphs.Edge>();
		for (jebl.evolution.graphs.Node node : getInternalNodes()) {
			if (node != getRootNode()) {
			    edges.add(((MutableRootedNode)node).getEdge());
			}
		}
		return edges;
	}

    /**
     * @param degree the number of edges connected to a node
     * @return a set containing all nodes in this graph of the given degree.
     */
    public Set<jebl.evolution.graphs.Node> getNodes(int degree) {
        Set<jebl.evolution.graphs.Node> nodes = new LinkedHashSet<jebl.evolution.graphs.Node>();
        for (jebl.evolution.graphs.Node node : getNodes()) {
            // Account for no anncesstor of root, assumed by default in getDegree
            final int deg = node.getDegree();
            if (deg == degree) nodes.add(node);
        }
        return nodes;
    }

    /**
     * Set the node heights from the current branch lengths.
     */
    private void calculateNodeHeights() {

        if (!lengthsKnown) {
            throw new IllegalArgumentException("Can't calculate node heights because branch lengths not known");
        }

        nodeLengthsToHeights(rootNode, 0.0);

        double maxHeight = 0.0;
        for (jebl.evolution.graphs.Node externalNode : getExternalNodes()) {
            if (((MutableRootedNode)externalNode).getHeight() > maxHeight) {
                maxHeight = ((MutableRootedNode)externalNode).getHeight();
            }
        }

        for (jebl.evolution.graphs.Node node : getNodes()) {
            ((MutableRootedNode)node).setHeight(maxHeight - ((MutableRootedNode)node).getHeight());
        }

        heightsKnown = true;
    }

    /**
     * Set the node heights from the current node branch lengths. Actually
     * sets distance from root so the heights then need to be reversed.
     */
    private void nodeLengthsToHeights(MutableRootedNode node, double height) {

        double newHeight = height;

        if (node.getLength() > 0.0) {
            newHeight += node.getLength();
        }

        node.setHeight(newHeight);

        for (jebl.evolution.graphs.Node child : node.getChildren()) {
            nodeLengthsToHeights((MutableRootedNode)child, newHeight);
        }
    }

    /**
     * Calculate branch lengths from the current node heights.
     */
    protected void calculateBranchLengths() {

        if (!hasLengths) {
            throw new IllegalArgumentException("Can't calculate branch lengths because node heights not known");
        }

        nodeHeightsToLengths(rootNode, getHeight(rootNode));

        lengthsKnown = true;
    }

    /**
     * Calculate branch lengths from the current node heights.
     */
    private void nodeHeightsToLengths(MutableRootedNode node, double height) {
        final double h = node.getHeight();
        node.setLength(h >= 0 ? height - h : 1);

        for (jebl.evolution.graphs.Node child : node.getChildren()) {
            nodeHeightsToLengths((MutableRootedNode)child, node.getHeight());
        }

    }

    public void setConceptuallyUnrooted(boolean intent) {
        conceptuallyUnrooted = intent;
    }

    public boolean conceptuallyUnrooted() {
        return conceptuallyUnrooted;
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

    protected MutableRootedNode rootNode = null;
    protected final Set<jebl.evolution.graphs.Node> internalNodes = new LinkedHashSet<jebl.evolution.graphs.Node>();
    private final Map<jebl.evolution.taxa.Taxon, jebl.evolution.graphs.Node> externalNodes = new LinkedHashMap<jebl.evolution.taxa.Taxon, jebl.evolution.graphs.Node>();

    private boolean heightsKnown = false;
    private boolean lengthsKnown = false;

    private boolean hasHeights = false;
    private boolean hasLengths = false;

    private boolean conceptuallyUnrooted = false;

    private class MutableRootedNode extends BaseNode {
        public MutableRootedNode(jebl.evolution.taxa.Taxon taxon) {
            this.children = Collections.unmodifiableList(new ArrayList<jebl.evolution.graphs.Node>());
            this.taxon = taxon;
        }

        public MutableRootedNode(List<? extends jebl.evolution.graphs.Node> children) {
            this.children = Collections.unmodifiableList(new ArrayList<jebl.evolution.graphs.Node>(children));
            this.taxon = null;
        }


        public void removeChild(jebl.evolution.graphs.Node node) {
            List<jebl.evolution.graphs.Node> c = new ArrayList<jebl.evolution.graphs.Node>(children);
            c.remove(node);
            children = Collections.unmodifiableList(c);
        }

        public void addChild(MutableRootedNode node) {
            List<jebl.evolution.graphs.Node> c = new ArrayList<jebl.evolution.graphs.Node>(children);
            c.add(node);
            node.setParent(this);
            children = Collections.unmodifiableList(c);
        }

        public void replaceChildren(List<MutableRootedNode> nodes) {
            for( MutableRootedNode n : nodes ) {
                n.setParent(this);
            }
            children = Collections.unmodifiableList(new ArrayList<jebl.evolution.graphs.Node>(nodes));
        }


        public jebl.evolution.graphs.Node getParent() {
            return parent;
        }

        public void setParent(jebl.evolution.graphs.Node parent) {
            this.parent = parent;
        }

        public List<jebl.evolution.graphs.Node> getChildren() {
            return children;
        }

        public double getHeight() {
            return height;
        }

        // height above latest tip
        public void setHeight(double height) {
            this.height = height;
        }

        // length of branch to parent
        public double getLength() {
            return length;
        }

        public void setLength(double length) {
            this.length = length;
        }

        public int getDegree() {
            return children.size() + (this==rootNode?0:1);
        }

        /**
         * returns the edge connecting this node to the parent node
         * @return the edge
         */
        public jebl.evolution.graphs.Edge getEdge() {
            if (edge == null) {
                edge = new BaseEdge() {
                    public double getLength() {
                        return length;
                    }
                };
            }

            return edge;
        }

        /**
         * For a rooted tree, getting the adjacencies is not the most efficient
         * operation as it makes a new set containing the children and the parent.
         * @return the adjacaencies
         */
        public List<jebl.evolution.graphs.Node> getAdjacencies() {
            List<jebl.evolution.graphs.Node> adjacencies = new ArrayList<jebl.evolution.graphs.Node>();
            if (children != null) adjacencies.addAll(children);
            if (parent != null) adjacencies.add(parent);
            return adjacencies;
        }

        public jebl.evolution.taxa.Taxon getTaxon() {
            return taxon;
        }

        public void setTaxa(jebl.evolution.taxa.Taxon to) {
            taxon = to;
        }

        private List<jebl.evolution.graphs.Node> children;
        private jebl.evolution.taxa.Taxon taxon;

        private jebl.evolution.graphs.Node parent;
        private double height;
        private double length;

        private jebl.evolution.graphs.Edge edge = null;

    }
}
