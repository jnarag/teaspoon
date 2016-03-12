package jebl.evolution.trees;

import jebl.evolution.taxa.Taxon;
import jebl.util.AttributableHelper;
import jebl.evolution.graphs.Graph;

import java.util.*;

/**
 * A simple, and initially immutable rooted tree implementation. All returned collections
 * are defensively copied. The implementation of Node is private. A number of methods are
 * provided that can be used to construct a tree (createExternalNode & createInternalNode).
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: SimpleRootedTree.java 935 2008-07-22 16:52:04Z rambaut $
 */
final public class SubtreeRootedTree implements jebl.evolution.trees.RootedTree {

    /**
     * Make a copy of the given rooted tree
     * @param tree a rooted tree
     * @param includedTaxa
     */
    public SubtreeRootedTree(jebl.evolution.trees.RootedTree tree, Set<Taxon> includedTaxa) {
        createNodes(tree, tree.getRootNode(), includedTaxa);
    }

    /**
     * Clones the entire tree structure from the given RootedTree.
     * @param tree
     * @param node
     * @param includedTaxa may be empty
     * @return
     */
    private SubtreeRootedNode createNodes(jebl.evolution.trees.RootedTree tree, jebl.evolution.graphs.Node node, Set<Taxon> includedTaxa) {

        SubtreeRootedNode newNode;

        if (tree.isExternal(node)) {                                                
            //System.out.println(node);
            newNode = createExternalNode(tree.getTaxon(node), includedTaxa);

        } else {
            List<SubtreeRootedNode> children = new ArrayList<SubtreeRootedNode>();
            for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
                SubtreeRootedNode newChild = createNodes(tree, child, includedTaxa);
                if (newChild != null) {
                    children.add(newChild);
                }
            }
            if (children.size() >= 2) {
                newNode = createInternalNode(children);
            } else if (children.size() == 1) {

                newNode = children.get(0);
            } else {

                newNode = null;
            }
        }

        if (newNode != null) {
//        final Map<String, Object> map = node.getAttributeMap();
//        if( ! map.isEmpty() ) {
            for( Map.Entry<String, Object> e : node.getAttributeMap().entrySet() ) {
                newNode.setAttribute(e.getKey(), e.getValue());
            }
            // }
            setHeight(newNode, tree.getHeight(node));
        }


        return newNode;
    }

    /**
     * Creates a new external node with the given taxon. See createInternalNode
     * for a description of how to use these methods.
     * @param taxon the taxon associated with this node
     * @return the created node reference
     */
    private SubtreeRootedNode createExternalNode(Taxon taxon, Set<Taxon> includedTaxa) {
        if( getTaxa().contains(taxon) ) {
            throw new IllegalArgumentException("duplicate taxon "+taxon.getName());
        }

        if (!includedTaxa.contains(taxon)) {
            return null;
        }

        SubtreeRootedNode node = new SubtreeRootedNode(taxon);
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
    private SubtreeRootedNode createInternalNode(List<? extends jebl.evolution.graphs.Node> children) {
        SubtreeRootedNode node = new SubtreeRootedNode(children);

        for (jebl.evolution.graphs.Node child : children) {
            ((SubtreeRootedNode)child).setParent(node);
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

        ((SubtreeRootedNode)node).setHeight(height);
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

        ((SubtreeRootedNode)node).setLength(length);
    }

    /**
     * @param node the node whose children are being requested.
     * @return the list of nodes that are the children of the given node.
     *         The list may be empty for a terminal node (a tip).
     */
    public List<jebl.evolution.graphs.Node> getChildren(jebl.evolution.graphs.Node node) {
        return new ArrayList<jebl.evolution.graphs.Node>(((SubtreeRootedNode)node).getChildren());
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
        return ((SubtreeRootedNode)node).getHeight();
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
        return ((SubtreeRootedNode)node).getLength();
    }

    /**
     * @param node the node whose parent is requested
     * @return the parent node of the given node, or null
     *         if the node is the root node.
     */
    public jebl.evolution.graphs.Node getParent(jebl.evolution.graphs.Node node) {
        if (!(node instanceof SubtreeRootedNode)) {
            throw new IllegalArgumentException("Node, " + node.toString() + " is not an instance of SubtreeRootedNode");
        }
        return ((SubtreeRootedNode)node).getParent();
    }

    public jebl.evolution.graphs.Edge getParentEdge(jebl.evolution.graphs.Node node) {
        if (!(node instanceof SubtreeRootedNode)) {
            throw new IllegalArgumentException("Node, " + node.toString() + " is not an instance of SubtreeRootedNode");
        }
        return ((SubtreeRootedNode)node).getEdge();
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
    public Set<Taxon> getTaxa() {
        return new LinkedHashSet<Taxon>(externalNodes.keySet());
    }

    /**
     * @param node the node whose associated taxon is being requested.
     * @return the taxon object associated with the given node, or null
     *         if the node is an internal node.
     */
    public Taxon getTaxon(jebl.evolution.graphs.Node node) {
        if (!(node instanceof SubtreeRootedNode)) {
            throw new IllegalArgumentException("Node, " + node.toString() + " is not an instance of SubtreeRootedNode");
        }
        return ((SubtreeRootedNode)node).getTaxon();
    }

    /**
     * @param node the node
     * @return true if the node is of degree 1.
     */
    public boolean isExternal(jebl.evolution.graphs.Node node) {
        if (!(node instanceof SubtreeRootedNode)) {
            throw new IllegalArgumentException("Node, " + node.toString() + " is not an instance of SubtreeRootedNode");
        }
        boolean result= ((SubtreeRootedNode)node).getChildren().size() == 0;
        return  result;//((SubtreeRootedNode)node).getChildren().size() == 0;
    }

    /**
     * @param taxon the taxon
     * @return the external node associated with the given taxon, or null
     *         if the taxon is not a member of the taxa set associated with this tree.
     */
    public jebl.evolution.graphs.Node getNode(Taxon taxon) {
        return externalNodes.get(taxon);
    }

    public void renameTaxa(Taxon from, Taxon to) {
        SubtreeRootedNode node = (SubtreeRootedNode)externalNodes.get(from);

        // TT: The javadoc doesn't specify whether renameTaxa() should fail or silently do nothing
        // if Taxon from doesn't exist. But the code already threw a NullPointerException before (bug 4824),
        // so it's probably ok to throw a more informative IllegalArgumentException instead.
        if (node == null) {
            throw new IllegalArgumentException("Unknown taxon " + from + "; can't rename to " + to);
        }

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
            edges.add(((SubtreeRootedNode)adjNode).getEdge());

        }
        return edges;
    }

    /**
     * @param node
     * @return the set of nodes that are attached by edges to the given node.
     */
    public List<jebl.evolution.graphs.Node> getAdjacencies(jebl.evolution.graphs.Node node) {
        return ((SubtreeRootedNode)node).getAdjacencies();
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
    public jebl.evolution.graphs.Edge getEdge(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) throws NoEdgeException {
        if (((SubtreeRootedNode)node1).getParent() == node2) {
            return ((SubtreeRootedNode)node1).getEdge();
        } else if (((SubtreeRootedNode)node2).getParent() == node1) {
            return ((SubtreeRootedNode)node2).getEdge();
        } else {
            throw new NoEdgeException();
        }
    }

    /**
     * @param node1
     * @param node2
     * @return the length of the edge connecting node1 and node2.
     * @throws Graph.NoEdgeException
     *          if the nodes are not directly connected by an edge.
     */
    public double getEdgeLength(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) throws NoEdgeException {
        if (((SubtreeRootedNode)node1).getParent() == node2) {
            if (heightsKnown) {
                return ((SubtreeRootedNode)node2).getHeight() - ((SubtreeRootedNode)node1).getHeight();
            } else {
                return ((SubtreeRootedNode)node1).getLength();
            }
        } else if (((SubtreeRootedNode)node2).getParent() == node1) {
            if (heightsKnown) {
                return ((SubtreeRootedNode)node1).getHeight() - ((SubtreeRootedNode)node2).getHeight();
            } else {
                return ((SubtreeRootedNode)node2).getLength();
            }
        } else {
            throw new NoEdgeException();
        }
    }

    /**
     * Returns an array of 2 nodes which are the nodes at either end of the edge.
     *
     * @param edge
     * @return an array of 2 edges
     */
    public jebl.evolution.graphs.Node[] getNodes(jebl.evolution.graphs.Edge edge) {
        for (jebl.evolution.graphs.Node node : getNodes()) {
            if (((SubtreeRootedNode)node).getEdge() == edge) {
                return new jebl.evolution.graphs.Node[] { node, ((SubtreeRootedNode)node).getParent() };
            }
        }
        return null;
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
                edges.add(((SubtreeRootedNode)node).getEdge());
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
            edges.add(((SubtreeRootedNode)node).getEdge());
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
                edges.add(((SubtreeRootedNode)node).getEdge());
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
            final int deg = node.getDegree() ;
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
            if (((SubtreeRootedNode)externalNode).getHeight() > maxHeight) {
                maxHeight = ((SubtreeRootedNode)externalNode).getHeight();
            }
        }

        for (jebl.evolution.graphs.Node node : getNodes()) {
            ((SubtreeRootedNode)node).setHeight(maxHeight - ((SubtreeRootedNode)node).getHeight());
        }

        heightsKnown = true;
    }

    /**
     * Set the node heights from the current node branch lengths. Actually
     * sets distance from root so the heights then need to be reversed.
     */
    private void nodeLengthsToHeights(SubtreeRootedNode node, double height) {

        double newHeight = height;

        if (node.getLength() > 0.0) {
            newHeight += node.getLength();
        }

        node.setHeight(newHeight);

        for (jebl.evolution.graphs.Node child : node.getChildren()) {
            nodeLengthsToHeights((SubtreeRootedNode)child, newHeight);
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
    private void nodeHeightsToLengths(SubtreeRootedNode node, double height) {
        final double h = node.getHeight();
        node.setLength(h >= 0 ? height - h : 1);

        for (jebl.evolution.graphs.Node child : node.getChildren()) {
            nodeHeightsToLengths((SubtreeRootedNode)child, node.getHeight());
        }

    }

    public void setConceptuallyUnrooted(boolean intent) {
        conceptuallyUnrooted = intent;
    }

    public boolean conceptuallyUnrooted() {
        return conceptuallyUnrooted;
    }

    public boolean isRoot(jebl.evolution.graphs.Node node) {
        return node == rootNode;
    }

    // Attributable IMPLEMENTATION

    public void setAttribute(String name, Object value) {
        if (helper == null) {
            helper = new AttributableHelper();
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

    private AttributableHelper helper = null;

    protected SubtreeRootedNode rootNode = null;
    protected final Set<jebl.evolution.graphs.Node> internalNodes = new LinkedHashSet<jebl.evolution.graphs.Node>();
    private final Map<Taxon, jebl.evolution.graphs.Node> externalNodes = new LinkedHashMap<Taxon, jebl.evolution.graphs.Node>();

    private boolean heightsKnown = false;
    private boolean lengthsKnown = false;

    private boolean hasHeights = false;
    private boolean hasLengths = false;

    private boolean conceptuallyUnrooted = false;

    private class SubtreeRootedNode extends BaseNode {
        public SubtreeRootedNode(Taxon taxon) {
            this.children = Collections.unmodifiableList(new ArrayList<jebl.evolution.graphs.Node>());
            this.taxon = taxon;
        }

        public SubtreeRootedNode(List<? extends jebl.evolution.graphs.Node> children) {
            this.children = Collections.unmodifiableList(new ArrayList<jebl.evolution.graphs.Node>(children));
            this.taxon = null;
        }

        public void removeChild(jebl.evolution.graphs.Node node) {
            List<jebl.evolution.graphs.Node> c = new ArrayList<jebl.evolution.graphs.Node>(children);
            c.remove(node);
            children = Collections.unmodifiableList(c);
        }

        public void addChild(SubtreeRootedNode node) {
            List<jebl.evolution.graphs.Node> c = new ArrayList<jebl.evolution.graphs.Node>(children);
            c.add(node);
            node.setParent(this);
            children = Collections.unmodifiableList(c);
        }

        public void replaceChildren(List<SubtreeRootedNode> nodes) {
            for( SubtreeRootedNode n : nodes ) {
                n.setParent(this);
            }
            children = Collections.unmodifiableList(new ArrayList<jebl.evolution.graphs.Node>(nodes));
        }

        void swapChildren(int i0, int i1) {
            ArrayList<jebl.evolution.graphs.Node> nc = new ArrayList<jebl.evolution.graphs.Node>(children);
            //there was a user reported crash where i0 was > size of the array of children nodes
            if (i0 < 0 || i0 >= nc.size() || i1 < 0 || i1 >= nc.size()) {
                throw new IllegalArgumentException("Tried to swap children ("+i0+","+i1+") on node with " + nc.size() + " children");
            }
            final jebl.evolution.graphs.Node ni0 = nc.get(i0);
            nc.set(i0, nc.get(i1));
            nc.set(i1, ni0);
            children = Collections.unmodifiableList(nc);
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
            return children.size() +(this==rootNode?0:1);
        }

        public void setTaxa(Taxon to) {
            taxon = to;
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

        public Taxon getTaxon() {
            return taxon;
        }

        private List<jebl.evolution.graphs.Node> children;
        private Taxon taxon;

        private jebl.evolution.graphs.Node parent;
        private double height;
        private double length;

        private jebl.evolution.graphs.Edge edge = null;
    }
}