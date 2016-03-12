package jebl.evolution.trees;

import jebl.evolution.graphs.Graph;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: FilteredRootedTree.java 936 2008-08-06 14:12:07Z rambaut $
 */
public abstract class FilteredRootedTree implements RootedTree {

    public FilteredRootedTree(final RootedTree source) {
        this.source = source;
    }

	public RootedTree getSource() {
		return source;
	}

    public boolean conceptuallyUnrooted() {
        return source.conceptuallyUnrooted();
    }

    public List<jebl.evolution.graphs.Node> getChildren(jebl.evolution.graphs.Node node) {
	    return source.getChildren(node);
    }

    public boolean hasHeights() {
        return source.hasHeights();
    }

    public double getHeight(jebl.evolution.graphs.Node node) {
        return source.getHeight(node);
    }

    public boolean hasLengths() {
        return source.hasLengths();
    }

    public double getLength(jebl.evolution.graphs.Node node) {
        return source.getLength(node);
    }

    public jebl.evolution.graphs.Node getParent(jebl.evolution.graphs.Node node) {
        return source.getParent(node);
    }

    public jebl.evolution.graphs.Node getRootNode() {
        return source.getRootNode();
    }

    public Set<jebl.evolution.graphs.Node> getExternalNodes() {
        return source.getExternalNodes();
    }

    public Set<jebl.evolution.graphs.Node> getInternalNodes() {
        return source.getInternalNodes();
    }

	public Set<jebl.evolution.graphs.Edge> getExternalEdges() {
		return source.getExternalEdges();
	}

	public Set<jebl.evolution.graphs.Edge> getInternalEdges() {
		return source.getInternalEdges();
	}

    public jebl.evolution.graphs.Node getNode(jebl.evolution.taxa.Taxon taxon) {
        return source.getNode(taxon);
    }

    public Set<jebl.evolution.taxa.Taxon> getTaxa() {
        return source.getTaxa();
    }

    public jebl.evolution.taxa.Taxon getTaxon(jebl.evolution.graphs.Node node) {
        return source.getTaxon(node);
    }

    public boolean isExternal(jebl.evolution.graphs.Node node) {
        return source.isExternal(node);
    }

    public List<jebl.evolution.graphs.Node> getAdjacencies(jebl.evolution.graphs.Node node) {
        return source.getAdjacencies(node);
    }

    public List<jebl.evolution.graphs.Edge> getEdges(jebl.evolution.graphs.Node node) {
        return source.getEdges(node);
    }

    public Set<jebl.evolution.graphs.Edge> getEdges() {
        return source.getEdges();
    }

	public jebl.evolution.graphs.Node[] getNodes(jebl.evolution.graphs.Edge edge) {
	    return source.getNodes(edge);
	}

    public jebl.evolution.graphs.Edge getEdge(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) throws Graph.NoEdgeException {
        return source.getEdge(node1, node2);
    }

    public double getEdgeLength(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) throws Graph.NoEdgeException {
        return source.getEdgeLength(node1, node2);
    }

    public Set<jebl.evolution.graphs.Node> getNodes() {
        return source.getNodes();
    }

    public Set<jebl.evolution.graphs.Node> getNodes(int degree) {
        return source.getNodes(degree);
    }

	public boolean isRoot(jebl.evolution.graphs.Node node) {
		return source.isRoot(node);
	}

    public void renameTaxa(jebl.evolution.taxa.Taxon from, jebl.evolution.taxa.Taxon to) {
        source.renameTaxa(from, to);
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

	// PRIVATE members

	protected final RootedTree source;
}