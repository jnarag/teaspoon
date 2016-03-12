package jebl.evolution.graphs;

/**
 * Represents a node in a graph or tree. In general it is
 * used only as a handle to traverse a graph or tree structure and
 * it has no methods or instance variables.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: Edge.java 295 2006-04-14 14:59:10Z rambaut $
 */
public interface Edge extends jebl.util.Attributable {

    double getLength();
}
