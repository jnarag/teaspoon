package jebl.evolution.io;

import jebl.evolution.graphs.Node;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.trees.SimpleRootedTree;
import jebl.evolution.trees.Tree;

import java.io.EOFException;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: NewickImporter.java 931 2008-07-01 01:31:28Z richardmoir $
 */
public class NewickImporter implements TreeImporter {
    private final ImportHelper helper;
    private boolean unquotedLabels;

    /**
     * Constructor
     * @param reader  tree text
     * @param unquotedLabels if true, try to read unqouted lables containing spaces
     */
    public NewickImporter(Reader reader, boolean unquotedLabels) {
        helper = new ImportHelper(reader);
        this.unquotedLabels = unquotedLabels;
    }

    /**
     * Returns an iterator over a set of elements of type T.
     *
     * @return an Iterator.
     */
    public Iterator<Tree> iterator() {
        return new Iterator<Tree>() {

            public boolean hasNext() {
                boolean hasNext = false;
                try {
                    hasNext = hasTree();
                } catch (IOException e) {
                    // deal with errors by stopping the iteration
                } catch (ImportException e) {
                    // deal with errors by stopping the iteration
                }
                return hasNext;
            }

            public Tree next() {
                Tree tree = null;
                try {
                    tree = importNextTree();
                } catch (IOException e) {
                    // deal with errors by stopping the iteration
                } catch (ImportException e) {
                    // deal with errors by stopping the iteration
                }
                if (tree == null) throw new NoSuchElementException("No more trees in this file");
                return tree;
            }

            public void remove() {
                throw new UnsupportedOperationException("operation is not supported by this Iterator");
            }
        };
    }

    public boolean hasTree() throws IOException, ImportException {
        try {
            helper.skipUntil("(");
            helper.unreadCharacter('(');
        } catch (EOFException e) {
            return false;
        }

        return true;
    }

    public Tree importNextTree() throws IOException, ImportException {

        try {
            helper.skipUntil("(");
            helper.unreadCharacter('(');

            return readTree();
        } catch (EOFException e) {
            //
            throw new ImportException("error");
        }
    }

    public List<Tree> importTrees() throws IOException, ImportException {
        List<Tree> trees = new ArrayList<Tree>();

        while (hasTree()) {
            final Tree t = importNextTree();
            if( t != null ) {
              trees.add(t);
            }
        }

        return trees;
    }

    private RootedTree readTree() throws IOException, ImportException {
        SimpleRootedTree tree = new SimpleRootedTree();

        readInternalNode(tree);

        return tree;
    }

    /**
     * Reads a internal or external node along with its incoming branch and
     * its complete subtree from helper and adds it to the tree. If the branch
     * length is not specified in the newick file, it is assumed to be 1.0.
     * This method calls {@link #readInternalNode(jebl.evolution.trees.SimpleRootedTree)}
     * or {@link #readExternalNode(jebl.evolution.trees.SimpleRootedTree)} to read
     * the internal or external node and add it to the tree.
     *
     * @param tree Tree to which to add the new node along with is incoming branch
     * @return the internal or external node read from helper and added to the tree.
     */
    private Node readBranch(SimpleRootedTree tree) throws IOException, ImportException
    {
        Node branch;

        if (helper.nextCharacter() == '(') { // it's an internal node
            branch = readInternalNode(tree);
        } else { // it's an external node
            branch = readExternalNode(tree);
        }

        if (helper.getLastDelimiter() == ':') {
            double length = helper.readDouble(",():;");
            tree.setLength(branch, length);
        } else {
        	tree.setLength(branch, 1.0);
        }

        return branch;
    }

    /**
     * Reads a node in. This could be a polytomy. Calls readBranch on each branch
     * in the node.
     */
    private Node readInternalNode(SimpleRootedTree tree) throws IOException, ImportException
    {
        List<Node> children = new ArrayList<Node>();

        // read the opening '('
        helper.readCharacter();

        // read the first child
        children.add(readBranch(tree));

        // We used to require that an internal node has at least 2 children, however this
        // caused bug 4303 and I'm not sure if any other code depends on the outdegree 2
        /*if (helper.getLastDelimiter() != ',') {
            throw new ImportException.BadFormatException("Missing ',' in tree");
        } */

        // read subsequent children (recursively)
        while (helper.getLastDelimiter() == ',') {
            children.add(readBranch(tree));
        }

        // should have had a closing ')'
        if (helper.getLastDelimiter() != ')') {
            throw new ImportException.BadFormatException("Missing closing ')' in tree");
        }

        final Node node = tree.createInternalNode(children);

        try {
          // find the next delimiter
          String token = helper.readToken(":(),;");

            if (token.length() > 0) {
                node.setAttribute("label", NexusImporter.parseValue(token));
            }

            // If there is a metacomment before the branch length indicator (:), then it is a node attribute
            NexusImporter.parseAndClearMetaComments(node, helper);

        } catch( EOFException e) {
            // Ok if we just finished
        }

        return node;
    }

    /**
     * Reads an external node from helper and adds it to the specified tree. As a side effect,
     * the new node is added to the tree and helper's read position is advanced.
     * @return the external node read from helper
     */
    private Node readExternalNode(SimpleRootedTree tree) throws IOException, ImportException
    {
        String label = helper.readToken(":(),;");
        while( unquotedLabels && helper.getLastDelimiter() == ' ' ) {
            label = label + " " + helper.readToken(":(),;");
        }
        try {
            return tree.createExternalNode(Taxon.getTaxon(label));
        } catch (IllegalArgumentException e) {
           throw new ImportException.DuplicateTaxaException(e.getMessage());
        }
    }
}
