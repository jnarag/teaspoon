package jebl.evolution.trees;

import jebl.evolution.taxa.Taxon;
import jebl.util.ProgressListener;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 *
 * A {@link TreeBuilder} that builds a consensus tree for a set of trees on identical leaf sets.
 * This abstract base class is designed to be extended.
 *
 * @author Joseph Heled
 * @version $Id: ConsensusTreeBuilder.java 850 2007-12-06 04:51:28Z twobeers $
 */

public abstract class ConsensusTreeBuilder<T extends Tree> implements TreeBuilder<T> {

    /** Name of attribute specifing amount of support for branch */
    final static public String DEFAULT_SUPPORT_ATTRIBUTE_NAME = "Consensus support(%)";

	private String supportAttributeName;
	private boolean supportAsPercent;

    /**
     * Supported consensus methods.
     */
    public static enum Method {
        GREEDY("Greedy"),
        MRCAC("MRCA Clustering");

        private final String name;

        Method(String name) {
           this.name = name;
        }

        @Override
        public String toString() {
            return getName();
        }

        public String getName() {
            return name;
        }
    }

    /** Number of external nodes/taxa */
    protected final int nExternalNodes;

    /** List of common taxa in all trees */
    protected final List<Taxon> taxons;

	/**
	 * Check for consistancy and establish the common taxa
	 * @param trees
	 */
	ConsensusTreeBuilder(Tree[] trees) {
		this(trees, DEFAULT_SUPPORT_ATTRIBUTE_NAME, true);
	}

    /**
     * Checks that all trees have the same taxa.
     * @param trees with identical taxa sets for which the consensus tree is to be built
     * @param supportAttributeName name of attribute (see {@link jebl.evolution.trees.Tree#getAttribute(String)}) describing a tree's amount of support
     * @param supportInPercent when true, support is in percent (0 - 100), otherwise in number of trees from the set.
     * @throws IllegalArgumentException if the trees don't have identical taxa sets, or if trees is empty
     */
    ConsensusTreeBuilder(Tree[] trees, String supportAttributeName, boolean supportInPercent) throws IllegalArgumentException {
        if (trees.length == 0) {
            throw new IllegalArgumentException("Expected at least one tree, but got none");
        }
        this.supportAttributeName = supportAttributeName;
	    this.supportAsPercent = supportInPercent;
        this.nExternalNodes = trees[0].getExternalNodes().size();
        this.taxons = Collections.unmodifiableList(new ArrayList<Taxon>(trees[0].getTaxa()));
        Utils.assertAllTreesHaveTheSameTaxa(Arrays.asList(trees));
    }

    /**
     * Returns a human readable name of this consensus tree building method
     * @return A human readable name of this consensus tree building method
     */
    abstract public String getMethodDescription(); 

    protected String getSupportDescription(double supportThreshold) {
        String supporDescription;
        if (supportThreshold == 1.0) {
           supporDescription = "Strict";
        } else if (supportThreshold == .5) {
           supporDescription = "Majority";
        } else {
            supporDescription = "Above " + (100*supportThreshold) + "% support";
        }
        return supporDescription;
    }

    public String getSupportAttributeName() {
		return supportAttributeName;
	}

	public boolean isSupportAsPercent() {
		return supportAsPercent;
	}

    public void addProgressListener(ProgressListener listener) {
        synchronized(listeners) {
            listeners.add(listener);
        }
    }

    /**
     * Stops a ProgressListener from receiving progress. Note: Due to threading, the
     * the listener may keep receiving progress for a short while after this method returns.
     * @param listener The ProgressListener for which to no longer report progress.
     */
    public void removeProgressListener(ProgressListener listener) {
        synchronized(listeners) {
            listeners.remove(listener);
        }
    }

    /**
     * Informs all {@link #addProgressListener(jebl.util.ProgressListener) added} ProgressListeners of the progress
     * @param fractionCompleted progress fraction (between 0 and 1)
     * @return true if operation canceled
     */
    protected boolean fireSetProgress(double fractionCompleted) {
	    boolean canceled = false;
        List<ProgressListener> listenersCopy;
        synchronized(listeners) { // create a copy because we don't want to hold the lock while calling the listeners
            listenersCopy = new ArrayList<ProgressListener>(listeners);
        }
        for (ProgressListener listener : listenersCopy) {
            if (listener.setProgress(fractionCompleted)) {
                canceled = true;
            }
        }
        return canceled;
    }

    private final List<ProgressListener> listeners = new ArrayList<ProgressListener>();
}