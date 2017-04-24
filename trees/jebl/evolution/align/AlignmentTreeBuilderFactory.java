package jebl.evolution.align;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import jebl.evolution.alignments.Alignment;
import jebl.evolution.trees.Tree;
import jebl.evolution.trees.TreeBuilder;
import jebl.evolution.trees.TreeBuilderFactory;
import jebl.util.CompositeProgressListener;
import jebl.util.ProgressListener;

/**
 * @author Joseph Heled
 * @version $Id: AlignmentTreeBuilderFactory.java 842 2007-11-12 22:33:48Z twobeers $
 */
public class AlignmentTreeBuilderFactory {
    private final static Logger logger = Logger.getLogger(AlignmentTreeBuilderFactory.class.getName());

    private static interface DistanceMatrixBuilder {
        jebl.evolution.distances.DistanceMatrix buildDistanceMatrix(final ProgressListener progressListener) throws CannotBuildDistanceMatrixException;
    }

    static public class Result {
        public final Tree tree;
        public final jebl.evolution.distances.DistanceMatrix distance;

        Result(Tree tree, jebl.evolution.distances.DistanceMatrix distance) {
            this.tree = tree;
            this.distance = distance;
        }
    }    

    /**
     * private utility method to get rid of former code duplication in the two build methods from
     * alignment and list of sequences with pairwise aligner.
     * @param distanceMatrixBuilder Encapsulation of all information required to build distance matrix
     * @param method the tree building method to use
     * @param _progressListener must not be null
     * @return A tree building result (containing a tree and a distance matrix)
     */
    private static Result build(DistanceMatrixBuilder distanceMatrixBuilder, TreeBuilderFactory.Method method, ProgressListener _progressListener) throws CannotBuildDistanceMatrixException {
        // The requirement that progress isn't null has only been added on 2006-12-29.
        // For a grace period, we check whether it is null and only warn.
        if (_progressListener == null) {
            logger.warning("ProgressListener is null");
            _progressListener = ProgressListener.EMPTY;
        }
        CompositeProgressListener progressListener = new CompositeProgressListener(_progressListener, .5, .5);

        progressListener.beginSubtask("Computing genetic distance for all pairs");
        long start = System.currentTimeMillis();
        jebl.evolution.distances.DistanceMatrix distanceMatrix = distanceMatrixBuilder.buildDistanceMatrix(progressListener);
        logger.fine("took " +(System.currentTimeMillis() - start) + " to build distance matrix");
        progressListener.beginSubtask("Building tree");
        TreeBuilder treeBuilder = TreeBuilderFactory.getBuilder(method, distanceMatrix);
        treeBuilder.addProgressListener(progressListener);
        Result result = new Result(treeBuilder.build(), distanceMatrix);
        treeBuilder.removeProgressListener(progressListener);
        return result;
    }

    /**
     * @param alignment Alignment to calculate distance matrix from
     * @param method the tree building method to use
     * @param model substitution model for distance matrix: JukesCantor, TamuraNei, HKY or F84.
     * @param progressListener must not be null. If you are not interested in progress, pass in ProgressListener.EMPTY
     * @param useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable If the alignment contains pairs of sequences that are not overlapping it is not possible to build
     *          a tree from the alignment. If this parameter is false, then CannotBuildDistanceMatrixException will be thrown. If this parameter
     *          is true, then twice the maximum distance between all other pairs will be used non-overlapping pairs.
     * @return A tree building result (containing a tree and a distance matrix)
     * @throws CannotBuildDistanceMatrixException only if useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable is false
     */
    static public Result build(final Alignment alignment, TreeBuilderFactory.Method method, final TreeBuilderFactory.DistanceModel model, ProgressListener progressListener, final boolean useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable)
            throws CannotBuildDistanceMatrixException
    {
        DistanceMatrixBuilder matrixBuilder = new DistanceMatrixBuilder() {
            public jebl.evolution.distances.DistanceMatrix buildDistanceMatrix(final ProgressListener progressListener) throws CannotBuildDistanceMatrixException {
                switch( model ) {
                    case F84:
                        return new jebl.evolution.distances.F84DistanceMatrix(alignment, progressListener);
                    case HKY:
                        return new jebl.evolution.distances.HKYDistanceMatrix(alignment, progressListener,useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable);
                    case TamuraNei:
                        return new jebl.evolution.distances.TamuraNeiDistanceMatrix(alignment, progressListener,useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable);
                    case JukesCantor:
                    default:
                        return new jebl.evolution.distances.JukesCantorDistanceMatrix(alignment, progressListener,useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable);
                }
            }
        };
        return build(matrixBuilder, method, progressListener);
    }

    /**
     * @param alignment Alignment to calculate distance matrix from
     * @param method the tree building method to use
     * @param model substitution model for distance matrix: JukesCantor, TamuraNei, HKY or F84.
     * @param progressListener must not be null. If you are not interested in progress, pass in ProgressListener.EMPTY
     * @return A tree building result (containing a tree and a distance matrix)
     */
    static public Result build(final Alignment alignment, TreeBuilderFactory.Method method, final TreeBuilderFactory.DistanceModel model, ProgressListener progressListener)
            throws CannotBuildDistanceMatrixException
    {
        return build(alignment,method,model,progressListener,false);
    }

    /**
     *
     * @param seqs Sequences to build distance matrix from
     * @param method method the tree building method to use
     * @param aligner pairwise aligner which will be used to calculate a pairwise distance
     * @param progressListener must not be null. If you are not interested in progress, pass in ProgressListener.EMPTY
     * @return A tree building result (containing a tree and a distance matrix)
     */
    static public Result build(final List<jebl.evolution.sequences.Sequence> seqs, TreeBuilderFactory.Method method, final PairwiseAligner aligner,
                               ProgressListener progressListener)
            throws CannotBuildDistanceMatrixException
    {
        DistanceMatrixBuilder matrixBuilder = new DistanceMatrixBuilder() {
            public jebl.evolution.distances.DistanceMatrix buildDistanceMatrix(final ProgressListener progressListener) throws CannotBuildDistanceMatrixException {
                return new SequenceAlignmentsDistanceMatrix(seqs, aligner, progressListener);
            }
        };
        return build(matrixBuilder, method, progressListener);
    }

    static public Result build(List<jebl.evolution.sequences.Sequence> seqs, TreeBuilderFactory.Method method, MultipleAligner aligner,
                                /*boolean needDistances, */ProgressListener progress, final boolean useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable)
            throws CannotBuildDistanceMatrixException
    {
       // needDistances = false;

        jebl.evolution.trees.SimpleRootedTree gtree = new jebl.evolution.trees.SimpleRootedTree();
        List<jebl.evolution.graphs.Node> nodes = new ArrayList<jebl.evolution.graphs.Node>();
        for(jebl.evolution.sequences.Sequence s : seqs) {
            jebl.evolution.graphs.Node tip = gtree.createExternalNode(s.getTaxon());
            nodes.add(tip);
        }

        int nnodes = nodes.size();
        while( nnodes > 1 ) {
           List<jebl.evolution.graphs.Node> upnodes = new ArrayList<jebl.evolution.graphs.Node>();
            for(int k = 0; k < nnodes/2; ++k) {
                upnodes.add(gtree.createInternalNode(nodes.subList(2*k,2*k+2)));
            }
            if( (nnodes & 1) != 0 ) {
                upnodes.add(nodes.get(nnodes - 1));
            }
            nodes = upnodes;
            nnodes = nodes.size();
        }

        final int alignWork = seqs.size()-1;
        final int treeWork = 1;
        final int matrixWork = 0; // needDistances ? 1 : 0;

        CompoundAlignmentProgressListener p = new CompoundAlignmentProgressListener(progress,
                                                                                    alignWork + treeWork + matrixWork);
        final ProgressListener minorProgress = p.getMinorProgress();

        progress.setMessage("Building alignment for guide");
        p.setSectionSize(alignWork);
        final Alignment alignment = aligner.doAlign(seqs, gtree, minorProgress);
        if (p.isCanceled()) {
            return null;
        }
        p.incrementSectionsCompleted(alignWork);

        final boolean isProtein = seqs.get(0).getSequenceType().getCanonicalStateCount() > 4;

        final TreeBuilderFactory.DistanceModel distanceModel =
                isProtein ? TreeBuilderFactory.DistanceModel.JukesCantor : TreeBuilderFactory.DistanceModel.HKY;

        p.setSectionSize(treeWork);
        progress.setMessage("Building guide tree from alignment");
        final Result result = build(alignment, method, distanceModel, minorProgress, useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable);
        //final Tree guideTree = result.tree;
        p.incrementSectionsCompleted(treeWork);
        return result;
        /*
        DistanceMatrix distanceMatrix = null;
        if( needDistances ) {
            p.setSectionSize(matrixWork);
            progress.setMessage("Computing genetic distance for all pairs");

            if( distanceModel == TreeBuilderFactory.DistanceModel.HKY ) {
                distanceMatrix = new HKYDistanceMatrix(alignment, minorProgress);
            } else {
                distanceMatrix = new JukesCantorDistanceMatrix(alignment, minorProgress);
            }
            p.incrementSectionsCompleted(matrixWork);
        }
        return new Result(guideTree, distanceMatrix);
        */
    }
}
