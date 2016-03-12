package jebl.evolution.align;

import jebl.evolution.distances.CannotBuildDistanceMatrixException;
import jebl.evolution.sequences.BasicSequence;
import jebl.evolution.sequences.SequenceType;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Matt Kearse
 * @version $Id: BartonSternberg.java 842 2007-11-12 22:33:48Z twobeers $
 *
 * Implements the BartonSternberg multiple sequence alignment algorithm.
 *
 * Note: this is not yet complete, it does not create an initial ordering
 * in which to add sequences to the profile.
 *
 * Also, after creating the profile, it just removes and adds each sequence back into
 * the profile a fixed number of times(currently two).
 */
public class BartonSternberg implements MultipleAligner {

    jebl.evolution.align.scores.Scores scores;
    NeedlemanWunschLinearSpaceAffine aligner;
    private int refinementIterations;
    private float gapOpen,gapExtend;
    private boolean freeGapsAtEnds;
    private boolean fastGuide;
    // if not null, scores are from estimate
    private jebl.evolution.align.scores.Scores origScores = null;

    private void establishScores(jebl.evolution.align.scores.Scores scores) {
        this.scores = scores;
        this.scores = jebl.evolution.align.scores.Scores.includeGaps(scores, -gapExtend, 0);
        aligner = new NeedlemanWunschLinearSpaceAffine(this.scores, gapOpen, gapExtend, freeGapsAtEnds);
    }

    public jebl.evolution.align.scores.Scores getEstimatedScores() {
        return origScores != null ? scores : null;
    }

    public BartonSternberg(jebl.evolution.align.scores.Scores scores, float gapOpen, float gapExtend, int refinementIterations,
                           boolean freeGapsAtEnds, boolean fastGuide) {
//        if (true) throw new RuntimeException("testing");
       this.gapOpen = gapOpen;
       this.gapExtend = gapExtend;
       this.freeGapsAtEnds = freeGapsAtEnds;

        this.fastGuide = fastGuide;

        this.refinementIterations = refinementIterations;
        establishScores(scores);
    }

    CompoundAlignmentProgressListener compoundProgress;
    /*private ProgressListener progress;
    private boolean cancelled = false;
    private int sectionsCompleted;
    private int totalSections;
    ProgressListener minorProgress = new ProgressListener() {
        public boolean setProgress(double fractionCompleted) {
            double totalProgress = (sectionsCompleted + fractionCompleted)/totalSections;
            if(progress.setProgress(totalProgress)) cancelled = true;
            return cancelled;
        }
    };*/

    // on entry from top (non recursive), compoundProgress should have allocated #tips - 1 utins of work

    private Profile align(jebl.evolution.trees.RootedTree tree, jebl.evolution.graphs.Node node, List<jebl.evolution.sequences.Sequence> seqs,
                          CompoundAlignmentProgressListener compoundProgress) {
        if( tree.isExternal(node) ) {
            final jebl.evolution.taxa.Taxon tax = tree.getTaxon(node);
            final int iSeq = Integer.parseInt(tax.getName());

            final Profile profile = new Profile(scores.getAlphabet().length());
            profile.addSequence(iSeq, seqs.get(iSeq).getString());
            return profile;
        }

        List<jebl.evolution.graphs.Node> children = tree.getChildren(node);                            assert( children.size() == 2 );
        final Profile left = align(tree, children.get(0), seqs, compoundProgress);
        if( compoundProgress.isCanceled() ) return null;

        final Profile right = align(tree, children.get(1), seqs, compoundProgress);
        if( compoundProgress.isCanceled() ) return null;

        compoundProgress.setSectionSize(1);
        AlignmentResult results[] = aligner.doAlignment(left, right, compoundProgress.getMinorProgress(), false);
        compoundProgress.incrementSectionsCompleted(1);
        if(compoundProgress.isCanceled()) return null;
        return Profile.combine(left, right, results[0], results[1]);
    }


    /**
     *
     * @param sourceSequences
     * @param progress
     * @param refineOnly if specified, then the input sequences are assumed to be aligned already,
     * and this function will only refine the alignment.
     */
    public final String[] align(List<jebl.evolution.sequences.Sequence> sourceSequences, jebl.util.ProgressListener progress, boolean refineOnly,
                          boolean estimateMatchMismatchCosts)
            throws CannotBuildDistanceMatrixException
    {
        if( origScores != null ) {
            establishScores(origScores);
        }

        final int numSequences = sourceSequences.size();

        Profile[] sequenceProfilesWithoutGaps = new Profile[numSequences];
        String[] sequencesWithoutGaps = new String[numSequences];
        for (int i = 0; i < numSequences; i++) {
            sequencesWithoutGaps[i] = Align.stripIllegalCharacters(sourceSequences.get(i).getString(), scores.getAlphabet(), false);
            sequenceProfilesWithoutGaps[i] = new Profile(i, sequencesWithoutGaps[i]);
        }

        int treeWork = refineOnly ? 0 : (fastGuide ? numSequences : numSequences*(numSequences - 1)/2);
        int alignmentWork = refineOnly ? 0 : numSequences - 1;
        int refinementWork = numSequences * refinementIterations;

        compoundProgress = new CompoundAlignmentProgressListener(progress,treeWork + refinementWork + alignmentWork);

        Profile profile;
        if (refineOnly) {
            String[] sequencesWithGaps = new String[numSequences];
            for (int i = 0; i < numSequences; i++) {
                sequencesWithGaps[i] = Align.stripIllegalCharacters(sourceSequences.get(i).getString(), scores.getAlphabet(), true);
            }
            profile = new Profile(Profile.calculateAlphabetSize(sequencesWithGaps));
            for (int i = 0; i < numSequences; i++) {
                assert(sequencesWithGaps[i].length() == sequencesWithGaps[0].length());
                profile.addSequence(i, sequencesWithGaps[i]);
            }
        } else {
            List<jebl.evolution.sequences.Sequence> sequencesForGuideTree = new ArrayList<jebl.evolution.sequences.Sequence>(sourceSequences.size());
            for (int i = 0; i < numSequences; i++) {
                jebl.evolution.sequences.Sequence s = sourceSequences.get(i);
                sequencesForGuideTree.add(new BasicSequence(s.getSequenceType(), jebl.evolution.taxa.Taxon.getTaxon("" + i), sequencesWithoutGaps[i]));
            }
            compoundProgress.setSectionSize(treeWork);
            // We want a binary rooted tree

            //long start = System.currentTimeMillis();
            final boolean estimateMatchCost = estimateMatchMismatchCosts && scores instanceof jebl.evolution.align.scores.NucleotideScores;

            final AlignmentTreeBuilderFactory.Result unrootedGuideTree =
                    fastGuide ?
                            AlignmentTreeBuilderFactory.build(sequencesForGuideTree, jebl.evolution.trees.TreeBuilderFactory.Method.NEIGHBOR_JOINING,
                                    this, compoundProgress.getMinorProgress(),true) :
                            AlignmentTreeBuilderFactory.build(sequencesForGuideTree, jebl.evolution.trees.TreeBuilderFactory.Method.NEIGHBOR_JOINING,
                                    aligner, compoundProgress.getMinorProgress());
            if (compoundProgress.isCanceled()) return null;
            //long duration = System.currentTimeMillis() - start;
            //System.out.println("took " + duration +  " for " + (fastGuide ? " fast" : "normal") + " guide tree");

            jebl.evolution.trees.RootedTree guideTree = jebl.evolution.trees.Utils.rootTreeAtCenter(unrootedGuideTree.tree);
            compoundProgress.incrementSectionsCompleted(treeWork);

            if( estimateMatchCost ) {
                final jebl.evolution.distances.DistanceMatrix distanceMat = unrootedGuideTree.distance;
                final double[][] distances = distanceMat.getDistances();
                double sum = 0.0;
                final int n = distances.length;
                for(int k = 0; k < n; ++k) {
                    for(int j = k+1; j < n; ++j) {
                        // ignore infinity and high values
                        sum += Math.min(5.0, distances[k][j]);
                    }
                }
                final double avg = sum / ((n * (n - 1)/2));

                final double percentmatches = 1 - (3.0/4.0) * (1 - Math.exp(-4.0 * avg / 3.0));

                origScores = scores;
                final jebl.evolution.align.scores.NucleotideScores nucleotideScores = new jebl.evolution.align.scores.NucleotideScores(scores, percentmatches);
                establishScores(nucleotideScores);
            }
            progress.setMessage("Building alignment");
            profile = align(guideTree, guideTree.getRootNode(), sequencesForGuideTree, compoundProgress);
            if (compoundProgress.isCanceled()) return null;
        }

        //now remove a single sequence, and we
        for (int j = 0; j < refinementIterations; j++) {
            String message = "Refining alignment";
            if(refinementIterations> 1) {
                message = message + " (iteration " +(j+1) + " of " + refinementIterations+ ")";
            }
            progress.setMessage(message);
            for (int i = 0; i < numSequences; ++i) {
//                if(j> 0&& i!= 8) continue;
//                Profile sequenceProfile = sequenceProfiles[i];
                boolean display = false;

                String sequence = profile.getSequence(i);
                if(j>= 0 && i== 8) {
//                    display = true;
                }
                if(display) {
                    System.out.println("remove sequence =" + sequence);
                    profile.print (true);
                }
                Profile sequenceProfile = new Profile(i, sequence);
                profile.remove(sequenceProfile);
//                aligner.setDebug(display);

                AlignmentResult results[] = aligner.doAlignment(profile, sequenceProfilesWithoutGaps[i], compoundProgress.getMinorProgress(), false);
//                aligner.setDebug(false);
                if (compoundProgress.isCanceled()) return null;
                compoundProgress.incrementSectionsCompleted(1);
                if(display){
                    profile.print(false);

                    System.out.println("result =" + results[0].size + "," + results[1].size + " from " + profile.length() + "," + sequenceProfile.length());
                }
                profile = Profile.combine(profile, sequenceProfilesWithoutGaps[i], results[0], results[1]);
                if(display) {
                    profile.print(true);
                }
            }
        }

        String[] results = new String[numSequences];
        for (int i = 0; i < numSequences; i++) {
            results[i]= profile.getSequence(i);
        }
        return results;
    }

    public static void main(String[] arguments) throws IOException, jebl.evolution.io.ImportException {
        File file = new File(arguments[0]);
        SequenceType sequenceType = SequenceType.AMINO_ACID;

        jebl.evolution.io.FastaImporter importer = new jebl.evolution.io.FastaImporter(file, sequenceType);
        List<jebl.evolution.sequences.Sequence> xsequences = importer.importSequences();
        List<String> sequenceStrings = new ArrayList<String>();
        int count = 0;
        int maximum = 10;
        for (jebl.evolution.sequences.Sequence sequence : xsequences) {
            BasicSequence basic = (BasicSequence) sequence;
            String string = basic.getCleanString();
            sequenceStrings.add(string);
            System.out.println(string);
            if(count++ >= maximum) break;
        }
        System.out.println ();
        count = 0;
        for (jebl.evolution.sequences.Sequence sequence : xsequences) {
            BasicSequence basic = (BasicSequence) sequence;
            String string = basic.getString();
            System.out.println(string);
            if (count++ >= maximum) break;
        }
        long start = System.currentTimeMillis();
        BartonSternberg alignment = new BartonSternberg( new jebl.evolution.align.scores.Blosum60(), 20, 1, 2, true, false);
        String[] sequences = sequenceStrings.toArray(new String[0]);
        System.out.println("aligning " + sequences.length);
        try {
            String results[] = alignment.align(xsequences, null, false, false);
            for (String result : results) {
                System.out.println(result);
            }
            System.out.println("took " +(System.currentTimeMillis() - start) + " milliseconds");
        } catch (CannotBuildDistanceMatrixException e) {
            e.printStackTrace();
        }
    }

    public jebl.evolution.alignments.Alignment doAlign(List<jebl.evolution.sequences.Sequence> seqs, jebl.evolution.trees.RootedTree guideTree, jebl.util.ProgressListener progress) {
        final int count = seqs.size();
        final CompoundAlignmentProgressListener p = new CompoundAlignmentProgressListener(progress, count - 1);

        Profile profile = align(guideTree, guideTree.getRootNode(), seqs, p);
        if (p.isCanceled()) {
            return null;
        }

        List<jebl.evolution.sequences.Sequence> aSeqs = new ArrayList<jebl.evolution.sequences.Sequence>(count);
        for (int i = 0; i < count; i++) {
            String seq = profile.getSequence(i);
            final jebl.evolution.sequences.Sequence s = seqs.get(i);
            aSeqs.add(new BasicSequence(s.getSequenceType(), s.getTaxon(), seq));
        }
        return new jebl.evolution.alignments.BasicAlignment(aSeqs);
    }


    public jebl.evolution.alignments.Alignment doAlign(jebl.evolution.alignments.Alignment a1, jebl.evolution.alignments.Alignment a2, jebl.util.ProgressListener progress) {
        List<jebl.evolution.sequences.Sequence> seqs1 = a1.getSequenceList();
        List<jebl.evolution.sequences.Sequence> seqs2 = a2.getSequenceList();

        final int size1 = seqs1.size();
        final int size2 = seqs2.size();

        final Profile profile1 = new Profile(a1, scores.getAlphabet().length());
        final Profile profile2 = new Profile(a2, scores.getAlphabet().length(), size1);

        AlignmentResult results[] = aligner.doAlignment(profile1, profile2, progress, false);
        if (progress.isCanceled()) {
            return null;
        }
        Profile profile = Profile.combine(profile1, profile2, results[0], results[1]);

        final int count = size1 + size2;
        List<jebl.evolution.sequences.Sequence> aSeqs = new ArrayList<jebl.evolution.sequences.Sequence>(count);
        for (int i = 0; i < count; i++) {
            final String seq = profile.getSequence(i);
            final jebl.evolution.sequences.Sequence s = (i < size1) ? seqs1.get(i) : seqs2.get(i - size1);
            aSeqs.add(new BasicSequence(s.getSequenceType(), s.getTaxon(), seq));
        }
        return new jebl.evolution.alignments.BasicAlignment(aSeqs);
    }

    public jebl.evolution.alignments.Alignment doAlign(jebl.evolution.alignments.Alignment alignment, jebl.evolution.sequences.Sequence sequence, jebl.util.ProgressListener progress) {

        for (jebl.evolution.sequences.Sequence seq : alignment.getSequenceList()) {
            if (seq.getTaxon().getName().equals(sequence.getTaxon().getName())) {
                throw new IllegalArgumentException("Sequence taxon " + sequence.getTaxon().getName() + " appears in alignment and sequence.");
            }
        }

        final Profile aprofile = new Profile(alignment, scores.getAlphabet().length(),1);

        final Profile sprofile = new Profile(scores.getAlphabet().length());
        sprofile.addSequence(0, sequence.getString());

        AlignmentResult results[] = aligner.doAlignment(aprofile, sprofile, progress, false);
        Profile profile = Profile.combine(aprofile, sprofile, results[0], results[1]);

        List<jebl.evolution.sequences.Sequence> seqs1 = alignment.getSequenceList();
        List<jebl.evolution.sequences.Sequence> seqs2 = new ArrayList<jebl.evolution.sequences.Sequence>(); seqs2.add(sequence);

        final int size1 = seqs1.size();
        final int size2 = seqs2.size();

        final int count = size1 + size2;
        List<jebl.evolution.sequences.Sequence> aSeqs = new ArrayList<jebl.evolution.sequences.Sequence>(count);
        for (int i = 0; i < count; i++) {
            final String seq = profile.getSequence(i);
            final jebl.evolution.sequences.Sequence s = (i < count-1) ? seqs1.get(i) : seqs2.get(0);
            aSeqs.add(new BasicSequence(s.getSequenceType(), s.getTaxon(), seq));
        }
        return new jebl.evolution.alignments.BasicAlignment(aSeqs);
    }

    public double getScore() {
        return aligner.getScore();
    }
}
