package jebl.evolution.distances;

import java.util.List;

import jebl.evolution.sequences.SequenceType;

/**
 *
 * @author Joseph Heled
 * @version $Id: ModelBasedDistanceMatrix.java 786 2007-09-18 22:45:40Z matt_kearse $
 *
 */
public class ModelBasedDistanceMatrix  {
    protected static final double MAX_DISTANCE = 1000.0;

    protected double freqR, freqY;

    /**
     * @param sequences
     * @return array holding the count of each canonical state in the sequences.
     *         E.g. for nucleotide sequences, this array will have length 4.
     */
    private int[] countStates(List<jebl.evolution.sequences.Sequence> sequences) {
        if (sequences.isEmpty()) {
            throw new IllegalArgumentException("No sequences passed in - unable to determine sequence type");
        }
        SequenceType sequenceType = sequences.get(0).getSequenceType();
        final int canonicalStateCount = sequenceType.getCanonicalStateCount();
        int[] counts = new int[canonicalStateCount];
        for( jebl.evolution.sequences.Sequence sequence : sequences ) {
            if (!sequence.getSequenceType().equals(sequenceType)) {
                throw new IllegalArgumentException("Sequences of mixed type");
            }
            for( int i : sequence.getStateIndices() ) {
                // ignore non definite states (ask alexei)
                if( i < canonicalStateCount ) {
                    ++counts[i];
                }
            }
        }
        return counts;
    }

    /**
     * Same as countStates, but if any of the counts would normally be 0, this
     * method adds 1 to each count to avoid counts of 0.
     *
     * @param sequences
     * @return approximation of state counts, each guaranteed to be > 0.
     */
    private int[] countStatesSafe(List<jebl.evolution.sequences.Sequence> sequences) {
        int[] counts = countStates(sequences);
        int numSequences = counts.length;

        boolean anyZero = false;

        for (int i=0; i < numSequences; i++) {
            anyZero |= (counts[i] == 0);
        }

        // if any of the counts are 0, adjust all of them by 1 to avoid
        // division by 0 in extreme cases
        if (anyZero) {
            for (int i = 0; i < numSequences; ++i) {
                counts[i]++;
            }
        }
        return counts;
    }

    private double[] getFrequenciesMaybeSafe(List<jebl.evolution.sequences.Sequence> sequences, boolean safe) {
        SequenceType sequenceType = sequences.get(0).getSequenceType();
        int[] counts = (safe ? countStatesSafe(sequences) : countStates(sequences));
        int canonicalStateCount = counts.length;
        double[] freqs = new double[canonicalStateCount];

        // calculate total number of residues
        long count = 0;
        for (int i=0; i < canonicalStateCount; i++) {
            count += counts[i];
        }
        for (int i=0; i < canonicalStateCount; i++) {
            freqs[i] = (double) counts[i] / (double) count;
        }

        if (sequenceType.equals(SequenceType.NUCLEOTIDE)) {
           freqR = freqs[jebl.evolution.sequences.Nucleotides.A_STATE.getIndex()] + freqs[jebl.evolution.sequences.Nucleotides.G_STATE.getIndex()];
           freqY = freqs[jebl.evolution.sequences.Nucleotides.C_STATE.getIndex()] + freqs[jebl.evolution.sequences.Nucleotides.T_STATE.getIndex()];
        }
        return freqs;
    }


    /**
     *
     * As a side effect, this method sets freqR and freqY if called on
     * nucleotide sequences.
     *
     * @param sequences A list of sequences of the same type
     * @return Approximation of the relative canonical state
     *         frequencies in the sequences; Each frequency
     *         is guaranteed to be > 0 (and therefore it can
     *         only be an approximation).
     */
    protected double[] getFrequenciesSafe(List<jebl.evolution.sequences.Sequence> sequences) {
        return getFrequenciesMaybeSafe(sequences, true);
     }

    protected double[] getFrequenciesSafe(jebl.evolution.alignments.Alignment alignment) {
        return getFrequenciesSafe(alignment.getSequenceList());
    }

    /**
     * As a side effect, this method sets freqR and freqY if called on
     * nucleotide sequences.
     *
     * @param sequences A list of sequences of the same type
     * @return Relative canonical state frequencies in the sequences;
     */
    protected double[] getFrequencies(List<jebl.evolution.sequences.Sequence> sequences) {
        return getFrequenciesMaybeSafe(sequences, false);
     }

    protected double[] getFrequencies(jebl.evolution.alignments.Alignment alignment) {
        return getFrequencies(alignment.getSequenceList());
    }
}