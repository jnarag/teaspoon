package jebl.evolution.align.scores;

import jebl.evolution.sequences.NucleotideState;
import jebl.evolution.sequences.Nucleotides;
import jebl.evolution.sequences.State;

import java.util.ArrayList;
import java.util.List;

/**
 * 
 * @author Richard Moir
 * @author Alexei Drummond
 *
 * @version $Id: NucleotideScores.java 916 2008-05-26 03:33:27Z matt_kearse $
 * 
 */
public class NucleotideScores extends Scores {

    float match = 5;
    float mismatchTransition = -4;
    float mismatchTransversion = -4;
    String name = "";
    private boolean includeAmbiguities;
    private String alphabet =
            Nucleotides.CANONICAL_STATES[0].getCode() +
                    Nucleotides.CANONICAL_STATES[1].getCode() +
                    Nucleotides.CANONICAL_STATES[2].getCode() +
                    Nucleotides.CANONICAL_STATES[3].getCode() +"U";

    public static final NucleotideScores IUB = new NucleotideScores(1.0f, -0.9f);
    public static final NucleotideScores CLUSTALW = new NucleotideScores(1.0f, 0.0f);

    protected NucleotideScores() {
    }

    public NucleotideScores(NucleotideScores scores) {
        name = scores.name;
        alphabet = scores.getAlphabet();
        match = scores.match;
        mismatchTransition = scores.mismatchTransition;
        mismatchTransversion = scores.mismatchTransversion;
    }

    /**
     * @param match match score
     * @param misMatch mismatch score
     */
    public NucleotideScores(float match, float misMatch) {
        this("", match, misMatch);
    }

    public NucleotideScores(float match, float misMatch, float ambiguousMatch) {
        this("", match, misMatch, misMatch, ambiguousMatch, false);
    }

    public NucleotideScores(String name, float match, float misMatch) {
        this(name, match, misMatch, misMatch, 0, true);
    }

    public NucleotideScores(String name, float match, float mismatchTransition, float mismatchTransversion) {
        this.name = name;
        buildScores(match, mismatchTransition, mismatchTransversion, 0, false);
    }

    public NucleotideScores(String name, float match, float mismatchTransition, float mismatchTransversion, float ambiguousMatch, boolean useWeightedAmbigousMatches) {
        this.name = name;
        buildScores(match, mismatchTransition, mismatchTransversion, ambiguousMatch, true, useWeightedAmbigousMatches);
    }

    public NucleotideScores(Scores scores, double percentmatches) {
        double match = Math.log(percentmatches/(4 *.25 *.25));
        double mismatch = Math.log((1-percentmatches)/(12 * .25 *.25));

        // normalize match from scores
        float ma = scores.score['A']['A'];
        float mm = (float)(mismatch * (ma/match));

        name = ((int)Math.round(100*percentmatches)) + "% similarity";
        buildScores(ma, mm, mm, 0, true);
        includeAdditionalCharacters(this, scores.getExtraResidues());
    }

    void buildScores(float match, float mismatchTransition, float mismatchTransversion, float ambiguousMatch, boolean includeAmbiguities) {
        buildScores(match, mismatchTransition, mismatchTransversion, ambiguousMatch, includeAmbiguities,false);
    }

    /**
     *
     * @param useWeightedAmbigousMatches true so that abiguities are converted to a fraction between 0 and 1 representing the fractional number of
     * matches of all canonical states represented by them. the score of such matches = mismatchScore + (match-mismatchScore)*fraction. For
     * example, if misamtch=0 and match = 1, then score(A,R)=0.5, score(T,R)=0, score(R,R)=0.5, score(B,B)=0.33.
     */
    void buildScores(float match, float mismatchTransition, float mismatchTransversion, float ambiguousMatch, boolean includeAmbiguities, boolean useWeightedAmbigousMatches) {

        this.match = match;
        this.mismatchTransition = mismatchTransition;
        this.mismatchTransversion = mismatchTransversion;
        this.includeAmbiguities = includeAmbiguities;

//        final int states = includeAmbiguities? Nucleotides.getStateCount():Nucleotides.getCanonicalStateCount();
        List<NucleotideState> states = new ArrayList<NucleotideState>();
        StringBuilder builder = new StringBuilder();
        for (NucleotideState state : Nucleotides.STATES) {
            if (state.isGap()) continue;
            if (state.isAmbiguous()&& !includeAmbiguities) continue;
            states.add (state);
            builder.append (state.getCode ());
        }
        // Add RNA "U" and the corresponding canonical state which is T_STATE to the list:
        alphabet = builder.toString() + "U";
        states.add(Nucleotides.T_STATE);

        int statesCount = states.size();
        float[][] scores = new float[statesCount][statesCount];
        for (int i = 0; i < statesCount; i++) {
            NucleotideState state1 = states.get(i);
            for (int j = 0; j < statesCount; j++) {
                NucleotideState state2 = states.get(j);
                float value;
                if (state1.equals(state2) && !useWeightedAmbigousMatches) {
                    value = match;
                }
                else if (state1.possiblyEqual(state2)) {
                    if (useWeightedAmbigousMatches) {
                        float min=Math.min(mismatchTransition, mismatchTransversion);                        
                        value = (float) (min + state1.fractionEqual(state2)*(match-min));
                    }
                    else {
                        value = ambiguousMatch;
                    }
                }
                else if (
                    (Nucleotides.isPurine(state1) && Nucleotides.isPurine(state2)) ||
                    (Nucleotides.isPyrimidine(state1) && Nucleotides.isPyrimidine(state2)) ) {
                        value = mismatchTransition;
                    } else {
                    value = mismatchTransversion;
                }

               /* float val = (i == j) ? match :
                        ((isPurine(i) == isPurine(j)) ? mismatchTransition : mismatchTransversion);
               */
                scores[i][j] = value;
            }
        }
        buildScores(scores);
    }

   /*
    private boolean isPurine(int state) {
        return Nucleotides.isPurine(Nucleotides.CANONICAL_STATES[state]);
    }
    */

    public String getName() {
        return name;
    }

    public final String getAlphabet() {
        return alphabet + getExtraResidues ();
    }

    public String toString() {
        String result = match + "/" + mismatchTransition;
        if(mismatchTransversion != mismatchTransition) {
            result = result + "/" + mismatchTransversion;
        }
        if(name.length()>  0){
            result = name + " (" + result + ")";
        }
        return result;
    }
}

