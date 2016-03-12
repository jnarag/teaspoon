/*
 * NucleotideState.java
 *
 * (c) 2002-2005 JEBL Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package jebl.evolution.sequences;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: NucleotideState.java 924 2008-06-18 05:58:28Z matt_kearse $
 */
public final class NucleotideState extends State {

    NucleotideState(String name, String stateCode, int index) {
        super(name, stateCode, index);
    }

    NucleotideState(String name, String stateCode, int index, NucleotideState[] ambiguities) {
        super(name, stateCode, index, ambiguities);
    }

    @Override
    public int compareTo(Object o) {
        // throws ClassCastException on across-class comparison
        NucleotideState that = (NucleotideState) o;
        return super.compareTo(that);
    }

    public boolean isGap() {
		return this == Nucleotides.GAP_STATE;
	}



}
