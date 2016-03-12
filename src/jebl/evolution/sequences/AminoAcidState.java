/*
 * AminoAcidState.java
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
 * @version $Id: AminoAcidState.java 810 2007-10-12 00:40:45Z twobeers $
 */
public final class AminoAcidState extends State {
    AminoAcidState(String name, String stateCode, int index) {
        super(name, stateCode, index);
    }

    AminoAcidState(String name, String stateCode, int index, AminoAcidState[] ambiguities) {
        super(name, stateCode, index, ambiguities);
    }

    @Override
    public int compareTo(Object o) {
        // throws ClassCastException on across-class comparison
        AminoAcidState that = (AminoAcidState) o;
        return super.compareTo(that);
    }

    // we do not need to override equals() and hashCode() because there is only one
    // unique instance of each state

    public boolean isGap() {
		return this == AminoAcids.GAP_STATE;
	}

	public boolean isStop() {
		return this == AminoAcids.STOP_STATE;
	}
}
