package teaspoon.adaptation;

public class SequenceInfo {
	String Name = new String("");
	String Accession = new String("");
	String Genotype = new String("");
	String Strain = new String("");
	String Gene = new String("");	
	String Species = new String("");
	String GenomeAccession = new String("");
	String Country = new String("");
	String Continent = new String("");
	private int[] Sequence;
	private String Taxon;
	double DecimalDate = 0.0;
	double Year = 0.0;
	double Gap = 0.0;
	public SequenceInfo(){
		this.Name = null;
		this.Accession = null;
		this.Gene = null;
		this.GenomeAccession = null;
		this.Genotype = null;
		this.Species = null;
		this.Continent = null;
		this.Country = null;
		this.setSequence(null);
		this.setTaxon(null);
		this.Year = 0.0;
		this.DecimalDate = 0.0;
		this.Strain = null;
		this.Gap=0.0;
	}
	/**
	 * @return the taxon
	 */
	public String getTaxon() {
		return Taxon;
	}
	/**
	 * @param taxon the taxon to set
	 */
	public void setTaxon(String taxon) {
		Taxon = taxon;
	}
	/**
	 * @return the sequence
	 */
	public int[] getSequence() {
		return Sequence;
	}
	/**
	 * @param sequence the sequence to set
	 */
	public void setSequence(int[] sequence) {
		Sequence = sequence;
	}
	
	


	
}
