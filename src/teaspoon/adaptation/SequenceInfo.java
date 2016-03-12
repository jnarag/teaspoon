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
	int[] Sequence;
	String Taxon;
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
		this.Sequence = null;
		this.Taxon = null;
		this.Year = 0.0;
		this.DecimalDate = 0.0;
		this.Strain = null;
		this.Gap=0.0;
	}
	
	


	
}
