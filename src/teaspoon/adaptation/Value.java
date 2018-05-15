package teaspoon.adaptation;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/jnarag">@jnarag</a>
 * @version 0.1
 * 
 * Class Value is a data holding (utility) object for TEASPOON analyses.
 * TODO @jnarag please label fields
 * 
 */
public class Value {
	/* Analysis data e.g. counts and estimates */
	// TODO redundant - no longer used in the code (from Sam's old code)
	double Nadapt;
	// TODO replacement polymorphisms, low-frequency class?
    double[] r_low;
	// TODO silent polymorphisms, low-frequency class
    double[] s_low;
    
	// TODO replacement polymorphisms, medium-frequency class
	double[] r_mid;
	// TODO silent polymorphisms, medium-frequency class
    double[] s_mid;
    
	// TODO neutral ratio, replacement:silent polymorphisms, mid-frequency class (alternative name rs_ratio_mid)
    double[] neutral_ratio;
    
	// TODO replacement polymorphisms, high-frequency class
    double[] r_high;
	// TODO silent polymorphisms, high-frequency class
    double[] s_high;
    
	// TODO ratio of silent:replacement polymorphisms, high-frequency class
    double[] sr_ratio_high;
	// TODO ratio of replacement:silent polymorphisms, high-frequency class
    //
    double[] rs_ratio_high;
	// TODO ratio of replacement:silent substitutions, low-frequency class
    double[] rs_low;
    
	// TODO total adaptations e.g. replacement in all classes? Or multinomial estimate?
    double[] adaptations;
	// TODO originally the bootstrap index but not currently used in the code
    double[] Bstrap;
	// TODO estimated theta (average pairwise genetic diversity)
    double[] theta;
	// TODO waterson's S (another method for estimating genetic diversity, unbiased pairwise differences estimator)
    double[] watterson_S;
	// TODO watterson's pi (another method for estimating genetic diversity, unbiased segregating sites estimator)
    double[] watterson_pi;
	// TODO tajima's D - number of segregating sites divided by number of pairwise differences
    double[] tajimas_D;
	
    /* Book-keeping stuff - track position in the dataset */
	String row;			// TODO timepoint
	String column;		// TODO dataset name
	String dataset;		// TODO redundant - no longer used
    int codons;			// TODO Number of codons
    double[] no_windows;// TODO Number of sliding windows (for deepGenome analysis, it will be no of codons minus site columns with no sequences
	
	/**
	 * Default no-arg constructor.
	 * Initialises data arrays with size N and other sensible defaults.
	 * @param N
	 */
    public Value(int N){
		Nadapt=0.0;
		row = "NA";
		column = "NA";
		dataset = "NA";
        codons = 0;
        Bstrap = new double[N];
		r_mid = new double[N];
        s_mid = new double[N];
        neutral_ratio = new double[N];
        adaptations = new double[N];
        r_high = new double[N];
        s_high = new double[N];
        r_low = new double[N];
        s_low = new double[N];
        sr_ratio_high = new double[N];
        rs_ratio_high = new double[N];
        rs_low = new double[N];
        theta = new double[N];
        watterson_pi = new double[N];
        watterson_S = new double[N];
        tajimas_D = new double[N];
        no_windows = new double[N];
	}
}
