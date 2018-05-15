package teaspoon.adaptation;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * Class Value is a data holding (utility) object for TEASPOON analyses. 
 * 
 * @author <a href="http://github.com/jnarag">@jnarag</a>
 * @version 0.1
 */
public class Value {

	/* Analysis data e.g. counts */
	
	double Nadapt;			// redundant - no longer used in the code (from Sam's old code)
    double[] r_low;			// replacement polymorphisms, low-frequency class	
    double[] s_low;			// silent polymorphisms, low-frequency class  
	
	double[] r_mid;			// replacement polymorphisms, medium-frequency class
    double[] s_mid;			// silent polymorphisms, medium-frequency class
    double[] neutral_ratio;	// neutral ratio, replacement:silent polymorphisms, mid-frequency class (alternative name rs_ratio_mid)
    
    double[] r_high;		// replacement polymorphisms, high-frequency class
    double[] s_high;		// silent polymorphisms, high-frequency class
    	
    double[] sr_ratio_high;	// ratio of silent:replacement polymorphisms, high-frequency class
    double[] rs_ratio_high;	// ratio of replacement:silent polymorphisms, high-frequency class	
    double[] rs_low;		// ratio of replacement:silent substitutions, low-frequency class
    
    /* Analysed data - estimates etc */
	
    double[] adaptations;	// total adaptations e.g. replacement in all classes? Or multinomial estimate?
    double[] Bstrap;		// originally the bootstrap index but not currently used in the code	
    double[] theta;			// estimated theta (average pairwise genetic diversity)
    double[] watterson_S;	// waterson's S (another method for estimating genetic diversity, unbiased pairwise differences estimator)
    double[] watterson_pi;	// watterson's pi (another method for estimating genetic diversity, unbiased segregating sites estimator)
    double[] tajimas_D;		// tajima's D - number of segregating sites divided by number of pairwise differences
	
    /* Book-keeping stuff - track position in the dataset */

    String row;				// timepoint
	String column;			// dataset name
	String dataset;			// redundant - no longer used
    int codons;				// Number of codons
    double[] no_windows;	// Number of sliding windows (for deepGenome analysis, it will be no of codons minus site columns with no sequences
	
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
