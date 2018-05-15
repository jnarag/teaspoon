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
	// TODO number of adaptive substitutions?
	double Nadapt;
	// TODO replacement substitutions, low-frequency class?
    double[] rl;
	// TODO silent substitutions, low-frequency class?
    double[] sl;
    
	// TODO replacement substitutions, medium-frequency class?
	double[] rm;
	// TODO silent substitutions, medium-frequency class?
    double[] sm;
    
	// TODO number replacement substitutions, total??
    double[] nr;
    
	// TODO replacement substitutions, high-frequency class?
    double[] rh;
	// TODO silent substitutions, high-frequency class?
    double[] sh;
    
	// TODO ratio of silent:replacement substitutions, high-frequency class?
    double[] sh_rh;
	// TODO ratio of replacement:silent substitutions, high-frequency class?
    // FIXME SHOULDN'T THIS BE MEDIUM FREQUENCY CLASS??
    double[] rh_sh;
	// TODO ratio of replacement:silent substitutions, low-frequency class?
    double[] rl_sl;
    
	// TODO total adaptations e.g. replacement in all classes? Or multinomial estimate?
    double[] adaptations;
	// TODO number of bootstraps or bootstrap replicate index?
    double[] Bstrap;
	// TODO estimated theta?
    double[] theta;
	// TODO waterson's S?
    double[] watterson_S;
	// TODO watterson's pi?
    double[] watterson_pi;
	// TODO tajima's D?
    double[] tajimas_D;
	
    /* Book-keeping stuff - track position in the dataset */
	String row;			// TODO Alignment row, e.g. sequence?
	String column;		// TODO Alignment col, e.g. site?
	String dataset;		// TODO Alignment dataset, e.g. sample?
    int codons;			// TODO Number of (? whole?) codons?
    double[] no_windows;// TODO Number of sliding windows (not window size?)
	
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
		rm = new double[N];
        sm = new double[N];
        nr = new double[N];
        adaptations = new double[N];
        rh = new double[N];
        sh = new double[N];
        rl = new double[N];
        sl = new double[N];
        sh_rh = new double[N];
        rh_sh = new double[N];
        rl_sl = new double[N];
        theta = new double[N];
        watterson_pi = new double[N];
        watterson_S = new double[N];
        tajimas_D = new double[N];
        no_windows = new double[N];
	}
}
