package teaspoon.app.utils;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * Class TeaspoonValues is a siteData holding (utility) object for TEASPOON analyses. 
 * 
 * @author <a href="http://github.com/jnarag">@jnarag</a>
 * @version 0.1.2
 */
public class TeaspoonValues {

	/* Analysis siteData e.g. counts of polymorphisms and their ratios */
	
	private double Nadapt;			// redundant - no longer used in the code (from Sam's old code)
    private double[] replacementsLowFreq;	// replacement polymorphisms, low-frequency class	
    private double[] silentsLowFreq;		// silentProb polymorphisms, low-frequency class  
	
	private double[] replacementsMidFreq;	// replacement polymorphisms, medium-frequency class
    private double[] silentsMidFreq;		// silentProb polymorphisms, medium-frequency class
    
    private double[] replacementsHighFreq;	// replacement polymorphisms, high-frequency class
    private double[] silentsHighFreq;		// silentProb polymorphisms, high-frequency class
    	
    private double[] neutralRatio;		// neutral ratio, replacement:silentProb polymorphisms, mid-frequency class (alternative name rs_ratio_mid)
    private double[] srRatioHighFreq;	// ratio of silentProb:replacement polymorphisms, high-frequency class
    private double[] rsRatioHighFreq;	// ratio of replacement:silentProb polymorphisms, high-frequency class	
    private double[] rsRatioLowFreq;	// ratio of replacement:silentProb substitutions, low-frequency class
    
    /* Analysed siteData - estimates etc */
	
    private double[] adaptations;	// total adaptations e.g. replacement in all classes? Or multinomial estimate?
    private double[] bootstrapIndex;		// originally the bootstrap index but not currently used in the code	
    private double[] theta;			// estimated theta (average pairwise genetic diversity)
    private double[] wattersonS;	// waterson's S (another method for estimating genetic diversity, unbiased pairwise differences estimator)
    private double[] wattersonPi;	// watterson's pi (another method for estimating genetic diversity, unbiased segregating sites estimator)
    private double[] tajimasD;		// tajima's D - number of segregating sites divided by number of pairwise differences
	
    /* Book-keeping stuff - track position in the dataset */

    private String row;		// timepoint
	private String column;			// dataset name
	private String dataset;			// redundant - no longer used
    private int numCodons;			// Number of num codons
    private double[] numSlidingWindows;	// Number of sliding windows (for deepGenome analysis, it will be no of numCodons minus site columns with no sequences
	
	/**
	 * Default no-arg constructor.
	 * Initialises siteData arrays with size numReplicates and other sensible defaults.
	 * @param numReplicates
	 */
    public TeaspoonValues(int N){
		setNadapt(0.0);
		setRow("NA");
		setColumn("NA");
		setDatasetLabel("NA");
        setNumCodons(0);
        setBootstrapIndex(new double[N]);
		setReplacementPolymorphsMidClass(new double[N]);
        setSilentPolymorphsMidClass(new double[N]);
        setNeutralRatio(new double[N]);
        setAdaptations(new double[N]);
        setReplacementSubsHigh(new double[N]);
        setSilentSubsHigh(new double[N]);
        setReplacementPolymorphsLowClass(new double[N]);
        setSilentPolymorphsLowClass(new double[N]);
        setSilentToReplacementPolymorphsHighClassRatio(new double[N]);
        setReplacementToSilentPolymorphsHighClassRatio(new double[N]);
        setReplacementToSilentPolymorphsLowClassRatio(new double[N]);
        setTheta(new double[N]);
        setWattersonPi(new double[N]);
        setWattersonS(new double[N]);
        setTajimasD(new double[N]);
        setNumWindows(new double[N]);
	}

	/**
	 * @return the numSlidingWindows
	 */
	public double[] getNumWindows() {
		return numSlidingWindows;
	}

	/**
	 * @param numSlidingWindows the numSlidingWindows to set
	 */
	protected void setNumWindows(double[] no_windows) {
		this.numSlidingWindows = no_windows;
	}

	/**
	 * @return the row
	 */
	public String getRow() {
		return row;
	}

	/**
	 * @param row the row to set
	 */
	public void setRow(String row) {
		this.row = row;
	}

	/**
	 * @return the silentsHighFreq
	 */
	public double[] getSilentSubsHigh() {
		return silentsHighFreq;
	}

	/**
	 * @param silentsHighFreq the silentsHighFreq to set
	 */
	public void setSilentSubsHigh(double[] s_high) {
		this.silentsHighFreq = s_high;
	}

	/**
	 * @return the replacementsHighFreq
	 */
	public double[] getReplacementSubsHigh() {
		return replacementsHighFreq;
	}

	/**
	 * @param replacementsHighFreq the replacementsHighFreq to set
	 */
	public void setReplacementSubsHigh(double[] r_high) {
		this.replacementsHighFreq = r_high;
	}

	/**
	 * @return the adaptations
	 */
	public double[] getAdaptations() {
		return adaptations;
	}

	/**
	 * @param adaptations the adaptations to set
	 */
	public void setAdaptations(double[] adaptations) {
		this.adaptations = adaptations;
	}

	/**
	 * @return the silentsMidFreq
	 */
	public double[] getSilentPolymorphsMidClass() {
		return silentsMidFreq;
	}

	/**
	 * @param silentsMidFreq the silentsMidFreq to set
	 */
	public void setSilentPolymorphsMidClass(double[] s_mid) {
		this.silentsMidFreq = s_mid;
	}

	/**
	 * @return the replacementsMidFreq
	 */
	public double[] getReplacementPolymorphsMidClass() {
		return replacementsMidFreq;
	}

	/**
	 * @param replacementsMidFreq the replacementsMidFreq to set
	 */
	public void setReplacementPolymorphsMidClass(double[] r_mid) {
		this.replacementsMidFreq = r_mid;
	}

	/**
	 * @return the silentsLowFreq
	 */
	public double[] getSilentPolymorphsLowClass() {
		return silentsLowFreq;
	}

	/**
	 * @param silentsLowFreq the silentsLowFreq to set
	 */
	public void setSilentPolymorphsLowClass(double[] s_low) {
		this.silentsLowFreq = s_low;
	}

	/**
	 * @return the replacementsLowFreq
	 */
	public double[] getReplacementPolymorphsLowClass() {
		return replacementsLowFreq;
	}

	/**
	 * @param replacementsLowFreq the replacementsLowFreq to set
	 */
	public void setReplacementPolymorphsLowClass(double[] r_low) {
		this.replacementsLowFreq = r_low;
	}

	/**
	 * @return the srRatioHighFreq
	 */
	public double[] getSilentToReplacementPolymorphsHighClassRatio() {
		return srRatioHighFreq;
	}

	/**
	 * @param srRatioHighFreq the srRatioHighFreq to set
	 */
	public void setSilentToReplacementPolymorphsHighClassRatio(
			double[] sr_ratio_high) {
		this.srRatioHighFreq = sr_ratio_high;
	}

	/**
	 * @return the rsRatioHighFreq
	 */
	public double[] getReplacementToSilentPolymorphsHighClassRatio() {
		return rsRatioHighFreq;
	}

	/**
	 * @param rsRatioHighFreq the rsRatioHighFreq to set
	 */
	public void setReplacementToSilentPolymorphsHighClassRatio(
			double[] rs_ratio_high) {
		this.rsRatioHighFreq = rs_ratio_high;
	}

	/**
	 * @return the neutralRatio e.g. replacement:silentProb polymorphisms in the mid-frequency class.
	 */
	public double[] getNeutralRatio() {
		return neutralRatio;
	}

	/**
	 * @param newNR the neutralRatio to set
	 */
	public void setNeutralRatio(double[] newNR) {
		this.neutralRatio = newNR;
	}

	/**
	 * @return the rs_low
	 */
	public double[] getReplacementToSilentPolymorphsLowClassRatio() {
		return rsRatioLowFreq;
	}

	/**
	 * @param rs_low the rs_low to set
	 */
	public void setReplacementToSilentPolymorphsLowClassRatio(double[] rs_low) {
		this.rsRatioLowFreq = rs_low;
	}

	/**
	 * @return the theta - pairwise genetic diversity estimate
	 */
	public double[] getTheta() {
		return theta;
	}

	/**
	 * @param theta the theta to set
	 */
	public void setTheta(double[] theta) {
		this.theta = theta;
	}

	/**
	 * @return the tajimasD - ratio of (count of segregating sites / number of pairwise differences)
	 */
	public double[] getTajimasD() {
		return tajimasD;
	}

	/**
	 * @param tajimasD the tajimasD to set
	 */
	public void setTajimasD(double[] tajimas_D) {
		this.tajimasD = tajimas_D;
	}

	/**
	 * @return the wattersonS - e.g. unbiased pairwise differences estimator
	 */
	public double[] getWattersonS() {
		return wattersonS;
	}

	/**
	 * @param wattersonS the wattersonS to set
	 */
	public void setWattersonS(double[] watterson_S) {
		this.wattersonS = watterson_S;
	}

	/**
	 * @return the wattersonPi - e.g. unbiased segregating sites estimator
	 */
	public double[] getWattersonPi() {
		return wattersonPi;
	}

	/**
	 * @param wattersonPi the wattersonPi to set
	 */
	public void setWattersonPi(double[] watterson_pi) {
		this.wattersonPi = watterson_pi;
	}

	/**
	 * @return the column
	 */
	public String getColumn() {
		return column;
	}

	/**
	 * @param column the column to set
	 */
	public void setColumn(String column) {
		this.column = column;
	}

	/**
	 * @return the numCodons
	 */
	public int getNumCodons() {
		return numCodons;
	}

	/**
	 * @param numCodons the numCodons to set
	 */
	public void setNumCodons(int numCodons) {
		this.numCodons = numCodons;
	}

	/**
	 * @return the nadapt
	 */
	public double getNadapt() {
		return Nadapt;
	}

	/**
	 * @param nadapt the nadapt to set
	 */
	public void setNadapt(double nadapt) {
		Nadapt = nadapt;
	}

	/**
	 * @return the dataset
	 */
	public String getDatasetLabel() {
		return dataset;
	}

	/**
	 * @param dataset the dataset to set
	 */
	public void setDatasetLabel(String dataset) {
		this.dataset = dataset;
	}

	/**
	 * @return the bootstrapIndex
	 */
	public double[] getBootstrapIndex() {
		return bootstrapIndex;
	}

	/**
	 * @param bootstrapIndex the bootstrapIndex to set
	 */
	public void setBootstrapIndex(double[] bootstrapIndex) {
		this.bootstrapIndex = bootstrapIndex;
	}
}
