/**
 * 
 */
package teaspoon.app.GUI.models;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 5 Dec 2017
 * @version 0.0.1
 */
public class TeaspoonModel {

	/*
	 * Data models
	 * @TODO checks and other sensible things
	 */
	private Object[][] data; 
	private Object[][] filesDataModel;
	private Object[][] maskTracksModel;
	
	/* 
	 * Teaspoon analysis variables
	 */
	private File analysisDirectory;				// Working directory to trawl for input alignments
	private int alignmentLength;				// Maximum (equal!) alignment length
	private int bootstraps;						// Number of bootstrap replicates
	private int slidingWindowSize;				// Size of sliding window
	private ArrayList<int[]> maskTracks;		// Array of masking tracks for alignment, each an int tuple of nucleotide positions 5'->3' e.g. [10,31]
	private Object[] sequenceAlignments;		// Sequence alignments in FASTA FIXME decide type...
	private Date[] samplingDates;				// Dates to associate with each alignment
	private int minReadCoverageDepth;			// Minimum read depth (# of reads) per nucleotide position required for analysis
	
	/*
	 * GUI variables
	 */
	private boolean areAlignmentsValidated;		// Have sequence alignments in wd been parsed with equal lengths?
	private boolean doSlidingWindowAnalysis;	// Enable sliding window analysis?
	private boolean doBootstrappedAnalysis;		// Conduct bootstrap sampling?
	
	/**
	 * Default no-arg constructor
	 */
	public TeaspoonModel(){
		// Set up some sensible defaults
		this.setData(new Object[1][1]);
		this.bootstraps = 100;
		this.slidingWindowSize = 50;
		this.minReadCoverageDepth = 95;
		this.areAlignmentsValidated = false;
		this.doSlidingWindowAnalysis = false;
		this.doBootstrappedAnalysis = true;
		
		
	}

	/**
	 * @return the siteData
	 */
	public Object[][] getData() {
		return data;
	}

	/**
	 * @param siteData the siteData to set
	 */
	private void setData(Object[][] data) {
		this.data = data;
	}

	/**
	 * @return the analysisDirectory
	 */
	public File getAnalysisDirectory() {
		return analysisDirectory;
	}

	/**
	 * @return the alignmentLength
	 */
	public int getAlignmentLength() {
		return alignmentLength;
	}

	/**
	 * @return the bootstraps
	 */
	public int getBootstraps() {
		return bootstraps;
	}

	/**
	 * @return the slidingWindowSize
	 */
	public int getSlidingWindowSize() {
		return slidingWindowSize;
	}

	/**
	 * @return the maskTracks
	 */
	public ArrayList<int[]> getMaskTracks() {
		return maskTracks;
	}

	/**
	 * @return the sequenceAlignments
	 */
	public Object[] getSequenceAlignments() {
		return sequenceAlignments;
	}

	/**
	 * @return the samplingDates
	 */
	public Date[] getSamplingDates() {
		return samplingDates;
	}

	/**
	 * @return the minReadCoverageDepth
	 */
	public int getMinReadCoverageDepth() {
		return minReadCoverageDepth;
	}

	/**
	 * @return the areAlignmentsValidated
	 */
	public boolean isAreAlignmentsValidated() {
		return areAlignmentsValidated;
	}

	/**
	 * @return the doSlidingWindowAnalysis
	 */
	public boolean isDoSlidingWindowAnalysis() {
		return doSlidingWindowAnalysis;
	}

	/**
	 * @return the doBootstrappedAnalysis
	 */
	public boolean isDoBootstrappedAnalysis() {
		return doBootstrappedAnalysis;
	}

	/**
	 * @param analysisDirectory the analysisDirectory to set
	 */
	public void setAnalysisDirectory(File analysisDirectory) {
		this.analysisDirectory = analysisDirectory;
	}

	/**
	 * @param alignmentLength the alignmentLength to set
	 */
	public void setAlignmentLength(int alignmentLength) {
		this.alignmentLength = alignmentLength;
	}

	/**
	 * @param bootstraps the bootstraps to set
	 */
	public void setBootstraps(int bootstraps) {
		this.bootstraps = bootstraps;
	}

	/**
	 * @param slidingWindowSize the slidingWindowSize to set
	 */
	public void setSlidingWindowSize(int slidingWindowSize) {
		this.slidingWindowSize = slidingWindowSize;
	}

	/**
	 * @param maskTracks the maskTracks to set
	 */
	public void setMaskTracks(ArrayList<int[]> maskTracks) {
		this.maskTracks = maskTracks;
	}

	/**
	 * @param sequenceAlignments the sequenceAlignments to set
	 */
	public void setSequenceAlignments(Object[] sequenceAlignments) {
		this.sequenceAlignments = sequenceAlignments;
	}

	/**
	 * @param samplingDates the samplingDates to set
	 */
	public void setSamplingDates(Date[] samplingDates) {
		this.samplingDates = samplingDates;
	}

	/**
	 * @param minReadCoverageDepth the minReadCoverageDepth to set
	 */
	public void setMinReadCoverageDepth(int minReadCoverageDepth) {
		this.minReadCoverageDepth = minReadCoverageDepth;
	}

	/**
	 * @param areAlignmentsValidated the areAlignmentsValidated to set
	 */
	public void setAreAlignmentsValidated(boolean areAlignmentsValidated) {
		this.areAlignmentsValidated = areAlignmentsValidated;
	}

	/**
	 * @param doSlidingWindowAnalysis the doSlidingWindowAnalysis to set
	 */
	public void setDoSlidingWindowAnalysis(boolean doSlidingWindowAnalysis) {
		this.doSlidingWindowAnalysis = doSlidingWindowAnalysis;
	}

	/**
	 * @param doBootstrappedAnalysis the doBootstrappedAnalysis to set
	 */
	public void setDoBootstrappedAnalysis(boolean doBootstrappedAnalysis) {
		this.doBootstrappedAnalysis = doBootstrappedAnalysis;
	}

	/**
	 * @return the filesDataModel
	 */
	public Object[][] getFilesDataModel() {
		return filesDataModel;
	}

	/**
	 * @return the maskTracksModel
	 */
	public Object[][] getMaskTracksModel() {
		return maskTracksModel;
	}

	/**
	 * @param filesDataModel the filesDataModel to set
	 */
	public void setFilesDataModel(Object[][] filesDataModel) {
		this.filesDataModel = filesDataModel;
	}

	/**
	 * @param maskTracksModel the maskTracksModel to set
	 */
	public void setMaskTracksModel(Object[][] maskTracksModel) {
		this.maskTracksModel = maskTracksModel;
	}
}
