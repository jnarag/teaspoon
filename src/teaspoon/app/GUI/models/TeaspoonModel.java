/**
 * 
 */
package teaspoon.app.GUI.models;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import teaspoon.app.TeaspoonMask;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.MainAlignmentParser;

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
public class TeaspoonModel extends AbstractTableModel{

	// JTable needed stuff - column names
	private final static String[] columnNames = {			
							//column index
		"Filename",			//0		File
		"Number Of Taxa",	//1		int
		"Number Of Sites",	//2		int
		"Decimal Date",		//3		float
		"is Ancestral?"		//4		boolean
	};				
	// JTable needed studd - column definitions
	private final static String[] columnDefinitions = {	
													//column index
		"Alignment = Alignment file",				//0		File
		"# taxa = Number of taxa",					//1		int
		"# sites (NT) = Number of nucleotides (NT)",//2		int
		"Date as decimal from 0 BCE",				//3		float
		"which alignment is considered ancestral"	//4		bool
	};
	// Column indices of Integers
	private final static Integer[] integerIndices = new Integer[]{1,2};
	// Column indices of Floats
	private final static Integer[] floatIndices = new Integer[]{3};
	// Default values for (hopefully) sizing the table, etc
	public final Object[] longValues = {
		null, 
		new Integer(0), 
		new Integer(0), 
		new Float(0),
		new Boolean(false)
		};
	
	/*
	 * Data models
	 * @TODO checks and other sensible things
	 */
	private Object[][] data; 
	private Object[][] filesDataModel;
	private TeaspoonMaskModel maskTracksModel;
	
	/* 
	 * Teaspoon analysis variables
	 */
	private File analysisDirectory;				// Working directory to trawl for input alignments
	private int alignmentLength;				// Maximum (equal!) alignment length
	private int bootstraps;						// Number of bootstrap replicates
	private int slidingWindowSize;				// Size of sliding window
	private ArrayList<TeaspoonMask> maskTracks;	// Array of masking tracks for alignment, each an int tuple of nucleotide positions 5'->3' e.g. [10,31]
	private File ancestralSequenceAlignment;	// Sequence alignment in FASTA of outgroup/ancestral
	private File[] mainSequenceAlignments;		// Sequence alignments in FASTA of main/focal files
	private Date[] samplingDates;				// Dates to associate with each alignment
	private int minReadCoverageDepth;			// Minimum read depth (# of reads) per nucleotide position required for analysis
	private double[][] binIntervals = {			// Intervals for the site frequency bins
			{0.0,0.15,0.75},
			{0.15,0.75,1.0}
		};			

	/*
	 * The main Parameters for any analysis - coulb be written as JSON etc eventually too
	 */
	private BhattAdaptationParameters parameters;
	
	/*
	 * GUI variables
	 */
	private boolean areAlignmentsValidated;		// Have sequence alignments in wd been parsed with equal lengths?
	private boolean doSlidingWindowAnalysis;	// Enable sliding window analysis?
	private boolean doBootstrappedAnalysis;		// Conduct bootstrap sampling?
	private boolean ancestralSequenceAlignmentSet; // Has an ancestral sequence alignment been specified?
	
	/**
	 * Default no-arg constructor
	 */
	public TeaspoonModel(){
		// Init parameters object
		this.parameters = new BhattAdaptationParameters();

		// Set up some sensible defaults
		Object[][] defaultData = {
				{new File("<filename>"),-1,-1,0.0,false}
		};
		this.setData(defaultData);
		this.bootstraps = 10;
		this.slidingWindowSize = 50;
		this.minReadCoverageDepth = 95;
		this.areAlignmentsValidated = false;
		this.doSlidingWindowAnalysis = false;
		this.doBootstrappedAnalysis = true;
		this.ancestralSequenceAlignmentSet = false;
		this.alignmentLength = -1;
		
		// pass those to parameters
		this.parameters.setBootstrapReplicates(bootstraps);
		this.parameters.setCustomBinSettings(binIntervals);
		
		// up the subsidiary masks table model
		this.maskTracksModel = new TeaspoonMaskModel();
		
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
	public ArrayList<TeaspoonMask> getMaskTracks() {
		return this.maskTracksModel.getMasks();
	}

	/**
	 * @return the sequenceAlignments
	 */
	public File[] getMainSequenceAlignments() {
		return mainSequenceAlignments;
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
	 * @return the double[2][3] representing site-frequency bin intervals
	 */
	public double[][] getCustomBinIntervals(){
		return binIntervals;
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
		this.parameters.setBootstrapReplicates(bootstraps);
	}

	/**
	 * @param slidingWindowSize the slidingWindowSize to set
	 */
	public void setSlidingWindowSize(int slidingWindowSize) {
		this.slidingWindowSize = slidingWindowSize;
	}


	/**
	 * Deprecated.
	 * Mask additions/removals should be handled by calls to the 
	 * TeaspoonMaskModel.addRow() and .removeRow() methods.
	 * @param sequenceAlignments the sequenceAlignments to set
	 */
	@Deprecated
	public void setMainSequenceAlignments(File[] sequenceAlignments) {
		//FIXME deprecated.
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
	public TeaspoonMaskModel getMaskTracksModel() {
		return this.maskTracksModel;
	}

	/**
	 * @param filesDataModel the filesDataModel to set
	 */
	public void setFilesDataModel(Object[][] filesDataModel) {
		this.filesDataModel = filesDataModel;
	}

	/**
	 * Returns the parameters as a snapshot, e.g the current 
	 * validated state of each, formatted as a BhattAdaptationParameters
	 * @return
	 */
	public BhattAdaptationParameters getParametersSnapshot() {
		// TODO Auto-generated method stub
		return this.parameters;
	}

	/* (non-Javadoc)
	 * @see javax.swing.table.TableModel#getRowCount()
	 */
	@Override
	public int getRowCount() {
		if(data != null){
			return getData().length;
		}else{
			return 0;
		}
	}

	/* (non-Javadoc)
	 * @see javax.swing.table.TableModel#getColumnCount()
	 */
	@Override
	public int getColumnCount() {
		return columnNames.length;
	}

	/* (non-Javadoc)
	 * @see javax.swing.table.TableModel#getValueAt(int, int)
	 */
	@Override
	public Object getValueAt(int rowIndex, int colIndex) {
		if(data != null){
			if(colIndex==0){
				// this is the first col e.g. file. toString() is too long..
				return ((File)data[rowIndex][colIndex]).getName();
			}else{
				return data[rowIndex][colIndex];
			}
		}else{
			return null;
		}
	}

	/*
	 * JTable uses this method to determine the default renderer/
	 * editor for each cell.  If we didn't implement this method,
	 * then the last column would contain text ("true"/"false"),
	 * rather than a check box.
	 */
	public Class getColumnClass(int c) {
		return getValueAt(0, c).getClass();
	}

	/**
	 * Get the definitions corresponding to each column in the data
	 * @return String[] of table definitions
	 */
	public String[] getTableColumnDefinitions() {
		// TODO Auto-generated method stub
		return this.columnDefinitions;
	}

	/*
	 * Ensures the column names are available.
	 * (non-Javadoc)
	 * @see javax.swing.table.AbstractTableModel#getColumnName(int)
	 */
	@Override
	public String getColumnName(int col) {
		return columnNames[col];
	}

	/*
	 * Don't need to implement this method unless your table's
	 * editable.
	 */
	public boolean isCellEditable(int row, int col) {
		//Note that the data/cell address is constant,
		//no matter where the cell appears onscreen.
		if (col < 2) {
			return false;
		} else {
			return true;
		}
	}

	/*
	 * Don't need to implement this method unless your table's
	 * data can change.
	 */
	public void setValueAt(Object value, int row, int col) {
	
		getData()[row][col] = value;
		fireTableCellUpdated(row, col);
	}

	/**
	 * @param newFileForRowData
	 */
	public void addFileRow(File newFileForRowData) {
		// first parse the new row
		BhattAdaptationFullSiteMatrix alignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(newFileForRowData).read());
		Object[] newRow = {
				(File)	newFileForRowData,
				(int) 	alignment.numberOfSequences(),
				(int) 	alignment.alignmentLength(),
				(float) 0.0,
				false};
		
		// now add to data table
		//if(data == null || ((String)data[0][0]).equals("<filename>") ){
		String filename = ((File)data[0][0]).getName();
		if(data == null || filename.equals("<filename>") ){
			// no data exists yet, create new
			data = new Object[1][TeaspoonModel.columnNames.length];
			data[0] = newRow;
		}else{
			// existing data, create a new matrix n+1 rows, k cols
			Object[][] newData = new Object[data.length+1][TeaspoonModel.columnNames.length];
			for(int i=0;i<data.length;i++){
				newData[i] = data[i];
			}
			newData[data.length] = newRow;
			data = newData;
			
		}
		this.fireTableDataChanged();
	}

	/**
	 * Updates the value of the ancestral alignment,
	 * and sets the main files list to be the remaining alignments
	 * @param selectedRow
	 * @throws FileNotFoundException 
	 */
	public void updateAncestralWithSelectedRow(int selectedRow) throws FileNotFoundException {
		// TODO Auto-generated method stub
		for(int i=0;i<data.length;i++){
			if(i==selectedRow){
				data[i][TeaspoonModel.columnNames.length-1] = true;
				this.ancestralSequenceAlignment = (File)data[i][0];
				this.parameters.setAncestralFile((File)data[i][0]);
				this.alignmentLength = (int) data[i][2];
				this.ancestralSequenceAlignmentSet = true;
			}else{
				data[i][TeaspoonModel.columnNames.length-1] = false;
			}
		}
		
		/* now update mainfiles list. 
		 * the mainfiles list will be every other file apart from
		 * the one currently set to ancestral
		 */
		File[] mainfiles = new File[this.data.length-1];
		int mainfileCount = 0;
		while(mainfileCount<mainfiles.length){
			for(Object[] row:data){
				// check whether this row is ancestral
				if(!(boolean)row[TeaspoonModel.columnNames.length-1]){
					mainfiles[mainfileCount] = (File) row[0];
					mainfileCount++;
				}
			}
		}
		this.mainSequenceAlignments = mainfiles;
		this.parameters.setInputFileList(mainfiles);
		
		// update table data
		this.fireTableDataChanged();
		
	}

	/**
	 * This returns the ancestral sequence alignment, assuming one has been set.
	 * @return the ancestralSequenceAlignment
	 */
	public File getAncestralSequenceAlignment() {
		return ancestralSequenceAlignment;
	}

	/**
	 * @param ancestralSequenceAlignment the ancestralSequenceAlignment to set
	 */
	public void setAncestralSequenceAlignment(File ancestralSequenceAlignment) {
		this.ancestralSequenceAlignment = ancestralSequenceAlignment;
	}

	/**
	 * @return true if an ancestral sequence alignment has been specified
	 */
	public boolean hasAncestralSequenceAlignmentBeenSet() {
		return this.ancestralSequenceAlignmentSet;
	}

	/**
	 * Checks to see whether all aligments in the datamodel are of identical length.
	 * Note this doesn't check for internal frameshifts/indels/other biological alignment
	 * errors - just that each has literally the same number of nucleotides.
	 * @return
	 */
	public boolean areAlignmentsOfEqualLength() {
		boolean equalLength = true;
		for(Object[] row:data){
			// get length of this row's alignment
			int length = (int)row[2];
			// true if and only if current equallength state true AND this row's length matches global alignmentLength
			equalLength = (equalLength && (length == this.alignmentLength) );
		}
		// TODO Auto-generated method stub
		return equalLength;
	}

	/**
	 * Adds a new TeaspoonMask to the internal TeaspoonMaskModel
	 * @param newMaskTrack
	 */
	public void addMaskRow(TeaspoonMask newMaskTrack) {
		// add the mask track
		this.maskTracksModel.addMaskRow(newMaskTrack);	
		// update my mask list
		this.maskTracks = this.maskTracksModel.getMasks();
	}

	/**
	 * @param selectedRow
	 */
	public void removeMaskWithSelectedRow(int selectedRow) {
		this.maskTracksModel.removeSelectedRow(selectedRow);
		
	}

	/**
	 * Deletes the selected row
	 * @param selectedRow
	 */
	public void removeFileWithSelectedRow(int selectedRow) {
		Object[][] newData = new Object[this.data.length-1][TeaspoonModel.columnNames.length];
		int newRowCounter = 0;
		for(int rowIndex=0;rowIndex<this.data.length;rowIndex++){
			if(rowIndex != selectedRow){
				newData[newRowCounter] = data[rowIndex];
				newRowCounter++;
			}
		}
		this.data = newData;
		this.fireTableDataChanged();
		
	}

	/**
	 * @param customBins
	 */
	public void setCustomBinIntervals(double[][] customBins) {
		this.parameters.setCustomBinSettings(customBins);
		this.binIntervals = customBins;
	}

}
