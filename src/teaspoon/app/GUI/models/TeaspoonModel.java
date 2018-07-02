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
	 * The main Parameters for any analysis - coulb be written as JSON etc eventually too
	 */
	private BhattAdaptationParameters parameters;
	
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
		// Init parameters object
		this.parameters = new BhattAdaptationParameters();

		// Set up some sensible defaults
		Object[][] defaultData = {
				{new File("<filename>"),0,0,0.0,false}
		};
		this.setData(defaultData);
		this.bootstraps = 10;
		this.slidingWindowSize = 50;
		this.minReadCoverageDepth = 95;
		this.areAlignmentsValidated = false;
		this.doSlidingWindowAnalysis = false;
		this.doBootstrappedAnalysis = true;
		
		// pass those to parameters
		this.parameters.setBootstrapReplicates(bootstraps);
		
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
		this.parameters.setBootstrapReplicates(bootstraps);
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
	 * Updates the value of the ancestral alignment
	 * @param selectedRow
	 * @throws FileNotFoundException 
	 */
	public void updateAncestralWithSelectedRow(int selectedRow) throws FileNotFoundException {
		// TODO Auto-generated method stub
		for(int i=0;i<data.length;i++){
			if(i==selectedRow){
				data[i][TeaspoonModel.columnNames.length-1] = true;
				this.parameters.setAncestralFile((File)data[i][0]);
			}else{
				data[i][TeaspoonModel.columnNames.length-1] = false;
			}
		}
		this.fireTableDataChanged();
		
	}

}
