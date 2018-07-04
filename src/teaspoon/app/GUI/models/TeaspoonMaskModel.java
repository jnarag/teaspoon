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

import teaspoon.app.TeaspoonMask;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.MainAlignmentParser;
import teaspoon.app.utils.RateEstimationBehaviour;

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
public class TeaspoonMaskModel extends AbstractTableModel{

	// JTable needed stuff - column names
	private final static String[] columnNames = {			
							//column index
		"Mask",				//0	TeaspoonMask
		"start",			//1	int
		"end",				//2	int
		"length",			//3	int
		"behaviour",		//4 RateEstimationBehaviour
		"neutral ratio",	//5	Double
		"do sliding-window?"//6	boolean
	};				
	// JTable needed studd - column definitions
	private final static String[] columnDefinitions = {	
													//column index
		"Mask hashcode",							//0	TeaspoonMask
		"mask start (index to 0)",					//1	int
		"mask end (index to 0)",					//2	int
		"length (nt)",								//3	int
		"ratio estimation method",					//4 RateEstimationBehaviour
		"neutral ratio",							//5	Double
		"which alignment is considered ancestral"	//6	bool
	};
	// Column indices of Integers
	private final static Integer[] integerIndices = new Integer[]{1,2,3};
	// Column indices of Floats
	private final static Integer[] doubleIndices = new Integer[]{5};
	// Default values for (hopefully) sizing the table, etc
	public final Object[] longValues = {
		null,
		new Integer(0), 
		new Integer(0), 
		new Integer(0), 
		RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED, 
		new Double(0),
		new Boolean(false)
		};
	
	/*
	 * Data models
	 * @TODO checks and other sensible things
	 */
	private Object[][] data; 
	
	/* 
	 * Teaspoon analysis variables
	 */
	private int alignmentLength;				// Maximum (equal!) alignment length
	private int slidingWindowSize;				// Size of sliding window
	
	/*
	 * GUI variables
	 */
	private boolean doSlidingWindowAnalysis;	// Enable sliding window analysis?
	
	/**
	 * Default no-arg constructor
	 */
	public TeaspoonMaskModel(){

		// Set up some sensible defaults
		Object[][] defaultData = {
				{
					new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED, null, 1),
					new Integer(0), 
					new Integer(0), 
					new Integer(0), 
					RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED, 
					new Double(0),
					new Boolean(false)
					}
		};
		this.setData(defaultData);
		this.slidingWindowSize = 50;
		this.doSlidingWindowAnalysis = false;
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
	 * @return the alignmentLength
	 */
	public int getAlignmentLength() {
		return alignmentLength;
	}


	/**
	 * @return the slidingWindowSize
	 */
	public int getSlidingWindowSize() {
		return slidingWindowSize;
	}

	/**
	 * @return the doSlidingWindowAnalysis
	 */
	public boolean isDoSlidingWindowAnalysis() {
		return doSlidingWindowAnalysis;
	}

	/**
	 * @param alignmentLength the alignmentLength to set
	 */
	public void setAlignmentLength(int alignmentLength) {
		this.alignmentLength = alignmentLength;
	}


	/**
	 * @param slidingWindowSize the slidingWindowSize to set
	 */
	public void setSlidingWindowSize(int slidingWindowSize) {
		this.slidingWindowSize = slidingWindowSize;
	}


	/**
	 * @param doSlidingWindowAnalysis the doSlidingWindowAnalysis to set
	 */
	public void setDoSlidingWindowAnalysis(boolean doSlidingWindowAnalysis) {
		this.doSlidingWindowAnalysis = doSlidingWindowAnalysis;
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
				return ((TeaspoonMask)data[rowIndex][colIndex]).hashCode();
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

	public ArrayList<TeaspoonMask> getMasks(){
		ArrayList<TeaspoonMask> masksList = new ArrayList<TeaspoonMask>();
		for(Object[] row:data){
			masksList.add((TeaspoonMask)row[0]);
		}
		return masksList;
	}
	/*
	 * Don't need to implement this method unless your table's
	 * editable.
	 */
	public boolean isCellEditable(int row, int col) {
		//Note that the data/cell address is constant,
		//no matter where the cell appears onscreen.
		//if (col < 2) {
			return false;
		//} else {
		//	return true;
		//}
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
	 * @param newMaskTrack
	 */
	public void addMaskRow(TeaspoonMask newMaskTrack) {
		// first set up the new row
		Object[] newRow = {
				newMaskTrack,
				(int) 	newMaskTrack.getFirstStart(),
				(int) 	newMaskTrack.getLastEnd(),
				(int) 	newMaskTrack.getLength(),
				(RateEstimationBehaviour) newMaskTrack.estimationBehaviour,
				(double)newMaskTrack.getNeutralRatio(),
				false};
		
		// now add to data table
		//if(data == null || ((String)data[0][0]).equals("<filename>") ){
		TeaspoonMask dummy = (TeaspoonMask)data[0][0];
		if(data == null || dummy.equals(new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED, null, 1)) ){
			// no data exists yet, create new
			data = new Object[1][TeaspoonMaskModel.columnNames.length];
			data[0] = newRow;
		}else{
			// existing data, create a new matrix n+1 rows, k cols
			Object[][] newData = new Object[data.length+1][TeaspoonMaskModel.columnNames.length];
			for(int i=0;i<data.length;i++){
				newData[i] = data[i];
			}
			newData[data.length] = newRow;
			data = newData;
			
		}
		this.fireTableDataChanged();
	}

	/**
	 * Deletes the selected mask row
	 * @param selectedRow
	 */
	public void removeSelectedRow(int selectedRow) {
		Object[][] newData = new Object[this.data.length-1][TeaspoonMaskModel.columnNames.length];
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

}
