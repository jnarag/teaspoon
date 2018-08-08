/**
 * 
 */
package teaspoon.app.GUI.views;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.XYSeries.XYSeriesRenderStyle;
import org.knowm.xchart.style.Styler.ChartTheme;
import org.knowm.xchart.style.Styler.LegendPosition;
import org.knowm.xchart.style.markers.SeriesMarkers;

import teaspoon.adaptation.DataSet;


/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.

 * A JFrame holding plotting (scatter, histogram) info
 * leans a lot on https://github.com/lonelyjoeparker/qmul-genome-convergence-pipeline/blob/master/src/uk/ac/qmul/sbcs/evolution/convergence/gui/views/AnalysesView.java
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 24 Jul 2018
 * @version 0.1
 */
public class SimpleRegressionPlottingFrame extends JFrame{
	/**
	 * 
	 */
	private static final long serialVersionUID = 4331637424478573163L;
	JPanel mainPanel, infoPanel, statsTablePanel, optionsPanel, chartPanel;
	JLabel label;
	JScrollPane statsTableScrollPane, wholeViewScrollPane;
	String internalText = "Some plotting data.";
    String currentScatterSeriesName = null; // keep track of scatter data to add/remove series
    String currentRegressionName = null;	// keep track of regression to add/remove series
    XYChart scatterChart;
	XChartPanel scatterChartPanel;
    JCheckBox plotLogX, plotLogY, collectOverlay;
    boolean doPlotLogX, doPlotLogY, doCollectOverlay;
    
    /**
     * No-arg constructor.
     */
	public SimpleRegressionPlottingFrame(){
		super("Data plotting");
		mainPanel = new JPanel();
		mainPanel.setLayout(new BoxLayout(mainPanel,BoxLayout.Y_AXIS));
		infoPanel = new JPanel(new FlowLayout());
		statsTablePanel = new JPanel();
		chartPanel = new JPanel();
		label = new JLabel("<html><center>"+internalText+"</html>");
		optionsPanel = new JPanel(new FlowLayout());
		plotLogX = new JCheckBox("Log-plot X-axis");
		plotLogY = new JCheckBox("Log-plot Y-axis");
		collectOverlay = new JCheckBox("Overlay existing plots");
		optionsPanel.add(plotLogX);
		optionsPanel.add(plotLogY);
		//optionsPanel.add(collectOverlay);
		scatterChart = this.getScatterChart();
		scatterChartPanel = new XChartPanel<XYChart>(scatterChart);
		scatterChartPanel.setSize(650,250);
		scatterChartPanel.setPreferredSize(new Dimension(650, 250));
		chartPanel.add(scatterChartPanel);
		chartPanel.setPreferredSize(new Dimension(650,500));
		//infoPanel.add(label);
		infoPanel.add(optionsPanel);
		mainPanel.add(statsTablePanel);
		mainPanel.add(infoPanel);
		mainPanel.add(chartPanel);
		wholeViewScrollPane = new JScrollPane(mainPanel);
		wholeViewScrollPane.setPreferredSize(new Dimension(650,600));
		add(wholeViewScrollPane);
		setSize(650,600);
		setLocationRelativeTo(null);
		//setVisible(true);
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		// quickly create/add action listeners for each checkbox
		class PlotLogXActionListener implements ActionListener{
			@Override
			public void actionPerformed(ActionEvent arg0) {
				doPlotLogX = plotLogX.isSelected();	
				redrawScatterChartWithExistingData(scatterChart);
			}}
		plotLogX.addActionListener(new PlotLogXActionListener());
		
		class PlotLogYActionListener implements ActionListener{
			@Override
			public void actionPerformed(ActionEvent arg0) {
				doPlotLogY = plotLogY.isSelected();					
				redrawScatterChartWithExistingData(scatterChart);
			}}
		plotLogY.addActionListener(new PlotLogYActionListener());
		
		class CollectOverlayActionListener implements ActionListener{
			@Override
			public void actionPerformed(ActionEvent arg0) {
				doCollectOverlay = collectOverlay.isSelected();					
			}}
		collectOverlay.addActionListener(new CollectOverlayActionListener());
	}
	
	/**
	 * Update the label in the plotting window (mainly debug function)
	 * @param newContent
	 */
	public void updateLabelContent(String newContent){
		internalText = newContent;
		label.setText("<html><center>"+internalText+"</html>");
	}


	public void redrawScatterChartWithExistingData(XYChart existing){
		XYSeries some = existing.getSeriesMap().get(currentScatterSeriesName);
		// Use an iterator to find the min of X data
		Iterator<Float> xItr = (Iterator<Float>) some.getXData().iterator();
		Float xMin = xItr.next();
		while(xItr.hasNext()){
			xMin = Math.min(xMin, xItr.next());
		}
		// Use an iterator to find the min of Y data
		Iterator<Float> yItr = (Iterator<Float>) some.getYData().iterator();
		Float yMin = yItr.next();
		while(xItr.hasNext()){
			yMin = Math.min(yMin, yItr.next());
		}
		try {
			if(doPlotLogX && (xMin > 0)){
				scatterChart.getStyler().setXAxisLogarithmic(true);
			}else{
				scatterChart.getStyler().setXAxisLogarithmic(false);
			}
			if(doPlotLogY && (yMin > 0)){
				scatterChart.getStyler().setYAxisLogarithmic(true);
			}else{
				scatterChart.getStyler().setYAxisLogarithmic(false);
			}
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try{
			scatterChartPanel.repaint();
		}catch(Exception ex){
			scatterChart.getStyler().setXAxisLogarithmic(false);
			scatterChart.getStyler().setYAxisLogarithmic(false);
			scatterChartPanel.repaint();
			ex.printStackTrace();
		}
		
	}
	
	/**
	 * Update a chart with bivariate data
	 * @param name
	 * @param xyData
	 */
	public void updateScatterChart(String name, ArrayList<Float[]> xyData){
		// Create x and y lists for various purposes
		List<Float> xData = new ArrayList<Float>();
		List<Float> yData = new ArrayList<Float>();
		double[][] regressionData = new double[xyData.size()][2];
		DescriptiveStatistics xDataSeries = new DescriptiveStatistics();
		
		// laboriously add points (though we're unlikely to get more than 30 observations so who cares?)
		int validatedDataIndex = 0;
		Iterator<Float[]> dataIterator = xyData.iterator();
		while(dataIterator.hasNext()){
			Float[] dataPoint = dataIterator.next();
			// check for NaNs
			if(!dataPoint[0].isNaN() && !dataPoint[1].isNaN()){
				xData.add(dataPoint[0]);
				yData.add(dataPoint[1]);
				regressionData[validatedDataIndex] = new double[]{dataPoint[0],dataPoint[1]};
				xDataSeries.addValue(dataPoint[0]);
				// only increment destination data if !isNaN
				validatedDataIndex++;			
			}
		}

		// fit a linear least-sq regression
		SimpleRegression linearFit = new SimpleRegression();
		linearFit.addData(regressionData);

		// remove the old series and regression
		scatterChart.removeSeries(currentScatterSeriesName);
		scatterChart.removeSeries(currentRegressionName);
		
		// plot the new series
		currentScatterSeriesName = name;
		try {
			// scatter points
			scatterChart.addSeries(name, xData, yData);
			if(doPlotLogX && ((xDataSeries.getMin()) > 0)){
				scatterChart.getStyler().setXAxisLogarithmic(true);
			}else{
				scatterChart.getStyler().setXAxisLogarithmic(false);
			}
			if(doPlotLogY && (org.apache.commons.lang3.math.NumberUtils.min(ArrayUtils.toPrimitive(yData.toArray(new Float[yData.size()]))) > 0)){
				scatterChart.getStyler().setYAxisLogarithmic(true);
			}else{
				scatterChart.getStyler().setYAxisLogarithmic(false);
			}
			
			// regression line plotting, first get equation and min/max vals
			currentRegressionName = "y = "+linearFit.getIntercept()+" + "+linearFit.getSlope()+" x;\n"+
										"(N="+linearFit.getN()+", R-sq="+linearFit.getRSquare()+")";
			double predict_y_min = linearFit.predict(xDataSeries.getMin());
			double predict_y_max = linearFit.predict(xDataSeries.getMax());
			
			// add a regression line
			List<Double> lineXdata = Arrays.asList(new Double[] {xDataSeries.getMin(),xDataSeries.getMax()});
			List<Double> lineYdata = Arrays.asList(new Double[] {predict_y_min,predict_y_max});
			XYSeries regression = scatterChart.addSeries(currentRegressionName,lineXdata,lineYdata);
			regression.setXYSeriesRenderStyle(XYSeriesRenderStyle.Line);
			regression.setMarker(SeriesMarkers.NONE);

		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		scatterChartPanel.repaint();
	}
	
	/**
	 * Get an XChart scatter plot (XYChart)
	 * @return
	 */
	public XYChart getScatterChart(){

		// Create Chart
		currentScatterSeriesName = "Dummy data - run an analysis to plot.";
		currentRegressionName = "Regression";
		scatterChart = new XYChartBuilder().title("Scatterplot").theme(ChartTheme.GGPlot2).height(300).width(650).build();

		// Customize Chart
		scatterChart.getStyler().setDefaultSeriesRenderStyle(XYSeriesRenderStyle.Scatter);
		scatterChart.getStyler().setLegendPosition(LegendPosition.InsideN);

		// Series
		List<Float> xData = new ArrayList<Float>();
		List<Float> yData = new ArrayList<Float>();
	    List<Double> errorBars = new ArrayList<Double>();
		Random random = new Random();
		int size = 10;
		for (int i = 0; i < size; i++) {
			float nextRandom = random.nextFloat();
			xData.add((float) Math.pow(10, nextRandom * 4));
			float yVal = (float) (100.0 *( (nextRandom*10) + i));
			yData.add(yVal);
			errorBars.add(yVal * .3);
		}
		scatterChart.addSeries(currentScatterSeriesName, xData, yData, errorBars);
		
		// add a regression line
		List<Float> lineXdata = Arrays.asList(new Float[] {(float) 0.1,(float) 1000.0});
		List<Float> lineYdata = Arrays.asList(new Float[] {(float) 100,(float) 1000.0});
		XYSeries regression = scatterChart.addSeries(currentRegressionName,lineXdata,lineYdata);
		regression.setXYSeriesRenderStyle(XYSeriesRenderStyle.Line);
		regression.setMarker(SeriesMarkers.NONE);

		return scatterChart;
	}

}
