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
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.XYSeries.XYSeriesRenderStyle;
import org.knowm.xchart.style.Styler.ChartTheme;
import org.knowm.xchart.style.Styler.LegendPosition;


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
    String currentScatterSeriesName = null;
    XYChart scatterChart;
	XChartPanel scatterChartPanel;
    JCheckBox plotLogX, plotLogY, collectOverlay;
    boolean doPlotLogX, doPlotLogY, doCollectOverlay;
    
    
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
		setVisible(true);
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
		// Create an x and y list
		Float[] xData = new Float[xyData.size()];
		Float[] yData = new Float[xyData.size()];
		int dataIndex = 0;
		Iterator<Float[]> dataIterator = xyData.iterator();
		while(dataIterator.hasNext()){
			Float[] dataPoint = dataIterator.next();
			xData[dataIndex] = dataPoint[0];
			yData[dataIndex] = dataPoint[1];
			dataIndex++;
		}

		scatterChart.removeSeries(currentScatterSeriesName);
		currentScatterSeriesName = name;
		try {
			scatterChart.addSeries(name, Arrays.asList(xData), Arrays.asList(yData));
			if(doPlotLogX && (org.apache.commons.lang3.math.NumberUtils.min(ArrayUtils.toPrimitive(xData)) > 0)){
				scatterChart.getStyler().setXAxisLogarithmic(true);
			}else{
				scatterChart.getStyler().setXAxisLogarithmic(false);
			}
			if(doPlotLogY && (org.apache.commons.lang3.math.NumberUtils.min(ArrayUtils.toPrimitive(yData)) > 0)){
				scatterChart.getStyler().setYAxisLogarithmic(true);
			}else{
				scatterChart.getStyler().setYAxisLogarithmic(false);
			}
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
		scatterChart = new XYChartBuilder().title("Scatterplot").theme(ChartTheme.GGPlot2).height(300).width(650).build();

		// Customize Chart
		scatterChart.getStyler().setDefaultSeriesRenderStyle(XYSeriesRenderStyle.Scatter);
		scatterChart.getStyler().setLegendPosition(LegendPosition.InsideN);

		// Series
		List<Float> xData = new ArrayList<Float>();
		List<Float> yData = new ArrayList<Float>();
		Random random = new Random();
		int size = 10;
		for (int i = 0; i < size; i++) {
			float nextRandom = random.nextFloat();
			xData.add((float) Math.pow(10, nextRandom * 10));
			yData.add((float) (1000000000.0 + nextRandom));
		}
		scatterChart.addSeries(currentScatterSeriesName, xData, yData);

		return scatterChart;
	}

}
