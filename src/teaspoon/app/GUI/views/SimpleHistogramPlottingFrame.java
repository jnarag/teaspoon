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
import java.util.TreeMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;

import org.apache.commons.lang3.ArrayUtils;
import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategoryChartBuilder;
import org.knowm.xchart.Histogram;
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.style.Styler.ChartTheme;
import org.knowm.xchart.style.Styler.LegendPosition;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 25 Jul 2018
 * @version 0.1
 */
/**
 * A JFrame holding plotting (scatter, histogram) info
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 3 Jul 2017
 * @version 0.1
 */
public class SimpleHistogramPlottingFrame extends JFrame{
	/**
	 * 
	 */
	private static final long serialVersionUID = 4331637424478573163L;
	JPanel mainPanel, infoPanel, statsTablePanel, optionsPanel, chartPanel;
	JLabel label;
	JTable statsTable;
	JScrollPane statsTableScrollPane, wholeViewScrollPane;
	String internalText = "Some plotting data.";
    String currentHistogramSeriesName = null;
	CategoryChart histogramChart;
	XChartPanel histogramChartPanel;
    JCheckBox plotLogX, plotLogY, collectOverlay;
    public SimpleHistogramPlottingFrame(){
		super("site-frequency histogram plotting");
		mainPanel = new JPanel();
		mainPanel.setLayout(new BoxLayout(mainPanel,BoxLayout.Y_AXIS));
		//chartPanel = new JPanel(new GridLayout(2,1));
		chartPanel = new JPanel();
		label = new JLabel("<html><center>"+internalText+"</html>");
		histogramChart = this.getHistogramChart();
		histogramChartPanel = new XChartPanel<CategoryChart>(histogramChart);
		histogramChartPanel.setSize(650,250);
		histogramChartPanel.setPreferredSize(new Dimension(650, 250));
		chartPanel.add(histogramChartPanel);
		chartPanel.setPreferredSize(new Dimension(650,250));
		//infoPanel.add(label);
		mainPanel.add(chartPanel);
		wholeViewScrollPane = new JScrollPane(mainPanel);
		wholeViewScrollPane.setPreferredSize(new Dimension(650,300));
		add(wholeViewScrollPane);
		setSize(650,300);
		setLocationRelativeTo(null);
		setVisible(true);
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	}
	
	/**
	 * Update the label in the plotting window (mainly debug function)
	 * @param newContent
	 */
	public void updateLabelContent(String newContent){
		internalText = newContent;
		label.setText("<html><center>"+internalText+"</html>");
	}

	/**
	 * Update a chart with univariate data as histogram
	 * @param name
	 * @param spectrumPlottingData
	 */
	public void updateHistogram(String name, ArrayList<Float[]> spectrumPlottingData){
		histogramChart.removeSeries(currentHistogramSeriesName);
		currentHistogramSeriesName = name;
		//int maxXdataOOM = (int)Math.log10(org.apache.commons.lang3.math.NumberUtils.max(ArrayUtils.toPrimitive(spectrumPlottingData)));
		try {
			// arrgh need to make sure the bins are sorted
			TreeMap<Float,Float> sortedBins = new TreeMap<Float,Float>();
			for(Float[] binData:spectrumPlottingData){
				sortedBins.put(binData[0], binData[1]);
			}
			// values of sortedBins *should* be sorted in natural order. add to counts with iterator
			ArrayList<Float> counts = new ArrayList<Float>();
			Iterator<Float> binIterator = sortedBins.values().iterator();
			while(binIterator.hasNext()){
				counts.add(binIterator.next());
			}
			// OK now plot it
			Histogram histogram = new Histogram(counts, spectrumPlottingData.size(), 0.0, 1.0);
			histogramChart.addSeries(name, histogram.getxAxisData(),histogram.getyAxisData());
			
			// try and get the x-axis labels' precision right...
			histogramChart.getStyler().setXAxisDecimalPattern("#.#"); // just hardcode it; we know it's 0.0->1.0
			/*
			if(maxXdataOOM > 0){
				// bigger than 10^1 so safe to ignore d.p.
			    histogramChart.getStyler().setXAxisDecimalPattern("#");
			}else{
				// 10^0 or smaller, fiddle precision on x-axis labels
				switch(maxXdataOOM){
				case(0):{
					histogramChart.getStyler().setXAxisDecimalPattern("#.#");
				}
				case(-1):{
					histogramChart.getStyler().setXAxisDecimalPattern("#.#");
				}
				case(-2):{
					histogramChart.getStyler().setXAxisDecimalPattern("#.#");
				}
				case(-3):{
					histogramChart.getStyler().setXAxisDecimalPattern("#.##");
				}
				case(-4):{
					histogramChart.getStyler().setXAxisDecimalPattern("#.##");
				}
				default:{
					histogramChart.getStyler().setXAxisDecimalPattern("#.###");
				}
				}
				
			}
			*/
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		histogramChartPanel.repaint();
	}

	public CategoryChart getHistogramChart(){
		/**
		 * Internal class to generate data according to a normalish distribution
		 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
		 * @since 4 Jul 2017
		 * @version 0.1
		 */
		final class DataGenerator{
		  private List<Double> getGaussianData(int count) {

			    List<Double> data = new ArrayList<Double>(count);
			    Random r = new Random();
			    for (int i = 0; i < count; i++) {
			      data.add(r.nextGaussian() * 10);
			    }
			    return data;
			  }
		}
		
		// Create Chart
	    CategoryChart histogramPlaceHolder = new CategoryChartBuilder().title("Histogram").theme(ChartTheme.GGPlot2).height(300).width(650).build();
	    // Customize Chart
	    histogramPlaceHolder.getStyler().setLegendPosition(LegendPosition.InsideNW);
	    histogramPlaceHolder.getStyler().setAvailableSpaceFill(.96);
	    histogramPlaceHolder.getStyler().setOverlapped(true);
	    histogramPlaceHolder.getStyler().setXAxisDecimalPattern("#");
	    Histogram histogram = new Histogram(new DataGenerator().getGaussianData(100), 20, -20, 20);
	    String dummyName = "Dummy data - load alignments and select site-freq to plot.";
	    currentHistogramSeriesName = dummyName;
	    histogramPlaceHolder.addSeries(dummyName, histogram.getxAxisData(), histogram.getyAxisData());
	    return histogramPlaceHolder;			
	}
}
