/**
 * 
 */
package teaspoon.tests;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import teaspoon.adaptation.BhattMethod;
import teaspoon.app.BhattAdaptationAnalysis;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.BhattAdaptationResults;
import teaspoon.app.utils.NullNeutralRatioException;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 8 Jun 2018
 * @version 0.1
 */
public class BhattAdaptationAnalysisTest {
	BhattAdaptationParameters parameters;
	static File debugAncestralFile = new File("./HCV_data/sub_053/FP7_05301_0.fasta");
	static File debugMainFile = new File("./HCV_data/sub_053/FP7_05302_0.3644.fasta");	
	
	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		parameters = new BhattAdaptationParameters();
		// now populate the list
		try {
			parameters.setAncestralFile(debugAncestralFile);
			parameters.setInputFile(debugMainFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * @throws java.lang.Exception
	 */
	@After
	public void tearDown() throws Exception {
	}

	/**
	 * Test method for {@link teaspoon.app.BhattAdaptationAnalysis#runWithFixedNR()}.
	 */
	@Test
	public final void testRunWithFixedNRcustomBins() {
		// parameters
		parameters.setDebugFlag(true);
		parameters.setNeutralRate(0.7186788);
		double[][] customBins = {
				{0.0, 0.00015, 0.9995},
				{0.00015, 0.9995, 1.0}
		};
		parameters.setCustomBinSettings(customBins);
		/*
		 * Specs for 1-timepoint analysis on HCV first timepoint, lowFreq table data:
		 */
		double no_silent_sites_low	= 338.54969;
		double no_replacement_sites_low = 690.12431;
		double r_low_ratio_s_low	= 2.038472727;
		double no_of_noneutral_sites = 446.67;
		// tolerable |(observed-expected)| for these tests:
		double testPrecision = 1.0;
		
		// Run the thing and get the output
		BhattAdaptationResults output = new BhattAdaptationAnalysis(parameters).runWithFixedNR();

		// Output result as text 
		output.printToText();
		
		// Detailed inspection based on the BhattMethod output
		BhattMethod bmOutput = output.getBhattSiteCounter();

		// Get the debug (observed site) data to play with)
		double[][] observedData = bmOutput.getObservedSiteDebugData();
		// Get the correct test data
		double[][] correctTestData = HCValignmentObservedSiteCounts.more_sites;
		// compare our output with test standard, first matrix dimensions
		assertTrue(observedData[0].length == correctTestData[0].length);
		assertTrue(observedData.length == correctTestData.length);
		// now compare valuewise
		for(int siteRow=0; siteRow<observedData.length; siteRow++){
			for(int siteCol=0; siteCol<observedData[0].length; siteCol++){
				assertTrue(observedData[siteRow][siteCol] == correctTestData[siteRow][siteCol]);
			}	
		}
		
		// NB return array index==0 as low frequency table
		// Formally test
		assertFalse(
				Math.abs( no_silent_sites_low	- bmOutput.getSilentSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertFalse(
				Math.abs( no_replacement_sites_low - bmOutput.getReplacementSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( r_low_ratio_s_low	- bmOutput.getReplacementToSilentRatesRatio()[0] ) < testPrecision
					);
		assertFalse(
				Math.abs( no_of_noneutral_sites - bmOutput.getNonNeutralSubstitutions()[0] ) < testPrecision
				);
	}

	/**
	 * Test method for {@link teaspoon.app.BhattAdaptationAnalysis#runWithEstimatedNR()}.
	 */
	@Test
	public final void testRunWithEstimatedNR() {
		// parameters
		parameters.setDebugFlag(true);
		double[][] customBins = {
				{0.0, 0.15, 0.75},
				{0.15, 0.75, 1.0}
		};
		parameters.setCustomBinSettings(customBins);
		try {
			parameters.setInputFile(new File("./HCV_data/sub_053/FP7_05305_1.3699.fasta"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		/*
		 * Specs for 1-timepoint analysis on HCV fourth timepoint, lowFreq table data:
		 */
		double no_silent_sites_low	= 383.5975838035207;
		double no_replacement_sites_low = 743.4424161964808;
		double r_low_ratio_s_low	= 1.9380789858605398;
		double no_of_noneutral_sites = 335.32347336899505; // probably not a good idea to test with this as it jumps around more than the others
		double estimated_neutral_ratio = 1.07526;

		// tolerable |(observed-expected)| for these tests:
		double testPrecision = 1.0;
		
		// Run the thing and get the output
		BhattAdaptationResults output = new BhattAdaptationAnalysis(parameters).runWithEstimatedNR();

		// Output result as text 
		output.printToText();
		
		// Detailed inspection based on the BhattMethod output
		BhattMethod bmOutput = output.getBhattSiteCounter();

		// Get the debug (observed site) data to play with)
		double[][] observedData = bmOutput.getObservedSiteDebugData();
		// Get the correct test data
		double[][] correctTestData = HCValignmentObservedSiteCounts.sites;
		// compare our output with test standard, first matrix dimensions
		assertTrue(observedData[0].length == correctTestData[0].length);
		assertTrue(observedData.length == correctTestData.length);
		// now compare valuewise
		for(int siteRow=0; siteRow<observedData.length; siteRow++){
			for(int siteCol=0; siteCol<observedData[0].length; siteCol++){
				assertTrue(observedData[siteRow][siteCol] == correctTestData[siteRow][siteCol]);
			}	
		}
		
		// Output result as text 
		// NB return array index==0 as low frequency table
		System.out.println( 
				bmOutput.getSilentSubstitutionsCountArray()[(int) 0] 		+ "," + 
				bmOutput.getReplacementSubstitutionsCountArray()[(int) 0] + "," + 
				bmOutput.getReplacementToSilentRatesRatio()[(int) 0] 		+ "," + 
				bmOutput.getNeutralRatio()								+ "," +
				bmOutput.getNonNeutralSubstitutions()[(int) 0]
						);
		// Formally test
		assertTrue(
				Math.abs( no_silent_sites_low	- bmOutput.getSilentSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( no_replacement_sites_low - bmOutput.getReplacementSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( r_low_ratio_s_low	- bmOutput.getReplacementToSilentRatesRatio()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( estimated_neutral_ratio - bmOutput.getNeutralRatio() ) < testPrecision
				);
		assertTrue(
				Math.abs( no_of_noneutral_sites - bmOutput.getNonNeutralSubstitutions()[0] ) < (testPrecision * 10) // bigger precision allowed for this one.
				);
	}

	/**
	 * Test method for {@link teaspoon.app.BhattAdaptationAnalysis#runWithFixedNRandBootstrap()}.
	 */
	@Test
	public final void testRunWithFixedNRandBootstrap() {
		fail("Not yet implemented"); // TODO
	}

	/**
	 * Test method for {@link teaspoon.app.BhattAdaptationAnalysis#runWithEstimatedNRandBootstrap()}.
	 */
	@Test
	public final void testRunWithEstimatedNRandBootstrap() {
		fail("Not yet implemented"); // TODO
	}

	/**
	 * Test method for {@link teaspoon.app.BhattAdaptationAnalysis#getAnalysisParameters()}.
	 */
	@Test
	public final void testGetAnalysisParameters() {
		fail("Not yet implemented"); // TODO
	}

	/**
	 * Test method for {@link teaspoon.app.BhattAdaptationAnalysis#setAnalysisParameters(teaspoon.app.utils.BhattAdaptationParameters)}.
	 */
	@Test
	public final void testSetAnalysisParameters() {
		fail("Not yet implemented"); // TODO
	}

	/**
	 * Test method for {@link teaspoon.app.BhattAdaptationAnalysis#runWithEstimatedNR()}.
	 */
	@Test
	public final void testRunWithEstimatedNRcustomBins() {
		// parameters
		parameters.setDebugFlag(true);
		double[][] customBins = {
				{0.0, 0.0015, 0.995},
				{0.0015, 0.995, 1.0}
		};
		parameters.setCustomBinSettings(customBins);
		try {
			parameters.setInputFile(new File("./HCV_data/sub_053/FP7_05305_1.3699.fasta"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
		/*
		 * Specs for 1-timepoint analysis on HCV fourth timepoint, lowFreq table data:
		 */
		double no_silent_sites_low	= 383.5975838035207;
		double no_replacement_sites_low = 743.4424161964808;
		double r_low_ratio_s_low	= 1.9380789858605398;
		double no_of_noneutral_sites = 335.32347336899505; // probably not a good idea to test with this as it jumps around more than the others
		double estimated_neutral_ratio = 1.07526;
	
		// tolerable |(observed-expected)| for these tests:
		double testPrecision = 1.0;
		
		// Run the thing and get the output
		BhattAdaptationResults output = new BhattAdaptationAnalysis(parameters).runWithEstimatedNR();
	
		// Output result as text 
		output.printToText();
		
		// Detailed inspection based on the BhattMethod output
		BhattMethod bmOutput = output.getBhattSiteCounter();
	
		// Get the debug (observed site) data to play with)
		double[][] observedData = bmOutput.getObservedSiteDebugData();
		// Get the correct test data
		double[][] correctTestData = HCValignmentObservedSiteCounts.sites;
		// compare our output with test standard, first matrix dimensions
		assertTrue(observedData[0].length == correctTestData[0].length);
		assertTrue(observedData.length == correctTestData.length);
		// now compare valuewise
		for(int siteRow=0; siteRow<observedData.length; siteRow++){
			for(int siteCol=0; siteCol<observedData[0].length; siteCol++){
				assertTrue(observedData[siteRow][siteCol] == correctTestData[siteRow][siteCol]);
			}	
		}
		
		// Output result as text 
		// NB return array index==0 as low frequency table
		System.out.println( 
				bmOutput.getSilentSubstitutionsCountArray()[(int) 0] 		+ "," + 
				bmOutput.getReplacementSubstitutionsCountArray()[(int) 0] + "," + 
				bmOutput.getReplacementToSilentRatesRatio()[(int) 0] 		+ "," + 
				bmOutput.getNeutralRatio()								+ "," +
				bmOutput.getNonNeutralSubstitutions()[(int) 0]
						);
		// Formally test
		assertFalse(
				Math.abs( no_silent_sites_low	- bmOutput.getSilentSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertFalse(
				Math.abs( no_replacement_sites_low - bmOutput.getReplacementSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( r_low_ratio_s_low	- bmOutput.getReplacementToSilentRatesRatio()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( estimated_neutral_ratio - bmOutput.getNeutralRatio() ) < testPrecision
				);
		assertFalse(
				Math.abs( no_of_noneutral_sites - bmOutput.getNonNeutralSubstitutions()[0] ) < (testPrecision * 10) // bigger precision allowed for this one.
				);
	}

	/**
	 * Test method for {@link teaspoon.app.BhattAdaptationAnalysis#runWithFixedNR()}.
	 */
	@Test
	public final void testRunWithFixedNR() {
		// parameters
		parameters.setDebugFlag(true);
		parameters.setNeutralRate(0.7186788);
		double[][] customBins = {
				{0.0, 0.15, 0.75},
				{0.15, 0.75, 1.0}
		};
		parameters.setCustomBinSettings(customBins);
		/*
		 * Specs for 1-timepoint analysis on HCV first timepoint, lowFreq table data:
		 */
		double no_silent_sites_low	= 338.54969;
		double no_replacement_sites_low = 690.12431;
		double r_low_ratio_s_low	= 2.038472727;
		double no_of_noneutral_sites = 446.67;
		// tolerable |(observed-expected)| for these tests:
		double testPrecision = 1.0;
		
		// Run the thing and get the output
		BhattAdaptationResults output = new BhattAdaptationAnalysis(parameters).runWithFixedNR();
	
		// Output result as text 
		output.printToText();
		
		// Detailed inspection based on the BhattMethod output
		BhattMethod bmOutput = output.getBhattSiteCounter();
	
		// Get the debug (observed site) data to play with)
		double[][] observedData = bmOutput.getObservedSiteDebugData();
		// Get the correct test data
		double[][] correctTestData = HCValignmentObservedSiteCounts.more_sites;
		// compare our output with test standard, first matrix dimensions
		assertTrue(observedData[0].length == correctTestData[0].length);
		assertTrue(observedData.length == correctTestData.length);
		// now compare valuewise
		for(int siteRow=0; siteRow<observedData.length; siteRow++){
			for(int siteCol=0; siteCol<observedData[0].length; siteCol++){
				assertTrue(observedData[siteRow][siteCol] == correctTestData[siteRow][siteCol]);
			}	
		}
		
		// NB return array index==0 as low frequency table
		// Formally test
		assertTrue(
				Math.abs( no_silent_sites_low	- bmOutput.getSilentSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( no_replacement_sites_low - bmOutput.getReplacementSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( r_low_ratio_s_low	- bmOutput.getReplacementToSilentRatesRatio()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( no_of_noneutral_sites - bmOutput.getNonNeutralSubstitutions()[0] ) < testPrecision
				);
	}

}
