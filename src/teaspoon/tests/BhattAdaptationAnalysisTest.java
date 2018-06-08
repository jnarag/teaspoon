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

import teaspoon.app.BhattAdaptationAnalysis;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.BhattAdaptationResults;

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
	static File debugAncestralFile = new File("./H7_1stWave.fasta");
	static File debugMainFile = new File("./PRD_waves_year_W2.fasta");	
	
	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		BhattAdaptationParameters parameters = new BhattAdaptationParameters();
		// now populate the list
		try {
			parameters.setAncestralFile(debugAncestralFile);
			parameters.setInputFile(debugMainFile);
			parameters.setNeutralRate(0.7186788);
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
	public final void testRunWithFixedNR() {
		// run the thing
		new BhattAdaptationAnalysis(parameters).runWithFixedNR();
	}

	/**
	 * Test method for {@link teaspoon.app.BhattAdaptationAnalysis#runWithEstimatedNR()}.
	 */
	@Test
	public final void testRunWithEstimatedNR() {
		fail("Not yet implemented"); // TODO
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

}
