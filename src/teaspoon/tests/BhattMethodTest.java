/**
 * 
 */
package teaspoon.tests;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import teaspoon.adaptation.BhattMethod;
import teaspoon.app.utils.MainAlignmentParser;
import teaspoon.app.utils.NullNeutralRatioException;
import teaspoon.app.utils.TeaspoonMethods;

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
public class BhattMethodTest {

	BhattMethod bm;
	String mainFile 		= "./HCV_data/main_HCVpacbio_filelist.JUnit.txt";	// most tests
	String fixedFile		= "./HCV_data/sub_053/FP7_05302_0.3644.fasta"; 	// for fixed NR test only
	String ancestralFile 	= "./HCV_data/ancestral_HCVpacbio_filelist.edited.txt";
	int [][] main;
	int [] ans_tmp;
    double[] nr = {0.7186788};
    double[] L = {0.0, 0.15, 0.75};
    double[] H = {0.15, 0.75, 1.0};
    boolean[] Nvec = {false,true,false};
    double[] prior = {1.0,1.0,1.0,1.0};
	private double[][] bins;
	
	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		// Stuff for I/O
		BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
		BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));

		if (reader1.ready()) {

			String ancesfilename = reader1.readLine().trim();
			MainAlignmentParser ra = new MainAlignmentParser(ancesfilename);

			int[][] ancestralMatrix = ra.readFASTA();

			if (reader2.ready()) {
				String mainFile = reader2.readLine().trim();
				MainAlignmentParser r = new MainAlignmentParser(mainFile);
				main = r.readFASTA();
				while((main == null || main.length < 10) && reader2.ready()) {
					mainFile = reader2.readLine().trim();
					r = new MainAlignmentParser(mainFile);
					main = r.readFASTA();
				}
				ans_tmp = ra.consensusArray(ancestralMatrix);
			}
		}
		
		// Stuff for GeneAnalysis
        nr = new double[1];
        bins = new double[2][L.length];
 
        for(int i=0;i<L.length;i++){
            bins[0][i]=L[i];
            bins[1][i]=H[i];
        }


	}

	/**
	 * @throws java.lang.Exception
	 */
	@After
	public void tearDown() throws Exception {
	}

	/**
	 * Test method for {@link teaspoon.adaptation.BhattMethod#BhattMethod()}.
	 */
	@Test(expected=RuntimeException.class)
	public final void testBhattMethod() {
		bm = new BhattMethod();
	}

	/**
	 * Test method for {@link teaspoon.adaptation.BhattMethod#BhattMethod(int[][], int[])}.
	 */
	@Test
	public final void testBhattMethodIntArrayArrayIntArray() {
		bm = new BhattMethod(main, ans_tmp);
		//fail("Not yet implemented"); // TODO
		//TODO a more sensible test
	}

	/**
	 * Test method for {@link teaspoon.adaptation.BhattMethod#setNeutralRatio(double)}.
	 */
	@Test
	public final void testSetNeutralRatio() {
		bm = new BhattMethod(main, ans_tmp);
		bm.setNeutralRatio(0.4);
		//fail("Not yet implemented"); // TODO
	}

	/**
	 * Test method for {@link teaspoon.adaptation.BhattMethod#inferCountsEstimatedNR(double[][], double[], boolean, boolean[])}.
	 */
	@Test
	public final void testInferCountsEstimatedNR() {
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
		
		System.out.println("Test BhattMethod; expected counts tested to "+testPrecision+" precision:");
		// test with overloaded constructor which takes debug flag
		bm = new BhattMethod(main, ans_tmp, true);
		bm.inferCountsEstimatedNR(bins, prior, true, Nvec);
				
		// Get the debug (observed site) data to play with)
		double[][] observedData = bm.getObservedSiteDebugData();
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
				bm.getSilentSubstitutionsCountArray()[(int) 0] 		+ "," + 
				bm.getReplacementSubstitutionsCountArray()[(int) 0] + "," + 
				bm.getReplacementToSilentRatesRatio()[(int) 0] 		+ "," + 
				bm.getNeutralRatio()								+ "," +
				bm.getNonNeutralSubstitutions()[(int) 0]
						);
		// Formally test
		assertTrue(
				Math.abs( no_silent_sites_low	- bm.getSilentSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( no_replacement_sites_low - bm.getReplacementSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( r_low_ratio_s_low	- bm.getReplacementToSilentRatesRatio()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( estimated_neutral_ratio - bm.getNeutralRatio() ) < testPrecision
				);
		assertTrue(
				Math.abs( no_of_noneutral_sites - bm.getNonNeutralSubstitutions()[0] ) < (testPrecision * 10) // bigger precision allowed for this one.
				);
	}

	/**
	 * Test method for {@link teaspoon.adaptation.BhattMethod#inferCountsFixedNR(double[][], double[], boolean, boolean[], double)}.
	 */
	@Test
	public final void testInferCountsFixedNR() {
		/*
		 * Specs for 1-timepoint analysis on HCV first timepoint, lowFreq table data:
		 */
		double no_silent_sites_low	= 338.54969;
		double no_replacement_sites_low = 690.12431;
		double r_low_ratio_s_low	= 2.038472727;
		double no_of_noneutral_sites = 690.815825;
		// tolerable |(observed-expected)| for these tests:
		double testPrecision = 1.0;
		
		// fudge to use the original input file (./HCV_data/sub_053/FP7_05302_0.3644.fasta) for this test only
		main = new MainAlignmentParser(this.fixedFile).readFASTA();
		
		System.out.println("Test BhattMethod; expected counts tested to "+testPrecision+" precision:");
		// test with overloaded constructor which takes debug flag
		bm = new BhattMethod(main, ans_tmp, true);
		try {
			bm.inferCountsFixedNR(bins, prior, true, Nvec, nr[0]);
		} catch (NullNeutralRatioException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// Get the debug (observed site) data to play with)
		double[][] observedData = bm.getObservedSiteDebugData();
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
		
		// Output result as text 
		// NB return array index==0 as low frequency table
		System.out.println( 
				bm.getSilentSubstitutionsCountArray()[(int) 0] 		+ "," + 
				bm.getReplacementSubstitutionsCountArray()[(int) 0] + "," + 
				bm.getReplacementToSilentRatesRatio()[(int) 0] 		+ "," + 
				bm.getNonNeutralSubstitutions()[(int) 0]
						);
		// Formally test
		assertTrue(
				Math.abs( no_silent_sites_low	- bm.getSilentSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( no_replacement_sites_low - bm.getReplacementSubstitutionsCountArray()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( r_low_ratio_s_low	- bm.getReplacementToSilentRatesRatio()[0] ) < testPrecision
					);
		assertTrue(
				Math.abs( no_of_noneutral_sites - bm.getNonNeutralSubstitutions()[0] ) < testPrecision
				);
	}
}
