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
	String mainFile 		= "./HCV_data/main_HCVpacbio_filelist.edited.txt";
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
		fail("Not yet implemented"); // TODO
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
		double no_of_noneutral_sites = 446.815825;
		double precision = 0.1;
		
		bm = new BhattMethod(main, ans_tmp);
		try {
			bm.inferCountsFixedNR(bins, prior, true, Nvec, nr[0]);
		} catch (NullNeutralRatioException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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
				Math.abs( no_silent_sites_low	- bm.getSilentSubstitutionsCountArray()[0] ) < precision
					);
		assertTrue(
				Math.abs( no_replacement_sites_low - bm.getReplacementSubstitutionsCountArray()[0] ) < precision
					);
		assertTrue(
				Math.abs( r_low_ratio_s_low	- bm.getReplacementToSilentRatesRatio()[0] ) < precision
					);
		assertTrue(
				Math.abs( no_of_noneutral_sites - bm.getNonNeutralSubstitutions()[0] ) < precision
				);
	}

}
