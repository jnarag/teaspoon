/**
 * 
 */
package teaspoon.tests;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import teaspoon.app.TeaspoonMask;
import teaspoon.app.standalone.TeaspoonCommandLineApp;
import teaspoon.app.standalone.TeaspoonMaskFactory;
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
 * @since 18 Jun 2018
 * @version 0.1
 */
public class TeaspoonCommandLineAppTest extends TeaspoonCommandLineApp {

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
	}

	/**
	 * @throws java.lang.Exception
	 */
	@After
	public void tearDown() throws Exception {
	}

	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonCommandLineApp#TeaspoonCommandLineApp(teaspoon.app.utils.BhattAdaptationParameters)}.
	 */
	@Test
	public final void testTeaspoonCommandLineAppBhattAdaptationParameters() {
		/*
		 * Set up an analysis to run as with the 1-timepoint analysis:
		 * (old options:
		 * ./HCV_data/ancestral_HCVpacbio_filelist.edited.txt ./HCV_data/main_HCVpacbio_filelist.edited.txt one)
		 * 
		 * ancestor:
		 * 	./HCV_data/sub_053/FP7_05301_0.fasta
		 * main:
		 * 	./HCV_data/sub_053/FP7_05302_0.3644.fasta
			./HCV_data/sub_053/FP7_05303_0.6137.fasta
			./HCV_data/sub_053/FP7_05304_0.8438.fasta
			./HCV_data/sub_053/FP7_05305_1.3699.fasta
			./HCV_data/sub_053/FP7_05306_1.7836.fasta
			./HCV_data/sub_053/FP7_05307_3.8986.fasta
			./HCV_data/sub_053/FP7_05308_6.8429.fasta
			./HCV_data/sub_053/FP7_05309_7.6849.fasta

		 * partition mask:
		 * 	./HCV_data/sub_053/mask
		 * rate:
		 * 	0.7186788
		 */
		double ratio = 0.7186788;
		File maskFile = new File("./HCV_data/sub_053/mask");
		File input = new File("./HCV_data/sub_053/FP7_05301_0.fasta");
		File[] inputList = new File[8];
		inputList[0] = new File("./HCV_data/sub_053/FP7_05302_0.3644.fasta");
		inputList[1] = new File("./HCV_data/sub_053/FP7_05303_0.6137.fasta");
		inputList[2] = new File("./HCV_data/sub_053/FP7_05304_0.8438.fasta");
		inputList[3] = new File("./HCV_data/sub_053/FP7_05305_1.3699.fasta");
		inputList[4] = new File("./HCV_data/sub_053/FP7_05306_1.7836.fasta");
		inputList[5] = new File("./HCV_data/sub_053/FP7_05307_3.8986.fasta");
		inputList[6] = new File("./HCV_data/sub_053/FP7_05308_6.8429.fasta");
		inputList[7] = new File("./HCV_data/sub_053/FP7_05309_7.6849.fasta");
		BhattAdaptationFullSiteMatrix alignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(input).readFASTA());
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> masks = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		int[] maskStartEnd = {0, alignment.alignmentLength()-1};
		ArrayList<int[]> maskRanges = new ArrayList<int[]>();
		maskRanges.add(maskStartEnd);
		masks.put(maskRanges,RateEstimationBehaviour.NEUTRAL_RATE_AVERAGED);
		masks.put(maskRanges,RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED);
		masks.put(maskRanges,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
	try {
			
			TeaspoonMaskFactory.writeMaskFileWithFixedRatio(maskFile, masks, alignment.alignmentLength(), ratio);
			/*
			 * TODO 	!!!
			 * FIXME 	!!!
			 * 
			 * Oh dear - because we're using the positions as the key we can't add more than one masking list at
			 * a time to the factory argument.
			 * 
			 * This is a pain in the arse.
			 * Either find another way to pass key-value or meh..
			 * 
			masks.put(maskRanges,RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED);
			TeaspoonMaskFactory.appendToMaskFileWithFixedRatio(maskFile, masks, alignment.alignmentLength(), ratio);
			masks.put(maskRanges,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
			TeaspoonMaskFactory.appendToMaskFileWithFixedRatio(maskFile, masks, alignment.alignmentLength(), ratio);
			 * 
			 */
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		BhattAdaptationParameters parameters = new BhattAdaptationParameters();
		try {
			parameters.setAncestralFile(input);
			parameters.setMaskFile(maskFile);
			parameters.setInputFileList(inputList);
			parameters.setBootstrapReplicates(10);
			parameters.setNeutralRate(ratio);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			new TeaspoonCommandLineApp(parameters);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		fail("Not yet implemented"); // TODO
	}

	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonCommandLineApp#main(java.lang.String[])}.
	 */
	@Test
	public final void testMain() {
		fail("Not yet implemented"); // TODO
	}

}
