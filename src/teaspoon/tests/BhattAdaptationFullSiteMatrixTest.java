/**
 * 
 */
package teaspoon.tests;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import teaspoon.app.TeaspoonBootstrap;
import teaspoon.app.TeaspoonMask;
import teaspoon.app.standalone.TeaspoonBootstrapFactory;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;
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
public class BhattAdaptationFullSiteMatrixTest {

	BhattAdaptationFullSiteMatrix mainPartition;
	File input;
	
	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		input = new File("./HCV_data/sub_053/FP7_05305_1.3699.fasta");
		mainPartition = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(input).readFASTA());
	}

	/**
	 * @throws java.lang.Exception
	 */
	@After
	public void tearDown() throws Exception {
	}


	/**
	 * Test method for {@link teaspoon.app.utils.BhattAdaptationFullSiteMatrix#loadAlignmentFile(java.io.File)}.
	 */
	@Test
	public final void testLoadAlignmentFile() {
		int[] testStart = {4,4,3,3,2,1,2,4,3,2,4,2,4,2,3,4,3,2,2,4,3,1,2,4,3,4,4};
		int[] result = new int[testStart.length];
		for(int i=0;i<testStart.length;i++){
			result[i] = mainPartition.getSiteMatrix()[0][i];
		}
		if(!TeaspoonBootstrapTest.compareIntArraysByValue(testStart, result)){
			fail("input as parsed doesn't match expect");
		}
	}

	/**
	 * Test method for {@link teaspoon.app.utils.BhattAdaptationFullSiteMatrix#subsampleByMask(teaspoon.app.TeaspoonMask)}.
	 */
	@Test
	public final void testSubsampleByMask() {
		/*
		 * Test specs
		 * 
		 * input:1	[4,4,3,3,2,1,2,4,3,2,4,2,4,2,3,4,3,2,2,4,3,1,2,4,3,4,4, ]
		 * input:n	[4,4,3,3,2,1,2,4,3,2,4,2,4,2,3,4,3,2,2,4,3,1,2,4,3,4,4, ]
		 * mask_mid		[T,f,f,T,T,T,f,f,f,T,f,f,T,T,T,f,f,f,T,f,f,T,T,T,f,f,f
		 * output:1	[4,    3,2,1,      2,    4,2,3,      2,    1,2,4,     ]
		 * 
		 */
		boolean [] maskPositions = {
				true,false,false,true,true,true,false,false,false,
				true,false,false,true,true,true,false,false,false,
				true,false,false,true,true,true,false,false,false};
		TeaspoonMask mask = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,maskPositions);
		int[][] subsampledPartition = mainPartition.subsampleByMask(mask).getSiteMatrix();
		int[] testResult = {4, 3, 2, 1, 2, 4, 2, 3, 2, 1, 2, 4};
		if(!TeaspoonBootstrapTest.compareIntArraysByValue(testResult, subsampledPartition[0])){
			fail("mask_mid positions do not match"); 
		}
		if(!TeaspoonBootstrapTest.compareIntArraysByValue(testResult, subsampledPartition[subsampledPartition.length-1])){
			fail("mask_mid positions do not match"); 
		}
	}

	/**
	 * Test method for {@link teaspoon.app.utils.BhattAdaptationFullSiteMatrix#deriveConsensus()}.
	 */
	@Test
	public final void testDeriveConsensus() {
		int[] ancestor = mainPartition.deriveConsensus();
		int[] testResult = {4, 4, 3, 3, 2, 1, 2, 4, 3, 2, 4};
		for(int i=0;i<10;i++){
			if(ancestor[i] != testResult[i]){
				fail("derived consensus does not match expected");
			}
		}
	}

	/**
	 * Test method for {@link teaspoon.app.utils.BhattAdaptationFullSiteMatrix#subsampleByBootstrap(teaspoon.app.TeaspoonBootstrap)}.
	 */
	@Test
	public final void testSubsampleByBoostrap() {
		/*
		 * Test specs
		 * 
		 * input:1	[4,4,3,3,2,1,2,4,3,2,4,2,4,2,3,4,3,2,2,4,3,1,2,4,3,4,4, 2,2,2,3,2, 2, 4, 2, 3, 3, 2, 2, 4, 1, 2, 1, 1, 2, 4, 1, 4, 2, 3, 2, 1, 1, 2, 1, 3, 2, 4, 2, 3, 3, 3, 4, 3, 4, 2, 4, 1, 2, 2, 1, 4, 3, 4, 2, 1, 2, 2, 1, 1, 2, 3, 1, 2, 4, 3, 2, 2, 2, 3, 1, 1, 2, 4, 2, 3, 1, 3, 2, 1, 4, 1, 3, 4, 2, 4, 1, 4, 3, 1, 1, 3, 2, 2, 3, 1, 4, 2, 1, 2, 2, 1, 2, 1, 4, 2, 4, 4, 3, 2, 1, 2, 2, 4, 2, 2, 2, 3, 3, 3, 1, 4, 3, 2, 3, 4, 4, 2, 2, 4, 4, 3, 2, 3, 4, 3, 1, 3, 3, 3, 4, 2, 3, 3, 3, 1, 1, 4, 1, 2, 3, 4, 2, 4, 1, 2, 1, 4, 3, 2, 4, 3, 3, 3, 4, 3, 4, 2, 4, 2, 4, 2, 1, 2, 4, 2, 2, 2, 1, 2, 4, 3, 4, 3, 3, 2, 4, 3, 2, 3, 2, 2, 1, 4, 1, 2, 2, 4, 3, 1, 1, 4, 3, 2, 4, 2, 2, 3, 2, 4, 4, 3, 1, 3, 4, 2, 4, 4, 4, 3, 1, 3, 1, 2, 3, 4, 2, 1, 2, 3, 4, 3, 3, 1, 2, 2, 4, 3, 1, 4, 3, 3, 4, 3, 3, 3, 2, 1, 3, 2, 3, 2, 2, 1, 2, 4, 2, 4, 4, 4, 3, 2, 4, 2, 2, 3, 2, 2, 2, 4, 4, 4, 1, 2, 1, 4, 2, 3, 3, 3, 3, 1, 2, 3, 4, 3, 4, 3, 4, 3, 3, 3, 3, 3, 4, 3, 2, 3, 4, 4, 4, 4, 4, 3, 3, 4, 2, 3, 3, 4, 2, 1, 1, 2, 4, 3, 4, 4, 2, 1, 2, 2, 4, 4, 2, 2, 3, 1, 2, 2, 4, 2, 3, 2, 2, 3, 4, 2, 1, 2, 4, 3, 3, 1, 2, 2, 1, 2, 2, 2, 1, 1, 3, 1, 2, 4, 3, 4, 1, 1, 4, 4, 3, 4, 4, 2, 2, 1, 4, 2, 4, 1, 2, 1, 2, 1, 3, 3, 1, 2, 1, 4, 1, 4, 2, 1, 2, 1, 3, 3, 1, 2, 1, 2, 1, 3, 1, 1, 4, 3, 3, 2, 4, 4, 3, 3, 3, 1, 2, 1, 4, 3, 1, 4, 3, 1, 4, 3, 1, 1, 4, 4, 3, 3, 1, 3, 2, 2, 2, 2, 1, 2, 2, 3, 2, 3, 1, 2, 3, 2, 4, 3, 3, 4, 2, 2, 4, 2, 3, 2, 2, 2, 1, 1, 2, 4, 2, 1, 4, 3, 1, 3, 3, 1, 4, 2, 2, 2, 3, 3, 3, 2, 3, 2, 2, 1, 4, 3, 3, 4, 2, 3, 1, 2, 2, 4, 3, 2, 4, 2, 3, 2, 1, 3, 3, 2, 3, 1, 2, 2, 1, 4, 4, 3, 3, 3, 3, 2, 1, 4, 2, 2, 4, 3, 3, 4, 4, 3, 3, 2, 1, 4, 1, 3, 2, 3, 4, 1, 2, 4, 4, 2, 1, 3, 2, 1, 4, 3, 2, 1, 1, 3, 2, 2, 1, 1, 4, 4, 3, 3, 3, 2, 2, 1, 1, 3, 3, 4, 4, 1, 4, 2, 2, 4, 3, 3, 4, 2, 4, 4, 3, 4, 4, 4, 2, 4, 2, 4, 4, 4, 3, 2, 4, 3, 3, 1, 3, 4, 2, 3, 1, 2, 3, 2, 4, 3, 1, 3, 1, 2, 4, 3, 4, 1, 3, 4, 2, 4, 2, 2, 3, 3, 3, 3, 5, 4, 2, 1, 3, 3, 2, 3, 3, 3, 1, 2, 3, 3, 1, 2, 2, 1, 2, 3, 2, 3, 2, 3, 3, 2, 2, 4, 4, 3, 2, 2, 3, 2, 1, 1, 4, 1, 4, 4, 2, 1, 3, 2, 2, 1, 3, 3, 3, 4, 4, 2, 4, 2, 1, 3, 2, 1, 3, 1, 1, 2, 3, 4, 1, 2, 1, 3, 2, 4, 2, 1, 4, 4, 1, 1, 2, 1, 2, 2, 1, 1, 2, 3, 3, 3, 1, 3, 2, 4, 3, 3, 2, 1, 2, 1, 4, 2, 1, 1, 2, 1, 3, 1, 1, 2, 4, 3, 2, 2, 2, 4, 2, 1, 1, 4, 4, 3, 4, 1, 1, 2, 3, 1, 2, 1, 3, 2, 2, 4, 3, 1, 1, 2, 1, 2, 1, 3, 3, 1, 4, 4, 2, 4, 4, 1, 3, 2, 2, 1, 3, 4, 2, 4, 3, 4, 4, 4, 4, 1, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 4, 4, 2, 1, 1, 2, 1, 3, 4, 4, 2, 2, 3, 3, 3, 4, 3, 4, 2, 2, 1, 3, 1, 3, 1, 3, 3, 2, 4, 1, 3, 2, 2, 1, 3, 2, 4, 3, 2, 1, 3, 2, 1, 3, 2, 2, 4, 4, 3, 1, 2, 1, 3, 2, 2, 4, 2, 2, 1, 1, 2, 1, 3, 3, 3, 2, 4, 3, 3, 3, 3, 2, 2, 2, 4, 2, 4, 4, 3, 3, 3, 1, 2, 4, 4, 1, 4, 2, 1, 3, 3, 2, 4, 1, 1, 4, 3, 1, 3, 4, 2, 2, 3, 1, 2, 1, 2, 1, 1, 3, 1, 2, 2, 2, 4, 1, 2, 4, 3, 2, 4, 3, 3, 1, 1, 2, 4, 1, 2, 1, 2, 3, 2, 2, 4, 2, 3, 1, 2, 2, 2, 4, 3, 2, 4, 3, 3, 1, 2, 3, 3, 4, 3, 2, 2, 4, 3, 2, 1, 4, 2, 2, 1, 2, 2, 3, 4, 3, 4, 3, 4, 3, 1, 2, 2, 2, 2, 3, 4, 1, 4, 1, 2, 4, 3, 4, 4, 4, 4, 1, 2, 1, 2, 2, 2, 1, 3, 2, 2, 2, 4, 3, 4, 3, 3, 4, 1, 3, 4, 2, 3, 3, 3, 1, 2, 1, 1, 2, 2, 3, 1, 4, 2, 3, 2, 2, 3, 4, 3, 3, 3, 3, 2, 4, 2, 2, 3, 1, 2, 1, 4, 1, 2, 1, 2, 3, 4, 3, 3, 3, 3, 2, 3, 1, 3, 1, 1, 4, 3, 1, 3, 1, 2, 2, 3, 1, 4, 3, 4, 3, 4, 4, 4, 2, 4, 3, 2, 4, 3, 1, 1, 2, 1, 2, 1, 1, 2, 2, 1, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 3, 2, 4, 4, 3, 3, 4, 4, 2, 3, 3, 2, 4, 3, 2, 1, 2, 3, 4, 3, 3, 1, 4, 3, 1, 1, 2, 1, 3, 2, 1, 2, 1, 3, 3, 3, 4, 4, 4, 1, 2, 1, 1, 1, 1, 4, 2, 2, 4, 3, 2, 3, 3, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 4, 3, 2, 1, 3, 4, 3, 4, 3, 1, 2, 2, 1, 2, 1, 1, 1, 4, 1, 1, 2, 3, 3, 2, 1, 2, 2, 4, 3, 3, 3, 3, 1, 4, 3, 4, 2, 2, 1, 1, 2, 4, 3, 1, 2, 4, 3, 2, 4, 4, 2, 1, 3, 1, 1, 1, 1, 2, 1, 4, 2, 2, 1, 3, 1, 1, 3, 2, 2, 1, 2, 2, 4, 1, 2, 1, 2, 2, 1, 1, 3, 4, 3, 4, 3, 3, 2, 4, 2, 1, 3, 3, 3, 2, 2, 4, 4, 3, 3, 2, 4, 4, 1, 2, 3, 2, 2, 1, 1, 3, 1, 4, 3, 2, 2, 4, 1, 3, 4, 2, 3, 1, 2, 4, 1, 4, 2, 2, 4, 4, 1, 2, 2, 3, 3, 2, 4, 3, 4, 3, 3, 2, 1, 2, 4, 1, 4, 2, 2, 3, 4, 3, 2, 1, 2, 1, 3, 4, 2, 1, 1, 2, 4, 1, 4, 1, 2, 2, 1, 4, 2, 4, 4, 4, 1, 1, 3, 3, 4, 2, 1, 3, 3, 1, 4, 3, 4, 1, 4, 3, 4, 2, 3, 3, 2, 3, 3, 3, 1, 4, 1, 3, 1, 3, 2, 1, 4, 2, 3, 3, 2, 4, 3, 3, 1, 2, 3, 2, 3, 3, 2, 2, 4, 3, 2, 1, 1, 2, 4, 3, 3, 1, 2, 4, 1, 3, 1, 3, 3, 3, 3, 1, 3, 2, 2, 2, 4, 3, 2, 1, 1, 2, 4, 4, 3, 3, 1, 1, 2, 1, 2, 1, 3, 3, 3, 1, 2, 2, 3, 2, 1, 2, 3, 3, 1, 3, 2, 4, 1, 4, 2, 4, 2, 2, 2, 2, 4, 3, 2, 4, 2, 2, 4, 2, 4, 2, 4, 1, 2, 2, 1, 2, 1, 2, 1, 3, 4, 3, 3, 2, 1, 1, 3, 4, 4, 2, 4, 2, 2, 2, 4, 4, 3, 2, 4, 2, 4, 4, 4, 2, 1, 2, 2, 1, 2, 3, 2, 4, 3, 2, 2, 2, 3, 2, 2, 2, 4, 2, 4, 2, 3, 1, 2, 4, 3, 3, 2, 4, 4, 3, 1, 4, 4, 2, 1, 2, 2, 4, 2, 2, 1, 2, 2, 1, 3, 1, 1, 2, 1, 4, 2, 3, 4, 4, 3, 1, 4, 3, 4, 2, 2, 1, 1, 4, 1, 2, 2, 4, 2, 4, 1, 4, 3, 3, 4, 3, 4, 4, 3, 3, 2, 4, 2, 1, 3, 2, 4, 3, 4, 3, 3, 4, 1, 4, 2, 4, 4, 3, 3, 3, 2, 2, 2, 4, 4, 1, 1, 3, 4, 3, 3, 3, 1, 3, 4, 1, 2, 1, 4, 4, 3, 4, 3, 2, 4, 4, 3, 2, 2, 4, 4, 2, 2, 4, 3, 2, 4, 3, 4, 4, 3, 3, 2, 3, 3, 1, 2, 3, 2, 2, 1, 3, 3, 2, 4, 2, 4, 3, 4, 3, 2, 3, 4, 3, 2, 4, 4, 1, 4, 3, 3, 1, 4, 3, 1, 4, 3, 2, 4, 3]
		 * input:n	[4,4,3,3,2,1,2,4,3,2,4,2,4,2,3,4,3,2,2,4,3,1,2,4,3,4,4, 2,2,2,3, 2, 2, 4, 2, 3, 3, 2, 2, 4, 1, 2, 1, 1, 2, 4, 1, 4, 2, 3, 2, 1, 1, 2, 1, 3, 2, 4, 2, 3, 3, 3, 4, 3, 4, 2, 4, 1, 2, 2, 1, 4, 3, 4, 2, 1, 2, 2, 1, 1, 2, 3, 1, 2, 4, 3, 2, 2, 2, 3, 1, 1, 2, 4, 2, 3, 1, 3, 2, 1, 4, 1, 3, 4, 2, 4, 1, 4, 3, 1, 1, 3, 2, 2, 3, 1, 4, 2, 1, 2, 2, 1, 2, 1, 4, 2, 4, 4, 3, 2, 1, 2, 2, 4, 2, 2, 2, 3, 3, 3, 1, 4, 3, 2, 3, 4, 4, 2, 2, 4, 4, 3, 2, 3, 4, 3, 1, 3, 3, 3, 4, 2, 3, 3, 3, 1, 1, 4, 1, 2, 3, 4, 2, 4, 1, 2, 1, 4, 3, 2, 4, 3, 3, 3, 4, 3, 4, 2, 4, 2, 4, 2, 1, 2, 4, 2, 2, 2, 1, 2, 4, 3, 4, 3, 3, 2, 4, 3, 2, 3, 2, 2, 1, 4, 1, 2, 2, 4, 3, 1, 1, 4, 3, 2, 4, 2, 2, 3, 2, 4, 4, 3, 1, 3, 4, 2, 4, 4, 4, 3, 1, 3, 1, 2, 3, 4, 2, 1, 2, 3, 4, 3, 3, 1, 2, 2, 4, 3, 1, 4, 3, 3, 4, 3, 3, 3, 2, 1, 3, 2, 3, 2, 2, 1, 2, 4, 2, 4, 4, 4, 3, 2, 4, 2, 2, 3, 2, 2, 2, 4, 4, 4, 1, 2, 1, 4, 2, 3, 3, 3, 3, 1, 2, 3, 4, 3, 4, 3, 4, 3, 3, 3, 3, 3, 4, 3, 2, 3, 4, 4, 4, 4, 4, 3, 3, 4, 2, 3, 3, 4, 2, 1, 1, 2, 4, 3, 4, 4, 2, 1, 2, 2, 4, 4, 2, 2, 3, 1, 2, 2, 4, 2, 3, 2, 2, 3, 4, 2, 1, 2, 4, 3, 3, 1, 2, 2, 1, 2, 2, 2, 1, 1, 3, 1, 2, 4, 3, 4, 1, 1, 4, 4, 3, 4, 4, 2, 2, 1, 4, 2, 4, 1, 2, 1, 2, 1, 3, 3, 1, 2, 1, 4, 1, 4, 2, 1, 2, 1, 3, 3, 1, 2, 1, 2, 1, 3, 1, 1, 4, 3, 3, 2, 4, 4, 3, 3, 3, 1, 2, 1, 4, 3, 1, 4, 3, 1, 4, 3, 1, 1, 4, 4, 3, 3, 1, 3, 2, 2, 2, 2, 1, 2, 2, 3, 2, 3, 1, 2, 3, 2, 4, 3, 3, 4, 2, 2, 4, 2, 3, 2, 2, 2, 1, 1, 2, 4, 2, 1, 4, 3, 1, 3, 3, 1, 4, 2, 2, 2, 3, 3, 3, 2, 3, 2, 2, 1, 4, 3, 3, 4, 2, 3, 1, 2, 2, 4, 3, 2, 4, 2, 3, 2, 1, 3, 3, 2, 3, 3, 2, 2, 1, 4, 4, 3, 3, 3, 3, 2, 1, 4, 2, 2, 4, 3, 3, 4, 4, 3, 3, 2, 1, 4, 1, 3, 2, 3, 4, 1, 2, 4, 4, 2, 1, 3, 2, 1, 4, 3, 2, 1, 1, 3, 2, 2, 1, 1, 4, 4, 3, 3, 3, 2, 2, 1, 1, 3, 3, 4, 4, 1, 4, 2, 2, 4, 3, 3, 4, 2, 4, 4, 3, 4, 4, 4, 2, 4, 2, 4, 4, 4, 3, 2, 4, 3, 3, 1, 3, 4, 2, 3, 1, 2, 3, 2, 4, 3, 1, 3, 1, 2, 4, 3, 4, 1, 3, 4, 2, 4, 2, 2, 3, 3, 3, 3, 3, 4, 2, 1, 3, 3, 2, 3, 3, 5, 1, 2, 3, 3, 1, 2, 2, 1, 2, 3, 2, 3, 2, 3, 3, 2, 2, 4, 4, 3, 2, 2, 3, 2, 1, 1, 4, 3, 4, 4, 2, 1, 3, 2, 2, 1, 3, 3, 3, 4, 4, 2, 4, 2, 1, 3, 2, 1, 3, 1, 1, 2, 3, 4, 1, 2, 1, 3, 2, 4, 2, 1, 4, 4, 1, 1, 2, 1, 2, 2, 1, 1, 2, 3, 3, 3, 1, 3, 2, 4, 3, 3, 2, 1, 2, 1, 4, 2, 1, 1, 2, 1, 3, 1, 1, 2, 4, 3, 2, 2, 2, 4, 2, 1, 1, 4, 4, 3, 4, 1, 1, 2, 3, 1, 2, 1, 3, 2, 2, 4, 3, 1, 1, 2, 1, 2, 1, 3, 3, 1, 4, 4, 2, 4, 4, 1, 3, 2, 2, 1, 3, 4, 2, 4, 3, 4, 4, 4, 4, 1, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 4, 4, 2, 1, 1, 2, 1, 3, 4, 4, 2, 2, 3, 3, 3, 4, 3, 4, 2, 2, 1, 3, 1, 3, 1, 3, 3, 2, 4, 1, 3, 2, 2, 1, 3, 2, 4, 3, 2, 1, 3, 2, 1, 3, 2, 2, 4, 4, 3, 1, 2, 1, 3, 2, 2, 4, 2, 2, 1, 1, 2, 1, 3, 3, 3, 2, 4, 3, 3, 3, 3, 2, 2, 2, 4, 2, 4, 4, 3, 3, 3, 1, 2, 4, 4, 1, 4, 2, 1, 3, 3, 2, 4, 1, 1, 4, 3, 1, 3, 4, 2, 2, 3, 1, 2, 1, 2, 1, 1, 3, 1, 2, 2, 2, 4, 1, 2, 4, 3, 2, 4, 3, 3, 1, 1, 2, 4, 1, 2, 1, 2, 3, 2, 2, 4, 2, 3, 1, 2, 2, 2, 4, 3, 2, 4, 3, 3, 1, 2, 3, 3, 4, 3, 2, 2, 4, 3, 2, 1, 4, 2, 2, 1, 2, 2, 3, 4, 3, 4, 3, 4, 3, 3, 2, 2, 2, 2, 3, 4, 1, 4, 1, 2, 4, 3, 4, 4, 4, 4, 1, 2, 1, 2, 2, 2, 1, 3, 2, 2, 2, 4, 3, 4, 3, 3, 4, 1, 3, 4, 2, 3, 3, 3, 1, 2, 1, 1, 2, 2, 3, 1, 4, 2, 3, 2, 2, 3, 4, 3, 3, 3, 3, 2, 4, 2, 2, 3, 1, 2, 1, 4, 1, 2, 1, 2, 3, 4, 3, 3, 3, 3, 2, 3, 1, 3, 1, 1, 4, 3, 1, 3, 1, 2, 2, 3, 1, 4, 3, 4, 3, 4, 4, 4, 2, 4, 3, 2, 4, 3, 1, 1, 2, 1, 2, 1, 1, 2, 2, 1, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 3, 2, 4, 4, 3, 3, 4, 4, 2, 3, 3, 2, 4, 3, 2, 1, 2, 3, 4, 3, 3, 1, 4, 3, 1, 1, 2, 1, 3, 2, 1, 2, 1, 3, 3, 3, 4, 4, 4, 1, 2, 1, 1, 1, 1, 4, 2, 2, 4, 3, 2, 3, 3, 2, 3, 3, 2, 2, 2, 2, 2, 2, 5, 4, 3, 2, 1, 3, 4, 3, 4, 3, 1, 2, 2, 1, 2, 1, 1, 1, 4, 1, 1, 2, 3, 3, 2, 1, 2, 2, 4, 3, 3, 3, 3, 1, 4, 3, 4, 2, 2, 1, 1, 2, 4, 3, 1, 2, 4, 3, 2, 4, 4, 2, 1, 3, 1, 1, 1, 1, 2, 1, 4, 2, 2, 1, 3, 1, 1, 3, 2, 2, 1, 2, 2, 4, 1, 2, 1, 2, 2, 1, 1, 3, 4, 3, 4, 3, 3, 2, 4, 2, 1, 3, 3, 3, 2, 2, 4, 4, 3, 3, 2, 4, 4, 1, 2, 3, 2, 2, 1, 1, 3, 1, 4, 3, 4, 2, 4, 1, 3, 4, 2, 3, 1, 2, 4, 1, 4, 2, 2, 4, 4, 1, 2, 2, 3, 3, 2, 4, 3, 4, 3, 3, 2, 1, 2, 4, 1, 4, 2, 2, 3, 4, 3, 2, 1, 2, 1, 3, 4, 2, 1, 1, 2, 4, 1, 4, 1, 2, 2, 1, 4, 2, 4, 4, 4, 1, 1, 3, 3, 4, 2, 1, 3, 3, 1, 4, 3, 4, 1, 4, 3, 4, 2, 3, 3, 2, 3, 3, 3, 1, 4, 1, 3, 1, 3, 2, 1, 4, 2, 3, 3, 2, 4, 3, 3, 1, 2, 3, 2, 3, 3, 2, 2, 4, 3, 2, 1, 1, 2, 4, 3, 3, 1, 2, 4, 1, 3, 1, 3, 3, 3, 3, 1, 3, 2, 2, 2, 4, 3, 2, 1, 1, 2, 4, 4, 3, 3, 1, 1, 2, 1, 2, 1, 3, 3, 3, 1, 2, 2, 3, 2, 1, 2, 3, 3, 1, 3, 2, 4, 1, 4, 2, 4, 2, 2, 2, 2, 4, 3, 2, 4, 2, 2, 4, 2, 4, 2, 4, 1, 2, 2, 1, 2, 1, 2, 1, 3, 4, 3, 3, 2, 1, 1, 3, 4, 4, 2, 4, 2, 2, 2, 4, 4, 3, 2, 4, 2, 4, 4, 4, 2, 1, 2, 2, 1, 2, 3, 2, 4, 3, 2, 2, 2, 3, 2, 2, 2, 4, 2, 4, 2, 3, 1, 2, 4, 3, 3, 2, 4, 4, 3, 1, 4, 4, 2, 1, 2, 2, 4, 2, 2, 1, 2, 2, 1, 3, 1, 1, 2, 1, 4, 2, 3, 4, 4, 3, 1, 4, 3, 4, 2, 2, 1, 1, 4, 1, 2, 2, 4, 2, 4, 1, 4, 3, 3, 4, 3, 4, 4, 3, 3, 2, 4, 2, 1, 3, 2, 4, 3, 4, 3, 3, 4, 1, 4, 2, 4, 4, 3, 3, 3, 2, 2, 2, 4, 4, 1, 1, 3, 4, 3, 3, 3, 1, 3, 4, 1, 2, 1, 4, 4, 3, 4, 3, 2, 4, 4, 3, 2, 4, 4, 4, 2, 2, 4, 3, 2, 4, 3, 4, 4, 3, 3, 2, 3, 3, 1, 2, 3, 2, 2, 1, 3, 3, 2, 4, 2, 4, 3, 4, 3, 2, 3, 4, 3, 2, 4, 4, 1, 4, 3, 3, 1, 4, 3, 1, 4, 3, 2, 4, 3]
		 * mask_mid		[2,0,0,0,0,1,0,0,0,0,0,0,0,2,2,0,0,0,1,0,0,1,1,1,0,0,0
		 * output:1	[4,4,1,2,2,3,3,2,1,2,4]
		 */
		int seed = 1;
		int replicates = 3;
		int[] testResult = {4,4,1,2,2,3,3,2,1,2,4};
		boolean [] maskPositions = {
				true,false,false,true,true,true,false,false,false,
				true,false,false,true,true,true,false,false,false,
				true,false,false,true,true,true,false,false,false};
		
		TeaspoonMask mask = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,maskPositions);
		TeaspoonBootstrap[] bootstraps = TeaspoonBootstrapFactory.generate(mask, mainPartition, replicates, seed);
		bootstraps[0].getBootstraps();
		for(int i=0;i<replicates;i++){
			int[][] subsampledPartition = mainPartition.subsampleByBootstrap(bootstraps[0]).getSiteMatrix();
			if(!TeaspoonBootstrapTest.compareIntArraysByValue(testResult, subsampledPartition[0])){
				fail("sampled arrays don't match");
			}
			
		}
	}

	/**
	 * Test method for {@link teaspoon.app.utils.BhattAdaptationFullSiteMatrix#appendTaxa(teaspoon.app.utils.BhattAdaptationFullSiteMatrix)}.
	 */
	@Test
	public final void testAppendTaxa() {
		int oldTaxaCount = mainPartition.numberOfSequences();
		int oldLength = mainPartition.alignmentLength();
		BhattAdaptationFullSiteMatrix newDoubledAlignment = mainPartition.appendTaxa(mainPartition);
		int newTaxaCount = newDoubledAlignment.numberOfSequences();
		int newLength = newDoubledAlignment.alignmentLength();
		// test alignment dimensions
		if(oldTaxaCount != mainPartition.numberOfSequences()){
			fail("old alignment taxa count has been modified");
		}
		if(oldLength != mainPartition.alignmentLength()){
			fail("old alignment length has been modified");
		}
		if(newTaxaCount != (oldTaxaCount*2)){
			fail("new alignment not twice as many seqs as old one");
		}
		if(newLength != oldLength){
			fail("new/old alignment lengths don't match");
		}
	}

	/**
	 * Test method for {@link teaspoon.app.utils.BhattAdaptationFullSiteMatrix#alignmentLength()}.
	 */
	@Test
	public final void testAlignmentLength() {
		int alignmentLength = mainPartition.alignmentLength();
		if(alignmentLength != 1680){
			fail("incorrect sequence length");
		}
	}

	/**
	 * Test method for {@link teaspoon.app.utils.BhattAdaptationFullSiteMatrix#numberOfSequences()}.
	 */
	@Test
	public final void testNumberOfSequences() {
		if(mainPartition.numberOfSequences() != 694){
			fail("incorrect number of taxa");
		}
	}

}
