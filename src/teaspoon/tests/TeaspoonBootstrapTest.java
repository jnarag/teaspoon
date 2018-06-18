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
import teaspoon.app.standalone.TeaspoonMaskFactory;
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
 * @since 15 Jun 2018
 * @version 0.1
 */
public class TeaspoonBootstrapTest {

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
	 * Test method for {@link teaspoon.app.TeaspoonBootstrap#TeaspoonBootstrap(int)}.
	 */
	@Test
	public final void testTeaspoonBootstrapInt() {
		TeaspoonBootstrap bootstrap = new TeaspoonBootstrap(30);
		bootstrap.getBootstraps();
	}

	/**
	 * Test method for {@link teaspoon.app.TeaspoonBootstrap#TeaspoonBootstrap(int, int)}.
	 */
	@Test
	public final void testTeaspoonBootstrapIntInt() {
		TeaspoonBootstrap bootstrap = new TeaspoonBootstrap(20,10);
		bootstrap.getBootstraps();
		TeaspoonBootstrap bootstrap2 = new TeaspoonBootstrap(20,10);
		bootstrap2.getBootstraps();
		TeaspoonBootstrap bootstrap3 = new TeaspoonBootstrap(20,1);
		bootstrap3.getBootstraps();
		if(!compareIntArraysByValue(bootstrap.getBootstraps(),bootstrap2.getBootstraps())){
			fail();
		}
		if(compareIntArraysByValue(bootstrap.getBootstraps(),bootstrap3.getBootstraps())){
			fail();
		}
		if(compareIntArraysByValue(bootstrap2.getBootstraps(),bootstrap3.getBootstraps())){
			fail();
		}
	}

	/**
	 * Test method for {@link teaspoon.app.TeaspoonBootstrap#TeaspoonBootstrap(int, teaspoon.app.TeaspoonMask)}.
	 */
	@Test
	public final void testTeaspoonBootstrapIntTeaspoonMask() {
		boolean [] maskPositions = {
				true,false,false,true,true,true,false,false,false,
				true,false,false,true,true,true,false,false,false,
				true,false,false,true,true,true,false,false,false};
		TeaspoonMask mask = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,maskPositions);
		TeaspoonBootstrap bootstrap = new TeaspoonBootstrap(27,mask);
		bootstrap.getBootstraps();
		TeaspoonBootstrap bootstrap2 = new TeaspoonBootstrap(27,mask);
		bootstrap2.getBootstraps();
		if(compareIntArraysByValue(bootstrap.getBootstraps(),bootstrap2.getBootstraps())){
			fail();
		}
		// test more heavily
		for(int testReplicate=0;testReplicate<100;testReplicate++){
			if(!compareIntArrayAndMask(new TeaspoonBootstrap(27,mask).getBootstraps(), maskPositions)){
				fail("test faliure: replicate "+testReplicate);
			}
		}
	}

	/**
	 * Test method for {@link teaspoon.app.TeaspoonBootstrap#TeaspoonBootstrap(int, teaspoon.app.TeaspoonMask, int)}.
	 */
	@Test
	public final void testTeaspoonBootstrapIntTeaspoonMaskInt() {
		boolean [] maskPositions = {true,false,false,true,true,true,false,false,false};
		boolean [] maskOtherPositions = {false,false,false,true,false,true,false,false,true};
		TeaspoonMask mask = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,maskPositions);
		TeaspoonMask otherMask = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,maskOtherPositions);
		int seed1 = 1;
		int seed2 = (int)System.currentTimeMillis()-1;
		TeaspoonBootstrap bootstrapOriginal = new TeaspoonBootstrap(9,mask,seed1);
		TeaspoonBootstrap bootstrapMatching = new TeaspoonBootstrap(9,mask,seed1);
		TeaspoonBootstrap bootstrapDifferentSeed = new TeaspoonBootstrap(9,mask,seed2);
		TeaspoonBootstrap bootstrapDifferentMask = new TeaspoonBootstrap(9,otherMask,seed1);
		if(!compareIntArraysByValue(bootstrapOriginal.getBootstraps(),bootstrapMatching.getBootstraps())){
			fail();
		}
		if(compareIntArraysByValue(bootstrapOriginal.getBootstraps(),bootstrapDifferentSeed.getBootstraps())){
			fail();
		}
		if(compareIntArraysByValue(bootstrapOriginal.getBootstraps(),bootstrapDifferentMask.getBootstraps())){
			fail();
		}
		if(compareIntArraysByValue(bootstrapMatching.getBootstraps(),bootstrapDifferentSeed.getBootstraps())){
			fail();
		}
		// test more heavily, mask directly tested
		for(int testReplicate=0;testReplicate<100;testReplicate++){
			if(!compareIntArrayAndMask(new TeaspoonBootstrap(9,mask).getBootstraps(), maskPositions)){
				fail("test faliure: replicate "+testReplicate);
			}
		}
		// test mask more heavily
		for(int testReplicate=0;testReplicate<100;testReplicate++){
			if(!compareIntArrayAndMask(new TeaspoonBootstrap(9,otherMask).getBootstraps(), maskOtherPositions)){
				fail("test faliure: replicate "+testReplicate);
			}
		}
		// test mask more heavily
		for(int testReplicate=0;testReplicate<100;testReplicate++){
			if(compareIntArrayAndMask(new TeaspoonBootstrap(9,mask).getBootstraps(), maskOtherPositions)){
				fail("test faliure: replicate "+testReplicate);
			}
		}
	}

	/**
	 * Test method for {@link teaspoon.app.TeaspoonBootstrap#getBootstraps()}.
	 */
	@Test
	public final void testGetBootstraps() {
		boolean [] maskPositions = {true,false,false,true,true,true,false,false,false};
		TeaspoonMask mask = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,maskPositions);
		int seed1 = 1;
		TeaspoonBootstrap bootstrapOriginal = new TeaspoonBootstrap(9,mask,seed1);
		int[] bootstraps = bootstrapOriginal.getBootstraps();
		int totalSampled = 0;
		for(int i=0;i<bootstraps.length;i++){totalSampled += bootstraps[i];}
		if(totalSampled != bootstraps.length){fail();}
	}

	/**
	 * Test method for {@link teaspoon.app.TeaspoonBootstrap#getSeed()}.
	 */
	@Test
	public final void testGetSeed() {
		boolean [] maskPositions = {true,false,false,true,true,true,false,false,false};
		TeaspoonMask mask = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,maskPositions);
		int seed1 = 1;
		TeaspoonBootstrap bootstrapOriginal = new TeaspoonBootstrap(9,mask,seed1);
		if(bootstrapOriginal.getSeed() != seed1){fail();}
	}
	
	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonBootstrapFactory#generate()}.
	 */
	@Test
	public final void testBootstrapFactory(){
		int seed = 1;
		int replicates = 3;
		boolean [] maskPositions = {
				true,false,false,true,true,true,false,false,false,
				true,false,false,true,true,true,false,false,false,
				true,false,false,true,true,true,false,false,false};
		File input = new File("./HCV_data/sub_053/FP7_05305_1.3699.fasta");
		TeaspoonMask mask = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,maskPositions);
		BhattAdaptationFullSiteMatrix mainPartition = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(input).readFASTA());
		TeaspoonBootstrap[] bootstraps = TeaspoonBootstrapFactory.generate(mask, mainPartition, replicates, seed);
		bootstraps[0].getBootstraps();
	}
	
	/**
	 * Utility method comparing integers in two arrays by value not hashcode
	 * @param array1
	 * @param array2
	 * @return
	 */
	public static boolean compareIntArraysByValue(int[] array1, int[] array2){
		boolean valid = true;
		if(array1.length != array2.length){
			valid = false;
		}else{
			for(int i=0;i<array1.length;i++){
				valid = valid && (array1[i] == array2[i]);
			}
		}
		return valid;
	}

	/**
	 * Utility method comparing a bootstrap against mask positions to check that
	 * only (mask==true) sites are selected.
	 * <p>NB not entirely guaranteed to be safe since it drops to false if sampled
	 * sites (bootstrap[i]>0) are in the mask's boolean array. Since even valid 
	 * sites might not get sampled, especially for big alignments, this test might 
	 * miss some. 
	 * @param array1
	 * @param array2
	 * @return
	 */
	public static boolean compareIntArrayAndMask(int[] array1, boolean[] array2){
		boolean valid = true;
		if(array1.length != array2.length){
			valid = false;
		}else{
			for(int i=0;i<array1.length;i++){
				if(array1[i]>0 && !array2[i]){
					valid = false;
				}
			}
		}
		return valid;
	}
}
