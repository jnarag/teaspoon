/**
 * 
 */
package teaspoon.tests;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import teaspoon.app.TeaspoonMask;
import teaspoon.app.standalone.TeaspoonMaskFactory;
import teaspoon.app.utils.RateEstimationBehaviour;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 14 Jun 2018
 * @version 0.1
 */
public class TeaspoonMaskFactoryTest extends TeaspoonMaskFactory {

	boolean [] mask = {false,false,false,true,true,true,false,false,false};
	boolean [] mask_end = {false,false,false,false,false,false,true,true,true};
	boolean [] mask_start = {true,true,true,false,false,false,false,false,false};

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
	 * Test method for {@link teaspoon.app.standalone.TeaspoonMaskFactory#parseFile(java.io.File)}.
	 */
	@Test
	public final void testParseFile() {
		TeaspoonMask[] masks = null;
		try {
			masks = TeaspoonMaskFactory.parseFile(new File("./HCV_data/mask"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		boolean validParse = true;
		if(masks[0].estimationBehaviour != RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED){
			validParse = false;
		}
		for(int position=0;position<mask.length;position++){
			if(mask[position] != masks[0].getPositions()[position]){
				validParse = false;
			}
		}
		if(!validParse){
			fail("Mask incorrect"); 
		}
	}

	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonMaskFactory#parseFile(java.io.File)}.
	 */
	@Test
	public final void testManyFilesParse(){
		// test a wider range of parse operations, against the maskfiles produced in other tests
		testWriteMaskFile();
		testWriteMaskFileWithFixedRatio();
		testParseFile(new File("test-mask-averaged"),new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AVERAGED,mask));
		testParseFile(new File("test-mask-fix-0.1"),new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,mask_start,0.1));
		testParseFile(new File("test-mask-fix-0.2"),new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,mask_start,0.2));
		testParseFile(new File("test-mask-fix-0.3"),new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,mask_start,0.3));
	}
	
	public final void testParseFile(File file, TeaspoonMask someMask){
		TeaspoonMask[] masks = null;
		try {
			masks = TeaspoonMaskFactory.parseFile(file);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		boolean validParse = true;
		if(masks[0].estimationBehaviour != someMask.estimationBehaviour){
			validParse = false;
		}
		for(int position=0;position<someMask.getLength();position++){
			if(someMask.getPositions()[position] != masks[0].getPositions()[position]){
				validParse = false;
			}
		}
		if(!validParse){
			fail("Mask incorrect"); 
		}
		
	}

	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonMaskFactory#writeMaskFileWithFixedRatio(java.io.File, java.util.HashMap, int, double)}.
	 */
	@Test
	public final void testWriteMaskFileWithFixedRatio() {
		ArrayList<int[]> positions = new ArrayList<int[]>();
		int[] region = {0,2};
		positions.add(region);
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map1 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map2 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map3 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		map1.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
		map2.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
		map3.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
		try {
			TeaspoonMaskFactory.writeMaskFileWithFixedRatio(new File("test-mask-fix-0.1"), map1, 9,0.1);
			TeaspoonMaskFactory.writeMaskFileWithFixedRatio(new File("test-mask-fix-0.2"), map2, 9,0.2);
			TeaspoonMaskFactory.writeMaskFileWithFixedRatio(new File("test-mask-fix-0.3"), map3, 9,0.3);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			fail("Not yet implemented"); // TODO
		}
	}

	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonMaskFactory#appendToMaskFile(java.io.File, java.util.HashMap, int)}.
	 */
	@Test
	public final void testAppendToMaskFile() {
		testWriteMaskFile();
		ArrayList<int[]> positions = new ArrayList<int[]>();
		int[] region = {6,8};
		positions.add(region);
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map1 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map2 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map3 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		map1.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED);
		map2.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_AVERAGED);
		map3.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
		try {
			TeaspoonMaskFactory.appendToMaskFile(new File("test-mask-aggregated"), map1, 9);
			TeaspoonMaskFactory.appendToMaskFile(new File("test-mask-averaged"), map2, 9);
			TeaspoonMaskFactory.appendToMaskFile(new File("test-mask-fixed"), map3, 9);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			fail("Not yet implemented"); // TODO
		}
	}

	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonMaskFactory#appendToMaskFileWithFixedRatio(java.io.File, java.util.HashMap, int, double)}.
	 */
	@Test
	public final void testAppendToMaskFileWithFixedRatio() {
		testWriteMaskFileWithFixedRatio();
		ArrayList<int[]> positions = new ArrayList<int[]>();
		int[] region = {6,8};
		positions.add(region);
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map1 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map2 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map3 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		map1.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
		map2.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
		map3.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
		try {
			TeaspoonMaskFactory.appendToMaskFileWithFixedRatio(new File("test-mask-fix-0.1"), map1, 9,0.1);
			TeaspoonMaskFactory.appendToMaskFileWithFixedRatio(new File("test-mask-fix-0.2"), map2, 9,0.2);
			TeaspoonMaskFactory.appendToMaskFileWithFixedRatio(new File("test-mask-fix-0.3"), map3, 9,0.3);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			fail("Not yet implemented"); // TODO
		}
	}

	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonMaskFactory#writeMaskFile(java.io.File, java.util.HashMap, int)}.
	 */
	@Test
	public final void testWriteMaskFile() {
		ArrayList<int[]> positions = new ArrayList<int[]>();
		int[] region = {3,5};
		positions.add(region);
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map1 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map2 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		HashMap<ArrayList<int[]>,RateEstimationBehaviour> map3 = new HashMap<ArrayList<int[]>,RateEstimationBehaviour>();
		map1.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED);
		map2.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_AVERAGED);
		map3.put(positions,RateEstimationBehaviour.NEUTRAL_RATE_FIXED);
		try {
			TeaspoonMaskFactory.writeMaskFile(new File("test-mask-aggregated"), map1, 9);
			TeaspoonMaskFactory.writeMaskFile(new File("test-mask-averaged"), map2, 9);
			TeaspoonMaskFactory.writeMaskFile(new File("test-mask-fixed"), map3, 9);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			fail("Not yet implemented"); // TODO
		}
	}

}
