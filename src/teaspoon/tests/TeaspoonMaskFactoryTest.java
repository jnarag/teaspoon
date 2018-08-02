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

	boolean [] mask_mid  = {false,false,false,true,true,true,false,false,false};
	boolean [] mask_all  = {true,true,true,true,true,true,true,true,true};
	boolean [] mask_none = {false,false,false,false,false,false,false,false,false};
	boolean [] mask_discontiguous = {true,true,true,false,false,false,true,true,true};
	boolean [] mask_at_end = {false,false,false,false,false,false,true,true,true};
	boolean [] mask_at_start = {true,true,true,false,false,false,false,false,false};

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
		for(int position=0;position<mask_mid.length;position++){
			if(mask_mid[position] != masks[0].getPositions()[position]){
				validParse = false;
			}
		}
		if(!validParse){
			fail("Mask incorrect"); 
		}
	}

	
	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonMaskFactory#combineMasksUnion(teaspoon.app.TeaspoonMask, teaspoon.app.TeaspoonMask)}.
	 */
	@Test
	public final void testCombineMasks(){
		TeaspoonMask mask_one;	// first
		TeaspoonMask mask_two;	// second
		TeaspoonMask combine;	// should be union
		
		mask_one = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,mask_at_start);
		mask_two = new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,mask_at_end);
		
		combine = TeaspoonMaskFactory.combineMasksUnion(mask_one,mask_two);
		for(int i=0;i<combine.getLength();i++){
			assertTrue(combine.getPositions()[i] == mask_discontiguous[i]);
		}
		
	}

	/**
	 * Test method for {@link teaspoon.app.standalone.TeaspoonMaskFactory#parseFile(java.io.File)}.
	 */
	@Test
	public final void testManyFilesParse(){
		// test a wider range of parse operations, against the maskfiles produced in other tests
		testWriteMaskFile();
		testParseFile(new File("test-mask-averaged"),new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AVERAGED,mask_mid));
		testParseFile(new File("test-mask-fix-0.1"),new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,mask_at_start,0.1));
		testParseFile(new File("test-mask-fix-0.2"),new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,mask_at_start,0.2));
		testParseFile(new File("test-mask-fix-0.3"),new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,mask_at_start,0.3));
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
	 * Test method for {@link teaspoon.app.standalone.TeaspoonMaskFactory#appendToMaskFile(java.io.File, java.util.HashMap, int)}.
	 */
	@Test
	public final void testAppendToMaskFile() {
		testWriteMaskFile();
		int[] region = {6,8};
		boolean[] positions = new boolean[10];
		for(int sequenceIndex=region[0];sequenceIndex<=region[1];sequenceIndex++){
			// flip positions in this range to true
			positions[sequenceIndex] = true;
		}
		ArrayList<TeaspoonMask> masksList = new ArrayList<TeaspoonMask>();
		masksList.add(new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,positions)) ;
		masksList.add(new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AVERAGED,positions))  ;
		masksList.add(new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,positions))  ;
		try {
			TeaspoonMaskFactory.writeMaskFile(new File("test-mask_mid-aggregated"), masksList);
			
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
		boolean[] positions = new boolean[10];
		int[] region = {3,5};
		for(int sequenceIndex=region[0];sequenceIndex<=region[1];sequenceIndex++){
			// flip positions in this range to true
			positions[sequenceIndex] = true;
		}
		ArrayList<TeaspoonMask> masksList = new ArrayList<TeaspoonMask>();
		masksList.add(new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,positions)) ;
		masksList.add(new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_AVERAGED,positions))  ;
		masksList.add(new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,positions))  ;
		try {
			TeaspoonMaskFactory.writeMaskFile(new File("test-mask_mid-aggregated"), masksList);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			fail("Not yet implemented"); // TODO
		}
	}

}
