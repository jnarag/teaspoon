/**
 * 
 */
package teaspoon.app.standalone;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.lang3.BooleanUtils;

import teaspoon.app.TeaspoonMask;
import teaspoon.app.utils.MainAlignmentParser;
import teaspoon.app.utils.RateEstimationBehaviour;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * TeaspoonMaskFactory provides simple methods to write sequence alignment
 * masks to and from files for adaptation analyses.
 * 
 * <p>Masks are comprised of a single ratio estimation specification and a 
 * list of booleans <i>exactly</i> as long as the sequence alignment.
 * 
 * <p><b>Note</b> that this is intended as a simple implementation, so very 
 * little safety checking is applied. In particular, there are no guarantees
 * that masks are as long as the sequence alignemt they refer to, or to each
 * other.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 10 May 2018
 * @version 0.1
 * @see TeaspoonMask
 * @see BhattAdaptationFullSitesMatrix
 * @see BhattAdaptationAnalysis
 * @see TeaspoonCommandLineApp
 * @see RateEstimationBehaviour
 */
public class TeaspoonMaskFactory {

	/**
	 * No-arg constructor deprecated at the moment, 
	 * methods are statics.
	 */
	@Deprecated
	public TeaspoonMaskFactory() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		/* 
		 * At the moment not likely to use this method
		 * but leaving open the possiblilty for user-generated
		 * maskfiles, e.g. 
		 * `java -jar MaskFactory.jar <alignment file> <mask_mid specification file>` etc
		 * 
		 */
		File ancestralAlignment, maskSpecification, maskOutput;
		try{
			ancestralAlignment = new File(args[0]);
			maskSpecification = new File(args[1]);
			maskOutput = new File(args[2]);
			/* sort input out*/
			int alignmentLength = (new MainAlignmentParser(ancestralAlignment).readFASTA())[0].length;
			System.out.println("found length "+alignmentLength+" positions");
			
			/* parse mask_mid spec */
			BufferedReader reader = new BufferedReader(new FileReader(maskSpecification));
			ArrayList<TeaspoonMask> masks = new ArrayList<TeaspoonMask>();
			double ratio = Double.NaN;

			while(reader.ready()){
				String[] tokens = reader.readLine().split("\\,");
				RateEstimationBehaviour behaviour;
				
				// check for the special characters
				int lower;
				int upper;
				if(tokens[1].equals("*")){
					lower = 0;
				}else{
					lower = Integer.parseInt(tokens[1]);
				}
				if(tokens[2].equals("*")){
					upper = alignmentLength-1;
				}else{
					upper = Integer.parseInt(tokens[2]);
				}

				// parse behaviour string
				// switch the first token to determine mask_mid behaviour...
				switch(tokens[0]){
				case("aggregated"):{
					behaviour = RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED;
					masks.add(initialiseMask(behaviour,lower,upper,alignmentLength));
					break;
					}
				case("averaged"):{
					behaviour = RateEstimationBehaviour.NEUTRAL_RATE_AVERAGED;
					masks.add(initialiseMask(behaviour,lower,upper,alignmentLength));
					break;
					}
				default:{
					// ... but check for a numeric to convert to static ratio
					try{
						if(Double.parseDouble(tokens[0])>=0){
							ratio = Double.parseDouble(tokens[0]);
							behaviour = RateEstimationBehaviour.NEUTRAL_RATE_FIXED;
							masks.add(initialiseMask(ratio,lower,upper,alignmentLength));
						}
					}catch (NumberFormatException ex){
						// we don't need to print the stack trace necessarily
						// just make sure sensible things happen
						// ex.printStackTrace();
						behaviour = RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED;
						masks.add(initialiseMask(behaviour,lower,upper,alignmentLength));
					}
					break;
					}
				}				
			}
			reader.close();
			
			/* write maskfile */
			TeaspoonMaskFactory.writeMaskFile(maskOutput, masks);
		}catch (Exception ex){
			System.err.println("Cannot parse arguments. Exiting.");
			ex.printStackTrace();
		}
	}

	/**
	 * Parses an alignment mask_mid file, returns a series of masks.
	 * @param maskFile
	 * @return
	 * @throws IOException 
	 */
	public static TeaspoonMask[] parseFile(File maskFile) throws IOException {
		BufferedReader maskBuffer = new BufferedReader(new FileReader(maskFile));
		ArrayList<TeaspoonMask> masks = new ArrayList<TeaspoonMask>();
		try{
			String line = null;
			while((line = maskBuffer.readLine()) != null){
				// parse each line in the maskfile
				RateEstimationBehaviour behaviour;
				double ratio = Double.NaN;
				boolean[] mask;
				String[] tokens = line.split(",");
				// switch the first token to determine mask_mid behaviour...
				switch(tokens[0]){
				case("aggregated"):{
					behaviour = RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED;
					break;
					}
				case("averaged"):{
					behaviour = RateEstimationBehaviour.NEUTRAL_RATE_AVERAGED;
					break;
					}
				default:{
					behaviour = RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED;
					break;
					}
				}
				// ... but check for a numeric to convert to static ratio
				try{
					if(Double.parseDouble(tokens[0])>=0){
						ratio = Double.parseDouble(tokens[0]);
						behaviour = RateEstimationBehaviour.NEUTRAL_RATE_FIXED;
					}
				}catch (NumberFormatException ex){
					// we don't need to print the stack trace necessarily
					// just make sure sensible things happen
					// ex.printStackTrace();
				}
				
				// Now parse the second token to get mask_mid positions themselves
				mask = new boolean[tokens[1].toCharArray().length];
				int maskIndex = 0;
				for(char c:tokens[1].toCharArray()){
					mask[maskIndex] = BooleanUtils.toBoolean(Integer.parseInt(c+""));
					maskIndex++;
				}
				
				// Should now have a mask_mid and at least some of the info we need, put into arraylist
				// First check ratio setting behaviour
				if(behaviour == RateEstimationBehaviour.NEUTRAL_RATE_FIXED && ratio >= 0){
					// we want a fixed-ratio mask_mid, use 3-arg constructor
					masks.add(new TeaspoonMask(behaviour,mask,ratio));
				}else{
					// we want an estimated ratio
					masks.add(new TeaspoonMask(behaviour,mask));
				}
			}
		}catch(Exception ex){
			ex.printStackTrace();
		}finally{
			maskBuffer.close();
		}
		
		// return the mask_mid list
		return masks.toArray(new TeaspoonMask[0]);
	}

	/**
	 * Writes a list of RateEstimationBehaviours and mask_mid specifications to a named file.
	 * Coordinates are indexed to 0 not 1, e.g. first codon reads {0,1,2} not {1,2,3}
	 * @param maskFile - File to write to.
	 * @param masks - hash of rate estimation behaviours (enum) and int[2] start, end positions for masks (one or more pairs per element)
	 * @param alignmentLength - maximum alignment length to pad to.
	 * @return boolean if write operation is relatively successful.
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static boolean appendToMaskFile(File maskFile,ArrayList<TeaspoonMask> masks) throws FileNotFoundException, IOException{
		// check the existing file does in fact exist
		if(!maskFile.canRead()){
			throw new FileNotFoundException("mask_mid file "+maskFile+" cannot be read");
		}
		
		// initialise buffer
		StringBuffer outputBuffer = new StringBuffer();
		
		// walk through masks
		Iterator<TeaspoonMask> maskItr = masks.iterator();
		while(maskItr.hasNext()){
			TeaspoonMask mask = maskItr.next();
			RateEstimationBehaviour behaviour = mask.estimationBehaviour;
			String behaviourString = "aggregated";
			switch(behaviour){
			case NEUTRAL_RATE_AGGREGATED:{
				behaviourString = "aggregated";
				break;
			}
			case NEUTRAL_RATE_AVERAGED:{
				behaviourString = "averaged";
				break;
			}
			case NEUTRAL_RATE_FIXED:{
				behaviourString = "0.0";
				if(mask.getNeutralRatio()>=0.0){
					behaviourString = ""+mask.getNeutralRatio();
				}
				break;
			}
			default:{
				behaviourString = "aggregated";
				break;
				
			}
			}
			outputBuffer.append(behaviourString+","+maskToString(mask.getPositions())+"\n");
		}

		// write buffer to file
		BufferedWriter writer = new BufferedWriter(new FileWriter(maskFile,true));
		writer.write(outputBuffer.toString());
		writer.close();
		return true;
	}

	/**
	 * Writes a list of RateEstimationBehaviours and mask_mid specifications to a named file.
	 * Coordinates are indexed to 0 not 1, e.g. first codon reads {0,1,2} not {1,2,3}
	 * @param maskFile - File to write to.
	 * @param masks - hash of rate estimation behaviours (enum) and int[2] start, end positions for masks (one or more pairs per element)
	 * @param alignmentLength - maximum alignment length to pad to.
	 * @return boolean if write operation is relatively successful.
	 * @throws IOException 
	 */
	public static boolean writeMaskFile(File maskFile,ArrayList<TeaspoonMask> masks) throws IOException{
		// initialise buffer
		StringBuffer outputBuffer = new StringBuffer();
		
		// walk through masks
		Iterator<TeaspoonMask> maskItr = masks.iterator();
		while(maskItr.hasNext()){
			TeaspoonMask mask = maskItr.next();
			RateEstimationBehaviour behaviour = mask.estimationBehaviour;
			String behaviourString = "aggregated";
			switch(behaviour){
			case NEUTRAL_RATE_AGGREGATED:{
				behaviourString = "aggregated";
				break;
			}
			case NEUTRAL_RATE_AVERAGED:{
				behaviourString = "averaged";
				break;
			}
			case NEUTRAL_RATE_FIXED:{
				behaviourString = "0.0";
				if(mask.getNeutralRatio()>=0.0){
					behaviourString = ""+mask.getNeutralRatio();
				}
				break;
			}
			default:{
				behaviourString = "aggregated";
				break;
				
			}
			}
			outputBuffer.append(behaviourString+","+maskToString(mask.getPositions())+"\n");
		}
		// write buffer to file
		BufferedWriter writer = new BufferedWriter(new FileWriter(maskFile));
		writer.write(outputBuffer.toString());
		writer.close();
		return true;
	}
	
	private static String maskToString(int[] mask){
		StringBuffer buffer = new StringBuffer();
		for(int i=0;i<mask.length;i++){
			buffer.append(Character.forDigit(mask[i],10));
		}
		return buffer.toString();
	}
	
	private static String maskToString(boolean[] mask){
		StringBuffer buffer = new StringBuffer();
		for(int i=0;i<mask.length;i++){
			if(mask[i]){
				buffer.append(1);
			}else{
				buffer.append(0);
			}
		}
		return buffer.toString();
	}

	/**
	 * Instantiates a single mask_mid with the given range and a fixed neutral ratio
	 * @param ratio
	 * @param start
	 * @param end
	 * @param length
	 * @return
	 */
	public static TeaspoonMask initialiseMask(double ratio, int start, int end, int length){
		/* check length >= end > start */
		if( (start > end) || (end > length) ){
			try {
				throw new Exception("Mask bounds are nonsensical; start and/or end overrun.");
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return null;
			}
		}else{
			boolean[]  mask = new boolean[length];
			for(int sequenceIndex=start;sequenceIndex<=end;sequenceIndex++){
				// flip positions in this range to true
				mask[sequenceIndex] = true;
			}
			return new TeaspoonMask(RateEstimationBehaviour.NEUTRAL_RATE_FIXED,mask,ratio);
		}
	}
	
	/**
	 * Instantiates a single mask_mid with the given range and rate estimation behaviour
	 * @param behaviour
	 * @param start
	 * @param end
	 * @param length
	 * @return
	 */
	public static TeaspoonMask initialiseMask(RateEstimationBehaviour behaviour, int start, int end, int length){
		/* check length >= end > start */
		if( (start > end) || (end > length) ){
			try {
				throw new Exception("Mask bounds are nonsensical; start and/or end overrun.");
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return null;
			}
		}else{
			boolean[]  mask = new boolean[length];
			for(int sequenceIndex=start;sequenceIndex<=end;sequenceIndex++){
				// flip positions in this range to true
				mask[sequenceIndex] = true;
			}
			return new TeaspoonMask(behaviour,mask);
		}
	}

	/**
	 * Performs the <b>union</b> (position-wise OR) of two sets of mask positions, e.g.:
	 * <pre>
	 * mask_1 = {1,1,0,0,0}
	 * mask_2 = {0,0,0,0,1}
	 * output = {1,1,0,0,1}
	 * <pre>
	 * 
	 * Note that mask behaviour is given by the first mask.
	 * 
	 * @param mask_one
	 * @param mask_two
	 * @return
	 */
	public static TeaspoonMask combineMasksUnion(TeaspoonMask mask_one,	TeaspoonMask mask_two) throws NullPointerException{
		if(mask_one.getLength() != mask_two.getLength()){
			throw new NullPointerException("Mask region lengths are not of equal length!");
		}else{
			TeaspoonMask combined;
			boolean[] positions = new boolean[mask_one.getLength()];
			for(int i=0;i<positions.length;i++){
				positions[i] = ( mask_one.getPositions()[i] || mask_two.getPositions()[i] );
			}
			if(mask_one.estimationBehaviour == RateEstimationBehaviour.NEUTRAL_RATE_FIXED){
				combined = new TeaspoonMask(mask_one.estimationBehaviour,positions,mask_one.getNeutralRatio());
			}else{
				combined = new TeaspoonMask(mask_one.estimationBehaviour,positions);
			}
			return combined;
		}
	}
}
