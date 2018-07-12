/**
 * 
 */
package teaspoon.app.utils;

import java.io.File;

import teaspoon.app.TeaspoonBootstrap;
import teaspoon.app.TeaspoonMask;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * The 'full' site matrix contains all positions in an input main- or ancestral-file
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattAdaptationFullSiteMatrix {

	// holds all the data for a file //
	int[][] siteMatrix;
	
	/**
	 * constructor takes an 
	 * @param readFASTA
	 */
	public BhattAdaptationFullSiteMatrix(int[][] fasta) {
		// TODO Auto-generated constructor stub
		siteMatrix = fasta;
	}

	/**
	 * @return the siteMatrix
	 */
	public int[][] getSiteMatrix() {
		return siteMatrix;
	}

	/**
	 * @param siteMatrix the siteMatrix to set
	 */
	public void setSiteMatrix(int[][] siteMatrix) {
		this.siteMatrix = siteMatrix;
	}

	/**
	 * Loads input to int matrix
	 */
	public void loadAlignmentFile(File input){}
	
	/**
	 * Given a list of mask_mid positions, returns submatrix.
	 * Mask may be shorter than alignment, but not longer.
	 * @param mask_mid
	 * @return
	 */
	public BhattAdaptationFullSiteMatrix subsampleByMask(TeaspoonMask mask){
		// check mask_mid isn't longer than the input alignment
		if(mask.getLength()>this.alignmentLength()){
			throw new ArrayIndexOutOfBoundsException("The alignment subsampling mask_mid ("+mask.getLength()+") is longer than the alignment("+this.alignmentLength()+")!");
		}
		boolean[] maskValues = mask.getPositions();
		int sampledSites = mask.getNumberOfValidPositions();
		int[][] newSiteMatrix = new int[this.numberOfSequences()][sampledSites];
		// walk through each taxon in the alignment picking only positions with mask_mid==true
		for(int taxon=0;taxon<this.numberOfSequences();taxon++){
			int[] subsampledTaxon = new int[sampledSites];
			int subsampleCounter = 0;
			// walk through positions
			for(int maskPosition=0;maskPosition<mask.getLength();maskPosition++){
				if(maskValues[maskPosition]){
					subsampledTaxon[subsampleCounter] = siteMatrix[taxon][maskPosition];
					subsampleCounter++;
				}
			}
			// put the subsampled sequence in the array
			newSiteMatrix[taxon] = subsampledTaxon;
		}
		return new BhattAdaptationFullSiteMatrix(newSiteMatrix);
	}
	
	/**
	 * Return a single bootstrap replicate.
	 * Bootstrap is assumed to be masked e.g. generated with TeaspoonBootstrap(mask_mid)
	 * @param bootstrap specifying which positions of the parent alignment to sample, and how many times
	 * @return
	 */
	public BhattAdaptationFullSiteMatrix subsampleByBootstrap(TeaspoonBootstrap bootstrap){
		// get the array telling us how much to sample by
		int[] samplingIntensities = bootstrap.getBootstraps();
		int lengthOfNewSequences = bootstrap.getNumberOfValidPositions();
		int[][] newMatrix = new int[this.numberOfSequences()][lengthOfNewSequences];
		// walk through each taxon
		for(int taxon = 0;taxon<this.numberOfSequences(); taxon++){
			// make a new array
			newMatrix[taxon] = new int[lengthOfNewSequences];
			// walk through the bootstrap copying sites in proportion to intensity
			int newSequencePositionIndex = 0;
			
			for(int oldSequencePosition = 0; oldSequencePosition<this.alignmentLength();oldSequencePosition++){
				int intensity = samplingIntensities[oldSequencePosition];				
				while(intensity>0){
					newMatrix[taxon][newSequencePositionIndex] = siteMatrix[taxon][oldSequencePosition];
					intensity --;
					newSequencePositionIndex++;
				}
			}
		}
		return new BhattAdaptationFullSiteMatrix(newMatrix);
	}
	
	
	/**
	 * Get a consensus
	 * @see {@link MainAlignmentParser#consensusArray()}
	 */
	public int[] deriveConsensus(){
		return MainAlignmentParser.consensusArray(siteMatrix);
	}
	
	/**
	 * Vertically joins two matrices e.g. appends taxa to an alignment.
	 * @param additionalSequences
	 * @throws ArrayIndexOutOfBoundsException if the two alignments do not have the same sequence lengths
	 * @return 
	 */
	public BhattAdaptationFullSiteMatrix appendTaxa(BhattAdaptationFullSiteMatrix additionalSequences) throws ArrayIndexOutOfBoundsException{
		int existingTaxa = this.numberOfSequences();
		int appendedTaxa = additionalSequences.numberOfSequences();
		int existingLength = this.alignmentLength();
		int appendedLength = additionalSequences.alignmentLength();
		if(appendedLength != existingLength){
			throw new ArrayIndexOutOfBoundsException("Alignments lengths to append do not match");
		}else{
			// should be safe to append the new sequences to the old matrix as rows starting at row 0.
			int[][] newAlignment = new int[existingTaxa+appendedTaxa][existingLength];
			// first copy existing sites
			for(int i=0;i<existingTaxa;i++){
				newAlignment[i] = siteMatrix[i];
			}
			// now copy new sequence rows, starting at row [existingTaxa]
			for(int i=0;i<appendedTaxa;i++){
				newAlignment[i+existingTaxa] = additionalSequences.siteMatrix[i];
			}
			// return new FSM
			return new BhattAdaptationFullSiteMatrix(newAlignment);
		}
	}

	/**
	 * @return length of the alignment (number of positions)
	 */
	public int alignmentLength() {
		return this.siteMatrix[0].length;
	}
	
	/**
	 * @return depth of the alignment (number of taxa)
	 */
	public int numberOfSequences(){
		return this.siteMatrix.length;
	}
}
