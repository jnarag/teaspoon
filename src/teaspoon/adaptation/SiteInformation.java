package teaspoon.adaptation;

/**
 * SiteInformation contains information about a site.
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/jnarag">@jnarag</a>
 * @version 0.1
 */
public class SiteInformation {

	int locusIndex;
	double 	numberOfDerived;
	int polymorphismCase;
	boolean hasAncestral;
	boolean isBadSite;
	SiteObservations[] siteData;
	double dirichletProb; // dirchlet probability / only for pybas-bhatt method
	double silentProb;
	double replacementProb;
	double totalNumBases;
	
	/**
	 * Default no-arg constructor.
	 */
	public SiteInformation(){
		polymorphismCase=9;
		hasAncestral=false;
		numberOfDerived=0.0;
		siteData = new SiteObservations[4];		// create array of teaspoon.adaptation.SiteObservations - an object that stores all info
		isBadSite=false;
	}
}
