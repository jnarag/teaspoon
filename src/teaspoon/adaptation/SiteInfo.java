package teaspoon.adaptation;

public class SiteInfo {

	int locus;
	double 	Numderived;
	int Case;
	boolean hasans;
	boolean badSite;
	Obs[] data;
	double Dprob; // dirchlet probability / only for pybas-bhatt method
	double Sil;
	double Rep;
	double totalNumBases;
	
	public SiteInfo(){
		Case=9;
		hasans=false;
		Numderived=0.0;
		data = new Obs[4];		// create array of teaspoon.adaptation.Obs - an object that stores all info
		badSite=false;
	}
}
