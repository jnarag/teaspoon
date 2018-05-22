package teaspoon.adaptation;

public class SiteObservations {

	int base;  //ans is 5
	double numObs;
	double rawNumObs;
	double prior;
	boolean inAncestral;

	public SiteObservations(){
		base=9;
		numObs=0.0;
		rawNumObs=0.0;
		prior=1.0; // change if need
		inAncestral=false;
	}
}
