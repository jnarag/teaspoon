package teaspoon.adaptation;

import java.util.*;
public class Mutation {
	

	ArrayList<Double> H = new ArrayList<Double>();
	ArrayList<Double> L = new ArrayList<Double>();
	int site;

	int base;

	ArrayList<Double> timeProb = new ArrayList<Double>();
	public Mutation(){

		
	}
	
	public void addProb(double val){
		timeProb.add(val);
	}
	
	public void addH(double val){
		H.add(val);
	}
	
	public void addL(double val){
		L.add(val);
	}
	

}
