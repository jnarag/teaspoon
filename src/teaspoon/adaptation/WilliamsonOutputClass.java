package teaspoon.adaptation;

import java.util.ArrayList;


public class WilliamsonOutputClass {
	ArrayList<Double> LowA;
	ArrayList<Double> HighA;
	ArrayList<Double> TotalA;
	ArrayList<Double> FixA;
	ArrayList<Double> MidA;
	ArrayList<Double> Low_S;
	ArrayList<Double> Low_R;
	ArrayList<Double> Mid_S;
	ArrayList<Double> Mid_R;
	ArrayList<Double> High_S;
	ArrayList<Double> High_R;
	ArrayList<Double> Fix_S;
	ArrayList<Double> Fix_R;
	ArrayList<double[]> BS;
	ArrayList<Double> LC;
	ArrayList<Double> HC;
	ArrayList<Double> dateavg;
	ArrayList<Integer> numSamples;
	double neut;
	double L;
	
	WilliamsonOutputClass(int X,int N){
		this.dateavg = new ArrayList<Double>();
		this.LC = new ArrayList<Double>();
		this.numSamples = new ArrayList<Integer>();
		this.HC = new ArrayList<Double>();
		this.LowA = new ArrayList<Double>();
		this.HighA = new ArrayList<Double>();
		this.TotalA = new ArrayList<Double>();
		this.FixA = new ArrayList<Double>();
		this.MidA = new ArrayList<Double>();
		this.Low_S = new ArrayList<Double>();
		this.Low_R = new ArrayList<Double>();
		this.Mid_S = new ArrayList<Double>();
		this.Mid_R = new ArrayList<Double>();
		this.High_S = new ArrayList<Double>();
		this.High_R = new ArrayList<Double>();
		this.Fix_S = new ArrayList<Double>();
		this.Fix_R = new ArrayList<Double>();
		BS = new ArrayList<double[]>();
		for(int i=0;i<X;i++){
			double[] tmp = new double[N];
			BS.add(tmp);
		}
		this.neut=0;
		this.L=0;
	}
	
	WilliamsonOutputClass(){
		this.dateavg = new ArrayList<Double>();
		this.LC = new ArrayList<Double>();
		this.numSamples = new ArrayList<Integer>();
		this.HC = new ArrayList<Double>();
		this.LowA = new ArrayList<Double>();
		this.HighA = new ArrayList<Double>();
		this.TotalA = new ArrayList<Double>();
		this.FixA = new ArrayList<Double>();
		this.MidA = new ArrayList<Double>();
		this.Low_S = new ArrayList<Double>();
		this.Low_R = new ArrayList<Double>();
		this.Mid_S = new ArrayList<Double>();
		this.Mid_R = new ArrayList<Double>();
		this.High_S = new ArrayList<Double>();
		this.High_R = new ArrayList<Double>();
		this.Fix_S = new ArrayList<Double>();
		this.Fix_R = new ArrayList<Double>();
		this.neut=0;
		this.L=0;
	}

	
	
}
