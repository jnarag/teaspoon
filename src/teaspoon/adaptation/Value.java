package teaspoon.adaptation;

public class Value {
	double Nadapt;
    double[] rl;
    double[] sl;
	double[] rm;
    double[] sm;
    double[] nr;
    double[] rh;
    double[] sh;
    double[] sh_rh;
    double[] rh_sh;
    double[] rl_sl;
    double[] adaptations;
    double[] Bstrap;
    double[] theta;
    double[] watterson_S;
    double[] watterson_pi;
    double[] tajimas_D;
	String row;
	String column;
	String dataset;
    int codons;
    double[] no_windows;
	
	public Value(int N){
		Nadapt=0.0;
		row = "NA";
		column = "NA";
		dataset = "NA";
        codons = 0;
        Bstrap = new double[N];
		rm = new double[N];
        sm = new double[N];
        nr = new double[N];
        adaptations = new double[N];
        rh = new double[N];
        sh = new double[N];
        rl = new double[N];
        sl = new double[N];
        sh_rh = new double[N];
        rh_sh = new double[N];
        rl_sl = new double[N];
        theta = new double[N];
        watterson_pi = new double[N];
        watterson_S = new double[N];
        tajimas_D = new double[N];
        no_windows = new double[N];
	}


	

}
