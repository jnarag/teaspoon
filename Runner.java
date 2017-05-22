package teaspoon;

import teaspoon.adaptation.Value;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class Runner {
	
	public static void main(String[] args) {
		
		String file = args[0];
		String datalocation = args[1];
		double low1 = Double.valueOf(args[2]);
		double low2 = Double.valueOf(args[3]);
		double mid1 = Double.valueOf(args[4]);
		double mid2 = Double.valueOf(args[5]);
		double high1 = Double.valueOf(args[6]);
		double high2 = Double.valueOf(args[7]);
		int N = Integer.valueOf(args[8]);
		double nur = Double.valueOf(args[9]);
		
		double[] low = {low1,low2};
		double[] mid = {mid1,mid2};
		double[] high = {high1,high2};
		
		System.out.println("The program is running");
		System.out.println("This is the file location:  " + datalocation);
		System.out.println("This is the output filename:  " + file);
		System.out.println("The choise of low neutral range is:  " + nur);
		System.out.println("The Number of bootraps is:  " + N);
		System.out.println("INTERVALS");
		System.out.println("low bound:  "+low[0]+"--->"+low[1]);
		System.out.println("mid bound:  "+mid[0]+"--->"+mid[1]);
		System.out.println("high bound:  "+high[0]+"--->"+high[1]);
		
		//h1n1
	//	double[] nr00={0.347958181,0.262833676,0.261582734,0.271944923,0.271944923,0.271944923,0.236294093,0.358890701,0.358890701,0.358890701,0.160869565,1.115241636,0.809906292};
	//	double[] nr05={0.189444501,0.15912031,0.136080332,0.125984252,0.125984252,0.125984252,0.120689655,0.152284264,0.152284264,0.152284264,0.07960199,1.130434783,0.573099415};
	//	double[] nr10={0.18256366,0.130558184,0.11713367,0.166089965,0.166089965,0.166089965,0.112903226,0.15862069,0.15862069,0.15862069,0.123076923,1.28,0.640776699};
	//	double[] nr15={0.175203447,0.09784324,0.10302409,0.19047619,0.19047619,0.19047619,0.125348189,0.153225806,0.153225806,0.153225806,0.105263158,1.16,0.846153846};
	//	double[] nr20={0.183238636,0.114858706,0.113293051,0.171122995,0.171122995,0.171122995,0.127659574,0.142857143,0.142857143,0.142857143,0.129032258,1.25,0.685185185};
	//	double[] nr25={0.196911197,0.126237624,0.107142857,0.154929577,0.154929577,0.154929577,0.152777778,0.15942029,0.15942029,0.15942029,0,1.25,0.608695652};
	//	double[] nr30={0.24,0.117437722,0.090778098,0.189473684,0.189473684,0.189473684,0.071428571,0.153846154,0.153846154,0.153846154,0,3,0.631578947};
	//	double[] nr35={0.266666667,0.117437722,0.090778098,0.105263158,0.105263158,0.105263158,0.076271186,0.10989011,0.10989011,0.10989011,0,3,0.631578947};
	//	double[] nr40={0.276497696,0.119771863,0.101286174,0.063157895,0.063157895,0.063157895,0.076271186,0.10989011,0.10989011,0.10989011,0,0,0.526315789};
	//	double[] nr45={0.262910798,0.125498008,0.10862069,0.063157895,0.063157895,0.063157895,0.076271186,0.120481928,0.120481928,0.120481928,0,0,0.526315789};
	
		
		// h3n2
		double[] nr00={0.337528771,0.365306046,0.362118466,0.431014151,0.431014151,0.431014151,0.278018168,0.314742237,0.314742237,0.314742237,0.199342105,1.047104247,1.09495549};
		double[] nr05={0.140006516,0.126206897,0.155296629,0.21875,0.21875,0.21875,0.273537755,0.171428571,0.171428571,0.171428571,0.193353474,0.382022472,0.8};
		double[] nr10={0.102505695,0.100700848,0.16549559,0.183908046,0.183908046,0.183908046,0.343782929,0.136170213,0.136170213,0.136170213,0.201005025,0.192307692,0.962962963};
		double[] nr15={0.095502779,0.098823529,0.175540327,0.172248804,0.172248804,0.172248804,0.341470588,0.112149533,0.112149533,0.112149533,0.208092486,0.227272727,1};
		double[] nr20={0.113621262,0.08988764,0.170897833,0.225,0.225,0.225,0.294510908,0.135135135,0.135135135,0.135135135,0.171428571,0.357142857,0.722222222};
		double[] nr25={0.071672355,0.062593145,0.158385093,0.258064516,0.258064516,0.258064516,0.365313653,0.036363636,0.036363636,0.036363636,0.142857143,0.25,0.818181818};
		double[] nr30={0.066176471,0.073043478,0.153681964,0.259259259,0.259259259,0.259259259,0.33496572,0.042553191,0.042553191,0.042553191,0.166666667,0.375,0.777777778};
		double[] nr35={0.069230769,0.067542214,0.15982242,0.269230769,0.269230769,0.269230769,0.272409779,0.044444444,0.044444444,0.044444444,0.176470588,0.5,0.75};
		double[] nr40={0.098092643,0.056470588,0.183079057,0.318181818,0.318181818,0.318181818,0.266846361,0.066666667,0.066666667,0.066666667,0.176470588,0.75,2};
		double[] nr45={0.102739726,0.054830287,0.170616114,0.358974359,0.358974359,0.358974359,0.275766017,0.074074074,0.074074074,0.074074074,0.214285714,1,1.5};
		
		double[] nr= new double[13];
		if(nur==0){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr00[i];
			}
		}
		if(nur==0.05){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr05[i];
			}
		}
		if(nur==0.1){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr10[i];
			}
		}
		if(nur==0.15){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr15[i];
			}
		}
		if(nur==0.2){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr20[i];
			}
		}
		if(nur==0.25){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr25[i];
			}
		}
		if(nur==0.3){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr30[i];
			}
		}
		if(nur==0.35){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr35[i];
			}
		}
		if(nur==0.4){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr40[i];
			}
		}
		if(nur==0.45){
			for(int i=0;i<nr.length;i++){
				nr[i]=nr45[i];
			}
		} 
		
		
		FluAnalysis2 ff = new FluAnalysis2(datalocation,low,mid,high,nr);
		Value[][] v = ff.clumpanalysish3n2BS(N);
//		teaspoon.adaptation.Value[][] v = ff.clumpanalysish1n1BS(N);
		// Send formatted output to the file.
	      try{
	    	    // Create file 
	    	    FileWriter fstream = new FileWriter(file);
	    	        BufferedWriter out = new BufferedWriter(fstream);
	    			for(int x=0;x<N;x++){
	    				for(int i=0;i<v.length;i++){
	    					for(int j=0;j<v[0].length;j++){
	    						out.write(new Double(v[i][j].Bstrap[x]).toString());

	    						out.write("	");
	    					}
	    					out.write("\n");
	    				}
	    			}

	    	    //Close the output stream
	    	    out.close();
	    	    }catch (Exception e){//Catch exception if any
	    	      System.err.println("Error: " + e.getMessage());
	    	    }
	}

}
