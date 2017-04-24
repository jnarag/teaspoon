package teaspoon.adaptation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;


public class CoalescentSim {

	int sample_size;
	double theta;
	Random generator = new Random();
	final int numResamples = 1000;
	double totalComparisons = 0;
	Double n;
	public CoalescentSim(){
		// default no-arg constructor
		throw new RuntimeException("ERROR: please input the number of samples");
	}

	public CoalescentSim(int NumberOfSamples, double thetavalue){
		sample_size=NumberOfSamples;
		theta = thetavalue;
		n = new Double(sample_size);
		totalComparisons = (n*(n-1.0))/2.0;
	}			
	// create a genealogy in newick format
	public ArrayList NewickGenealogy(){
		ArrayList<String> Vk = new ArrayList<String>();	// create vector of initial nodes
		ArrayList<String> PointStorage = new ArrayList<String>();
		for(int i= 0;i<sample_size;i++){	
			Integer tmp = new Integer(i+1);
			Vk.add(tmp.toString());
		}
		Iterator<String> initialit = Vk.iterator();
		StringBuffer initialsbuf = new StringBuffer();
		while (initialit.hasNext()){
			initialsbuf.append(initialit.next()).append("  ");
		}
		PointStorage.add(initialsbuf.toString());

//		start simulation
		for(int k=0;k<sample_size-1;k++){		
			Collections.shuffle(Vk);	// randomise
			ArrayList Vkshuffle = (ArrayList) Vk.clone();	// clone
			// pick two random numbers ik and jk from Vk
			Vk.remove(Vkshuffle.get(0));	// pick ik
			Vk.remove(Vkshuffle.get(1));	// pick jk
			// make subset of the two.
			StringBuffer sbuf = new StringBuffer();
			sbuf.append("(").append(Vkshuffle.get(0)).append(",").append(Vkshuffle.get(1)).append(")");
			Vk.add(sbuf.toString());	// concatanate
			Iterator<String> it = Vk.iterator();
			StringBuffer sbuf2 = new StringBuffer();
			while (it.hasNext()){
				sbuf2.append(it.next()).append("  ");
			}
			PointStorage.add(sbuf2.toString());
		}
		return PointStorage;
	}
	// create a genealogy
	public Node[] Genealogy(){
		Node[] tree = new Node[2*sample_size-1];
		Node[] list = new Node[sample_size];
		// initialise tree full of nodes
		for(int in=0;in<2*sample_size-1;in++){
			tree[in] = new Node(in);
		}
		for(int in=0;in<sample_size;in++){
			list[in] = tree[in];
		}
		double t=0;
		for(int in = sample_size; in >1;in--){
			t += -2.0 * Math.log(1.0-generator.nextDouble())/ ((double) in*((double) in-1.0));
			tree[2*sample_size-in].time = t;
		}

		for(int in = sample_size;in>1;in--){		
			double pick = in*generator.nextDouble();	// pick an index from list
			list[(int)pick].ancestor = tree[2*sample_size - in]; //set ancestor
			tree[2*sample_size - in].desc1 = list[(int) pick];	// add the random descendent
			list[(int) pick] = list[in-1];	// make node unpickable
			pick = (in-1)*generator.nextDouble();	// pick new 
			list[(int)pick].ancestor = tree[2*(int)sample_size - in];	//set ancestor
			tree[2*sample_size - in].desc2 = list[(int) pick]; // add the random descendent
			list[(int) pick] = tree[2*sample_size - in];  // update pickable nodes
		}
		// add mutation to genealogy
		for(int node=0;node<2*sample_size-2;node++){
			double time = 0;
			time = tree[node].ancestor.time - tree[node].time;
			tree[node].nmuts = generator.nextPoisson(time*(theta)/2.0);
		}

		return tree;	

	}

	public int[][] incidenceMatrix(){
		// creates and incidence matrix. each new column is a new mutation
		// each row is a sample. if there is a 1 present then that specific mutation is present, otherwise there is a zero
		Node[] tree = Genealogy();
		int k=0;
		int total_muts = 0;
		for(int node=0;node<2*sample_size-2;node++){
			total_muts += tree[node].nmuts;
		}
		int[][] incidence = new int[sample_size][total_muts];
		for(int in=0;in<tree.length;in++){
			if(tree[in].nmuts>0){	// if mutation has occured on this branch
				for(int i=0;i<tree[in].nmuts;i++){
					int[] vec = new int[sample_size];
					findDesc(tree[in],vec);
					for(int j=0;j<sample_size;j++){
						incidence[j][k] = vec[j];
					}
					k++;
				}
			}

		}
		return incidence;
	}


//	estimate quantities from coalescent	
	public double[] CoalescentEstimators(){
		// methods for finding estimates are the same as those in Tajima.java and teaspoon.adaptation.FuAndLi.java
		double pwd =0; double sing =0; 
		double[] est = new double[5];
		int[][] incidence = incidenceMatrix();
		est[0] = incidence[0].length;  //number of segregating sites
		double Es = 0.0;
		for(int i=1;i< sample_size;i++){
			Es += 1.0/(i);		// harmonic series estimated from coalescent
		}
		est[1] = est[0]/Es;		//theta-s
		int len = incidence.length;
		for(int j=0;j<incidence[0].length;j++){
			int temp=0;
			for(int i=0;i<incidence.length;i++){
				if(incidence[i][j]==1){
					temp++;
				}
			}
			if(temp ==1){sing++;}
			pwd += temp*(len-temp);
		}
		est[2] = pwd/totalComparisons;
		double thetan = ((sample_size-1.0)/(sample_size))*sing;
		est[3] = thetan;
		est[4] = sing;
		return est;
	}

	public double[] TajimaCS(){
		// stores the tails of the distributions 0.05% and 99.5%
		// in a size 2 array. 0 and 1 are the D tail values
		double[] tailvalues = new double[2];
		double[] store = new double[numResamples];
		int run = 0;

		if(theta == 0 ){
			tailvalues[0] = -1;  
			tailvalues[1] = 1;
		} else {

			while(run<numResamples){
				double D=0;
				double[] Est = CoalescentEstimators();
				if(Est[0]==0.0){
					store[run] = 0.0;
					run++;
				} else {
//					calculation is the same as in teaspoon.adaptation.DiversityStats.java
					double a1 = 0.0; double a2 = 0.0;
					int n=sample_size;
					for(double i=1;i<n;i++){
						a1 += (1.0/i);
						a2 += (1.0/(Math.pow(i,2)));
					}
					double b1 = (n+1.0)/(3.0*(n-1.0));
					double b2 = (2.0*(Math.pow(n,2) + n + 3.0))/(9.0*n*(n-1.0));
					double c1 = b1 - (1.0/a1);
					double c2 = b2 - ((n+2.0)/(a1*n)) + (a2/Math.pow(a1, 2));
					double e1 = c1/a1;
					double e2 = c2/(Math.pow(a1,2) + a2);
					D = (Est[2] - Est[1])/(Math.sqrt((e1*Est[0])+(e2*Est[0]*(Est[0]-1.0))));

					store[run] = D;
					run++;
				}
			}
			Arrays.sort(store);

			
			tailvalues[0] = store[(int)(0.025*numResamples)];
			tailvalues[1] = store[(int)(0.975*numResamples)];			
		}
		return tailvalues;
	}

//	FIX
	public double[] FuAndLiCS(){
		// stores the tails of the distributions 0.05% and 99.5%
		// in a size 4 array. 0 and 1 are the D tail values and 2 and 3 are the F tail values
		double[] tailvalues = new double[4];  // stores tail values
		double[] storeD = new double[numResamples];
		double[] storeF = new double[numResamples];
		int run = 0;

		if(theta == 0 ){
			tailvalues[0] = -1;
			tailvalues[1] = 1;
			tailvalues[2] = -1;
			tailvalues[3] = 1;
		} else {
			while(run<numResamples){ // run numResamples times
				double[] Est = CoalescentEstimators(); // find new estimators for every run
				if(Est[0]==0){
					storeD[run] = 0;
					storeF[run] = 0;
					run++;
				} else {
				// calculation is the same as in teaspoon.adaptation.FuAndLi.java
				double Dstar = 0.0; double Fstar=0.0;
				double an = 0.0; double bn = 0.0; double an1 = 0.0;
				for(double i=1;i<n;i++){
					an += 1.0/i;
					bn += 1.0/(i*i);
				}
				for(double i=1;i==n;i++){
					an1 += 1.0/i;
				}
				// Dstar terms **********
				double betaterm1 = 1.0/(Math.pow(an, 2.0)+bn);
				double betaterm2 = bn/Math.pow(an,2.0);
				double betaterm3 = (2.0/n)*(1.0+(1.0/an)-an+(an/n));
				double betaD = betaterm1*(betaterm2-betaterm3-(1.0/(n*n)));
				double alphaD = (1.0/an)*(((n-1.0)/n) - (1.0/an)) - betaD;
				// calcualtion of Dstar
				Dstar = (Est[1]-Est[3])/(Math.sqrt((alphaD*Est[0])+(betaD*Est[0]*Est[0])));
				// Fstar terms  ***********
				double betafterm1 = 1.0/((n*n)+bn);
				double betafterm2 = ((2.0*n*n*n) + ((n*n)*110.0) - (255.0*n) + 153.0)/(n*n*9.0 * (n-1.0));
				double betafterm3 = (2.0*(n-1.0)*an)/Math.pow(n,2.0);
				double betaF = betafterm1 * (betafterm2 + betafterm3 + ((8.0*bn)/n));
				double alphafterm1 = ((4.0*n*n) + (19.0*n) + 3.0 - ((12.0*(n+1.0))*an1))/(3.0*n*(n-1.0));
				double alphaF = ((1.0/an) * alphafterm1) - betaF;
				// calculation of Fstar
				Fstar = (Est[2]-Est[3])/(Math.sqrt((alphaF*Est[0])+(betaF*Est[0]*Est[0])));

				storeD[run] = Dstar;
				storeF[run] = Fstar;
				run++;
			}
			}

			Arrays.sort(storeD);
			Arrays.sort(storeF);
			tailvalues[0] = storeD[(int)(0.025*numResamples)];
			tailvalues[1] = storeD[(int)(0.975*numResamples)];
			tailvalues[2] = storeF[(int)(0.025*numResamples)];
			tailvalues[3] = storeF[(int)(0.975*numResamples)];
		}
		return tailvalues;
	}

	// takes actual values. the first is  D and the second H
	public double[] FayAndWuCS(double[] actualvalue){
		double[] counter = new double[sample_size];
		double countD = 0;
		double countH = 0;
		for(int run=0;run<numResamples;run++){
			int[][] incidence = incidenceMatrix();
			for(int j=0;j<incidence[0].length;j++){
				int num = 0;
				for(int i=0;i<n;i++){
					num += incidence[i][j];
				}
				counter[num]++;
			}
			double thetaw = 0; // number of segregating sites
			double part1 = 0;
			double part2 = 0;
			// loop goes from 1 -> n-1 where n is the sample size
			// counter goes from 0 -> n-2. which traverses all entries except the last - equvalent to n-1
			for(int k=1;k<n;k++){
				part1 += counter[k-1];
				part2 += 1.0/k;
			}
			thetaw = part1*(1.0/part2);
			double thetak = 0; // number of pairwise differences
			for(int k=1;k<n;k++){
				thetak += (2.0*counter[k-1]*k*(n-k))/(n*(n-1.0));
			}
			double thetah = 0; // weighted
			for(int k=1;k<n;k++){
				thetah += (2.0*counter[k-1]*k*k)/(n*(n-1.0));
			}
			double D = thetak - thetaw;
			double H = thetak - thetah;
			if(actualvalue[0]<D){
				countD++;
			}
			if(actualvalue[1]<H){
				countH++;
			}	
		}
		double Dvalue = countD/numResamples;
		double Hvalue = countH/numResamples;
		double[] values = {Dvalue,Hvalue};
		return values;
	}
	// recursive FUNCTIONs;
	// returns a vector showing which terminal nodes are connected to a node
	// recursive method.
	public int[] findDesc(Node node,int[] isTerminal){
		if(node.desc1 ==null || node.desc2 ==null){	// if node is terminal
			isTerminal[node.index] = 1;
			return isTerminal; 
		}
		// recusion to recuse down nodes
		findDesc(node.desc1,isTerminal);			
		findDesc(node.desc2,isTerminal);
		return isTerminal;
	}

//	estimate quantities from coalescent	
	public double[] CoalescentEstimatorsMean(){
		// methods for finding estimates are the same as those in Tajima.java and teaspoon.adaptation.FuAndLi.java

		double Es = 0.0;
		for(int i=1;i< sample_size;i++){
			Es += 1.0/(i);		// harmonic series estimated from coalescent
		}
		double[] est = new double[5];
		for(int run=0;run<20;run++){
			double pwd =0; double sing =0; 
			int[][] incidence = incidenceMatrix();
			double numS =  incidence[0].length;  
			est[0] += incidence[0].length;  //number of segregating sites
			est[1] += numS/Es;
			int len = incidence.length;
			for(int j=0;j<incidence[0].length;j++){
				int temp=0;
				for(int i=0;i<incidence.length;i++){
					if(incidence[i][j]==1){
						temp++;
					}
				}
				if(temp ==1){sing++;}
				pwd += temp*(len-temp);
			}
			est[2] += pwd/totalComparisons;
			double thetan = ((sample_size-1.0)/(sample_size))*sing;
			est[3] += thetan;
			est[4] += sing;
		}
		for(int i=0;i<est.length;i++){
			est[i] = est[i]/numResamples;
		}
		return est;
	}



	// knuth algorithm
	// generates poisson random variable with mean lambda
	public double poissonRV(double lambda){
		double L = Math.exp(-lambda); double k=0.0; double p=1.0;
		while(p >= L){k=k+1.0;p = p*generator.nextDouble();}
		return k-1.0;			
	}


}















