package teaspoon.adaptation;

import java.util.Calendar;
import java.util.Random;

import teaspoon.app.utils.TeaspoonMethods;
import cc.mallet.types.Dirichlet;

public class SiteEstMulti  {

	private static Random rng = new Random(Calendar.getInstance().getTimeInMillis() + Thread.currentThread().getId());

	double N = 5.0; //set number of replicates
	Double number_bins;
	boolean[] WhichBins;
	double[][] bins;
	int NumSample;
	public final  int[][] integer_matrix;
	public final int[] integer_ancestral;

	public  int[] codon_ancestral;
	public  int[][] codon_matrix;
	public  boolean[] bad_sites_list;
	TeaspoonMethods preprocess = new TeaspoonMethods();


	public final String[] AA =	{"K","numReplicates","K","numReplicates","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H","P","P","P","P",
			"R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V",
			"X","Y","X","Y","S","S","S","S","X","C","W","C","L","F","L","F","?","-","?" };

	public final double[] F3X4 ={0.00384245,0.03020593,0.02497593,0.01702419,0.0030501,0.02397721,0.01982568,0.01351366,0.00169545,
			0.01332811,0.01102042,0.00751178,0.00355278,0.02792876,0.02309304,0.01574077,0.00388819,0.03056552,0.02527326,0.01722686,
			0.00308642,0.02426265,0.0200617,0.01367453,0.00171563,0.01348678,0.01115161,0.00760121,0.00359507,0.02826125,0.02336796,
			0.01592816,0.00616393,0.04845534,0.04006555,0.02730964,0.00489288,0.03846344,0.03180369,0.02167816,0.00271978,0.02138052,
			0.01767859,0.01205015,0.00569924,0.04480239,0.03704508,0.02525082,0,0.0188787,0,0.01064012,0.00190632,0.01498576,0.01239105,
			0.00844604,0,0.00833007,0.00688776,0.00469486,0.00222048,0.01745548,0.01443315,0.00983798,0.0,0.0,0.0};

	public final double[] ENC = {0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,
			0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,
			0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,
			0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,
			0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0,0.016393443,0,0.016393443,0.016393443,0.016393443,0.016393443,
			0.016393443,0,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443,0.016393443};

	public SiteEstMulti(){
		// default no-arg constructor
		throw new RuntimeException("please input the raw intger matrix and the ancestral matrix");
	}
	public SiteEstMulti(int[][] m, int[] a){
		// creates good matricies which contain no gaps or sequencing errors
		TeaspoonMethods creatematrix = new TeaspoonMethods();
		this.integer_matrix = m;
		this.integer_ancestral = a;
		this.codon_ancestral = creatematrix.makeCodon(a);
		this.codon_matrix = creatematrix.make_codon(m);
		//	integerMatrix = creatematrix.get_object(m, a, "base");
		//	integerAncestralArray = creatematrix.get_ancestral_object(m,a,"base");
		//	codon_ancestral = creatematrix.get_ancestral_object(m, a, "codon");	
		bad_sites_list = creatematrix.InvalidSites(integer_matrix, integer_ancestral);	
		NumSample = integer_matrix.length;
		//		System.out.println("MCMC samples: "+numReplicates);
	}

	public SiteEstMulti(int[][] m, int[] a, boolean[] badsites){
		// creates good matricies which contain no gaps or sequencing errors
		TeaspoonMethods creatematrix = new TeaspoonMethods();
		this.integer_matrix = m;
		this.integer_ancestral = a;
		this.codon_ancestral = creatematrix.makeCodon(a);
		this.codon_matrix = creatematrix.make_codon(m);
		this.bad_sites_list = badsites;
		NumSample = integer_matrix.length;
	}



	// set number of binsMatrix
	public double SetNumBins(double value){
		number_bins=value;
		return number_bins;
	}

	private boolean[] SetBins(boolean[] NeutralVec){
		WhichBins = NeutralVec;
		return WhichBins;
	}

	private double[][] SetBinsValues(double[][] binVec){
		bins = binVec;
		return bins;
	}

	public String AminoFinder(int codon){
		int[] bases = preprocess.codonSplitter(codon);
		int codonnum = getcodonnumber(bases[0], bases[1], bases[2]);
		return AA[codonnum];
	}



	public SiteInformation SiteInformation(int site){
		SiteInformation SI = new SiteInformation();
		double numbase = 4;	// number of bases 

		SiteObservations[] data = new SiteObservations[(int)numbase];		// create array of teaspoon.adaptation.SiteObservations - an object that stores all info
		for(int i=0;i<numbase;i++){
			data[i] = new SiteObservations();			// initialises the teaspoon.adaptation.SiteObservations objects
		}
		for(int i=0;i<numbase;i++){
			data[i].base=i+1;
			data[i].rawNumObs = preprocess.num_of_base(integer_matrix, i+1, site);
			data[i].numObs = data[i].rawNumObs / (double) integer_matrix.length;
			data[integer_ancestral[site]-1].inAncestral=true;	
			//tests if base is ansestral
		}

		// store siteData in SI object
		SI.siteData = data;

		for(int i=0;i<numbase;i++){
			if(SI.siteData[integer_ancestral[site]-1].numObs!=0.0){
				SI.hasAncestral=true;		// test if site has ansestralbase
			}
			if(SI.siteData[i].numObs!=0.0 && data[i].inAncestral==false) {
				SI.numberOfDerived++; // site has no ancestral bases
			}
		}

		if (SI.numberOfDerived==0 && SI.hasAncestral==true){
			SI.polymorphismCase=1;// invariat
		}
		else if (SI.numberOfDerived==1 && SI.hasAncestral==false){
			SI.polymorphismCase=2;// fixed
		}
		else if (SI.numberOfDerived==1 && SI.hasAncestral==true){
			SI.polymorphismCase=3;// 1 state derived and ans
		}
		else if (SI.numberOfDerived==2 && SI.hasAncestral==false){
			SI.polymorphismCase=4;// 2 state derived no ans
		}
		else if (SI.numberOfDerived==2 && SI.hasAncestral==true){
			SI.polymorphismCase=5;// 2 state derived and ans
		}
		else if (SI.numberOfDerived==3 && SI.hasAncestral==false){
			SI.polymorphismCase=6;// 3 state derived no ans
		}
		else if (SI.numberOfDerived==3 && SI.hasAncestral==true){
			SI.polymorphismCase=7;// 3 state derived and ans
		}

		return SI;
	}

	public SiteInformation SiteInformation_IncludeMainGaps(int site){
		SiteInformation SI = new SiteInformation();
		double numbase = 4;	// number of bases 

		SiteObservations[] data = new SiteObservations[(int)numbase];		// create array of teaspoon.adaptation.SiteObservations - an object that stores all info
		for(int i=0;i<numbase;i++){
			data[i] = new SiteObservations();			// initialises the teaspoon.adaptation.SiteObservations objects
		}
		double TotalNumBases=0.0;
		for(int i=0;i<numbase;i++){
			data[i].base=i+1;
			data[i].rawNumObs = preprocess.num_of_base(integer_matrix, i+1, site);
			TotalNumBases+=data[i].rawNumObs;
			data[integer_ancestral[site]-1].inAncestral=true;	//tests if base is ansestral			
		}
		// calculates the number of observed ignoring gaps i.e as a sum of the total number of bases
		for(int i=0;i<numbase;i++){
			data[i].numObs = data[i].rawNumObs/TotalNumBases;
		}
		SI.totalNumBases=TotalNumBases;
		// store siteData in SI object
		SI.siteData = data;

		for(int i=0;i<numbase;i++){
			if(SI.siteData[integer_ancestral[site]-1].numObs!=0.0){
				SI.hasAncestral=true;		// test if site has ansestralbase
			}
			if(SI.siteData[i].numObs!=0.0 && data[i].inAncestral==false) {
				SI.numberOfDerived++; // site has no ancestral bases
			}
		}

		if (SI.numberOfDerived==0 && SI.hasAncestral==true){
			SI.polymorphismCase=1;// invariat
		}
		else if (SI.numberOfDerived==1 && SI.hasAncestral==false){
			SI.polymorphismCase=2;// fixed
		}
		else if (SI.numberOfDerived==1 && SI.hasAncestral==true){
			SI.polymorphismCase=3;// 1 state derived and ans
		}
		else if (SI.numberOfDerived==2 && SI.hasAncestral==false){
			SI.polymorphismCase=4;// 2 state derived no ans
		}
		else if (SI.numberOfDerived==2 && SI.hasAncestral==true){
			SI.polymorphismCase=5;// 2 state derived and ans
		}
		else if (SI.numberOfDerived==3 && SI.hasAncestral==false){
			SI.polymorphismCase=6;// 3 state derived no ans
		}
		else if (SI.numberOfDerived==3 && SI.hasAncestral==true){
			SI.polymorphismCase=7;// 3 state derived and ans
		}

		return SI;
	}


	public SiteInformation DirichletSRFreq(double u, double v,int site,double[] prior,boolean needprior) {
		SiteInformation Info = SiteInformation(site);
		double numbase = 4;	// number of bases 
		double[] observations = new double[(int)numbase];
		double[] sampler = new double[(int) N];

		if(needprior){
			for(int i=0;i<numbase;i++){		
				if(Info.siteData[i].inAncestral){
					prior[i]=1.0;
				}else{
					prior[i]=1.0/3.0;
				}
			}
		}		
		for(int i=0;i<numbase;i++){		
			Info.siteData[i].prior = prior[i];
			observations[i] = Info.siteData[i].rawNumObs+Info.siteData[i].prior;	// add observations into an array to find Dirichlet(01+valuesToSampleFrom....0k+valuesToSampleFrom) where valuesToSampleFrom is the prior
		}
//	Dirichlet D = new Dirichlet(observations);
		TeaspoonRandomSampler S = new TeaspoonRandomSampler(observations);

		double count = 0;
		for(int i=0;i<(int) N;i++){
	//		double[] point = D.nextDistribution();	
			double[] point = S.sampleDirichlet();
			for(int x=0;x<numbase;x++){
				if(Info.siteData[x].inAncestral==false){
					sampler[i]+=point[x];
				}
			}	
			if(sampler[i]>u && sampler[i]<=v){
				count++;
			}
		}

		Info.dirichletProb = (count/N);
		return Info;
	}

	public SiteInformation BetaSiteFreq(double u, double v,int site){
		SiteInformation Info = SiteInformation(site);
		double numbase = 4;	// number of bases 
		double derivedPrior = 1;
		double NonderivedPrior = 1;
		double[] observations=new double[2];
		observations[0]=0;
		observations[1]=0;
		if(Info.polymorphismCase==1){
			observations[0] = 0;	// add observations into an array to find Dirichlet(01+valuesToSampleFrom....0k+valuesToSampleFrom) where valuesToSampleFrom is the prior
		}
		else {
			for(int i=0;i<numbase;i++){			
				if(Info.siteData[i].inAncestral==false && Info.siteData[i].rawNumObs>0){
					observations[0] = Info.siteData[i].rawNumObs;	// add observations into an array to find Dirichlet(01+valuesToSampleFrom....0k+valuesToSampleFrom) where valuesToSampleFrom is the prior
				}
			}
		}
		observations[1] = (NumSample)-observations[0];
		observations[0] = observations[0]+derivedPrior;
		observations[1] = observations[1]+NonderivedPrior;
		Dirichlet D = new Dirichlet(observations);
		double count = 0;
		for(int i=0;i<N;i++){
			double[] point = D.nextDistribution();	
			if(point[0]>u && point[0]<v){
				count++;
			}

		}
		Info.dirichletProb = count/N;
		return Info;
	}


	public double[] BetaParameters(double mean, double n){
		double a = mean*n;
		double b = n-(mean*n);
		double[] val = {a,b};
		return val;

	}

	public double[] DirichletParameters(double[] means, double n){
		double[] val = new double[means.length];
		for(int i=0;i>means.length;i++){
			val[i] = n*means[i];
		}
		return val;
	}

	public double[] Degeneracy(){
		double count1fold=0;double count2fold=0;double count3fold=0;
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				double z=0;
				// pos 1 ********************************************************************
				SiteInformation Info = SiteInformation(site);
				if(Info.polymorphismCase!=1){z++;}
				Info = SiteInformation(site+1);
				if(Info.polymorphismCase!=1){z++;}
				Info = SiteInformation(site+2);
				if(Info.polymorphismCase!=1){z++;}
				if(z==1){
					count1fold++;
				} else if(z==2){
					count2fold++;
				} else if(z==3){
					count3fold++;
				}	
			}
		}
		double[] val = new double[3];
		val[0]=count1fold;
		val[1]=count2fold;
		val[2]=count3fold;
		return val;
	}

	//	*********************************************************************************
	//	nei gojobori method	
	public double[] NGpossible(){
		// main NG method for all sites
		double[][] possible = SRcodonInvariant();
		double[] identity = new double[2];
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				SiteInformation Info = SiteInformation(site);
				SiteInformation Info1 = SiteInformation(site+1);
				SiteInformation Info2 = SiteInformation(site+2);
				double sil=0;
				double rep=0;
				// find SR for pos1 ********************************************

				for(int i=0;i<integer_matrix.length;i++){
					int tmp = getcodonnumber(integer_matrix[i][site],integer_matrix[i][site+1],integer_matrix[i][site+2]);
					if(Double.isNaN(possible[tmp][0])==false){// && Info.Case==1){
						sil += possible[tmp][0];
						rep += 1.0-possible[tmp][0];
					}
				}
				identity[0] += sil/integer_matrix.length;
				identity[1] += rep/integer_matrix.length;

				sil=0;rep=0;
				// find SR for pos2 ********************************************

				for(int i=0;i<integer_matrix.length;i++){
					int tmp = getcodonnumber(integer_matrix[i][site],integer_matrix[i][site+1],integer_matrix[i][site+2]);
					if(Double.isNaN(possible[tmp][1])==false){// && Info1.Case==1){
						sil += possible[tmp][1];
						rep += 1.0-possible[tmp][1];
					}
				}
				identity[0] += sil/integer_matrix.length;
				identity[1] += rep/integer_matrix.length;

				sil=0;rep=0;
				// find SR for pos3 ********************************************

				for(int i=0;i<integer_matrix.length;i++){
					int tmp = getcodonnumber(integer_matrix[i][site],integer_matrix[i][site+1],integer_matrix[i][site+2]);
					if(Double.isNaN(possible[tmp][2])==false){// && Info2.Case==1){
						sil += possible[tmp][2];
						rep += 1.0-possible[tmp][2];
					}
				}
				identity[0] += sil/integer_matrix.length;
				identity[1] += rep/integer_matrix.length;

			}
		}
		return identity;
		// or add smoothing
	}

	public double[] NGoriginal(){
		// main NG method for all sites
		double[] identity = new double[2];
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				SiteInformation Info = SiteInformation(site);
				SiteInformation Info1 = SiteInformation(site+1);
				SiteInformation Info2 = SiteInformation(site+2);

				double sil=0;
				double rep=0;
				// find SR for pos1 ********************************************
				if(Info.polymorphismCase>1){
					for(int i=0;i<integer_matrix.length;i++){
						int[] tmp ={integer_matrix[i][site],integer_matrix[i][site+1],integer_matrix[i][site+2]};
						int[] tmpans ={integer_ancestral[site],integer_ancestral[site+1],integer_ancestral[site+2]};
						double[] sr = NGpathway(tmp,tmpans);
						if(sr[0]!=2){
							sil+=sr[0];
							rep+=1.0-sr[0];
						}
					}
					identity[0] += sil;
					identity[1] += rep;
				}
				// find SR for pos2 ********************************************
				sil=0;rep=0;
				if(Info1.polymorphismCase>1){
					for(int i=0;i<integer_matrix.length;i++){
						int[] tmp ={integer_matrix[i][site],integer_matrix[i][site+1],integer_matrix[i][site+2]};
						int[] tmpans ={integer_ancestral[site],integer_ancestral[site+1],integer_ancestral[site+2]};
						double[] sr = NGpathway(tmp,tmpans);
						if(sr[1]!=2){
							sil+=sr[1];
							rep+=1.0-sr[1];
						}
					}
					identity[0] += sil;
					identity[1] += rep;
				}
				// find SR for pos3 ********************************************
				sil=0;rep=0;
				if(Info2.polymorphismCase>1){
					for(int i=0;i<integer_matrix.length;i++){
						int[] tmp ={integer_matrix[i][site],integer_matrix[i][site+1],integer_matrix[i][site+2]};
						int[] tmpans ={integer_ancestral[site],integer_ancestral[site+1],integer_ancestral[site+2]};
						double[] sr = NGpathway(tmp,tmpans);
						if(sr[2]!=2){
							sil+=sr[2];
							rep+=1.0-sr[2];
						}
					}
					identity[0] += sil;
					identity[1] += rep;
				}
			}
		}

		return identity;
	}


	public double[][] NGmethod(int site,int codonsite,double[][] id){
		// HANDLES VARINT SITES ***************************************************
		// main NG method for all sites


		double[][] identity = new double[3][2];
		int[] ancestralbases = new int[3];
		ancestralbases[0]=integer_ancestral[site];
		ancestralbases[1]=integer_ancestral[site+1];
		ancestralbases[2]=integer_ancestral[site+2];
		int[] mainbases = new int[3];
		int[] count = new int[3];

		for(int i=0;i<integer_matrix.length;i++){
			// adds codon bases
			mainbases[0]=integer_matrix[i][site];
			mainbases[1]=integer_matrix[i][site+1];
			mainbases[2]=integer_matrix[i][site+2];
			double[] tmp = NGpathway(ancestralbases,mainbases);
			if(tmp[0]!=2.0){
				identity[0][0] += tmp[0];
				count[0]++;
			}
			if(tmp[1]!=2.0){
				identity[1][0] += tmp[1];
				count[1]++;
			}
			if(tmp[2]!=2.0){
				identity[2][0] += tmp[2];
				count[2]++;
			}
		}

		identity[0][0]=identity[0][0]/count[0];
		identity[1][0]=identity[1][0]/count[1];
		identity[2][0]=identity[2][0]/count[2];


		identity[0][1] = 1.0-identity[0][0];
		identity[1][1] = 1.0-identity[1][0];
		identity[2][1] = 1.0-identity[2][0];



		// HANDLES INVARIANT SITES ***************************************************
		// gets SR for invariant sites

		SiteInformation Info = SiteInformation(site);
		SiteInformation Info1 = SiteInformation(site+1);
		SiteInformation Info2 = SiteInformation(site+2);


		if(Info.polymorphismCase==1){	
			identity[0][0] = id[0][0];
			identity[0][1] = id[0][1];
		}
		// Position 2
		if(Info1.polymorphismCase==1){
			identity[1][0] = id[1][0];
			identity[1][1] = id[1][1];

		}
		// Position 3
		if(Info2.polymorphismCase==1){
			identity[2][0] = id[2][0];
			identity[2][1] = id[2][1];
		}

		return identity;
		// or add smoothing
	}

	public double[][] NGmethod(int site,int codonsite){
		// HANDLES VARINT SITES ***************************************************
		// main NG method for all sites


		double[][] identity = new double[3][2];
		int[] ancestralbases = new int[3];
		ancestralbases[0]=integer_ancestral[site];
		ancestralbases[1]=integer_ancestral[site+1];
		ancestralbases[2]=integer_ancestral[site+2];
		int[] mainbases = new int[3];
		int[] count = new int[3];

		for(int i=0;i<integer_matrix.length;i++){
			// adds codon bases
			mainbases[0]=integer_matrix[i][site];
			mainbases[1]=integer_matrix[i][site+1];
			mainbases[2]=integer_matrix[i][site+2];
			double[] tmp = NGpathway(ancestralbases,mainbases);
			if(tmp[0]!=2.0){
				identity[0][0] += tmp[0];
				count[0]++;
			}
			if(tmp[1]!=2.0){
				identity[1][0] += tmp[1];
				count[1]++;
			}
			if(tmp[2]!=2.0){
				identity[2][0] += tmp[2];
				count[2]++;
			}
		}

		identity[0][0]=identity[0][0]/count[0];
		identity[1][0]=identity[1][0]/count[1];
		identity[2][0]=identity[2][0]/count[2];


		identity[0][1] = 1.0-identity[0][0];
		identity[1][1] = 1.0-identity[1][0];
		identity[2][1] = 1.0-identity[2][0];
		return identity;
		// or add smoothing
	}

	public double[][] NGmethodNew(int site,int codonsite,double[][] CF){
		// HANDLES VARINT SITES ***************************************************
		// main NG method for all sites


		double[][] identity = new double[3][2];
		int[] ancestralbases = new int[3];
		ancestralbases[0]=integer_ancestral[site];
		ancestralbases[1]=integer_ancestral[site+1];
		ancestralbases[2]=integer_ancestral[site+2];
		int[] mainbases = new int[3];
		int[] count = new int[3];

		for(int i=0;i<integer_matrix.length;i++){
			// adds codon bases
			mainbases[0]=integer_matrix[i][site];
			mainbases[1]=integer_matrix[i][site+1];
			mainbases[2]=integer_matrix[i][site+2];
			double[] tmp = NGpathway(ancestralbases,mainbases);
			if(tmp[0]!=2.0){
				identity[0][0] += tmp[0];
				count[0]++;
			}
			if(tmp[1]!=2.0){
				identity[1][0] += tmp[1];
				count[1]++;
			}
			if(tmp[2]!=2.0){
				identity[2][0] += tmp[2];
				count[2]++;
			}
		}

		identity[0][0]=identity[0][0]/count[0];
		identity[1][0]=identity[1][0]/count[1];
		identity[2][0]=identity[2][0]/count[2];


		identity[0][1] = 1.0-identity[0][0];
		identity[1][1] = 1.0-identity[1][0];
		identity[2][1] = 1.0-identity[2][0];



		// HANDLES INVARIANT SITES ***************************************************
		// gets SR for invariant sites

		SiteInformation Info = SiteInformation(site);
		SiteInformation Info1 = SiteInformation(site+1);
		SiteInformation Info2 = SiteInformation(site+2);
		count = new int[3];
		double[][] id = new double[3][2];
		if(Info.polymorphismCase==1){
			for(int i=0;i<integer_matrix.length;i++){
				int tmp = getcodonnumber(integer_matrix[i][site],integer_matrix[i][site+1],integer_matrix[i][site+2]);

				id[0][0] += CF[tmp][0];
				count[0]++;
			}
			id[0][0]=id[0][0]/count[0];
			identity[0][0]=id[0][0];identity[0][1]=1.0-id[0][0];
		}
		// Position 2
		if(Info1.polymorphismCase==1){
			for(int i=0;i<integer_matrix.length;i++){
				int tmp = getcodonnumber(integer_matrix[i][site],integer_matrix[i][site+1],integer_matrix[i][site+2]);
				id[1][0] += CF[tmp][1];
				count[1]++;
			}
			id[1][0]=id[1][0]/count[1];
			identity[1][0]=id[1][0];identity[1][1]=1.0-id[1][0];
		}
		// Position 3
		if(Info2.polymorphismCase==1){
			for(int i=0;i<integer_matrix.length;i++){
				int tmp = getcodonnumber(integer_matrix[i][site],integer_matrix[i][site+1],integer_matrix[i][site+2]);
				id[2][0] += CF[tmp][2];
				count[2]++;
			}
			id[2][0]=id[2][0]/count[2];
			identity[2][0]=id[2][0];identity[2][1]=1.0-id[2][0];
		}



		return identity;
		// or add smoothing
	}

	public double[][] NGinvariant(){
		int[] mainbases = new int[3];
		double L = (double) integer_matrix.length;
		double N = (double) integer_matrix[0].length;
		double[][] identity = new double[3][2];
		double[][] Avgidentity = new double[3][2];

		// Most sequences do not have even base frequencies. 
		// this bit of code scales for this case by calculated ATGC frequencies for the whole alignment

		double[] cf = EmpiricalCodonFreq();

		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			identity = new double[3][2]; // set new identity vector
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) {  // check for bad sites
				for(int sequence=0;sequence<integer_matrix.length;sequence++){		// loop through sequences

					double cf1Total = 0; //reset base frequency total for each codon
					double cf2Total = 0;
					double cf3Total = 0;
					mainbases[0]=integer_matrix[sequence][site]; //main bases
					mainbases[1]=integer_matrix[sequence][site+1];
					mainbases[2]=integer_matrix[sequence][site+2];
					double sil = 0; //silentProb count
					double rep=0; //rep count
					//********************************************************************************************************************************************
					//pos 1
					sil=0;
					// Get total of probabilities of possible bases
					// i.e if an invariant sites exists at pos 1, and is A, then sum the base frequencies of CGT (a form of rescaling)
					for(int i=1;i<5;i++){
						int codonNumber =  getcodonnumber(i,mainbases[1],mainbases[2]);
						if(mainbases[0]!=i){ // && codonNumber!=48 && codonNumber!=50 && codonNumber!=56){
							cf1Total += cf[codonNumber];
						}			
					}

					for(int i=1;i<5;i++){  // loop through all bases
						if(i != mainbases[0]){ // make sure actual base is not counted
							//	int actualCodonNumber =  getcodonnumber(mainbases[0],mainbases[1],mainbases[2]); //ancestral codon
							int intermediateCodonNumber =  getcodonnumber(i,mainbases[1],mainbases[2]); //psedo codon
							//	sil += SilentOrReplacement(actualCodonNumber,intermediateCodonNumber)*(cf[intermediateCodonNumber]/cf1Total); //check silentProb change or replacement change
							int[] actualCodon = {mainbases[0],mainbases[1],mainbases[2]};
							int[] IntermediateCodon = {i,mainbases[1],mainbases[2]};
							double[] SR = NGpathway(actualCodon,IntermediateCodon);
							sil += SR[0]*(cf[intermediateCodonNumber]/cf1Total);							
						}
					}
					rep = 1.0-sil;
					identity[0][0]+=sil;identity[0][1]+=rep;
					//********************************************************************************************************************************************
					//pos 2
					sil=0;
					for(int i=1;i<5;i++){
						int codonNumber =  getcodonnumber(mainbases[0],i,mainbases[2]);
						if(mainbases[1]!=i){ // && codonNumber!=48 && codonNumber!=50 && codonNumber!=56){
							cf2Total += cf[codonNumber];
						}			
					}


					for(int i=1;i<5;i++){  // loop through all bases
						if(i != mainbases[1]){ // make sure actual base is not counted
							//	int actualCodonNumber =  getcodonnumber(mainbases[0],mainbases[1],mainbases[2]); //ancestral codon
							int intermediateCodonNumber =  getcodonnumber(mainbases[0],i,mainbases[2]); //psedo codon
							//	sil += SilentOrReplacement(actualCodonNumber,intermediateCodonNumber)*(cf[intermediateCodonNumber]/cf2Total);  //check silentProb change or replacement change
							int[] actualCodon = {mainbases[0],mainbases[1],mainbases[2]};
							int[] IntermediateCodon = {mainbases[0],i,mainbases[2]};
							double[] SR = NGpathway(actualCodon,IntermediateCodon);
							sil += SR[1]*(cf[intermediateCodonNumber]/cf2Total);	
						}
					}
					rep = 1.0-sil;
					identity[1][0]+=sil;identity[1][1]+=rep;
					//********************************************************************************************************************************************
					//pos 3
					sil=0;
					for(int i=1;i<5;i++){
						int codonNumber =  getcodonnumber(mainbases[0],mainbases[1],i);
						if(mainbases[2]!=i){ // && codonNumber!=48 && codonNumber!=50 && codonNumber!=56){
							cf3Total += cf[codonNumber];
						}			
					}


					for(int i=1;i<5;i++){  // loop through all bases
						if(i != mainbases[2]){ // make sure actual base is not counted
							//	int actualCodonNumber =  getcodonnumber(mainbases[0],mainbases[1],mainbases[2]); //ancestral codon
							int intermediateCodonNumber =  getcodonnumber(mainbases[0],mainbases[1],i); //psedo codon
							//	sil += SilentOrReplacement(actualCodonNumber,intermediateCodonNumber)*(cf[intermediateCodonNumber]/cf3Total);  //check silentProb change or replacement change
							int[] actualCodon = {mainbases[0],mainbases[1],mainbases[2]};
							int[] IntermediateCodon = {mainbases[0],mainbases[1],i};
							double[] SR = NGpathway(actualCodon,IntermediateCodon);
							sil += SR[2]*(cf[intermediateCodonNumber]/cf3Total);	
						}
					}
					rep = 1.0-sil;
					identity[2][0]+=sil;identity[2][1]+=rep;
					//********************************************************************************************************************************************
				}
			}
			Avgidentity[0][0]+=identity[0][0]/L;
			Avgidentity[1][0]+=identity[1][0]/L;
			Avgidentity[2][0]+=identity[2][0]/L;
			Avgidentity[0][1]+=identity[0][1]/L;
			Avgidentity[1][1]+=identity[1][1]/L;
			Avgidentity[2][1]+=identity[2][1]/L;
		}
		Avgidentity[0][0]=Avgidentity[0][0]/N;
		Avgidentity[1][0]=Avgidentity[1][0]/N;
		Avgidentity[2][0]=Avgidentity[2][0]/N;
		Avgidentity[0][1]=Avgidentity[0][1]/N;
		Avgidentity[1][1]=Avgidentity[1][1]/N;
		Avgidentity[2][1]=Avgidentity[2][1]/N;
		return Avgidentity;

	}

	// calculates the possible changes at each site
	public double[] NGpossibleChanges(int[] seqcodon,double[] bf1,double[] bf2, double[] bf3){
		// A,C,G,T = 1,2,3,4
		//remember to change
		int actualCodonNumber =  getcodonnumber(seqcodon[0],seqcodon[1],seqcodon[2]);
		double[] possibleSR = new double[3];
		double bf1Total = 0;
		double bf2Total = 0;
		double bf3Total = 0;
		// Get total of probabilities of possible bases
		// i.e if an invariant sites exists at pos 1, and is A, then sum the base frequencies of CGT (a form of rescaling)
		for(int i=1;i<5;i++){
			if(seqcodon[0]!=i){
				bf1Total += bf1[i-1];
			}

		}
		for(int i=1;i<5;i++){
			if(seqcodon[1]!=i){
				bf2Total += bf2[i-1];
			}

		}
		for(int i=1;i<5;i++){
			if(seqcodon[2]!=i){
				bf3Total += bf3[i-1];
			}

		}
		// Position 1
		for(int i=1;i<5;i++){
			int codonNumber =  getcodonnumber(i,seqcodon[1],seqcodon[2]);
			if(seqcodon[0]!=i){
				possibleSR[0] +=  SilentOrReplacement(actualCodonNumber,codonNumber)*(bf1[i-1]/bf1Total);
			}
		}
		// Position 2
		for(int i=1;i<5;i++){
			int codonNumber =  getcodonnumber(seqcodon[0],i,seqcodon[2]);
			if(seqcodon[1]!=i){
				possibleSR[1] +=  SilentOrReplacement(actualCodonNumber,codonNumber)*(bf2[i-1]/bf2Total);
			}
		}
		// Position 3
		for(int i=1;i<5;i++){
			int codonNumber =  getcodonnumber(seqcodon[0],seqcodon[1],i);
			if(seqcodon[2]!=i){
				possibleSR[2] +=  SilentOrReplacement(actualCodonNumber,codonNumber)*(bf3[i-1]/bf3Total);
			}
		}
		return possibleSR;
	}


	public double[] NGpathway(int[] a, int[] b){ //a = ancestral b=main
		// 1 silentProb, 0 replacement, 2 invariant
		int degen = 0;

		int AnscodonNumber =  getcodonnumber(a[0],a[1],a[2]);
		int codonNumber =  getcodonnumber(b[0],b[1],b[2]);
		boolean[] diff = whichDiff(a,b);


		double[] identity = new double[3];
		identity[0]=2;identity[1]=2;identity[2]=2;
		for(int i=0;i<3;i++){if(diff[i]){degen++;}} // find how many sites are degererate
		// 1 site different
		if(degen==1){
			if(diff[0]){identity[0]=SilentOrReplacement(AnscodonNumber,codonNumber);}
			if(diff[1]){identity[1]=SilentOrReplacement(AnscodonNumber,codonNumber);}
			if(diff[2]){identity[2]=SilentOrReplacement(AnscodonNumber,codonNumber);}
		}
		// 2 sites different
		if(degen==2){ 
			//			Pathway 1  (a) X11-X21-X22 and  (b) X11-X12-X22	  *Both pathways occur with equal probability*
			if(diff[0]==false){
				// (a) X11-X21-X22
				int a_codonNumber =  getcodonnumber(a[0],b[1],a[2]);
				identity[1]=SilentOrReplacement(AnscodonNumber,a_codonNumber);
				identity[2]=SilentOrReplacement(a_codonNumber,codonNumber);
				// (b) X11-X12-X22
				int b_codonNumber =  getcodonnumber(a[0],a[1],b[2]);
				identity[2]+=SilentOrReplacement(AnscodonNumber,b_codonNumber);
				identity[1]+=SilentOrReplacement(b_codonNumber,codonNumber);
				identity[1]=identity[1]/2;identity[2]=identity[2]/2; //average over all pathways
			}
			//			Pathway 2  (a) 1X1-2X1-2X2 and (b) 1X1-1X2-2X2		
			if(diff[1]==false){
				// (a) 1X1-2X1-2X2
				int a_codonNumber =  getcodonnumber(b[0],a[1],a[2]);
				identity[0]=SilentOrReplacement(AnscodonNumber,a_codonNumber);
				identity[2]=SilentOrReplacement(a_codonNumber,codonNumber);

				// (b) 1X1-1X2-2X2	
				int b_codonNumber =  getcodonnumber(a[0],a[1],b[2]);
				identity[2]+=SilentOrReplacement(AnscodonNumber,b_codonNumber);
				identity[0]+=SilentOrReplacement(b_codonNumber,codonNumber);
				identity[0]=identity[0]/2;identity[2]=identity[2]/2; //average over all pathways
			}
			//			Pathway 3   (a) 11X-21X-22X and  (b) 11X-12X-22X	
			if(diff[2]==false){
				// (a) 11X-21X-22X 
				int a_codonNumber =  getcodonnumber(b[0],a[1],a[2]);
				identity[0]=SilentOrReplacement(AnscodonNumber,a_codonNumber);
				identity[1]=SilentOrReplacement(a_codonNumber,codonNumber);
				// (b) 11X-12X-22X
				int b_codonNumber =  getcodonnumber(a[0],b[1],a[2]);
				identity[1]+=SilentOrReplacement(AnscodonNumber,b_codonNumber);
				identity[0]+=SilentOrReplacement(b_codonNumber,codonNumber);
				identity[0]=identity[0]/2;identity[1]=identity[1]/2; //average over all pathways
			}
		}

		// 3 sites different
		if(degen==3){

			// Pathway 1 111-211-221-222	*x -(a)-(b)- x*        
			int a_codonNumber = getcodonnumber(b[0],a[1],a[2]);
			int b_codonNumber = getcodonnumber(b[0],b[1],a[2]);
			identity[0]=SilentOrReplacement(AnscodonNumber,a_codonNumber);
			identity[1]=SilentOrReplacement(a_codonNumber,b_codonNumber);
			identity[2]=SilentOrReplacement(b_codonNumber,codonNumber);
			// Pathway 2 111-211-212-222       
			a_codonNumber = getcodonnumber(b[0],a[1],a[2]);
			b_codonNumber = getcodonnumber(b[0],a[1],b[2]);
			identity[0]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
			identity[2]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
			identity[1]+=SilentOrReplacement(b_codonNumber,codonNumber);				
			// Pathway 3 111-121-221-222    
			a_codonNumber = getcodonnumber(a[0],b[1],a[2]);
			b_codonNumber = getcodonnumber(b[0],b[1],a[2]);
			identity[1]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
			identity[0]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
			identity[2]+=SilentOrReplacement(b_codonNumber,codonNumber);		
			// Pathway 4  111-121-122-222     
			a_codonNumber = getcodonnumber(a[0],b[1],a[2]);
			b_codonNumber = getcodonnumber(a[0],b[1],b[2]);
			identity[1]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
			identity[2]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
			identity[0]+=SilentOrReplacement(b_codonNumber,codonNumber);		
			// Pathway 5 111-112-212-222      
			a_codonNumber = getcodonnumber(a[0],a[1],b[2]);
			b_codonNumber = getcodonnumber(b[0],a[1],b[2]);
			identity[2]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
			identity[0]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
			identity[1]+=SilentOrReplacement(b_codonNumber,codonNumber);
			// Pathway 6 111-112-122-222       
			a_codonNumber = getcodonnumber(a[0],a[1],b[2]);
			b_codonNumber = getcodonnumber(a[0],b[1],b[2]);
			identity[2]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
			identity[1]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
			identity[0]+=SilentOrReplacement(b_codonNumber,codonNumber);

			identity[0]=identity[0]/6.0;identity[1]=identity[1]/6.0;identity[2]=identity[2]/6.0;
		}


		return identity;
	}

	// evaluats wheter a codon is silentProb (1) or replacement (0)
	public double SilentOrReplacement(int ansnumber, int number){
		double identity = 0.0;
		if(AA[ansnumber].equals(AA[number]) ){identity=1.0;} else {identity=0.0;}
		return identity;
	}
	// finds which bases are different from ancestral codon
	// identifies which bases are different between anscestral codon and main codon - for use in nei gojobori pathways
	public boolean[] whichDiff(int[] anscodon,int[] seqcodon){
		boolean[] flag = new boolean[3];
		flag[0]=true;flag[1]=true;flag[2]=true;
		if(anscodon[0]==seqcodon[0]){
			flag[0]=false;
		}
		if(anscodon[1]==seqcodon[1]){
			flag[1]=false;
		}
		if(anscodon[2]==seqcodon[2]){
			flag[2]=false;
		}
		return flag;
	}
	//	*********************************************************************************

	public double[][] srSmoothing(int site,int codon){
		double effectiveS = 0.0;
		double[][] identity = new double[3][2];
		int[] ancestralbases = preprocess.codonSplitter(codon_ancestral[codon]);
		int codonNumber =  getcodonnumber(ancestralbases[0],ancestralbases[1],ancestralbases[2]);
		// equal codon probabilities
		double pr1 = 0.0;double pr2 = 0.0;double pr3 = 0.0;	
		pr1=0.071;pr2=0.0;pr3=0.653;    //F3X4 model probabilities 

		// pos 1 ********************************************************************
		int silcount = 0;
		int repcount = 0;
		SiteInformation Info = SiteInformation(site);
		if(Info.polymorphismCase==1){				
			identity[0][0]=pr1;
			identity[0][1]=1.0-pr1;		
		} else {
			int num = 0;
			for(int sample=0;sample<integer_matrix.length;sample++){	// loop through all samples				
				if(integer_matrix[sample][site]!=integer_ancestral[site]){
					num++;
					int derivedcodonNumber = 0;
					double pos1_identity = 0.0;	//position 1

					derivedcodonNumber = getcodonnumber(integer_matrix[sample][site],ancestralbases[1],ancestralbases[2]);
					if(AA[codonNumber].equals(AA[derivedcodonNumber]) ){
						pos1_identity=1.0;
					} else {
						pos1_identity=0.0;
					}

					if(pos1_identity==1.0){
						silcount++;
					}
					else if(pos1_identity==0.0){
						repcount++;
					} 
				}
			}
			//smoothing		
			double[] observations = new double[2];
			observations[0]=silcount;
			observations[1]=repcount;
			// add priors
			double pos1mean = pr1;
			double[] priors = BetaParameters(pos1mean,effectiveS);
			observations[0]+=priors[0];
			observations[1]+=priors[1];
			double silentcount = 0;
			double replacementcount = 0;
			silentcount = observations[0]/(observations[0]+observations[1]);
			replacementcount = 1.0-silentcount;

			identity[0][0] = silentcount;
			identity[0][1] = replacementcount;
		}
		// pos 2 *******************************************************************
		silcount = 0;
		repcount = 0;
		Info = SiteInformation(site+1);
		if(Info.polymorphismCase==1){
			identity[1][0]=pr2;
			identity[1][1]=1.0-pr2;			
		} else {
			int num = 0;
			for(int sample=0;sample<integer_matrix.length;sample++){	// loop through all samples
				if(integer_matrix[sample][site+1]!=integer_ancestral[site+1]){	
					num++;
					int derivedcodonNumber = 0;
					double pos2_identity = 0.0;	//position 1
					derivedcodonNumber = getcodonnumber(ancestralbases[0],integer_matrix[sample][site+1],ancestralbases[2]);
					if(AA[codonNumber].equals(AA[derivedcodonNumber]) ){
						pos2_identity=1.0;
					} else {
						pos2_identity=0.0;
					}

					if(pos2_identity==1.0){
						silcount++;
					}	
					if(pos2_identity==0.0){
						repcount++;
					}
				}
			}
			//smoothing
			double[] observations = new double[2];
			observations = new double[2];
			observations[0]=silcount;
			observations[1]=repcount;
			// add priors
			double pos2mean = pr2;
			double[] priors = BetaParameters(pos2mean,effectiveS);
			observations[0]+=priors[0];
			observations[1]+=priors[1];
			double silentcount = 0;
			double replacementcount = 0;
			silentcount = observations[0]/(observations[0]+observations[1]);
			replacementcount = 1.0-silentcount;
			identity[1][0] = silentcount;
			identity[1][1] = replacementcount; 

		}

		// pos 3  *******************************************************************
		silcount = 0;
		repcount = 0;
		Info = SiteInformation(site+2);
		if(Info.polymorphismCase==1){
			identity[2][0]=pr3;
			identity[2][1]=1.0-pr3;
		} else {
			int num = 0;
			for(int sample=0;sample<integer_matrix.length;sample++){	// loop through all samples
				if(integer_matrix[sample][site+2]!=integer_ancestral[site+2]){	
					num++;
					int derivedcodonNumber = 0;
					double pos3_identity = 0.0;	//position 1
					derivedcodonNumber = getcodonnumber(ancestralbases[0],ancestralbases[1],integer_matrix[sample][site+2]);
					if(AA[codonNumber].equals(AA[derivedcodonNumber]) ){
						pos3_identity=1.0;
					} else {
						pos3_identity=0.0;
					}

					if(pos3_identity==1.0){
						silcount++;
					}
					if(pos3_identity==0.0){
						repcount++;
					}
				}
			}
			//smoothing
			double[] observations = new double[2];
			observations = new double[2];
			observations[0]=silcount;
			observations[1]=repcount;
			// add priors
			double pos3mean = pr3;
			double[] priors = BetaParameters(pos3mean,effectiveS);
			observations[0]+=priors[0];
			observations[1]+=priors[1];
			double silentcount = 0;
			double replacementcount = 0;
			silentcount = observations[0]/(observations[0]+observations[1]);
			replacementcount = 1.0-silentcount;
			identity[2][0] = silentcount;
			identity[2][1] = replacementcount; 

		}
		return identity;
	}

	public double[] SiteFreq(double u,double v, double[] prior, boolean needPrior){
		double[][] temp = new double[3][2];
		double rho = 0.0;
		double sigma = 0.0;
		double[] finalans = new double[3];
		double total = 0.0;	

		//		double[][] Invidentity = NGinvariant();
		//	double[][] Invidentity = SRcodonInvariant();
		double[][] Invidentity= InvariantNew();
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				//		temp = NGmethodNew(site,codon,Invidentity);
				temp = NGmethod(site,codon,Invidentity);
				SiteInformation inf = DirichletSRFreq(u, v, site, prior, needPrior);//DirichletSiteFreq(u,v,site,prior,needPrior);		// find site freq for pos1
				//	if(inf.Case>1){
				if(inf.polymorphismCase>1){	
					total += inf.dirichletProb;
					sigma += inf.dirichletProb*temp[0][0];
					rho += inf.dirichletProb*temp[0][1];
				}
				SiteInformation inf2 = DirichletSRFreq(u, v, site + 1, prior, needPrior); //DirichletSiteFreq(u,v,site+1,prior,needPrior);		// find site freq for pos2
				//		if(inf2.Case>1){
				if(inf2.polymorphismCase>1){		
					total += inf2.dirichletProb;
					sigma += inf2.dirichletProb*temp[1][0];
					rho += inf2.dirichletProb*temp[1][1];
				}
				SiteInformation inf3 = DirichletSRFreq(u, v, site + 2, prior, needPrior); // DirichletSiteFreq(u,v,site+2,prior,needPrior);		// find site freq for pos3
				//	if(inf3.Case>1){
				if(inf3.polymorphismCase>1){		
					total += inf3.dirichletProb;
					sigma += inf3.dirichletProb*temp[2][0];
					rho +=  inf3.dirichletProb*temp[2][1];
				}
			} 		
		}
		
		
		
		finalans[0] = sigma;
		finalans[1] = rho;
		finalans[2] = total;
		return finalans;		

	}
	
	
	public double[] SiteFreqNew(double u,double v, double[] prior, boolean needPrior){
		double[][] temp = new double[3][2];
		double rho = 0.0;
		double sigma = 0.0;
		double[] finalans = new double[3];
		double total = 0.0;	

		//		double[][] Invidentity = NGinvariant();
		//	double[][] Invidentity = SRcodonInvariant();
		double[][] Invidentity= InvariantNew();
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				//		temp = NGmethodNew(site,codon,Invidentity);
				temp = NGmethod(site,codon,Invidentity);
				SiteInformation inf = DirichletSRFreq(u, v, site, prior, needPrior); //DirichletSiteFreq(u,v,site,prior,needPrior);		// find site freq for pos1
				//	if(inf.Case>1){
				if(inf.polymorphismCase>1){	
					total += inf.dirichletProb;
					sigma += inf.dirichletProb*temp[0][0];
					rho += inf.dirichletProb*temp[0][1];
				}
				SiteInformation inf2 = DirichletSRFreq(u, v, site + 1, prior, needPrior); // DirichletSiteFreq(u,v,site+1,prior,needPrior);		// find site freq for pos2
				//		if(inf2.Case>1){
				if(inf2.polymorphismCase>1){		
					total += inf2.dirichletProb;
					sigma += inf2.dirichletProb*temp[1][0];
					rho += inf2.dirichletProb*temp[1][1];
				}
				SiteInformation inf3 = DirichletSRFreq(u, v, site + 2, prior, needPrior); //DirichletSiteFreq(u,v,site+2,prior,needPrior);		// find site freq for pos3
				//	if(inf3.Case>1){
				if(inf3.polymorphismCase>1){		
					total += inf3.dirichletProb;
					sigma += inf3.dirichletProb*temp[2][0];
					rho +=  inf3.dirichletProb*temp[2][1];
				}
			} 		
		}
		
		
		
		finalans[0] = sigma;
		finalans[1] = rho;
		finalans[2] = total;
		return finalans;		

	}
	//	******************	

	public double[][] intervals(){
		double n = number_bins.doubleValue();
		double bin_size = 1/n;
		double[][] intervals = new double[2][number_bins.intValue()]; // stores upper and lower bounds
		intervals[0][0] = 0;
		intervals[1][0] = (bin_size);
		for(int i=1; i< n; i++){
			intervals[0][i] = (intervals[1][i-1]); 
			intervals[1][i] = (intervals[1][i-1]+bin_size);
		}
		return intervals;	
	}

	public double totalNoAdapt(boolean[] NeutralVec , double[][] binsvalues, double NeutralRatio,double[] prior, boolean needPrior){
		SetBins(NeutralVec);
		SetBinsValues(binsvalues);
		double totaladapt=0;
		double[][] finalmat = new double[3][bins[0].length];
		for(int i=0;i<finalmat[0].length;i++){
			double[] temp = SiteFreq(bins[0][i], bins[1][i],prior,needPrior);
			finalmat[0][i] = temp[0];   // number silentProb
			finalmat[1][i] = temp[1];	// number replacement	
			finalmat[2][i] = temp[2];	// total number
		}
		NeutralRatio=0;
		double nsil=0;
		double nrep=0;
		for(int i=0;i<finalmat[0].length;i++){
			if(WhichBins[i]){
				nsil += finalmat[0][i];
				nrep += finalmat[1][i];
			}
		}
		if(nsil!=0){
			NeutralRatio = nrep/nsil;
		}

		for(int i=0;i<finalmat[0].length;i++){
			if(WhichBins[i]==false){
				totaladapt += finalmat[1][i]-finalmat[0][i]*NeutralRatio;
			}
		}
		if(totaladapt<0){
			totaladapt=0;
		}

		return totaladapt;
	}

	public double[][] ValueMat(boolean[] NeutralVec , double[][] binsvalues, double NeutralRatio,double[] prior, boolean needPrior){
		SetBins(NeutralVec);
		SetBinsValues(binsvalues);
		SetNumBins(binsvalues[0].length);
		double totaladapt=0;
		double[][] finalmat = new double[3][bins[0].length];
		for(int i=0;i<finalmat[0].length;i++){
			double[] temp = SiteFreq(bins[0][i], bins[1][i],prior,needPrior);
			finalmat[0][i] = temp[0];   // number silentProb
			finalmat[1][i] = temp[1];	// number replacement	
			finalmat[2][i] = temp[2];	// total number
		}
		for(int i=0;i<finalmat[0].length;i++){
			if(WhichBins[i]==false){
				totaladapt += finalmat[1][i]-finalmat[0][i]*NeutralRatio;
			}
		}


		return finalmat;

	}


	//	evenly spaced binsMatrix
	public double[][] Value_Matrix(boolean[] NeutralVec, double[][] binsvalues,double[] prior, boolean needPrior){
		SetBins(NeutralVec);
		SetBinsValues(binsvalues);
		SetNumBins(binsvalues[0].length);
		double[][] finalmat = new double[6][binsvalues[0].length]; 

		for(int i=0;i<number_bins.intValue();i++){
			double[] temp = SiteFreq(binsvalues[0][i], binsvalues[1][i],prior,needPrior);
			finalmat[0][i] = temp[0];   // number silentProb
			finalmat[1][i] = temp[1];	// number replacement	
			finalmat[2][i] = temp[2];	// total number

			if(temp[1]!=0){
				finalmat[3][i] = temp[1]/temp[0];	// Silent/replacement ratio
			} else {
				finalmat[3][i] = Double.NaN;
			}
		}
		// calcualte neutral ratio. Vector NeutralVec provides boolean true false if a given bin is neutral
		// user has to specify this bin
		double counter = 0;
		double neutralratio=0;
		for(int i=0;i<number_bins.intValue();i++){
			if(WhichBins[i]){
				if(finalmat[0][i]!=0){
					neutralratio += finalmat[1][i]/finalmat[0][i];
					counter	++;
				}
			}
		}
		neutralratio = neutralratio/counter; //average
		//	neutralRatio=1.75533276;

		// update neutral ratios for the output matrix
		for(int t=0;t<number_bins.intValue();t++){
			finalmat[4][t] = neutralratio;
		}

		// number non neutral
		for(int i=0;i<number_bins.intValue();i++){
			finalmat[5][i] = finalmat[1][i]*(1.0 - ((finalmat[0][i]/finalmat[1][i])*(neutralratio)));	// number non neutral
			if(finalmat[5][i]<0){
				finalmat[5][i]=0.0;
			}
		}

		return finalmat;

	}




	public double[][] find_identityMK(int site, int codon){
		double[][] identity = new double[3][2];
		int[] ancestralbases = preprocess.codonSplitter(codon_ancestral[codon]);
		int codonNumber =  getcodonnumber(ancestralbases[0],ancestralbases[1],ancestralbases[2]);
		int silcount = 0;
		// pos 1
		for(int i=0;i<integer_matrix.length;i++){	// loop through all samples
			int derivedcodonNumber = 0;
			double pos1_identity = 0.0;	//position 1
			if (codonNumber==8 || codonNumber==10 || codonNumber==24 ||
					codonNumber==25 || codonNumber==26 || codonNumber==27 ||
					codonNumber==28 || codonNumber==29 || codonNumber==30 ||
					codonNumber==31 || codonNumber==60 || codonNumber==62) {
				pos1_identity = 0.5;
			}
			if(pos1_identity==0.5){
				derivedcodonNumber = getcodonnumber(integer_matrix[i][site],ancestralbases[1],ancestralbases[2]);
				if(AA[codonNumber].equals(AA[derivedcodonNumber]) ){
					pos1_identity=1.0;
				} else {
					pos1_identity=0.0;
				}
			}
			if(pos1_identity==1.0){
				silcount++;
			}
		}
		identity[0][0] = silcount/integer_matrix.length;
		identity[0][1] = 1.0 - identity[0][0]; 
		// pos 2
		silcount = 0;
		for(int i=0;i<integer_matrix.length;i++){	// loop through all samples
			int derivedcodonNumber = 0;
			double pos2_identity = 0.0;	//position 1


			if (codonNumber==48 || codonNumber==50 || codonNumber==56) {
				pos2_identity = 0.5;
			}
			if(pos2_identity==0.5){
				derivedcodonNumber = getcodonnumber(ancestralbases[0],integer_matrix[i][site+1],ancestralbases[2]);
				if(AA[codonNumber].equals(AA[derivedcodonNumber]) ){
					pos2_identity=1.0;
				} else {
					pos2_identity=0.0;
				}
			}
			if(pos2_identity==1.0){
				silcount++;
			}	
		}
		identity[1][0] = silcount/integer_matrix.length;
		identity[1][1] = 1.0 - identity[1][0];
		// pos 3	
		silcount = 0;
		for(int i=0;i<integer_matrix.length;i++){	// loop through all samples
			int derivedcodonNumber = 0;
			double pos3_identity = 1.0;	//position 1

			if (codonNumber==14 || codonNumber==58) {
				pos3_identity= 0.0;
			} else if (codonNumber==0 || codonNumber==1 || codonNumber==2 || codonNumber==3 ||
					codonNumber==9 || codonNumber==11 || codonNumber==12 || codonNumber==13 ||
					codonNumber==15 || codonNumber==16 || codonNumber==17 || codonNumber==18 ||
					codonNumber==19 || codonNumber==32 || codonNumber==33 || codonNumber==34 ||
					codonNumber==35 || codonNumber==48 || codonNumber==49 || codonNumber==50 ||
					codonNumber==51 || codonNumber==56 || codonNumber==57 || codonNumber==59 ||
					codonNumber==61 || codonNumber==63) {
				pos3_identity = 0.5;
			}
			if(pos3_identity==0.5){
				derivedcodonNumber = getcodonnumber(ancestralbases[0],ancestralbases[1],integer_matrix[i][site+2]);
				if(AA[codonNumber].equals(AA[derivedcodonNumber]) ){
					pos3_identity=1.0;
				} else {
					pos3_identity=0.0;
				}
			}
			if(pos3_identity==1.0){
				silcount++;
			}
		}
		identity[2][0] = silcount/integer_matrix.length;
		identity[2][1] = 1.0 - identity[2][0];
		return identity;

	}

	public int getcodonnumber(int pos1, int pos2, int pos3){
		if (pos1<5 && pos2<5 && pos3<5) { // normal codon
			return (16*(pos1-1)) + (4*(pos2-1)) + (1*(pos3-1));
		}
		else if (pos1==5 && pos2==5 && pos3==5) { // codon is NNN
			return 64;
		}
		else if (pos1==5 && pos2==5 && pos3==5) { // codon is ---
			return 65;
		}
		return 66; // codon is unrecognised
	}

	public double Distance(){
		double p=0;
		double d=0;
		double count1=0;
		double count2=0;
		for(int i=0;i<integer_matrix.length;i++){
			count1=0;
			count2=0;
			p=0;
			for(int j=0;j<integer_matrix[0].length;j++){
				if(bad_sites_list[j]==false){
					if(integer_matrix[i][j]!=integer_ancestral[j]){count1++;} //counting differences
					count2++;
				}
			}
			p = count1/count2;

			//	d += count1/count2;  //normal uncorrected PWD
			d += (-3.0/4.0)*Math.log(1.0-((4.0/3.0)*p));	//  jukes cantor correction
		}
		return (d/integer_matrix.length);
	}

	public double K2P(){
		double P=0;
		double Q=0;
		double d=0;
		double count2=0;
		for(int i=0;i<integer_matrix.length;i++){

			count2=0;
			P=0;Q=0;
			for(int j=0;j<integer_matrix[0].length;j++){
				if(bad_sites_list[j]==false){
					if(integer_matrix[i][j]==integer_ancestral[j]){/*do nothing*/} 
					else if(integer_matrix[i][j]==1 && integer_ancestral[j] == 3 || integer_matrix[i][j]==3 && integer_ancestral[j] == 1 || integer_matrix[i][j]==2 && integer_ancestral[j] == 4 || integer_matrix[i][j]==4 && integer_ancestral[j] == 2){P++;} 
					else {Q++;} 
					count2++;
				}
			}
			P=P/count2;
			Q=Q/count2;
			// jukes cantor correction
			d += ((-1.0/2.0)*Math.log(1.0 - (2.0*P) - Q)) - ((-1/4)*Math.log(1.0-(2.0*Q)));	
		}
		return (d/integer_matrix.length);
	}


	public double howManyMultiPoly(){
		double count=0;
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				SiteInformation s = SiteInformation(site);			
				if(s.polymorphismCase==4 ||s.polymorphismCase==5 ||s.polymorphismCase==6 ||s.polymorphismCase==7 ){
					count++;
				}
				s = SiteInformation(site+1);			
				if(s.polymorphismCase==4 ||s.polymorphismCase==5 ||s.polymorphismCase==6 ||s.polymorphismCase==7 ){
					count++;
				}
				s = SiteInformation(site+2);			
				if(s.polymorphismCase==4 ||s.polymorphismCase==5 ||s.polymorphismCase==6 ||s.polymorphismCase==7 ){
					count++;
				}
			}
		}
		return count;
	}

	public double howManySingPoly(){
		double count=0;
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				SiteInformation s = SiteInformation(site);			
				if(s.polymorphismCase==3){
					count++;
				}
				s = SiteInformation(site+1);			
				if(s.polymorphismCase==3){
					count++;
				}
				s = SiteInformation(site+2);			
				if(s.polymorphismCase==3){
					count++;
				}
			}
		}
		return count;
	}

	public double howManyFixed(){
		double count=0;
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				SiteInformation s = SiteInformation(site);			
				if(s.polymorphismCase==2){
					count++;
				}
				s = SiteInformation(site+1);			
				if(s.polymorphismCase==2){
					count++;
				}
				s = SiteInformation(site+2);			
				if(s.polymorphismCase==2){
					count++;
				}
			}
		}
		return count;
	}

	public double howManyInvariant(){
		double count=0;
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				SiteInformation s = SiteInformation(site);			
				if(s.polymorphismCase==1){
					count++;
				}
				s = SiteInformation(site+1);			
				if(s.polymorphismCase==1){
					count++;
				}
				s = SiteInformation(site+2);			
				if(s.polymorphismCase==1){
					count++;
				}
			}
		}
		return count;
	}

	public double howManyBad(){
		double count=0;
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				// do nothing	
			} else {
				count=count+3;
			}
		}
		return count;
	}





	public static double sampleGamma(double k, double theta) {
		boolean accept = false;
		if (k < 1) {
			// Weibull algorithm
			double c = (1 / k);
			double d = ((1 - k) * Math.pow(k, (k / (1 - k))));
			double u, v, z, e, x;
			do {
				u = rng.nextDouble();
				v = rng.nextDouble();
				z = -Math.log(u);
				e = -Math.log(v);
				x = Math.pow(z, c);
				if ((z + e) >= (d + x)) {
					accept = true;
				}
			} while (!accept);
			return (x * theta);
		} else {
			// Cheng's algorithm
			double b = (k - Math.log(4));
			double c = (k + Math.sqrt(2 * k - 1));
			double lam = Math.sqrt(2 * k - 1);
			double cheng = (1 + Math.log(4.5));
			double u, v, x, y, z, r;
			do {
				u = rng.nextDouble();
				v = rng.nextDouble();
				y = ((1 / lam) * Math.log(v / (1 - v)));
				x = (k * Math.exp(y));
				z = (u * v * v);
				r = (b + (c * y) - x);
				if ((r >= ((4.5 * z) - cheng)) ||
						(r >= Math.log(z))) {
					accept = true;
				}
			} while (!accept);
			return (x * theta);
		}
	}

	public double[] unfoldedSiteFreq(){
		double[] sf = new double[10];
		//	double sil[][] = new double[3][2];
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				//	sil = srSmoothing(site,codon);
				SiteInformation Info = SiteInformation(site);
				if(Info.polymorphismCase>0){
					double numbase = 4;	// number of bases 
					double[] observations = new double[(int)numbase];
					for(int i=0;i<numbase;i++){		
						observations[i] = Info.siteData[i].numObs;	// add observations into an array to find Dirichlet(01+valuesToSampleFrom....0k+valuesToSampleFrom) where valuesToSampleFrom is the prior
					}
					double[] temp = unfoldedCount(observations,Info);
					for(int i=0;i<temp.length;i++){
						sf[i]+=temp[i];
					}
				}
				Info = SiteInformation(site+1);
				if(Info.polymorphismCase>0){
					double numbase = 4;	// number of bases 
					double[] observations = new double[(int)numbase];
					for(int i=0;i<numbase;i++){		
						observations[i] = Info.siteData[i].numObs;	// add observations into an array to find Dirichlet(01+valuesToSampleFrom....0k+valuesToSampleFrom) where valuesToSampleFrom is the prior
					}
					double[] temp = unfoldedCount(observations,Info);
					for(int i=0;i<temp.length;i++){
						sf[i]+=temp[i];;
					}
				}
				Info = SiteInformation(site+2);
				if(Info.polymorphismCase>0){
					double numbase = 4;	// number of bases 
					double[] observations = new double[(int)numbase];
					for(int i=0;i<numbase;i++){		
						observations[i] = Info.siteData[i].numObs;	// add observations into an array to find Dirichlet(01+valuesToSampleFrom....0k+valuesToSampleFrom) where valuesToSampleFrom is the prior
					}
					double[] temp = unfoldedCount(observations,Info);
					for(int i=0;i<temp.length;i++){
						sf[i]+=temp[i];;
					}
				}
			}
		}
		return sf;
	}

	public double[] unfoldedCount(double[] observations, SiteInformation Info){
		int numbase=4;
		double[] sf = new double[10];
		double freq=0;
		for(int x=0;x<numbase;x++){
			if(Info.siteData[x].inAncestral==false){
				freq+=observations[x];
			}
		}
		if(freq>=0 && freq <=0.1){
			sf[0]++;
		}
		if(freq>0.1 && freq <=0.2){
			sf[1]++;
		}
		if(freq>0.2 && freq <=0.3){
			sf[2]++;
		}
		if(freq>0.3 && freq <=0.4){
			sf[3]++;
		}
		if(freq>0.4 && freq <=0.5){
			sf[4]++;
		}
		if(freq>0.5 && freq <=0.6){
			sf[5]++;
		}
		if(freq>0.6 && freq <=0.7){
			sf[6]++;
		}
		if(freq>0.7 && freq <=0.8){
			sf[7]++;
		}
		if(freq>0.8 && freq <=0.9){
			sf[8]++;
		}
		if(freq>0.9 && freq <=1.0){
			sf[9]++;
		}
		return sf;
	}

	public double[] EmpiricalCodonFreq(){
		double[] CodonFreq = new double[AA.length];
		int[] mainbases = new int[3];

		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				for(int sequence=0;sequence<integer_matrix.length;sequence++){	
					mainbases[0]=integer_matrix[sequence][site]; //main bases
					mainbases[1]=integer_matrix[sequence][site+1];
					mainbases[2]=integer_matrix[sequence][site+2];
					int codonNumber =  getcodonnumber(mainbases[0],mainbases[1],mainbases[2]);
					CodonFreq[codonNumber]++;				
				}
			}
		}
		double sum=0;
		for(int i=0;i<CodonFreq.length;i++){
			sum+=CodonFreq[i];
		}
		for(int i=0;i<CodonFreq.length;i++){
			CodonFreq[i]=CodonFreq[i]/sum;
		}

		return CodonFreq;

	}

	public double[][] BaseComposition(){
		double[][] pconfig = new double [3][4];
		int counter1 = 0;
		int counter2 = 0;
		int counter3 = 0;
		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { 
				for(int j=0;j<integer_matrix.length;j++){
					// site 1
					if(integer_matrix[j][site]==1){
						pconfig[0][0]++;
						counter1++;
					}
					if(integer_matrix[j][site]==2){
						pconfig[0][1]++;
						counter1++;
					}
					if(integer_matrix[j][site]==3){
						pconfig[0][2]++;
						counter1++;
					}
					if(integer_matrix[j][site]==4){
						pconfig[0][3]++;
						counter1++;
					}

					// site 2
					if(integer_matrix[j][site+1]==1){
						pconfig[1][0]++;
						counter2++;
					}
					if(integer_matrix[j][site+1]==2){
						pconfig[1][1]++;
						counter2++;
					}
					if(integer_matrix[j][site+1]==3){
						pconfig[1][2]++;
						counter2++;
					}
					if(integer_matrix[j][site+1]==4){
						pconfig[1][3]++;
						counter2++;
					}

					// site 3
					if(integer_matrix[j][site+2]==1){
						pconfig[2][0]++;
						counter3++;
					}
					if(integer_matrix[j][site+2]==2){
						pconfig[2][1]++;
						counter3++;
					}
					if(integer_matrix[j][site+2]==3){
						pconfig[2][2]++;
						counter3++;
					}
					if(integer_matrix[j][site+2]==4){
						pconfig[2][3]++;
						counter3++;
					}

				}
			}
		}

		for(int j=0;j<pconfig[0].length;j++){
			pconfig[0][j] = pconfig[0][j]/counter1;
		}

		for(int j=0;j<pconfig[0].length;j++){
			pconfig[1][j] = pconfig[1][j]/counter2;
		}

		for(int j=0;j<pconfig[0].length;j++){
			pconfig[2][j] = pconfig[2][j]/counter3;
		}
		return pconfig;
	}

	public double[][] TransitionProbabilityMatrix(){
		double[][] CodonTrans = new double[64][64];

		double count=0;
		for(int i=0;i<integer_matrix[0].length-2;i=i+3){
			int anscodonNumber = getcodonnumber(integer_ancestral[i],integer_ancestral[i+1],integer_ancestral[i+2]);
			for(int j=0;j<integer_matrix.length;j++){
				int codonNumber = getcodonnumber(integer_matrix[j][i],integer_matrix[j][i+1],integer_matrix[j][i+2]);
				if(codonNumber != anscodonNumber){
					CodonTrans[anscodonNumber][codonNumber]++;
					count++;
				}
			}
		}

		for(int i=0;i<CodonTrans.length;i++){
			for(int j=0;j<CodonTrans.length;j++){
				CodonTrans[i][j] = CodonTrans[i][j]/count;
			}
		}
		return CodonTrans;

	}

	public double[][] InvariantNew(){
		double count1=0;double count2=0;double count3=0;
		double sil1=0; double rep1=0;
		double sil2=0; double rep2=0;
		double sil3=0; double rep3=0;

		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) {  // check for bad sites
				SiteInformation Info1 = SiteInformation(site);			
				/// position 1
				if(Info1.polymorphismCase==3){
					for(int i=0;i<4;i++){
						if(Info1.siteData[i].rawNumObs<100 && Info1.siteData[i].inAncestral==false){
							double[][] temp =  NGmethod(site,codon);
							sil1 += temp[0][0];
							rep1 += temp[0][1];
							count1++;
						}
					}
				}
				SiteInformation  Info2 = SiteInformation(site+1);			
				/// position 2
				if(Info2.polymorphismCase==3){
					for(int i=0;i<4;i++){
						if(Info2.siteData[i].rawNumObs<100 && Info2.siteData[i].inAncestral==false){
							double[][] temp =  NGmethod(site,codon);
							sil2 += temp[1][0];
							rep2 += temp[1][1];
							count1++;
						}
					}
				}
				SiteInformation Info3 = SiteInformation(site+2);		
				/// position 3
				if(Info3.polymorphismCase==3){
					for(int i=0;i<4;i++){
						if(Info3.siteData[i].rawNumObs<100 && Info3.siteData[i].inAncestral==false){
							double[][] temp =  NGmethod(site,codon);
							sil3 += temp[2][0];
							rep3 += temp[2][1];
							count1++;
						}
					}
				}
			}
		}
		System.out.println(count1+"\t"+count2+"\t"+count3);
		double[][] id = new double[3][2];

		id[0][0] = sil1/(sil1+sil2+sil3);
		id[1][0] = sil2/(sil1+sil2+sil3);
		id[2][0] = sil3/(sil1+sil2+sil3);

		id[0][1] = rep1/(rep1+rep2+rep3);
		id[1][1] = rep2/(rep1+rep2+rep3);
		id[2][1] = rep3/(rep1+rep2+rep3);

		return id;
	}

	public double[][] SRcodonInvariant(){
		double[][] SRcodonStore = new double[AA.length][3];
		double[] freq = F3X4;
		double[][] transitions = TransitionProbabilityMatrix();

		double sil = 0; //silentProb count
		for(int i=1;i<5;i++){ //position 1
			for(int j=1;j<5;j++){
				for(int k=1;k<5;k++){
					int[] codon = {i,j,k};
					int codonNumber = getcodonnumber(i,j,k);
					// position 1 ****************************************************************************************		
					sil=0;	
					double cf=0;
					for(int x=1;x<5;x++){  // loop through all bases
						if(x != codon[0]){ // make sure actual base is not counted
							int IntermediateCodonNumber = getcodonnumber(x,codon[1],codon[2]);
							//		cf+=freq[IntermediateCodonNumber];
							cf+=transitions[codonNumber][IntermediateCodonNumber];
						}
					}
					for(int x=1;x<5;x++){  // loop through all bases
						if(x != codon[0]){ // make sure actual base is not counted
							int[] IntermediateCodon = {x,codon[1],codon[2]};
							int IntermediateCodonNumber = getcodonnumber(x,codon[1],codon[2]);
							double[] SR = NGpathway(codon,IntermediateCodon);
							//		sil += SR[0]*(freq[IntermediateCodonNumber])/cf;
							sil += SR[0]*(transitions[codonNumber][IntermediateCodonNumber]/cf);
						}
					}
					SRcodonStore[codonNumber][0] = sil;
					// position 2 ****************************************************************************************
					sil=0;	
					cf=0;
					for(int x=1;x<5;x++){  // loop through all bases
						if(x != codon[1]){ // make sure actual base is not counted
							int IntermediateCodonNumber = getcodonnumber(codon[0],x,codon[2]);
							//		cf+=freq[IntermediateCodonNumber];
							cf+=transitions[codonNumber][IntermediateCodonNumber];
						}
					}

					for(int x=1;x<5;x++){  // loop through all bases
						if(x != codon[1]){ // make sure actual base is not counted
							int[] IntermediateCodon = {codon[0],x,codon[2]};
							int IntermediateCodonNumber = getcodonnumber(codon[0],x,codon[2]);
							double[] SR = NGpathway(codon,IntermediateCodon);
							//		sil += SR[1]*(freq[IntermediateCodonNumber])/cf;
							sil += SR[1]*(transitions[codonNumber][IntermediateCodonNumber]/cf);
						}
					}
					SRcodonStore[codonNumber][1] = sil;
					// position 3 ****************************************************************************************
					sil=0;	
					cf=0;
					for(int x=1;x<5;x++){  // loop through all bases
						if(x != codon[2]){ // make sure actual base is not counted
							int IntermediateCodonNumber = getcodonnumber(codon[0],codon[1],x);
							//		cf+=freq[IntermediateCodonNumber];
							cf+=transitions[codonNumber][IntermediateCodonNumber];
						}
					}


					for(int x=1;x<5;x++){  // loop through all bases
						if(x != codon[2]){ // make sure actual base is not counted
							int[] IntermediateCodon = {codon[0],codon[1],x};
							int IntermediateCodonNumber = getcodonnumber(codon[0],codon[1],x);
							double[] SR = NGpathway(codon,IntermediateCodon);
							//		sil += SR[2]*(freq[IntermediateCodonNumber])/cf;
							sil += SR[2]*(transitions[codonNumber][IntermediateCodonNumber]/cf);
						}
					}
					SRcodonStore[codonNumber][2] = sil;
				}
			}		
		}

		return SRcodonStore;
	}

	public void PrintMatrix(double[][] mat){
		for(int i=0;i<mat[0].length;i++){
			for(int j=0;j<mat.length;j++){
				System.out.print(mat[i][j]+"\t");
			}
			System.out.println();
		}
	}


}




