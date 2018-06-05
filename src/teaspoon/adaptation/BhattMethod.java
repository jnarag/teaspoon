package teaspoon.adaptation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import teaspoon.app.utils.NullNeutralRatioException;
import teaspoon.app.utils.TeaspoonMethods;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * <p>BhattMethod is a class for counting of synonymous (silent) and nonsynonymous (replacement) 
 * polymorphisms to obtain site-frequency class information to calculate neutral rates of evolution
 * and numbers of adaptive substitutions.
 * <p>It includes methods to obtain observed changes by fractional counting, expected changes by 
 * probabalistic sampling, and utility methods.
 * @see doi:10.1371/journal.pcbi.1004694.g001 Raghwani et al (2016) for details
 * 
 * @author <a href="http://github.com/jnarag">@jnarag</a>
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 5 Jun 2018
 * @version 0.1
 */
public class BhattMethod {

    double numReplicates = 500.0; //set number of replicates
    double numBins;
    boolean[] whichBins;
    double[][] binsMatrix;
    int numSample;
    public final  int[][] integerMatrix;
    public final int[] integerAncestralArray;
    double neutralRatio;
    double neutralBins;
    private double[] silentCountArray;
    private double[] replacementCountArray;
    double[] totalCountArray;
    double[] totalCountArrayInvariantSitesExcluded;
    double[] totalCountArrayInvariantSitesOnly;
    private double[] replacementToSilentRatio;
    private double[] nonNeutralSubstitutions;
    double adaptation;
    double deleteriousLoad;
    public  int[][] codonMatrix;
    public  boolean[] badSitesList;

    public final String[] AA =	
    	{"K","N","K","N","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H","P","P","P","P",
            "R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V",
            "X","Y","X","Y","S","S","S","S","X","C","W","C","L","F","L","F","?","-","?" };

    @Deprecated
    /**
     * Deprecated no-arg constructor.
     * @throws RuntimeException
     */
    public BhattMethod(){
        throw new RuntimeException("please input the raw integer matrix and the ancestral matrix");
    }

    /**
     * 
     * @param m
     * @param a
     */
    public BhattMethod(int[][] m,int[] a){
        this.integerMatrix = m;
        this.integerAncestralArray = a;
        badSitesList = findInvalidSites(integerMatrix, integerAncestralArray);
    }

    /**
     * Set the neutral ratio.
     * @param newNR
     */
    public void setNeutralRatio(double newNR){
        this.neutralRatio = newNR;
    }

    /**
     * 
     * @param FullStore
     * @param prior
     * @param needprior
     * @return
     */
    public ArrayList<Mutation> tracking(ArrayList<Mutation> FullStore,double[] prior,boolean needprior){
    	int c=0;
    	for(int site=0;site<integerMatrix[0].length;site++){
    		//************************************************************************************************************************
    		// dirichlet site frequency spectrum
    		SiteInformation Info = calculateSiteInformation(site);
    		double numbase = 4;	// number of bases
    		double[] observations = new double[(int)numbase];

    		if(needprior){
    			for(int i=0;i<numbase;i++){
    				if(Info.siteData[i].inAncestral){
    					prior[i]=1.0;
    				}else{
    					prior[i]=1.0/3.0; // change later
    				}
    			}
    		}

    		for(int i=0;i<numbase;i++){
    			Info.siteData[i].prior = prior[i];
    			/*
    			 *  add observations into an array to find 
    			 *  Dirichlet(01+valuesToSampleFrom....0k+valuesToSampleFrom) 
    			 *  where valuesToSampleFrom is the prior
    			 */
    			observations[i] = Info.siteData[i].rawNumObs+Info.siteData[i].prior;	
    		}

    		TeaspoonRandomSampler S = new TeaspoonRandomSampler(observations);
    		double[][] dist = new double[(int)numReplicates][4];
    		double[] mean = new double[4];

    		for(int i=0;i<(int) numReplicates;i++){ // Sampler loop
    			double[] point = S.sampleDirichlet();
    			dist[i]=point;
    			for(int x=0;x<numbase;x++){
    				mean[x]+=point[x];
    			}
    		}	// end loop

    		mean[0]=mean[0]/numReplicates;
    		mean[1]=mean[1]/numReplicates;
    		mean[2]=mean[2]/numReplicates;
    		mean[3]=mean[3]/numReplicates;

    		double[] A = new double[(int)numReplicates];
    		double[] C = new double[(int)numReplicates];
    		double[] G = new double[(int)numReplicates];
    		double[] T = new double[(int)numReplicates];

    		for(int s=0;s<numReplicates;s++){
    			A[s] = dist[s][0];
    			C[s] = dist[s][1];
    			G[s] = dist[s][2];
    			T[s] = dist[s][3];
    		}

    		Arrays.sort(A);
    		Arrays.sort(C);
    		Arrays.sort(G);
    		Arrays.sort(T);

    		for(int s=0;s<numReplicates;s++){
    			dist[s][0]=A[s];
    			dist[s][1]=C[s];
    			dist[s][2]=G[s];
    			dist[s][3]=T[s];
    		}

    		double high = numReplicates*0.975;double low = numReplicates*0.025;

    		for(int x=0;x<numbase;x++){
    			if(Info.siteData[x].inAncestral==false){
    				FullStore.get(c).addProb(mean[x]);
    				FullStore.get(c).addH(dist[(int) high][x]);
    				FullStore.get(c).addL(dist[(int) low][x]);
    				c++;
    			}
    		}
    	}
    	return FullStore;
    }

    /**
     * Estimate expected site freqs from Dirichlet
     * @param u - site frequency range lower
     * @param v - site frequency range upper
     * @param site
     * @param prior
     * @param needprior
     * @return
     */
    public SiteInformation dirichletSiteFreq(double u, double v,int site,double[] prior,boolean needprior) {
        SiteInformation Info = calculateSiteInformation(site);
        double numbase = 4;	// number of bases
        double[] observations = new double[(int)numbase];
        double[] sampler = new double[(int) numReplicates];

        //set the prior - otherwise it will be a flat prior?
        if(needprior){
            for(int i=0;i<numbase;i++){
                if(Info.siteData[i].inAncestral){
                    prior[i]=1.0;       //I think this should be a class variable
                }else{                  // Also I don't understand how this gives you a uniform distribution?
                    prior[i]=1.0/3.0;
                }
            }
        }

        // no. of observations is equal to numbases
        for(int i=0;i<numbase;i++){
            Info.siteData[i].prior = prior[i]; // is prior a prior probability or a  parameter of a distirbution?
            observations[i] = Info.siteData[i].rawNumObs + Info.siteData[i].prior;	// add observations into an array to find Dirichlet(01+valuesToSampleFrom....0k+valuesToSampleFrom) where valuesToSampleFrom is the prior
        }

        TeaspoonRandomSampler S = new TeaspoonRandomSampler(observations);

        double count = 0;
        //bootstrap replicate based on observations (which are essential prior plus observations)
        for(int i=0;i<(int) numReplicates;i++){
            double[] point = S.sampleDirichlet();      // randomly draw the base frequencies from a dirichlet distribution which is based on the observed+prior
 
            for(int x=0;x<numbase;x++){
                if(Info.siteData[x].inAncestral==false){
                    sampler[i]+=point[x];    //derived freq from the Dirichlet distribution 500 X
                }
            }
            
            if(sampler[i]>u && sampler[i]<=v){
                count++;
            }
        }

        Info.dirichletProb = (count/numReplicates);    // average derived frequency (posterior probability)
        return Info;
    }

    /**
     * Estimate expected site freqs from beta-binomial
     * @param u - site frequency range lower
     * @param v - site frequency range upper
     * @param site
     * @param prior
     * @param needprior
     * @return
     */
    public SiteInformation betaSiteFreq(double u, double v, int site, double[] prior, boolean needprior) {

    	SiteInformation Info = calculateSiteInformation(site);
    	double[] observations = new double[2];
    	double[] sampler = new double[(int)numReplicates];

    	//numbase should be a class variable

    	if(needprior){
    		prior[0] = 1;       //using a uniform prior...
            prior[1] = 1;
        }

        for(int i = 0; i < 4; i++) {
            if(Info.siteData[i].inAncestral==false) {
                observations[0] += (Info.siteData[i].rawNumObs);
            }
            else{
                observations[1] = Info.siteData[i].rawNumObs ;//+ prior[0];
            }
        }

        observations[0] = observations[0]+prior[0];   //derived frequency
        observations[1] = observations[1]+prior[1];   // ancestral frequency
        TeaspoonRandomSampler S = new TeaspoonRandomSampler(observations);
        double count = 0;

        for(int i=0; i < (int)numReplicates; i++) {
            double point = S.sampleBeta();
            sampler[i] = point;

            if(sampler[i]> u && sampler[i] < v){
                 count++;
            }
        }

        Info.dirichletProb = (count/numReplicates);

        return Info;
    }

    /**
     * Estimate expected site freqs from Dirichlet, silent/replacement
     * @param u - site frequency range lower
     * @param v - site frequency range upper
     * @param site
     * @param prior
     * @param needprior
     * @param pos
     * @return
     */
    public SiteInformation dirichletSRFreq(double u, double v,int site,double[] prior,boolean needprior,int pos) {
        SiteInformation siteInformation = calculateSiteInformation(site);
        double numbase = 4;	// number of bases
        double[] observations = new double[(int)numbase];
        double[] sampler = new double[(int) numReplicates];

        if(needprior){
            for(int i=0;i<numbase;i++){
                if(siteInformation.siteData[i].inAncestral){
                    prior[i]=1.0;
                }else{
                    prior[i]=1.0/3.0;
                }
            }
        }
 
        for(int i=0;i<numbase;i++){
            siteInformation.siteData[i].prior = prior[i];
            observations[i] = siteInformation.siteData[i].rawNumObs+siteInformation.siteData[i].prior;	// add observations into an array to find Dirichlet(01+valuesToSampleFrom....0k+valuesToSampleFrom) where valuesToSampleFrom is the prior
        }

        TeaspoonRandomSampler randomSampler = new TeaspoonRandomSampler(observations);
        double[] isSil = new double[4];

        //checking if site is silentProb or replacement
        for(int i=0;i<numbase;i++){
            if(pos==1){
                int tmp = getcodonnumber(i+1,integerAncestralArray[site+1],integerAncestralArray[site+2]);
                int tmp2 = getcodonnumber(integerAncestralArray[site],integerAncestralArray[site+1],integerAncestralArray[site+2]);
                isSil[i]=determineSilentOrReplacement(tmp,tmp2);
            }

            if(pos==2){
                int tmp = getcodonnumber(integerAncestralArray[site-1],i+1,integerAncestralArray[site+1]);
                int tmp2 = getcodonnumber(integerAncestralArray[site-1],integerAncestralArray[site],integerAncestralArray[site+1]);
                isSil[i]=determineSilentOrReplacement(tmp,tmp2);
            }

            if(pos==3){
                int tmp = getcodonnumber(integerAncestralArray[site-2],integerAncestralArray[site-1],i+1);
                int tmp2 = getcodonnumber(integerAncestralArray[site-2],integerAncestralArray[site-1],integerAncestralArray[site]);
                isSil[i]=determineSilentOrReplacement(tmp,tmp2);
            }
        }

        double[] silentSampler = new double[(int) numReplicates];
        double[] replacementSampler = new double[(int) numReplicates];
        double silent=0;
        double replacement=0;
        double count = 0;

        for(int i=0;i<(int) numReplicates;i++){
            double[] point = randomSampler.sampleDirichlet();
 
            for(int x=0;x<numbase;x++){
                if(siteInformation.siteData[x].inAncestral==false){  //site is invariant or fixed  (no ancestral site?)
                    sampler[i]+=point[x];
                    silentSampler[i]+=point[x]*isSil[x];
                    replacementSampler[i]+=point[x]*(1.0-isSil[x]);
                }
            }

            if(sampler[i]>u && sampler[i]<=v){
                count++;
            }
            
            if(silentSampler[i]>u && silentSampler[i]<=v){
                silent++;
            }
            
            if(replacementSampler[i]>u && replacementSampler[i]<=v){
                replacement++;
            }
        }

        siteInformation.dirichletProb = (count/numReplicates);
        siteInformation.silentProb = (silent/numReplicates);
        siteInformation.replacementProb = (replacement/numReplicates);
        return siteInformation;
    }

    /**
     * Probably Nei-Gojobori (1986) stuff 
     * @param site
     * @param codonsite
     * @return
     */
    public double[][] nGmethod(int site,int codonsite){
        // HANDLES VARINT SITES ***************************************************
        // main NG method for all sites
        double[][] identity = new double[3][2];
        int[] ancestralbases = new int[3];
        ancestralbases[0]=integerAncestralArray[site];
        ancestralbases[1]=integerAncestralArray[site+1];
        ancestralbases[2]=integerAncestralArray[site+2];
        int[] mainbases = new int[3];
        int[] count = new int[3];

        for(int i=0;i<integerMatrix.length;i++){
            // adds codon bases
            mainbases[0]=integerMatrix[i][site];
            mainbases[1]=integerMatrix[i][site+1];
            mainbases[2]=integerMatrix[i][site+2];
            double[] tmp = nGpathway(ancestralbases,mainbases);  // per site
            if(tmp[0]!=2.0 && tmp[0]!=3.0 ){
                identity[0][0] += tmp[0];
                count[0]++;
            }
            if(tmp[1]!=2.0 && tmp[1]!=3.0 ){
                identity[1][0] += tmp[1];
                count[1]++;
            }
            if(tmp[2]!=2.0 && tmp[2]!=3.0 ){
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

    /**
     * 
     * Probably Nei-Gojobori (1986) stuff 
     * @param a
     * @param b
     * @return
     */
    public double[] nGpathway(int[] a, int[] b){
    //a = ancestral b=main
        // 1 silentProb, 0 replacement, 2 invariant
        int degen = 0;
        int AnscodonNumber =  getcodonnumber(a[0],a[1],a[2]);
        int codonNumber =  getcodonnumber(b[0],b[1],b[2]);
        boolean[] diff = whichDiff(a,b);


        double[] identity = new double[3];
        identity[0]=3;identity[1]=3;identity[2]=3;
        if(b[0]<5 || b[1]<5 || b[2]<5){  // Checks whether a site is a gap.
            identity[0]=2;identity[1]=2;identity[2]=2;  // if no gaps set identity to assume invariance i.e 2
            for(int i=0;i<3;i++){if(diff[i]){degen++;}} // find how many sites are degererate
            // 1 site different ****************************************************************************************
            if(degen==1){
                if(diff[0]){identity[0]=determineSilentOrReplacement(AnscodonNumber,codonNumber);}
                if(diff[1]){identity[1]=determineSilentOrReplacement(AnscodonNumber,codonNumber);}
                if(diff[2]){identity[2]=determineSilentOrReplacement(AnscodonNumber,codonNumber);}
            }
  
            // 2 sites different ****************************************************************************************
            if(degen==2){
                //			Pathway 1  (a) X11-X21-X22 and  (b) X11-X12-X22	  *Both pathways occur with equal probability*
                if(diff[0]==false){
                    // (a) X11-X21-X22
                    int a_codonNumber =  getcodonnumber(a[0],b[1],a[2]);
                    identity[1]=determineSilentOrReplacement(AnscodonNumber,a_codonNumber);
                    identity[2]=determineSilentOrReplacement(a_codonNumber,codonNumber);
                    // (b) X11-X12-X22
                    int b_codonNumber =  getcodonnumber(a[0],a[1],b[2]);
                    identity[2]+=determineSilentOrReplacement(AnscodonNumber,b_codonNumber);
                    identity[1]+=determineSilentOrReplacement(b_codonNumber,codonNumber);
                    identity[1]=identity[1]/2;identity[2]=identity[2]/2; //average over all pathways
                }
                //			Pathway 2  (a) 1X1-2X1-2X2 and (b) 1X1-1X2-2X2
                if(diff[1]==false){
                    // (a) 1X1-2X1-2X2
                    int a_codonNumber =  getcodonnumber(b[0],a[1],a[2]);
                    identity[0]=determineSilentOrReplacement(AnscodonNumber,a_codonNumber);
                    identity[2]=determineSilentOrReplacement(a_codonNumber,codonNumber);

                    // (b) 1X1-1X2-2X2
                    int b_codonNumber =  getcodonnumber(a[0],a[1],b[2]);
                    identity[2]+=determineSilentOrReplacement(AnscodonNumber,b_codonNumber);
                    identity[0]+=determineSilentOrReplacement(b_codonNumber,codonNumber);
                    identity[0]=identity[0]/2;identity[2]=identity[2]/2; //average over all pathways
                }
                //			Pathway 3   (a) 11X-21X-22X and  (b) 11X-12X-22X
                if(diff[2]==false){
                    // (a) 11X-21X-22X
                    int a_codonNumber =  getcodonnumber(b[0],a[1],a[2]);
                    identity[0]=determineSilentOrReplacement(AnscodonNumber,a_codonNumber);
                    identity[1]=determineSilentOrReplacement(a_codonNumber,codonNumber);
                    // (b) 11X-12X-22X
                    int b_codonNumber =  getcodonnumber(a[0],b[1],a[2]);
                    identity[1]+=determineSilentOrReplacement(AnscodonNumber,b_codonNumber);
                    identity[0]+=determineSilentOrReplacement(b_codonNumber,codonNumber);
                    identity[0]=identity[0]/2;identity[1]=identity[1]/2; //average over all pathways
                }
            }

            // 3 sites different ****************************************************************************************
            if(degen==3){
                // Pathway 1 111-211-221-222	*x -(a)-(b)- x*
                int a_codonNumber = getcodonnumber(b[0],a[1],a[2]);
                int b_codonNumber = getcodonnumber(b[0],b[1],a[2]);
                identity[0]=determineSilentOrReplacement(AnscodonNumber,a_codonNumber);
                identity[1]=determineSilentOrReplacement(a_codonNumber,b_codonNumber);
                identity[2]=determineSilentOrReplacement(b_codonNumber,codonNumber);
                // Pathway 2 111-211-212-222
                a_codonNumber = getcodonnumber(b[0],a[1],a[2]);
                b_codonNumber = getcodonnumber(b[0],a[1],b[2]);
                identity[0]+=determineSilentOrReplacement(AnscodonNumber,a_codonNumber);
                identity[2]+=determineSilentOrReplacement(a_codonNumber,b_codonNumber);
                identity[1]+=determineSilentOrReplacement(b_codonNumber,codonNumber);
                // Pathway 3 111-121-221-222
                a_codonNumber = getcodonnumber(a[0],b[1],a[2]);
                b_codonNumber = getcodonnumber(b[0],b[1],a[2]);
                identity[1]+=determineSilentOrReplacement(AnscodonNumber,a_codonNumber);
                identity[0]+=determineSilentOrReplacement(a_codonNumber,b_codonNumber);
                identity[2]+=determineSilentOrReplacement(b_codonNumber,codonNumber);
                // Pathway 4  111-121-122-222
                a_codonNumber = getcodonnumber(a[0],b[1],a[2]);
                b_codonNumber = getcodonnumber(a[0],b[1],b[2]);
                identity[1]+=determineSilentOrReplacement(AnscodonNumber,a_codonNumber);
                identity[2]+=determineSilentOrReplacement(a_codonNumber,b_codonNumber);
                identity[0]+=determineSilentOrReplacement(b_codonNumber,codonNumber);
                // Pathway 5 111-112-212-222
                a_codonNumber = getcodonnumber(a[0],a[1],b[2]);
                b_codonNumber = getcodonnumber(b[0],a[1],b[2]);
                identity[2]+=determineSilentOrReplacement(AnscodonNumber,a_codonNumber);
                identity[0]+=determineSilentOrReplacement(a_codonNumber,b_codonNumber);
                identity[1]+=determineSilentOrReplacement(b_codonNumber,codonNumber);
                // Pathway 6 111-112-122-222
                a_codonNumber = getcodonnumber(a[0],a[1],b[2]);
                b_codonNumber = getcodonnumber(a[0],b[1],b[2]);
                identity[2]+=determineSilentOrReplacement(AnscodonNumber,a_codonNumber);
                identity[1]+=determineSilentOrReplacement(a_codonNumber,b_codonNumber);
                identity[0]+=determineSilentOrReplacement(b_codonNumber,codonNumber);
                identity[0]=identity[0]/6.0;identity[1]=identity[1]/6.0;identity[2]=identity[2]/6.0;
            }
        }

        //System.out.println(identity);
        return identity;
    }

    /**
     * No idea what this does..
     * Probably Nei-Gojobori (1986) stuff 
     * @return
     */
    public double[][] nGpossible(){
        int[] mainbases = new int[3];
        double L = (double) integerMatrix.length;
        double N = (double) integerMatrix[0].length;
        double[][] identity = new double[3][2];
        double[][] FinalIdentity = new double[(int)N][2];
        for(int site=0,codon=0; site<integerMatrix[0].length-2;site=site+3,codon++){
            identity = new double[3][2]; // set new identity vector
            if (badSitesList[site] == false && badSitesList[site+1] == false && badSitesList[site+2] == false) {  // check for bad sites
                for(int sequence=0;sequence<integerMatrix.length;sequence++){		// loop through sequences

                    mainbases[0]=integerMatrix[sequence][site]; //main bases
                    mainbases[1]=integerMatrix[sequence][site+1];
                    mainbases[2]=integerMatrix[sequence][site+2];
                    double sil = 0; //silentProb count
                    double rep=0; //rep count
                    //********************************************************************************************************************************************
                    //pos 1
                    sil=0;
                    for(int i=1;i<5;i++){  // loop through all bases
                        if(i != mainbases[0]){ // make sure actual base is not counted
                            int[] actualCodon = {mainbases[0],mainbases[1],mainbases[2]};
                            int[] IntermediateCodon = {i,mainbases[1],mainbases[2]};
                            double[] SR = nGpathway(actualCodon,IntermediateCodon);
                            sil += SR[0]*(1.0/3.0);
                        }
                    }
                    rep = 1.0-sil;
                    identity[0][0]+=sil;identity[0][1]+=rep;
                    //********************************************************************************************************************************************
                    //pos 2
                    sil=0;
                    for(int i=1;i<5;i++){  // loop through all bases
                        if(i != mainbases[1]){ // make sure actual base is not counted
                            int[] actualCodon = {mainbases[0],mainbases[1],mainbases[2]};
                            int[] IntermediateCodon = {mainbases[0],i,mainbases[2]};
                            double[] SR = nGpathway(actualCodon,IntermediateCodon);
                            sil += SR[1]*(1.0/3.0);
                        }
                    }
                    rep = 1.0-sil;
                    identity[1][0]+=sil;identity[1][1]+=rep;
                    //********************************************************************************************************************************************
                    //pos 3
                    sil=0;
                    for(int i=1;i<5;i++){  // loop through all bases
                        if(i != mainbases[2]){ // make sure actual base is not counted
                            int[] actualCodon = {mainbases[0],mainbases[1],mainbases[2]};
                            int[] IntermediateCodon = {mainbases[0],mainbases[1],i};
                            double[] SR = nGpathway(actualCodon,IntermediateCodon);
                            sil += SR[2]*(1.0/3.0);
                        }
                    }
                    rep = 1.0-sil;
                    identity[2][0]+=sil;identity[2][1]+=rep;
                    //********************************************************************************************************************************************
                }
            }
            FinalIdentity[site][0]=identity[0][0]/L;
            FinalIdentity[site+1][0]=identity[1][0]/L;
            FinalIdentity[site+2][0]=identity[2][0]/L;
            FinalIdentity[site][1]=identity[0][1]/L;
            FinalIdentity[site+1][1]=identity[1][1]/L;
            FinalIdentity[site+2][1]=identity[2][1]/L;


        }
        return FinalIdentity;

    }

    /**
     * calculates a silentProb and replacement site frequency between ranges u and v
     * @param u
     * @param v
     * @param prior
     * @param needPrior
     * @return double[] list of data: 
     * 	silentProb count [0]; 
     * 	replacement count [1]; 
     * 	total count [2]; 
     * 	total variant [3]; 
     * 	total invariant [4]
     */
    public double[] calculateSiteFrequenciesWithinInterval(double u,double v, double[] prior, boolean needPrior){
        double[][] temp = new double[3][2];
        double rho1 = 0.0;
        double sigma1 = 0.0;
        double[] finalans = new double[5];
        double total1 = 0.0;
        double rho2 = 0.0;
        double sigma2 = 0.0;
        double total2 = 0.0;

        ArrayList<SiteInformation> Inv = new ArrayList<SiteInformation>();

        double old_rho = 0;
        for(int site=0,codon=0; site<integerMatrix[0].length-2;site+=3,codon++){
            //System.out.println(site+","+codon);

            if (badSitesList[site] == false && badSitesList[site+1] == false && badSitesList[site+2] == false) {
                temp = nGmethod(site,codon);
                // find site freq for pos1
                SiteInformation inf = betaSiteFreq(u,v,site,prior,needPrior);
                double raw_rho = 0;

                //System.out.println(site+","+inf.Dprob);
                //teaspoon.adaptation.SiteInformation inf = DirichletSRFreq(u,v,site,prior,needPrior, 1);
                if(inf.polymorphismCase==1){
                    Inv.add(inf);
                } else {
                    total1 += inf.dirichletProb;
                    sigma1 += inf.dirichletProb*temp[0][0];
                    rho1 += inf.dirichletProb*temp[0][1];
                    raw_rho += temp[0][1];

                }

                // find site freq for pos2
                SiteInformation inf2 = betaSiteFreq(u,v,site+1,prior,needPrior);
                //System.out.println(site+1);
                //teaspoon.adaptation.SiteInformation inf2 = DirichletSRFreq(u,v,site+1,prior,needPrior, 2);
                if(inf2.polymorphismCase==1){
                    Inv.add(inf2);
                } else {
                    total1 += inf2.dirichletProb;
                    sigma1 += inf2.dirichletProb*temp[1][0];
                    rho1 += inf2.dirichletProb*temp[1][1];
                    raw_rho += temp[1][1];
                }

                // find site freq for pos3
                SiteInformation inf3 = betaSiteFreq(u,v,site+2,prior,needPrior);
                //teaspoon.adaptation.SiteInformation inf3 = DirichletSRFreq(u,v,site+2,prior,needPrior, 3);
                if(inf3.polymorphismCase==1){
                    Inv.add(inf3);
                } else {
                    total1 += inf3.dirichletProb;
                    sigma1 += inf3.dirichletProb*temp[2][0];
                    rho1 += inf3.dirichletProb*temp[2][1];
                    raw_rho += temp[2][1];
                }
            }
        }

        //System.out.println(highFreqRSites);
        double S = sigma1/total1;   // proportion of variant sites that are silentProb
        double R = rho1/total1;     // proportion of variant sites that are replacement

        Iterator<SiteInformation> It =  Inv.iterator();

        // infers the frequencies for invariant sites probabilistically
        while(It.hasNext()){
            SiteInformation Element = It.next();
            total2 += Element.dirichletProb;     // posterior frequency of observed invariant site
            sigma2 += Element.dirichletProb*S;   // posterior frequency of observed invariant site being silentProb
            rho2 += Element.dirichletProb*R;	 // posterior frequency of observed invariant site being replacement
        }

        finalans[0] = sigma1+sigma2;  // silentProb count
        finalans[1] = rho1+rho2;      // replacement count
        finalans[2] = total1+total2;  // total count
        finalans[3] = total1;         // total variant
        finalans[4] = total2;         // total invariant

        return finalans;
    }
 
    /**
     * calculates a silentProb and replacement site frequency between ranges u and v
     * <p><i>Note: From Samir Bhatt's thesis pp. 145 (Discussion,  section 5.4); seems like this used the Nei & Gojobori (1986) counting method.</i>
     * @param u
     * @param v
     * @param prior
     * @param needPrior
     * @return
     * @see BhatMethod#MethodNG
     */
    public double[] calculateSiteFrequenciesWithinIntervalNG(double u,double v, double[] prior, boolean needPrior){
        double[][] tempMatrix = new double[3][2];
        double[] outputArray = new double[5];
        double[][] nGpossible = nGpossible();
        double rho 		= 0.0;
        double sigma 	= 0.0;
        double total 	= 0.0;
        double rhoP 	= 0.0;
        double sigmaP 	= 0.0;


        ArrayList<SiteInformation> invariants = new ArrayList<SiteInformation>();
        for(int siteIndex=0,codonIndex=0; siteIndex<integerMatrix[0].length-2;siteIndex=siteIndex+3,codonIndex++){
            if (badSitesList[siteIndex] == false && badSitesList[siteIndex+1] == false && badSitesList[siteIndex+2] == false) {
                tempMatrix = nGmethod(siteIndex,codonIndex);
                //inf - sitefreq for a specific nucleotide in the alignment

                //find site freq for pos1
                SiteInformation siteInformationPosOne = betaSiteFreq(u,v,siteIndex,prior,needPrior);
                if(siteInformationPosOne.polymorphismCase==1){
                    invariants.add(siteInformationPosOne);
                } else {
                    total += siteInformationPosOne.dirichletProb;
                    sigma += siteInformationPosOne.dirichletProb*tempMatrix[0][0];
                    rho += siteInformationPosOne.dirichletProb*tempMatrix[0][1];
                    sigmaP+=nGpossible[siteIndex][0];
                    rhoP+=nGpossible[siteIndex][1];
                }

                // find site freq for pos2
                SiteInformation siteInformationPosTwo = betaSiteFreq(u,v,siteIndex+1,prior,needPrior);
                if(siteInformationPosTwo.polymorphismCase==1){
                    invariants.add(siteInformationPosTwo);
                } else {
                    total += siteInformationPosTwo.dirichletProb;
                    sigma += siteInformationPosTwo.dirichletProb*tempMatrix[1][0];
                    rho += siteInformationPosTwo.dirichletProb*tempMatrix[1][1];
                    sigmaP+=nGpossible[siteIndex+1][0];
                    rhoP+=nGpossible[siteIndex+1][1];
                }

                // find site freq for pos3
                SiteInformation siteInformationPosThree = betaSiteFreq(u,v,siteIndex+2,prior,needPrior);
                if(siteInformationPosThree.polymorphismCase==1){
                    invariants.add(siteInformationPosThree);
                } else {
                    total += siteInformationPosThree.dirichletProb;
                    sigma += siteInformationPosThree.dirichletProb*tempMatrix[2][0];
                    rho +=  siteInformationPosThree.dirichletProb*tempMatrix[2][1];
                    sigmaP+=nGpossible[siteIndex+2][0];
                    rhoP+=nGpossible[siteIndex+2][1];
                }
            }
        }

        //So why is this is section different from the equivalent code in SiteFreq(...) method?
        //This does not contain the raw counts and has less entries in the final matrix
        double S = sigma/sigmaP;
        double R = rho/sigmaP;
        double Dn = (-3.0/4.0)*Math.log(1.0-(R*4.0/3.0));
        double Ds = (-3.0/4.0)*Math.log(1.0-(S*4.0/3.0));

        //finalans length = 5, so what are the values for 4th and 5th element
        outputArray[0] = Ds; outputArray[1]=Dn; outputArray[2] = total;

//        finalans[0] = sigma1+sigmaP;
//        finalans[1] = rho1+rhoP;
//        finalans[2] = total1;//+totalP;
//        finalans[3] = total1;
//        finalans[4] = total2;

        return outputArray;
    }

    /**
     * From Samir Bhatt's thesis pp. 145 (Discussion - section 5.4); seems like this used the Nei & Gojobori (1986) counting method.
     * @param binsvalues
     * @param prior
     * @param needPrior
     * @param which
     */
    public void inferCountsWithNeiGojoboriCounting(double[][] binsvalues,double[] prior, boolean needPrior,boolean[] which){
        this.numBins = binsvalues[0].length;
        double[][] finalmat = new double[6][binsvalues[0].length];
        double[][] totals = new double[3][binsvalues[0].length];
        for(int i=0;i< (int) numBins;i++){
            double[] temp = calculateSiteFrequenciesWithinIntervalNG(binsvalues[0][i], binsvalues[1][i],prior,needPrior);
            finalmat[0][i] = temp[0];   // number silentProb  - Ds
            finalmat[1][i] = temp[1];	// number replacement - Dn

            totals[0][i] = temp[2];
            totals[1][i] = temp[3];
            totals[2][i] = temp[4];
            if(temp[1]!=0){
                finalmat[3][i] = temp[1]/temp[0];	// Replacement/Silent ratio
            } else {
                finalmat[3][i] = Double.NaN;
            }
        }
        this.setSilentSubstitutionsCountArray(finalmat[0]);
        this.setReplacementSubstitutionsCountArray(finalmat[1]);
        this.totalCountArray = totals[0];
        this.totalCountArrayInvariantSitesExcluded = totals[1];
        this.totalCountArrayInvariantSitesOnly = totals[2];
        this.setReplacementToSilentRatesRatio(finalmat[3]);
        this.whichBins = which;
        // calcualte neutral ratio. Vector NeutralVec provides boolean true false if a given bin is neutral
        // user has to specify this bin
        double counter = 0;
        double NR=0;
        for(int i=0;i<(int) numBins;i++){
            if(whichBins[i]){
                if(finalmat[0][i]!=0){
                    NR += finalmat[1][i]/finalmat[0][i];
                    counter	++;
                }
            }
        }
        neutralRatio = NR/counter; //average

        // number non neutral
        for(int i=0;i<(int) numBins;i++){
            finalmat[5][i] = finalmat[1][i]*(1.0 - ((finalmat[0][i]/finalmat[1][i])*(neutralRatio)));	// number non neutral
            if(finalmat[5][i]<0){
                finalmat[5][i]=0.0;
            }
        }
        this.setNonNeutralSubstitutions(finalmat[5]);
        int flag=0;
        for(int i=0;i< (int) numBins;i++){
            if(whichBins[i]==false && flag==0){
                this.deleteriousLoad+=getNonNeutralSubstitutions()[i];
            }
            if(whichBins[i]){flag=1;}
            if(whichBins[i]==false && flag==1){
                this.adaptation+=getNonNeutralSubstitutions()[i];
            }
        }
    }

    /**
     * Does not seem to be used anywhere.
     * TODO consider deprecating this
     * @since 2018 May 22
     * @param binsvalues
     * @param prior
     * @param needPrior
     */
    public void inferCountsUnusedBhattMethod(double[][] binsvalues,double[] prior, boolean needPrior){
        this.numBins = binsvalues[0].length;
        double[][] finalmat = new double[6][binsvalues[0].length];
        for(int i=0;i< (int) numBins;i++){
            double[] temp = calculateSiteFrequenciesWithinInterval(binsvalues[0][i], binsvalues[1][i],prior,needPrior);
            finalmat[0][i] = temp[0];   // number silentProb
            finalmat[1][i] = temp[1];	// number replacement
            finalmat[2][i] = temp[2];	// total number

            if(temp[1]!=0){
                finalmat[3][i] = temp[1]/temp[0];	// Silent/replacement ratio
            } else {
                finalmat[3][i] = Double.NaN;
            }
        }
        this.setSilentSubstitutionsCountArray(finalmat[0]);
        this.setReplacementSubstitutionsCountArray(finalmat[1]);
        this.totalCountArray = finalmat[2];
        this.setReplacementToSilentRatesRatio(finalmat[3]);
        this.neutralRatio = getReplacementToSilentRatesRatio()[0];
        this.neutralBins = 0;

        for(int i=0;i< (int) numBins;i++){
            if(getReplacementToSilentRatesRatio()[i]<neutralRatio){
                this.neutralRatio=getReplacementToSilentRatesRatio()[i];
                this.neutralBins = i;
            }
        }



        // number non neutral
        for(int i=0;i<(int) numBins;i++){
            finalmat[5][i] = finalmat[1][i]*(1.0 - ((finalmat[0][i]/finalmat[1][i])*(neutralRatio)));	// number non neutral
            if(finalmat[5][i]<0){
                finalmat[5][i]=0.0;
            }
        }
        this.setNonNeutralSubstitutions(finalmat[5]);
        for(int i=0;i< (int) numBins;i++){
            if(i<neutralBins){
                this.deleteriousLoad+=getNonNeutralSubstitutions()[i];
            }
            if(i>neutralBins){
                this.adaptation+=getNonNeutralSubstitutions()[i];
            }
        }


    }

    /**
     * Runs the Bhatt counting method, neutral rate will be estimated.
     * Prior for the binomial is usually assumed to be beta[1,1].
     * See Raghwani et al (2016)
     * @param binsvalues - site frequency bin ranges
     * @param prior - parametise the beta distribution as prior for the binomial sampler
     * @param needPrior
     * @param which
     */
    public void inferCountsEstimatedNR(double[][] binsvalues,double[] prior, boolean needPrior,boolean[] which){
        this.numBins = binsvalues[0].length;
        double[][] finalmat = new double[6][binsvalues[0].length];
        double[][] totals = new double[3][binsvalues[0].length];
        for(int i=0;i< (int) numBins;i++){
            double[] temp = calculateSiteFrequenciesWithinInterval(binsvalues[0][i], binsvalues[1][i],prior,needPrior);
            finalmat[0][i] = temp[0];   // number silentProb
            finalmat[1][i] = temp[1];	// number replacement

            totals[0][i] = temp[2];    //total count
            totals[1][i] = temp[3];    //total variant
            totals[2][i] = temp[4];    //total invariant
            if(temp[1]!=0){
                finalmat[3][i] = temp[1]/temp[0];	// Silent/replacement ratio
            } else {
                finalmat[3][i] = Double.NaN;
            }
        }
        this.setSilentSubstitutionsCountArray(finalmat[0]);
        this.setReplacementSubstitutionsCountArray(finalmat[1]);
        this.totalCountArray = totals[0];
        this.totalCountArrayInvariantSitesExcluded = totals[1];
        this.totalCountArrayInvariantSitesOnly = totals[2];
        this.setReplacementToSilentRatesRatio(finalmat[3]);
        this.whichBins = which;

        // calcualate neutral ratio. Vector NeutralVec provides boolean true false if a given bin is neutral
        // user has to specify this bin
        double counter = 0;
        double NR=0;
        for(int i=0;i<(int) numBins;i++){
            if(whichBins[i]){
                if(finalmat[0][i]!=0){
                    NR += finalmat[1][i]/finalmat[0][i];
                    counter	++;
                }
            }
        }
        neutralRatio = NR/counter; //average NR

        // number non neutral
        for(int i=0;i<(int) numBins;i++){
            finalmat[5][i] = finalmat[1][i]*(1.0 - ((finalmat[0][i]/finalmat[1][i])*(neutralRatio)));	// number non neutral
            if(finalmat[5][i]<0){
                finalmat[5][i]=0.0;
            }
        }
        this.setNonNeutralSubstitutions(finalmat[5]);
        int flag=0;
        for(int i=0;i< (int) numBins;i++){
            if(whichBins[i]==false && flag==0){
                this.deleteriousLoad+=getNonNeutralSubstitutions()[i];
            }
            if(whichBins[i]){flag=1;}
            if(whichBins[i]==false && flag==1){
                this.adaptation+=getNonNeutralSubstitutions()[i];
            }
        }
    }

    /**
     * The Bhatt probabalistic counting method; neutral rate is fixed. 
     * Prior for the binomial is usually assumed to be beta[1,1].
     * See Raghwani et al (2016)
     * @param binsvalues - site frequency bin ranges
     * @param prior - parametise the beta distribution as prior for the binomial sampler
     * @param needPrior
     * @param which
     * @param NR - the fixed neutral (replacement:silent) ratio for the mid-frequency site class
     */
    public void inferCountsFixedNR(double[][] binsvalues,double[] prior, boolean needPrior,boolean[] which,double NR) throws NullNeutralRatioException{
        /*
         * It should not be possible to run this method unless a <b>sensible</b> neutral ratio is passed
         */
    	if(((Double)NR).equals(null)){
    		throw new NullNeutralRatioException();
    	}
    	this.numBins = binsvalues[0].length;
        double[][] finalmat = new double[6][binsvalues[0].length];
        double[][] totals = new double[3][binsvalues[0].length];
        for(int i=0;i< (int) numBins;i++){
            //System.out.println(prior);
            double[] temp = calculateSiteFrequenciesWithinInterval(binsvalues[0][i], binsvalues[1][i],prior,needPrior);
            finalmat[0][i] = temp[0];   // number silentProb, sigma
            finalmat[1][i] = temp[1];	// number replacement, rho

            totals[0][i] = temp[2];		// total count
            totals[1][i] = temp[3];		// total variant
            totals[2][i] = temp[4];		// total invariant
            if(temp[1]!=0){      // was temp[1] but should be temp[0]
                finalmat[3][i] = temp[1]/temp[0];	// Silent/replacement ratio
            } else {
                finalmat[3][i] = Double.NaN;
            }
        }
        this.setSilentSubstitutionsCountArray(finalmat[0]);
        this.setReplacementSubstitutionsCountArray(finalmat[1]);
        this.totalCountArray = totals[0];
        this.totalCountArrayInvariantSitesExcluded = totals[1];
        this.totalCountArrayInvariantSitesOnly = totals[2];
        this.setReplacementToSilentRatesRatio(finalmat[3]);
        this.whichBins = which;

        neutralRatio = NR; //average

        // number non neutral
        for(int i=0;i<(int) numBins;i++){
            finalmat[5][i] = finalmat[1][i]*(1.0 - ((finalmat[0][i]/finalmat[1][i])*(neutralRatio)));	// number non neutral
            if(finalmat[5][i]<0){
                finalmat[5][i]=0.0;
            }
        }
        this.setNonNeutralSubstitutions(finalmat[5]);
        int flag=0;
        for(int i=0;i< (int) numBins;i++){
            if(whichBins[i]==false && flag==0){
                this.deleteriousLoad+=getNonNeutralSubstitutions()[i];
            }
            if(whichBins[i]){flag=1;}
            if(whichBins[i]==false && flag==1){
                this.adaptation+=getNonNeutralSubstitutions()[i];
            }
        }
    }






    //************************* Utility Subroutines ************************************
 
    /**
     * Finds invalid sites. Invalid sites are sites which (i) have a gap in the main alignment (ii) have a gap in the main alignment
     * @param integerMatrix
     * @param integerArray
     * @return
     */
    public boolean[] findInvalidSites(int[][] integerMatrix, int[] integerArray){
        boolean isSiteBad = false;
        boolean[] badlist = new boolean[integerMatrix[0].length];
        for (int siteIndex = 0; siteIndex< integerMatrix[0].length; siteIndex++){	// for sites
            // flag any sites with gaps or invalid characters
            if(calculateBaseCount(integerMatrix,5,siteIndex)>0){
                isSiteBad = true;  // check sequence alignment
            }
         
            if(integerArray[siteIndex]>4){  // check anscetor does not have invalid sites
                isSiteBad = true; // check ancestral sequence
            }

            // if site flagged for any of the above reasons label as a bad site
            badlist[siteIndex]=isSiteBad;
            // reset flag status
            isSiteBad = false;							
        }
        return badlist;
    }
    
    /**
     * calculates number of bases
     * @param matrix
     * @param base
     * @param site
     * @return
     */
    public double calculateBaseCount(int[][] matrix, int base, int site){
        double count = 0.0;
        for (int i=0; i< matrix.length; i++){
            if (matrix[i][site] == base){
                count++;						// counter
            }
        }
        return count;
    }

    /**
     * Calculates site information.
     * @param site
     * @return
     */
    public SiteInformation calculateSiteInformation(int site){
        SiteInformation siteInformation = new SiteInformation();
        siteInformation.locusIndex=site;
        double numbase = 4;	// number of bases

        SiteObservations[] data = new SiteObservations[(int)numbase];		// create array of teaspoon.adaptation.SiteObservations - an object that stores all info
        for(int i=0;i<numbase;i++){
            data[i] = new SiteObservations();			// initialises the teaspoon.adaptation.SiteObservations objects
        }
        double TotalNumBases=0.0;
        for(int i=0;i<numbase;i++){
            data[i].base=i+1;
            data[i].rawNumObs = TeaspoonMethods.num_of_base(integerMatrix, i+1, site); //
            TotalNumBases+=data[i].rawNumObs;

            if(integerAncestralArray[site]-1 == -1) {
                System.out.println(site+" x");
            }

            data[integerAncestralArray[site]-1].inAncestral=true;	//tests if base is ansestral
        }

        // calculates the number of observed ignoring gaps i.e as a sum of the total number of bases
        for(int i=0;i<numbase;i++){
            data[i].numObs = data[i].rawNumObs/TotalNumBases;
        }
        
        siteInformation.totalNumBases=TotalNumBases;
        // store siteData in SI object
        siteInformation.siteData = data;

        for(int i=0;i<numbase;i++){
            if(siteInformation.siteData[integerAncestralArray[site]-1].numObs!=0.0){
                siteInformation.hasAncestral=true;		// test if site has ansestralbase
            }
            
            if(siteInformation.siteData[i].numObs!=0.0 && data[i].inAncestral==false) {
                siteInformation.numberOfDerived++; // site has no ancestral bases
            }
        }

        /* Determine which of the seven polymorphism types at this position */
        if (siteInformation.numberOfDerived==0 && siteInformation.hasAncestral==true){
            siteInformation.polymorphismCase=1;// invariant
        } else if (siteInformation.numberOfDerived==1 && siteInformation.hasAncestral==false){
            //System.out.println("fixed: "+site);
            siteInformation.polymorphismCase=2;// fixed
        } else if (siteInformation.numberOfDerived==1 && siteInformation.hasAncestral==true){
            siteInformation.polymorphismCase=3;// 1 state derived and ans
        } else if (siteInformation.numberOfDerived==2 && siteInformation.hasAncestral==false){
            siteInformation.polymorphismCase=4;// 2 state derived no ans
        } else if (siteInformation.numberOfDerived==2 && siteInformation.hasAncestral==true){
            siteInformation.polymorphismCase=5;// 2 state derived and ans
        } else if (siteInformation.numberOfDerived==3 && siteInformation.hasAncestral==false){
            siteInformation.polymorphismCase=6;// 3 state derived no ans
        } else if (siteInformation.numberOfDerived==3 && siteInformation.hasAncestral==true){
            siteInformation.polymorphismCase=7;// 3 state derived and ans
        }

        return siteInformation;
    }

    /**
     * calculates codon number
     * @param pos1
     * @param pos2
     * @param pos3
     * @return
     */
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

   /**
    * finds wheter a change is silentProb or replacement
    * @param ansnumber
    * @param number
    * @return double - 1.0 or 0.0
    */
    public double determineSilentOrReplacement(int ansnumber, int number){
        //System.out.println(ansnumber+","+ number);
        if(checkAA(ansnumber) && checkAA(number) ) {
            if (AA[ansnumber].equals(AA[number])) {
                return 1.0;
            } else {
                return 0.0;
            }
        }
        return 0.0;
    }

    /**
     * Is this a valid AA? If not (gap, NNN or unrecognised)
     * @param number
     * @return
     */
    private boolean checkAA(int number) {

        boolean validity = false;
        if(number >=0 && number <=64) {

            validity = true;
        }

        return validity;
    }

    // 
    /**
     * finds which bases are different from ancestral codon, 
     * identifies which bases are different between anscestral codon 
     * and main codon - for use in nei gojobori pathways
     * @param anscodon
     * @param seqcodon
     * @return
     */
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

    /**
     * Utility method, prints array elements as tab-delimited to STOUT
     * @param array
     */
    public static void printDoubleArrayElements(double[] array){
        for(int i=0;i<array.length;i++){
            System.out.print(array[i]+"\t");
        }
        System.out.println();
    }

    public Store createBlocks(int blocksize,int length, int[] Sampling){
        Store MatrixStore = new Store();
        MatrixStore.BlockMatrix = makeSeqBlocks(blocksize,length);
        MatrixStore.AnsBlockMatrix = makeAnsBlocks(blocksize,length);
        MatrixStore.RandomisedBlockMatrix = makeSeqBlocks(blocksize,length);
        MatrixStore.RandomisedAnsBlockMatrix = makeAnsBlocks(blocksize,length);
        MatrixStore.SamplingArray =  Sampling;

        for(int j=0;j<MatrixStore.BlockMatrix[0].length;j++){
            MatrixStore.RandomisedAnsBlockMatrix[j] = MatrixStore.AnsBlockMatrix[MatrixStore.SamplingArray[j]];
            for(int i=0;i<MatrixStore.BlockMatrix.length;i++){
                MatrixStore.RandomisedBlockMatrix[i][j] = MatrixStore.BlockMatrix[i][MatrixStore.SamplingArray[j]];
            }
        }

        MatrixStore.RandomisedIntegerMatrix = new int[MatrixStore.RandomisedBlockMatrix.length][MatrixStore.RandomisedBlockMatrix[0].length*blocksize];
        MatrixStore.RandomisedIntegerAncestral = new int[MatrixStore.RandomisedAnsBlockMatrix.length*blocksize];
        int k=0;
        for(int j=0;j<MatrixStore.RandomisedAnsBlockMatrix.length;j++){
            for(int x=0;x<MatrixStore.RandomisedAnsBlockMatrix[j].sites.length;x++){
                MatrixStore.RandomisedIntegerAncestral[k] =  MatrixStore.RandomisedAnsBlockMatrix[j].sites[x];
                k++;
            }
        }


        for(int i=0;i<MatrixStore.RandomisedBlockMatrix.length;i++){
            k=0;
            for(int j=0;j<MatrixStore.RandomisedBlockMatrix[0].length;j++){
                for(int x=0;x<MatrixStore.RandomisedBlockMatrix[i][j].sites.length;x++){
                    MatrixStore.RandomisedIntegerMatrix[i][k] =  MatrixStore.RandomisedBlockMatrix[i][j].sites[x];
                    k++;
                }
            }
        }

        return MatrixStore;
    }

    public BlockStruct[][] makeSeqBlocks(int blocksize,int length){
        double numblocks = length/blocksize;
        BlockStruct[][] blockmat = new BlockStruct[integerMatrix.length][(int) numblocks];
        int[] temp = new int[blocksize];
        for (int site = 0,x=0; site < length - (blocksize-1); site = site + blocksize,x++) {
            for(int i=0;i<integerMatrix.length;i++){
                int k=0;
                for(int j=site;j<site+blocksize;j++){
                    temp[k] = integerMatrix[i][j];
                    k++;
                }
                blockmat[i][x] = new BlockStruct(blocksize);
                for(int t=0;t<temp.length;t++){
                    blockmat[i][x].sites[t] = temp[t];
                }
            }
        }

        return blockmat;
    }

    public BlockStruct[] makeAnsBlocks(int blocksize,int length){
        double numblocks = length/blocksize;
        BlockStruct[] blockmat = new BlockStruct[(int) numblocks];

        int[] temp = new int[blocksize];
        for (int site = 0,x=0; site < length - (blocksize-1); site = site + blocksize,x++) {
            int k=0;
            for(int j=site;j<site+blocksize;j++){
                temp[k] = integerAncestralArray[j];
                k++;
            }
            blockmat[x] = new BlockStruct(blocksize);
            for(int t=0;t<temp.length;t++){
                blockmat[x].sites[t] = temp[t];
            }
        }
        return blockmat;
    }

    public BlockStruct[][] makeSeqBlocksOverlapping(int blocksize,int length){
        int numblocks = length-blocksize+1;
        BlockStruct[][] blockmat = new BlockStruct[integerMatrix.length][(int) numblocks];
        int[] temp = new int[blocksize];
        for (int site = 0,x=0; site < length - (blocksize-1); site=site+3,x++) {
            for(int i=0;i<integerMatrix.length;i++){
                int k=0;
                for(int j=site;j<site+blocksize;j++){
                    temp[k] = integerMatrix[i][j];
                    k++;
                }
                blockmat[i][x] = new BlockStruct(blocksize);
                for(int t=0;t<temp.length;t++){
                    blockmat[i][x].sites[t] = temp[t];
                }
            }
        }
        return blockmat;
    }

    public BlockStruct[] makeAnsBlocksOverlapping(int blocksize,int length){
        double numblocks = length-blocksize+1;
        BlockStruct[] blockmat = new BlockStruct[(int) numblocks];

        int[] temp = new int[blocksize];
        for (int site = 0,x=0; site < length - (blocksize-1); site = site+3,x++) {
            int k=0;
            for(int j=site;j<site+blocksize;j++){
                temp[k] = integerAncestralArray[j];
                k++;
            }
            blockmat[x] = new BlockStruct(blocksize);
            for(int t=0;t<temp.length;t++){
                blockmat[x].sites[t] = temp[t];
            }
        }
        return blockmat;
    }

	/**
	 * @return the replacementCountArray
	 */
	public double[] getReplacementSubstitutionsCountArray() {
		return replacementCountArray;
	}

	/**
	 * @param newReplacementCountArray the replacementCountArray to set
	 */
	public void setReplacementSubstitutionsCountArray(double[] newReplacementCountArray) {
		replacementCountArray = newReplacementCountArray;
	}

	/**
	 * @return the silentCountArray
	 */
	public double[] getSilentSubstitutionsCountArray() {
		return silentCountArray;
	}

	/**
	 * @param newSilentCountArray the silentCountArray to set
	 */
	public void setSilentSubstitutionsCountArray(double[] newSilentCountArray) {
		silentCountArray = newSilentCountArray;
	}

	/**
	 * @return the replacementSilentRatio
	 */
	public double[] getReplacementToSilentRatesRatio() {
		return replacementToSilentRatio;
	}

	/**
	 * @param newReplacementSilentRatio the replacementSilentRatio to set
	 */
	public void setReplacementToSilentRatesRatio(double[] newReplacementSilentRatio) {
		replacementToSilentRatio = newReplacementSilentRatio;
	}

	/**
	 * @return the nonNeutralSubstitutions
	 */
	public double[] getNonNeutralSubstitutions() {
		return nonNeutralSubstitutions;
	}

	/**
	 * @param newNonNeutralSubstitutions the nonNeutralSubstitutions to set
	 */
	public void setNonNeutralSubstitutions(double[] newNonNeutralSubstitutions) {
		nonNeutralSubstitutions = newNonNeutralSubstitutions;
	}
}
