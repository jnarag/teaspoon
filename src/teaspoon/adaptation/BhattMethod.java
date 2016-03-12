package teaspoon.adaptation;

import java.util.*;
public class BhattMethod {

    double N = 500.0; //set number of replicates
    double number_bins;
    boolean[] WhichBins;
    double[][] bins;
    int NumSample;
    public final  int[][] integer_matrix;
    public final int[] integer_ancestral;
    double neutralratio;
    double neutralbin;
    double[] SilentCountArray;
    double[] ReplacementCountArray;
    double[] TotalCountArray;
    double[] TotalCountArrayNoInvariant;
    double[] TotalCountArrayInvariantOnly;
    double[] ReplacementSilentRatio;
    double[] NonNeutralSubstitutions;
    double Adaptation;
    double DeleteriousLoad;

    public  int[][] codon_matrix;
    public  boolean[] bad_sites_list;
    Methods preprocess = new Methods();


    public final String[] AA =	{"K","N","K","N","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H","P","P","P","P",
            "R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V",
            "X","Y","X","Y","S","S","S","S","X","C","W","C","L","F","L","F","?","-","?" };

    public BhattMethod(){
        throw new RuntimeException("please input the raw integer matrix and the ancestral matrix");
    }

    public BhattMethod(int[][] m,int[] a){
        this.integer_matrix = m;
        this.integer_ancestral = a;
        bad_sites_list = InvalidSites(integer_matrix, integer_ancestral);
    }

    public void setNeutralRatio(double val){
        this.neutralratio=val;
    }

    public ArrayList<Mutation> Tracking(ArrayList<Mutation> FullStore,double[] prior,boolean needprior){
        int c=0;
        for(int site=0;site<integer_matrix[0].length;site++){
            //************************************************************************************************************************
            // dirichlet site frequency spectrum
            SiteInfo Info = SiteInformation(site);
            double numbase = 4;	// number of bases
            double[] observations = new double[(int)numbase];
            if(needprior){
                for(int i=0;i<numbase;i++){
                    if(Info.data[i].inans){
                        prior[i]=1.0;
                    }else{
                        prior[i]=1.0/3.0; // change later
                    }
                }
            }
            for(int i=0;i<numbase;i++){
                Info.data[i].Prior = prior[i];
                observations[i] = Info.data[i].rawNObs+Info.data[i].Prior;	// add observations into an array to find Dirichlet(01+p....0k+p) where p is the prior
            }

            Samplers S = new Samplers(observations);
            double[][] dist = new double[(int)N][4];
            double[] mean = new double[4];
            for(int i=0;i<(int) N;i++){ // Sampler loop
                double[] point = S.Dirichlet();
                dist[i]=point;
                for(int x=0;x<numbase;x++){
                    mean[x]+=point[x];
                }
            }	// end loop
            mean[0]=mean[0]/N;mean[1]=mean[1]/N;mean[2]=mean[2]/N;mean[3]=mean[3]/N;

            double[] A = new double[(int)N];
            double[] C = new double[(int)N];
            double[] G = new double[(int)N];
            double[] T = new double[(int)N];
            for(int s=0;s<N;s++){
                A[s] = dist[s][0];
                C[s] = dist[s][1];
                G[s] = dist[s][2];
                T[s] = dist[s][3];
            }

            Arrays.sort(A);
            Arrays.sort(C);
            Arrays.sort(G);
            Arrays.sort(T);

            for(int s=0;s<N;s++){
                dist[s][0]=A[s];
                dist[s][1]=C[s];
                dist[s][2]=G[s];
                dist[s][3]=T[s];
            }

            double high = N*0.975;double low = N*0.025;

            for(int x=0;x<numbase;x++){
                if(Info.data[x].inans==false){
                    FullStore.get(c).addProb(mean[x]);
                    FullStore.get(c).addH(dist[(int) high][x]);
                    FullStore.get(c).addL(dist[(int) low][x]);
                    c++;
                }
            }
            //************************************************************************************************************************

				
		/*		//************************************************************************************************************************
				// dirichlet site frequency spectrum
				Info = SiteInformation(site+1);
				 numbase = 4;	// number of bases 
				observations = new double[(int)numbase];
				if(needprior){
					for(int i=0;i<numbase;i++){		
						if(Info.data[i].inans){
							prior[i]=1.0;
						}else{
							prior[i]=1.0/3.0; // change later
						}
					}
				}		
				for(int i=0;i<numbase;i++){		
					Info.data[i].Prior = prior[i];
					observations[i] = Info.data[i].rawNObs+Info.data[i].Prior;	// add observations into an array to find Dirichlet(01+p....0k+p) where p is the prior
				}

				S = new teaspoon.adaptation.Samplers(observations);
	
				dist = new double[(int)N][4];
				 mean = new double[4];
				for(int i=0;i<(int) N;i++){ // Sampler loop
					double[] point = S.Dirichlet();
					dist[i]=point;
					for(int x=0;x<numbase;x++){					
						mean[x]+=point[x];
					}	
				}	// end loop
				mean[0]=mean[0]/N;mean[1]=mean[1]/N;mean[2]=mean[2]/N;mean[3]=mean[3]/N;
				M = new teaspoon.adaptation.Mutation[3];
				c=0;
				for(int x=0;x<numbase;x++){
					if(Info.data[x].inans==false){
						M[c].base=x+1;
						M[c].prob=mean[x];	
						int ans = getcodonnumber(integer_ancestral[site],integer_ancestral[site+1],integer_ancestral[site+2]);
						int main = getcodonnumber(integer_ancestral[site],x+1,integer_ancestral[site+2]);
						M[c].id = SilentOrReplacement(ans,main);
						M[c].site=site+1;
						M[c].addProb(mean[x]);
						teaspoon.adaptation.Store.add(M[c]);
						c++;
					}
				}
				//************************************************************************************************************************

				
				//************************************************************************************************************************
				// dirichlet site frequency spectrum
				Info = SiteInformation(site+2);
				 numbase = 4;	// number of bases 
				observations = new double[(int)numbase];
				if(needprior){
					for(int i=0;i<numbase;i++){		
						if(Info.data[i].inans){
							prior[i]=1.0;
						}else{
							prior[i]=1.0/3.0; // change later
						}
					}
				}		
				for(int i=0;i<numbase;i++){		
					Info.data[i].Prior = prior[i];
					observations[i] = Info.data[i].rawNObs+Info.data[i].Prior;	// add observations into an array to find Dirichlet(01+p....0k+p) where p is the prior
				}

				S = new teaspoon.adaptation.Samplers(observations);

				dist = new double[(int)N][4];
				 mean = new double[4];
				for(int i=0;i<(int) N;i++){ // Sampler loop
					double[] point = S.Dirichlet();
					dist[i]=point;
					for(int x=0;x<numbase;x++){					
						mean[x]+=point[x];
					}	
				}	// end loop
				mean[0]=mean[0]/N;mean[1]=mean[1]/N;mean[2]=mean[2]/N;mean[3]=mean[3]/N;
				M = new teaspoon.adaptation.Mutation[3];
				c=0;
				for(int x=0;x<numbase;x++){
					if(Info.data[x].inans==false){
						M[c].base=x+1;
						M[c].prob=mean[x];	
						int ans = getcodonnumber(integer_ancestral[site],integer_ancestral[site+1],integer_ancestral[site+2]);
						int main = getcodonnumber(integer_ancestral[site],integer_ancestral[site+1],x+1);
						M[c].id = SilentOrReplacement(ans,main);
						M[c].site=site+2;
						M[c].addProb(mean[x]);
						teaspoon.adaptation.Store.add(M[c]);
						c++;
					}
				}*/
            //************************************************************************************************************************
        }
        return FullStore;
    }

    public SiteInfo DirichletSiteFreq(double u, double v,int site,double[] prior,boolean needprior) {
        SiteInfo Info = SiteInformation(site);
        double numbase = 4;	// number of bases
        double[] observations = new double[(int)numbase];
        double[] sampler = new double[(int) N];

        //set the prior - otherwise it will be a flat prior?
        if(needprior){
            for(int i=0;i<numbase;i++){
                if(Info.data[i].inans){
                    prior[i]=1.0;       //I think this should be a class variable
                }else{                  // Also I don't understand how this gives you a uniform distribution?
                    prior[i]=1.0/3.0;
                }
            }
        }

        // no. of observations is equal to numbases
        for(int i=0;i<numbase;i++){
            Info.data[i].Prior = prior[i]; // is prior a prior probability or a  parameter of a distirbution?
            observations[i] = Info.data[i].rawNObs + Info.data[i].Prior;	// add observations into an array to find Dirichlet(01+p....0k+p) where p is the prior
        }

        Samplers S = new Samplers(observations);

        double count = 0;
        //bootstrap replicate based on observations (which are essential prior plus observations)
        for(int i=0;i<(int) N;i++){
            double[] point = S.Dirichlet();      // randomly draw the base frequencies from a dirichlet distribution which is based on the observed+prior
            for(int x=0;x<numbase;x++){
                if(Info.data[x].inans==false){
                    sampler[i]+=point[x];    //derived freq from the Dirichlet distribution 500 X
                }
            }
            if(sampler[i]>u && sampler[i]<=v){
                count++;
            }
        }

        Info.Dprob = (count/N);             // average derived frequency   (posterior probability)
        return Info;
    }

    public SiteInfo BetaSiteFreq(double u, double v, int site, double[] prior, boolean needprior) {

        SiteInfo Info = SiteInformation(site);

        double[] observations = new double[2];
        double[] sampler = new double[(int)N];

        //numbase should be a class variable

        if(needprior){

            prior[0] = 1;       //using a uniform prior...
            prior[1] = 1;

        }

        for(int i = 0; i < 4; i++) {
            if(Info.data[i].inans==false) {
                observations[0] += (Info.data[i].rawNObs);
            }
            else{
                observations[1] = Info.data[i].rawNObs ;//+ prior[0];
            }
        }

        observations[0] = observations[0]+prior[0];   //derived frequency
        observations[1] = observations[1]+prior[1];   // ancestral frequency

        Samplers S = new Samplers(observations);

        double count = 0;

        for(int i=0; i < (int)N; i++) {

            double point = S.Beta();
            sampler[i] = point;

            if(sampler[i]> u && sampler[i] < v){
                 count++;
            }
        }

        Info.Dprob = (count/N);



        return Info;
    }

    public SiteInfo DirichletSRFreq(double u, double v,int site,double[] prior,boolean needprior,int pos) {
        SiteInfo Info = SiteInformation(site);
        double numbase = 4;	// number of bases
        double[] observations = new double[(int)numbase];
        double[] sampler = new double[(int) N];

        if(needprior){
            for(int i=0;i<numbase;i++){
                if(Info.data[i].inans){
                    prior[i]=1.0;
                }else{
                    prior[i]=1.0/3.0;
                }
            }
        }
        for(int i=0;i<numbase;i++){
            Info.data[i].Prior = prior[i];
            observations[i] = Info.data[i].rawNObs+Info.data[i].Prior;	// add observations into an array to find Dirichlet(01+p....0k+p) where p is the prior
        }

        Samplers S = new Samplers(observations);
        double[] isSil = new double[4];

        //checking if site is silent or replacement
        for(int i=0;i<numbase;i++){
            if(pos==1){
                int tmp = getcodonnumber(i+1,integer_ancestral[site+1],integer_ancestral[site+2]);
                int tmp2 = getcodonnumber(integer_ancestral[site],integer_ancestral[site+1],integer_ancestral[site+2]);
                isSil[i]=SilentOrReplacement(tmp,tmp2);
            }
            if(pos==2){
                int tmp = getcodonnumber(integer_ancestral[site-1],i+1,integer_ancestral[site+1]);
                int tmp2 = getcodonnumber(integer_ancestral[site-1],integer_ancestral[site],integer_ancestral[site+1]);
                isSil[i]=SilentOrReplacement(tmp,tmp2);
            }
            if(pos==3){
                int tmp = getcodonnumber(integer_ancestral[site-2],integer_ancestral[site-1],i+1);
                int tmp2 = getcodonnumber(integer_ancestral[site-2],integer_ancestral[site-1],integer_ancestral[site]);
                isSil[i]=SilentOrReplacement(tmp,tmp2);
            }
        }

        double[] silsampler = new double[(int) N];
        double[] repsampler = new double[(int) N];
        double silent=0;
        double replacement=0;

        double count = 0;
        for(int i=0;i<(int) N;i++){
            double[] point = S.Dirichlet();
            for(int x=0;x<numbase;x++){
                if(Info.data[x].inans==false){  //site is invariant or fixed  (no ancestral site?)
                    sampler[i]+=point[x];
                    silsampler[i]+=point[x]*isSil[x];
                    repsampler[i]+=point[x]*(1.0-isSil[x]);
                }
            }
            if(sampler[i]>u && sampler[i]<=v){
                count++;
            }
            if(silsampler[i]>u && silsampler[i]<=v){
                silent++;
            }
            if(repsampler[i]>u && repsampler[i]<=v){
                replacement++;
            }

        }

        Info.Dprob = (count/N);
        Info.Sil = (silent/N);
        Info.Rep = (replacement/N);
        return Info;
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
            double[] tmp = NGpathway(ancestralbases,mainbases);  // per site
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

    public double[] NGpathway(int[] a, int[] b){
    //a = ancestral b=main
        // 1 silent, 0 replacement, 2 invariant
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
                if(diff[0]){identity[0]=SilentOrReplacement(AnscodonNumber,codonNumber);}
                if(diff[1]){identity[1]=SilentOrReplacement(AnscodonNumber,codonNumber);}
                if(diff[2]){identity[2]=SilentOrReplacement(AnscodonNumber,codonNumber);}
            }
            // 2 sites different ****************************************************************************************
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

            // 3 sites different ****************************************************************************************
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
        }

        //System.out.println(identity);
        return identity;
    }

    public double[][] NGpossible(){
        int[] mainbases = new int[3];
        double L = (double) integer_matrix.length;
        double N = (double) integer_matrix[0].length;
        double[][] identity = new double[3][2];
        double[][] FinalIdentity = new double[(int)N][2];
        for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
            identity = new double[3][2]; // set new identity vector
            if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) {  // check for bad sites
                for(int sequence=0;sequence<integer_matrix.length;sequence++){		// loop through sequences

                    mainbases[0]=integer_matrix[sequence][site]; //main bases
                    mainbases[1]=integer_matrix[sequence][site+1];
                    mainbases[2]=integer_matrix[sequence][site+2];
                    double sil = 0; //silent count
                    double rep=0; //rep count
                    //********************************************************************************************************************************************
                    //pos 1
                    sil=0;
                    for(int i=1;i<5;i++){  // loop through all bases
                        if(i != mainbases[0]){ // make sure actual base is not counted
                            int[] actualCodon = {mainbases[0],mainbases[1],mainbases[2]};
                            int[] IntermediateCodon = {i,mainbases[1],mainbases[2]};
                            double[] SR = NGpathway(actualCodon,IntermediateCodon);
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
                            double[] SR = NGpathway(actualCodon,IntermediateCodon);
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
                            double[] SR = NGpathway(actualCodon,IntermediateCodon);
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

    // calculates a silent and replacement site frequency between ranges u and v
    public double[] SiteFreq(double u,double v, double[] prior, boolean needPrior){
        double[][] temp = new double[3][2];
        double rho1 = 0.0;
        double sigma1 = 0.0;
        double[] finalans = new double[5];
        double total1 = 0.0;
        double rho2 = 0.0;
        double sigma2 = 0.0;
        double total2 = 0.0;

        ArrayList<SiteInfo> Inv = new ArrayList<SiteInfo>();


        double old_rho = 0;
        for(int site=0,codon=0; site<integer_matrix[0].length-2;site+=3,codon++){
            //System.out.println(site+","+codon);


            if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) {
                temp = NGmethod(site,codon);
                // find site freq for pos1
                SiteInfo inf = BetaSiteFreq(u,v,site,prior,needPrior);

                double raw_rho = 0;

                //System.out.println(site+","+inf.Dprob);
                //teaspoon.adaptation.SiteInfo inf = DirichletSRFreq(u,v,site,prior,needPrior, 1);
                if(inf.Case==1){
                    Inv.add(inf);
                } else {
                    total1 += inf.Dprob;
                    sigma1 += inf.Dprob*temp[0][0];
                    rho1 += inf.Dprob*temp[0][1];
                    raw_rho += temp[0][1];

                }

                // find site freq for pos2
                SiteInfo inf2 = BetaSiteFreq(u,v,site+1,prior,needPrior);
                //System.out.println(site+1);
                //teaspoon.adaptation.SiteInfo inf2 = DirichletSRFreq(u,v,site+1,prior,needPrior, 2);
                if(inf2.Case==1){
                    Inv.add(inf2);
                } else {
                    total1 += inf2.Dprob;
                    sigma1 += inf2.Dprob*temp[1][0];
                    rho1 += inf2.Dprob*temp[1][1];
                    raw_rho += temp[1][1];


                }

                // find site freq for pos3
                SiteInfo inf3 = BetaSiteFreq(u,v,site+2,prior,needPrior);
                //teaspoon.adaptation.SiteInfo inf3 = DirichletSRFreq(u,v,site+2,prior,needPrior, 3);
                if(inf3.Case==1){
                    Inv.add(inf3);
                } else {
                    total1 += inf3.Dprob;
                    sigma1 += inf3.Dprob*temp[2][0];
                    rho1 += inf3.Dprob*temp[2][1];
                    raw_rho += temp[2][1];


                }
            }
        }

         //System.out.println(highFreqRSites);
        double S = sigma1/total1;   // proportion of variant sites that are silent
        double R = rho1/total1;     // proportion of variant sites that are replacement

        Iterator<SiteInfo> It =  Inv.iterator();

        // infers the frequencies for invariant sites probabilistically
        while(It.hasNext()){
            SiteInfo Element = It.next();
            total2 += Element.Dprob;     // posterior frequency of observed invariant site
            sigma2 += Element.Dprob*S;   // posterior frequency of observed invariant site being silent
            rho2 += Element.Dprob*R;	 // posterior frequency of observed invariant site being replacement
        }

        finalans[0] = sigma1+sigma2;  // silent count
        finalans[1] = rho1+rho2;      // replacement count
        finalans[2] = total1+total2;  // total count
        finalans[3] = total1;         // total variant
        finalans[4] = total2;         // total invariant

        return finalans;

    }

    // calculates a silent and replacement site frequency between ranges u and v
    public double[] SiteFreqNG(double u,double v, double[] prior, boolean needPrior){
        double[][] temp = new double[3][2];
        double rho1 = 0.0;
        double sigma1 = 0.0;
        double[] finalans = new double[5];
        double total1 = 0.0;
        double rhoP = 0.0;
        double sigmaP = 0.0;

        double[][] NGpossible = NGpossible();

        ArrayList<SiteInfo> Inv = new ArrayList<SiteInfo>();
        for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
            if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) {
                temp = NGmethod(site,codon);
                //inf - sitefreq for a specific nucleotide in the alignment

                //find site freq for pos1
                SiteInfo inf = BetaSiteFreq(u,v,site,prior,needPrior);
                if(inf.Case==1){
                    Inv.add(inf);
                } else {
                    total1 += inf.Dprob;
                    sigma1 += inf.Dprob*temp[0][0];
                    rho1 += inf.Dprob*temp[0][1];
                    sigmaP+=NGpossible[site][0];
                    rhoP+=NGpossible[site][1];
                }

                // find site freq for pos2
                SiteInfo inf2 = BetaSiteFreq(u,v,site+1,prior,needPrior);
                if(inf2.Case==1){
                    Inv.add(inf2);
                } else {
                    total1 += inf2.Dprob;
                    sigma1 += inf2.Dprob*temp[1][0];
                    rho1 += inf2.Dprob*temp[1][1];
                    sigmaP+=NGpossible[site+1][0];
                    rhoP+=NGpossible[site+1][1];
                }

                // find site freq for pos3
                SiteInfo inf3 = BetaSiteFreq(u,v,site+2,prior,needPrior);
                if(inf3.Case==1){
                    Inv.add(inf3);
                } else {
                    total1 += inf3.Dprob;
                    sigma1 += inf3.Dprob*temp[2][0];
                    rho1 +=  inf3.Dprob*temp[2][1];
                    sigmaP+=NGpossible[site+2][0];
                    rhoP+=NGpossible[site+2][1];
                }
            }
        }

        //So why is this is section different from the equivalent code in SiteFreq(...) method?
        //This does not contain the raw counts and has less entries in the final matrix

        double S = sigma1/sigmaP;
        double R = rho1/sigmaP;
        double Dn = (-3.0/4.0)*Math.log(1.0-(R*4.0/3.0));
        double Ds = (-3.0/4.0)*Math.log(1.0-(S*4.0/3.0));


        //finalans length = 5, so what are the values for 4th and 5th element
        finalans[0] = Ds; finalans[1]=Dn; finalans[2] = total1;

//        finalans[0] = sigma1+sigmaP;
//        finalans[1] = rho1+rhoP;
//        finalans[2] = total1;//+totalP;
//        finalans[3] = total1;
//        finalans[4] = total2;

        return finalans;

    }

    public void MethodNG(double[][] binsvalues,double[] prior, boolean needPrior,boolean[] which){
        this.number_bins = binsvalues[0].length;
        double[][] finalmat = new double[6][binsvalues[0].length];
        double[][] totals = new double[3][binsvalues[0].length];
        for(int i=0;i< (int) number_bins;i++){
            double[] temp = SiteFreqNG(binsvalues[0][i], binsvalues[1][i],prior,needPrior);
            finalmat[0][i] = temp[0];   // number silent  - Ds
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
        this.SilentCountArray = finalmat[0];
        this.ReplacementCountArray = finalmat[1];
        this.TotalCountArray = totals[0];
        this.TotalCountArrayNoInvariant = totals[1];
        this.TotalCountArrayInvariantOnly = totals[2];
        this.ReplacementSilentRatio = finalmat[3];
        this.WhichBins = which;
        // calcualte neutral ratio. Vector NeutralVec provides boolean true false if a given bin is neutral
        // user has to specify this bin
        double counter = 0;
        double NR=0;
        for(int i=0;i<(int) number_bins;i++){
            if(WhichBins[i]){
                if(finalmat[0][i]!=0){
                    NR += finalmat[1][i]/finalmat[0][i];
                    counter	++;
                }
            }
        }
        neutralratio = NR/counter; //average

        // number non neutral
        for(int i=0;i<(int) number_bins;i++){
            finalmat[5][i] = finalmat[1][i]*(1.0 - ((finalmat[0][i]/finalmat[1][i])*(neutralratio)));	// number non neutral
            if(finalmat[5][i]<0){
                finalmat[5][i]=0.0;
            }
        }
        this.NonNeutralSubstitutions = finalmat[5];
        int flag=0;
        for(int i=0;i< (int) number_bins;i++){
            if(WhichBins[i]==false && flag==0){
                this.DeleteriousLoad+=NonNeutralSubstitutions[i];
            }
            if(WhichBins[i]){flag=1;}
            if(WhichBins[i]==false && flag==1){
                this.Adaptation+=NonNeutralSubstitutions[i];
            }
        }


    }

    public void Method(double[][] binsvalues,double[] prior, boolean needPrior){
        this.number_bins = binsvalues[0].length;
        double[][] finalmat = new double[6][binsvalues[0].length];
        for(int i=0;i< (int) number_bins;i++){
            double[] temp = SiteFreq(binsvalues[0][i], binsvalues[1][i],prior,needPrior);
            finalmat[0][i] = temp[0];   // number silent
            finalmat[1][i] = temp[1];	// number replacement
            finalmat[2][i] = temp[2];	// total number

            if(temp[1]!=0){
                finalmat[3][i] = temp[1]/temp[0];	// Silent/replacement ratio
            } else {
                finalmat[3][i] = Double.NaN;
            }
        }
        this.SilentCountArray = finalmat[0];
        this.ReplacementCountArray = finalmat[1];
        this.TotalCountArray = finalmat[2];
        this.ReplacementSilentRatio = finalmat[3];
        this.neutralratio = ReplacementSilentRatio[0];
        this.neutralbin = 0;

        for(int i=0;i< (int) number_bins;i++){
            if(ReplacementSilentRatio[i]<neutralratio){
                this.neutralratio=ReplacementSilentRatio[i];
                this.neutralbin = i;
            }
        }



        // number non neutral
        for(int i=0;i<(int) number_bins;i++){
            finalmat[5][i] = finalmat[1][i]*(1.0 - ((finalmat[0][i]/finalmat[1][i])*(neutralratio)));	// number non neutral
            if(finalmat[5][i]<0){
                finalmat[5][i]=0.0;
            }
        }
        this.NonNeutralSubstitutions = finalmat[5];
        for(int i=0;i< (int) number_bins;i++){
            if(i<neutralbin){
                this.DeleteriousLoad+=NonNeutralSubstitutions[i];
            }
            if(i>neutralbin){
                this.Adaptation+=NonNeutralSubstitutions[i];
            }
        }


    }

    public void Method(double[][] binsvalues,double[] prior, boolean needPrior,boolean[] which){
        this.number_bins = binsvalues[0].length;
        double[][] finalmat = new double[6][binsvalues[0].length];
        double[][] totals = new double[3][binsvalues[0].length];
        for(int i=0;i< (int) number_bins;i++){
            double[] temp = SiteFreq(binsvalues[0][i], binsvalues[1][i],prior,needPrior);
            finalmat[0][i] = temp[0];   // number silent
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
        this.SilentCountArray = finalmat[0];
        this.ReplacementCountArray = finalmat[1];
        this.TotalCountArray = totals[0];
        this.TotalCountArrayNoInvariant = totals[1];
        this.TotalCountArrayInvariantOnly = totals[2];
        this.ReplacementSilentRatio = finalmat[3];
        this.WhichBins = which;

        // calcualate neutral ratio. Vector NeutralVec provides boolean true false if a given bin is neutral
        // user has to specify this bin
        double counter = 0;
        double NR=0;
        for(int i=0;i<(int) number_bins;i++){
            if(WhichBins[i]){
                if(finalmat[0][i]!=0){
                    NR += finalmat[1][i]/finalmat[0][i];
                    counter	++;
                }
            }
        }
        neutralratio = NR/counter; //average

        // number non neutral
        for(int i=0;i<(int) number_bins;i++){
            finalmat[5][i] = finalmat[1][i]*(1.0 - ((finalmat[0][i]/finalmat[1][i])*(neutralratio)));	// number non neutral
            if(finalmat[5][i]<0){
                finalmat[5][i]=0.0;
            }
        }
        this.NonNeutralSubstitutions = finalmat[5];
        int flag=0;
        for(int i=0;i< (int) number_bins;i++){
            if(WhichBins[i]==false && flag==0){
                this.DeleteriousLoad+=NonNeutralSubstitutions[i];
            }
            if(WhichBins[i]){flag=1;}
            if(WhichBins[i]==false && flag==1){
                this.Adaptation+=NonNeutralSubstitutions[i];
            }
        }


    }

    public void Method(double[][] binsvalues,double[] prior, boolean needPrior,boolean[] which,double NR){
        this.number_bins = binsvalues[0].length;
        double[][] finalmat = new double[6][binsvalues[0].length];
        double[][] totals = new double[3][binsvalues[0].length];
        for(int i=0;i< (int) number_bins;i++){
            //System.out.println(prior);
            double[] temp = SiteFreq(binsvalues[0][i], binsvalues[1][i],prior,needPrior);
            finalmat[0][i] = temp[0];   // number silent
            finalmat[1][i] = temp[1];	// number replacement

            totals[0][i] = temp[2];
            totals[1][i] = temp[3];
            totals[2][i] = temp[4];
            if(temp[1]!=0){      // was temp[1] but should be temp[0]
                finalmat[3][i] = temp[1]/temp[0];	// Silent/replacement ratio
            } else {
                finalmat[3][i] = Double.NaN;
            }
        }
        this.SilentCountArray = finalmat[0];
        this.ReplacementCountArray = finalmat[1];
        this.TotalCountArray = totals[0];
        this.TotalCountArrayNoInvariant = totals[1];
        this.TotalCountArrayInvariantOnly = totals[2];
        this.ReplacementSilentRatio = finalmat[3];
        this.WhichBins = which;

        neutralratio = NR; //average

        // number non neutral
        for(int i=0;i<(int) number_bins;i++){
            finalmat[5][i] = finalmat[1][i]*(1.0 - ((finalmat[0][i]/finalmat[1][i])*(neutralratio)));	// number non neutral
            if(finalmat[5][i]<0){
                finalmat[5][i]=0.0;
            }
        }
        this.NonNeutralSubstitutions = finalmat[5];
        int flag=0;
        for(int i=0;i< (int) number_bins;i++){
            if(WhichBins[i]==false && flag==0){
                this.DeleteriousLoad+=NonNeutralSubstitutions[i];
            }
            if(WhichBins[i]){flag=1;}
            if(WhichBins[i]==false && flag==1){
                this.Adaptation+=NonNeutralSubstitutions[i];
            }
        }


    }






    //************************* Utility Subroutines ************************************
    // invalid sites are sites which (i) have a gap in the main alignment (ii) have a gap in the main alignment
    public boolean[] InvalidSites(int[][] integer_matrix, int[] integer_array){
        boolean flag = false;
        boolean[] badlist = new boolean[integer_matrix[0].length];
        for (int i = 0; i< integer_matrix[0].length; i++){	// for sites
            // flag any sites with gaps or invalid characters
            if(num_of_base(integer_matrix,5,i)>0){
                flag=true;  // check sequence alignment
            }
            if(integer_array[i]>4){  // check anscetor does not have invalid sites
                flag =true; // check ancestral sequence
            }
            //			if site flagged for any of the above reasons label as a bad site
            badlist[i]=flag;
            flag = false;							// reset flag status
        }
        return badlist;
    }
    //calculates number of bases
    public double num_of_base(int[][] matrix, int base, int site){
        double count = 0.0;
        for (int i=0; i< matrix.length; i++){
            if (matrix[i][site] == base){
                count++;						// counter
            }
        }
        return count;
    }

    public SiteInfo SiteInformation(int site){
        SiteInfo SI = new SiteInfo();
        SI.locus=site;
        double numbase = 4;	// number of bases

        Obs[] data = new Obs[(int)numbase];		// create array of teaspoon.adaptation.Obs - an object that stores all info
        for(int i=0;i<numbase;i++){
            data[i] = new Obs();			// initialises the teaspoon.adaptation.Obs objects
        }
        double TotalNumBases=0.0;
        for(int i=0;i<numbase;i++){
            data[i].base=i+1;
            data[i].rawNObs = preprocess.num_of_base(integer_matrix, i+1, site); //
            TotalNumBases+=data[i].rawNObs;

            if(integer_ancestral[site]-1 == -1) {

                System.out.println("x");
            }
            data[integer_ancestral[site]-1].inans=true;	//tests if base is ansestral
        }
        // calculates the number of observed ignoring gaps i.e as a sum of the total number of bases
        for(int i=0;i<numbase;i++){
            data[i].NObs = data[i].rawNObs/TotalNumBases;
        }
        SI.totalNumBases=TotalNumBases;
        // store data in SI object
        SI.data = data;

        for(int i=0;i<numbase;i++){
            if(SI.data[integer_ancestral[site]-1].NObs!=0.0){
                SI.hasans=true;		// test if site has ansestralbase
            }
            if(SI.data[i].NObs!=0.0 && data[i].inans==false) {
                SI.Numderived++; // site has no ancestral bases
            }
        }

        if (SI.Numderived==0 && SI.hasans==true){
            SI.Case=1;// invariant
        }
        else if (SI.Numderived==1 && SI.hasans==false){
            //System.out.println("fixed: "+site);
            SI.Case=2;// fixed
        }
        else if (SI.Numderived==1 && SI.hasans==true){
            SI.Case=3;// 1 state derived and ans
        }
        else if (SI.Numderived==2 && SI.hasans==false){
            SI.Case=4;// 2 state derived no ans
        }
        else if (SI.Numderived==2 && SI.hasans==true){
            SI.Case=5;// 2 state derived and ans
        }
        else if (SI.Numderived==3 && SI.hasans==false){
            SI.Case=6;// 3 state derived no ans
        }
        else if (SI.Numderived==3 && SI.hasans==true){
            SI.Case=7;// 3 state derived and ans
        }

        return SI;
    }

    // calculates codon number
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

    // finds wheter a change is silent or replacement
    public double SilentOrReplacement(int ansnumber, int number){
        double identity = 0.0;
        //System.out.println(ansnumber+","+ number);

        if(checkAA(ansnumber) && checkAA(number) ) {
            if (AA[ansnumber].equals(AA[number])) {
                identity = 1.0;
            } else {
                identity = 0.0;
            }
        }
        return identity;
    }

    private boolean checkAA(int number) {

        boolean validity = false;
        if(number >=0 && number <=64) {

            validity = true;
        }

        return validity;
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

    public void print(double[] array){
        for(int i=0;i<array.length;i++){
            System.out.print(array[i]+"\t");
        }
        System.out.println();
    }

    public Store CreateBlocks(int blocksize,int length, int[] Sampling){
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
        BlockStruct[][] blockmat = new BlockStruct[integer_matrix.length][(int) numblocks];
        int[] temp = new int[blocksize];
        for (int site = 0,x=0; site < length - (blocksize-1); site = site + blocksize,x++) {
            for(int i=0;i<integer_matrix.length;i++){
                int k=0;
                for(int j=site;j<site+blocksize;j++){
                    temp[k] = integer_matrix[i][j];
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
                temp[k] = integer_ancestral[j];
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
        BlockStruct[][] blockmat = new BlockStruct[integer_matrix.length][(int) numblocks];
        int[] temp = new int[blocksize];
        for (int site = 0,x=0; site < length - (blocksize-1); site=site+3,x++) {
            for(int i=0;i<integer_matrix.length;i++){
                int k=0;
                for(int j=site;j<site+blocksize;j++){
                    temp[k] = integer_matrix[i][j];
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
                temp[k] = integer_ancestral[j];
                k++;
            }
            blockmat[x] = new BlockStruct(blocksize);
            for(int t=0;t<temp.length;t++){
                blockmat[x].sites[t] = temp[t];
            }
        }
        return blockmat;
    }

}
