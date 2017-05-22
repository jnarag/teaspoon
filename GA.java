package teaspoon;
import java.util.Arrays;
import java.util.Random;

import flanagan.analysis.Regression;
import teaspoon.adaptation.WilliamsonOutputClass;

public class GA {
	final int tolerance = 100;
	final int minZeros = 116; //20% of length
	final int maxZeros = 450; //80% of length
	final int SeqLength = 566;
	final int InitOnes = 280;
	final int popsize = 100;
	final int mutationrate = 15;
	final int crossoverlength = 280;
	final int proportion = 85; // cut of for selecting individuals top 20 out of 100
	int[][] population;
	double neutralratio;
	Random generator = new Random();
	public GA(){
		population = new int[popsize][SeqLength];
		System.out.println("Number of runs: "+tolerance+" population size: "+ popsize+" mutation rate: " +mutationrate);
	}


	public int[] Initialise(){
		int[] randomised = new int[SeqLength/3];
		for(int i=0;i<InitOnes;i++){
			Random generator = new Random();
			int randint = generator.nextInt(randomised.length-1);	
			if(randomised[randint]!=1){
				randomised[randint]=1;
			}
		}
		return randomised;
	}

	public void InitialisePopulation(){	
		for(int i=0;i<popsize;i++){
			population[i]=Initialise();
		}
	}

	public int[] Mutation(int[] partition){
		int count = CheckZeros(partition);
		if(count<minZeros){
			partition=AddOne(partition);
		} else if(count>maxZeros){
			partition=AddZero(partition);
		} else {
			partition=Swap(partition);
		}
		return partition;
	}

	public void createNewPopulation(int[][] selectedIndv){
		int[][] Newpopulation = new int[popsize][SeqLength];
		for(int i=0;i<popsize;i=i+2){
			
			// cross over
			int randint = generator.nextInt(selectedIndv.length-1);	
			int randint2 = generator.nextInt(selectedIndv.length-1);	
			int[] tmp = new int[selectedIndv[0].length];
			int[] tmp2 = new int[selectedIndv[0].length];
			for(int x=0;x<selectedIndv[0].length;x++){
				if(x<crossoverlength){
					tmp[x] = selectedIndv[randint][x];
					tmp2[x] = selectedIndv[randint2][x];
				}else {
					tmp[x] = selectedIndv[randint2][x];
					tmp2[x] = selectedIndv[randint][x];
				}
			}
			// mutation
			Newpopulation[i]=Mutation(tmp);
			Newpopulation[i+1]=Mutation(tmp2);
		}
		this.population=Newpopulation;

		
	}
	
	public int[] Swap(int[] partition){
		for(int m=0;m<mutationrate;m++){
			Random generator = new Random();
			int randint = generator.nextInt(partition.length-1);	
			int randint2 = generator.nextInt(partition.length-1);	
			partition[randint]=1;
			partition[randint2]=0;		
		}
		return partition;
	}

	public int[] AddOne(int[] partition){
		for(int m=0;m<mutationrate;m++){
			Random generator = new Random();
			int randint = generator.nextInt(partition.length-1);	
			int randint2 = generator.nextInt(partition.length-1);	
			partition[randint]=1;
			partition[randint2]=1;		
		}
		return partition;
	}
	public int[] AddZero(int[] partition){
		for(int m=0;m<mutationrate;m++){
			Random generator = new Random();
			int randint = generator.nextInt(partition.length-1);	
			int randint2 = generator.nextInt(partition.length-1);	
			partition[randint]=1;
			partition[randint2]=1;		
		}
		return partition;
	}

	public int CheckZeros(int[] partition){
		int count=0;
		for(int i=0;i<partition.length;i++){
			if(partition[i]==0){
				count++;
			}
		}
		return count;
	}
	public int CheckOnes(int[] partition){
		int count=0;
		for(int i=0;i<partition.length;i++){
			if(partition[i]==1){
				count++;
			}
		}
		return count;
	}

	public double Calculate(int[] partition){
		FluAnalysis2 F= new FluAnalysis2();
		String Name = "H1";
		String Num = "#4"; //type
		String Geno = "B1"; //type
		int boot = 1;
		// create types
		String[] type = {Num+Geno,Num+Geno+"a",Num+Geno+"b",Num+Geno+"c",Num+Geno+"d",Num+Geno+"e",Num+Geno+"f",Num+Geno+"g",Num+Geno+"h"};

		String[] dates ={"1984","1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //ea


		WilliamsonOutputClass woc2 = F.RunnerBS("/Users/sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,partition,0);
		double S=0; double R=0;

		for(int i=0;i<woc2.Mid_S.size();i++){
			S+=woc2.Mid_S.get(i);
			R+=woc2.Mid_R.get(i);
		}
		neutralratio=R/S;
		WilliamsonOutputClass woc = F.RunnerBS("/Users/sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,partition,1);


		double[] y=new double [woc.TotalA.size()];
		double[] x=new double [woc.TotalA.size()];
		double initdate = woc.dateavg.get(0);
		for(int i=0;i<woc.TotalA.size();i++){
			y[i]=woc.TotalA.get(i);
			x[i]=woc.dateavg.get(i)-initdate;
		}
		Regression Reg = new Regression(x,y);
		Reg.polynomial(1,0);
		double[] est = Reg.getBestEstimates();
		return est[0];
	}

	


	public int[] Maximise(){
		// initialise
		System.out.println("Creating initial population");
		InitialisePopulation();
		double[] fitness = new double[popsize];
		System.out.println("Calculating fitness of init population");
		for(int i=0;i<population.length;i++){
			fitness[i]=Calculate(population[i]);
		}
		Arrays.sort(fitness);
		System.out.println("Taking top 20% of individuals of init population");
		System.out.print("These are the values \t");
		for(int i=proportion;i<popsize;i++){
			System.out.print(fitness[i]+"\t");
		}
		System.out.println();
		
		int[][] selectedindividuals = new int[popsize-proportion][SeqLength];
		int k=0;
		for(int i=proportion;i<popsize;i++){
			selectedindividuals[k] = population[i];
			k++;
		}
		
		// main loop
		for(int x=0;x<tolerance;x++){
			createNewPopulation(selectedindividuals);
			fitness = new double[popsize];
			for(int i=0;i<population.length;i++){
				fitness[i]=Calculate(population[i]);
			}
			Arrays.sort(fitness);
			System.out.print("These are the fitness values \t");
			for(int i=proportion;i<popsize;i++){
				System.out.print(fitness[i]+"\t");
			}
			System.out.println();
			
			selectedindividuals = new int[popsize-proportion][SeqLength];
			k=0;
			for(int i=proportion;i<popsize;i++){
				selectedindividuals[k] = population[i];
				k++;
			}
		}
		
		System.out.println("Maximisation finished...... \n\n Now I will print the best selected individuals\n\n");
		for(int i=0;i<selectedindividuals.length;i++){
			for(int j=0;j<selectedindividuals[0].length;j++){
				System.out.print(selectedindividuals[i][j]);
			}
			System.out.println();
		}
		

		return selectedindividuals[k-1];
	}



}
