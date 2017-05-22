package teaspoon;

import java.util.ArrayList;

import flanagan.analysis.Regression;

public class Predict {

	ArrayList<Double> xdata;
	ArrayList<Double> ydata;
	ArrayList<double[]> boots;



	public Predict(ArrayList<Double> x, ArrayList<Double> y,ArrayList<double[]> boot){
		this.xdata=x;
		this.ydata=y;
		this.boots=boot;
	}
	
	public Predict(ArrayList<Double> x, ArrayList<Double> y){
		this.xdata=x;
		this.ydata=y;
	}
	

	

	public double[] FitPoly(int or){
		int order=or;
		double[] x = new double[xdata.size()];
		double[] y = new double[ydata.size()];   
		double initdate = xdata.get(0);
		for(int i=0;i<xdata.size();i++){
			x[i]=xdata.get(i)-initdate;
			y[i]=ydata.get(i);
		}
		Regression reg = new Regression(x, y);
		reg.polynomial(order,0);
		double[] est = reg.getBestEstimates();
		double[] Beta = new double[order];
		for(int i=0;i<est.length;i++){
			Beta[i]=est[i];
		}
		return Beta;
	}
	
	public double[] bestModel(){
		
		double[] x = new double[xdata.size()];                  
		for(int i=0;i<xdata.size();i++){
			x[i]=xdata.get(i)-xdata.get(0);
		}
		double B[][] = new double[boots.get(0).length][boots.size()];
		double[] modelList = new double[boots.get(0).length];
		for(int i=0;i<boots.get(0).length;i++){
			for(int j=0;j<boots.size();j++){
				B[i][j] = boots.get(j)[i];
			}
		}
		for(int i=0;i<B.length;i++){
			double[] ystar = B[i];
			modelList[i]=WhichModel(x,ystar);		
		}
		double lin=0;double quad=0;double cube=0;
		
		for(double i: modelList){
			if(i==1.0){
				lin++;
			} else if(i==2.0){
				quad++;
			} else if(i==3.0){
				cube++;
			}
		}
		double[] p = {lin/modelList.length,quad/modelList.length,cube/modelList.length};
		return p;
	}
	
	public double[] bestModelAvg(){
		double[] x = new double[xdata.size()];                  
		for(int i=0;i<xdata.size();i++){
			x[i]=xdata.get(i)-xdata.get(0);
		}
		double B[][] = new double[boots.get(0).length][boots.size()];
		for(int i=0;i<boots.get(0).length;i++){
			for(int j=0;j<boots.size();j++){
				B[i][j] = boots.get(j)[i];
			}
		}
		double[] mses = new double[3];
		for(int i=0;i<B.length;i++){
			double[] ystar = B[i];
			double[] tmp=MSEval(x,ystar);	
			mses[0]+=tmp[0];mses[1]+=tmp[1];mses[2]+=tmp[2];	
		}
		
		for(int i=0;i<mses.length;i++){
			mses[i]=mses[i]/B.length;
		}
		
		return mses;
	}

	
	public double[][] bestModelDist(){
		double[] x = new double[xdata.size()];                  
		for(int i=0;i<xdata.size();i++){
			x[i]=xdata.get(i)-xdata.get(0);
		}
		double B[][] = new double[boots.get(0).length][boots.size()];
		for(int i=0;i<boots.get(0).length;i++){
			for(int j=0;j<boots.size();j++){
				B[i][j] = boots.get(j)[i];
			}
		}
		double[][] mses = new double[B.length][3];
		for(int i=0;i<B.length;i++){
			double[] ystar = B[i];
			double[] tmp=MSEval(x,ystar);	
			mses[i][0]=tmp[0];
			mses[i][1]=tmp[1];
			mses[i][2]=tmp[2];	
		}
		return mses;
	}
	
	public double[] FitLinear(){
		double[] x = new double[xdata.size()];                  
		for(int i=0;i<xdata.size();i++){
			x[i]=xdata.get(i)-xdata.get(0);
		}
		double[] ests = new double[boots.get(0).length];
		double B[][] = new double[boots.get(0).length][boots.size()];
		for(int i=0;i<boots.get(0).length;i++){
			for(int j=0;j<boots.size();j++){
				B[i][j] = boots.get(j)[i];
			}
		}
		
		
		for(int i=0;i<B.length;i++){
			double[] ystar = B[i];
			Regression reg = new Regression(x, ystar);
			reg.polynomial(1,0);
			double[] est = reg.getBestEstimates();
			ests[i]=est[0];
		}
		return ests;
	}
	
	
	public double[] FitPoly(double[] x, double[] y,int or){
		int order=or;
		Regression reg = new Regression(x, y);
		reg.polynomial(order,0);
		double[] est = reg.getBestEstimates();
		double[] Beta = new double[order];
		for(int i=0;i<est.length;i++){
			Beta[i]=est[i];
		}
		return Beta;
	}

	public double CrossValidation(int or){
		// cast to double
		int order=or;
		double[] x = new double[xdata.size()];
		double[] y = new double[ydata.size()];      
		double initdate = xdata.get(0);
		for(int i=0;i<xdata.size();i++){
			x[i]=xdata.get(i)-initdate;
			y[i]=ydata.get(i);
		}
		double MeanSqErr = 0;double r=0;
		for(int i=order+1;i<x.length;i++){
			double[] xtrain = new double[i];
			double[] ytrain = new double[i];

			for(int z=0;z<i;z++){
				xtrain[z] = x[z];
				ytrain[z] = y[z];
			}
			
			double[] est = FitPoly(xtrain,ytrain,order);

			double ypredict=0;
			if(order==1){
				ypredict = est[0]*x[i];
			} else if(order==2){
				ypredict = est[0]*x[i] + est[1]*x[i]*x[i] ;
			}else if(order==3){
				ypredict = est[0]*x[i] + est[1]*x[i]*x[i] + est[2]*x[i]*x[i]*x[i] ;
			}
			MeanSqErr += ((y[i]-ypredict)*(y[i]-ypredict));
			r++;
		}
		MeanSqErr = MeanSqErr/r;
		return MeanSqErr;		
	}

	public double CrossValidation(double[] x, double[] y,int or){
		double initdate = x[0];
		for(int i=0;i<xdata.size();i++){
			x[i]=x[i]-initdate;
		}
		
		// cast to double
		int order=or;
		double MeanSqErr = 0;double r=0;
		for(int i=order+1;i<x.length;i++){
			double[] xtrain = new double[i];
			double[] ytrain = new double[i];
			
			for(int z=0;z<i;z++){
				xtrain[z] = x[z];
				ytrain[z] = y[z];
			}
			double[] est = FitPoly(xtrain,ytrain,order);

			double ypredict=0;
			if(order==1){
				ypredict = est[0]*x[i];
			} else if(order==2){
				ypredict = est[0]*x[i] + est[1]*x[i]*x[i] ;
			}else if(order==3){
				ypredict = est[0]*x[i] + est[1]*x[i]*x[i] + est[2]*x[i]*x[i]*x[i] ;
			}
			MeanSqErr += ((y[i]-ypredict)*(y[i]-ypredict));
			r++;
		}
		MeanSqErr = MeanSqErr/r;
		return MeanSqErr;		
	}
	
	
	public double LOOCV(double[] x, double[] y,int or){
		double initdate = x[0];
		for(int i=0;i<xdata.size();i++){
			x[i]=x[i]-initdate;
		}
		int order=or;
		double MeanSqErr = 0;double r=0;
		int trainlength = x.length-1;
		for(int i=0;i<x.length;i++){ // number of times
			double[] xtrain = new double[trainlength];
			double[] ytrain = new double[trainlength];
			int z=0;
			for(int j=0;j<x.length;j++){
				if(j!=i){
					xtrain[z]=x[j];
					ytrain[z]=y[j];
					z++;
				}
			}
			double[] est = FitPoly(xtrain,ytrain,order);

			double ypredict=0;
			if(order==1){
				ypredict = est[0]*x[i];
			} else if(order==2){
				ypredict = est[0]*x[i] + est[1]*x[i]*x[i] ;
			}else if(order==3){
				ypredict = est[0]*x[i] + est[1]*x[i]*x[i] + est[2]*x[i]*x[i]*x[i] ;
			}
			MeanSqErr += ((y[i]-ypredict)*(y[i]-ypredict));
			r++;
			
		}
		MeanSqErr = MeanSqErr/r;
		return MeanSqErr;	
	}
	public double LOOCV(int or){
		double[] x = new double[xdata.size()];
		double[] y = new double[ydata.size()];      
		double initdate = xdata.get(0);
		for(int i=0;i<xdata.size();i++){
			x[i]=xdata.get(i)-initdate;
			y[i]=ydata.get(i);
		}
		int order=or;
		double MeanSqErr = 0;double r=0;
		int trainlength = x.length-1;
		for(int i=0;i<x.length;i++){ // number of times
			double[] xtrain = new double[trainlength];
			double[] ytrain = new double[trainlength];
			int z=0;
			for(int j=0;j<x.length;j++){
				if(j!=i){
					xtrain[z]=x[j];
					ytrain[z]=y[j];
					z++;
				}
			}
			double[] est = FitPoly(xtrain,ytrain,order);

			double ypredict=0;
			if(order==1){
				ypredict = est[0]*x[i];
			} else if(order==2){
				ypredict = est[0]*x[i] + est[1]*x[i]*x[i] ;
			}else if(order==3){
				ypredict = est[0]*x[i] + est[1]*x[i]*x[i] + est[2]*x[i]*x[i]*x[i] ;
			}
			MeanSqErr += ((y[i]-ypredict)*(y[i]-ypredict));
			r++;
			
		}
		MeanSqErr = MeanSqErr/r;
		return MeanSqErr;	
	}
	
	

	public double[] MSEval(){
		double[] mse = new double[3];
		mse[0] = CrossValidation(1);
		mse[1] = CrossValidation(2);
		mse[2] = CrossValidation(3);	
		return mse;
	}
	
	public double[] MSEval(double[] x, double[] y){
		double[] mse = new double[3];
		mse[0] = CrossValidation(x,y,1);
		mse[1] = CrossValidation(x,y,2);
		mse[2] = CrossValidation(x,y,3);	
		return mse;
	}
	
	
	public double WhichModel(double[] x, double[] y){
		double[] mse = new double[3];
		mse[0] = CrossValidation(x,y,1);
		mse[1] = CrossValidation(x,y,2);
		mse[2] = CrossValidation(x,y,3);
		
		double min =mse[0] ;
		double ind = 1;
		for ( int k=0; k<mse.length; k++ ){
			if ( mse[k] < min ) {
				min = mse[k];
				ind = k+1;
				}
		}
		return ind;
	}
	public double WhichModel(){
		double[] mse = new double[3];
		mse[0] = CrossValidation(1);
		mse[1] = CrossValidation(2);
		mse[2] = CrossValidation(3);
		double min =mse[0] ;
		double ind = 1;
		for ( int k=0; k<mse.length; k++ ){
			if ( mse[k] < min ) {
				min = mse[k];
				ind = k+1;
				}
		}
		return ind;
	}
	


	public void printxdata(){
		for(int i=0;i<xdata.size();i++){
			System.out.println(xdata.get(i));
		}
	}
	public void printxdata(double[] x){
		for(int i=0;i<x.length;i++){
			System.out.println(x[i]);
		}
	}
	public void printydata(){
		for(int i=0;i<ydata.size();i++){
			System.out.println(ydata.get(i));
		}
	}

}
