package teaspoon.adaptation;

public class TimeSeriesAligner {

	public TimeSeriesAligner(){

	}

	public double[][] Aligner(double [][] data,double[] newdays) {
		double[] days = data[0]; // input days vector
		double[] price = data[1]; // input price vector
		int L = days.length;	// old vector length
		int newL = newdays.length; // new vector length
		double[] newprice = new double[newL]; // new aligned price vector
		int index=0;
		for(int k=1;k<L;k++){
			while(newdays[index]<days[k]){
				if(newdays[index]<days[0]){ // first if statement to check for missing siteData and assign -1 to these
					newprice[index]=-1;
				} 
				if(newdays[index]>=days[k-1] && newdays[index]<days[k]){ // MAIN STATEMENT: assigns prices to the new time series based on the old time series
					newprice[index]=price[k-1];
				} 
				index++;
			}		
		}
		for(int i=index;i<newL;i++){  // this completes the new time series by adding the last value of the old series to every subsequent new one
			newprice[i] = price[L-1];
		}
		
		
		// make siteData into matrix for exporting
		double[][] aligned = new double[2][newdays.length];
		for(int i=0;i<aligned[0].length;i++){
			aligned[0][i]=newdays[i];
			aligned[1][i]=newprice[i];
		}
		
		return aligned;	
	}

// Output 
	public static void main(String[] args) {
		
		double[] days = {2.4,5,5.5,5.7};	//define original times
		double[] prices = {104.2,109.2,112.2,114}; // define original series
		double[][] data = new double[2][prices.length];
		for(int i=0;i<days.length;i++){
			data[0][i]=days[i];
			data[1][i]=prices[i];
		}
		double[] newdays = {1,2,3,4,5,7};	// define new time series
		TimeSeriesAligner A = new TimeSeriesAligner(); // Initialise
		double[][] results =A.Aligner(data,newdays); // Run function
		// Output everything
		System.out.println("Original Data");
		for(int i=0;i<data.length;i++){
			for(int j=0;j<data[0].length;j++){
				System.out.print(data[i][j]+"\t");
			}
			System.out.println();
		}
		System.out.println("\n");
		System.out.println("Aligned Data");
	
		
		
		for(int i=0;i<results.length;i++){
			for(int j=0;j<results[0].length;j++){
				System.out.print(results[i][j]+"\t");
			}
			System.out.println();
		}


	}

}
