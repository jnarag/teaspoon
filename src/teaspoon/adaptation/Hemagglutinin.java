package teaspoon.adaptation;

public class Hemagglutinin {
	public int[] HAH3 = new int[1698/3];
	public int[] HAH1 = new int[1698/3];
	public int[] HA = new int[1698/3];
	public int[] H1 = new int[1698/3];
	public int[] H3 = new int[1698/3];
	public int[] FluBMix = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,2,0,0,0,0,0,0,0,1,2,0,2,1,1,0,1,0,2,2,1,2,2,0,2,0,1,1,1,1,2,0,1,2,2,0,2,2,1,2,0,0,1,0,2,2,1,2,1,0,1,2,0,0,22,0,2,2,0,0,0,0,0,0,0,0,0,1,1,2,0,1,1,2,0,2,2,0,2,0,0,0,0,1,1,1,2,0,1,0,2,0,0,0,0,0,1,0,2,2,1,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,2,2,2,1,12,0,2,2,0,2,2,1,2,0,2,2,0,2,1,2,1,0,1,2,2,2,2,1,2,1,1,0,1,0,0,0,0,0,0,0,2,-1,-1,-1,2,2,2,2,1,2,2,1,2,0,2,0,1,1,2,0,2,2,2,0,1,0,0,0,0,0,0,0,0,00,2,2,2,1,0,2,2,1,0,1,1,1,2,2,0,2,0,1,0,1,0,2,2,2,2,1,2,1,2,0,1,1,1,2,1,2,2,1,2,1,1,2,2,2,1,2,2,0,0,0,0,0,0,0,0,1,0,2,2,2,2,2,0,1,0,1,0,0,00,0,0,0,0,0,1,0,0,0,0,1,2,2,0,1,0,0,2,0,2,1,2,1,2,2,2,0,1,0,0,0,1,2,1,0,1,1,2,1,2,1,1,0,0,0,2,2,1,1,0,0,0,1,2,1,0,1,1,0,2,2,2,0,1,0,0,1,0,01,2,0,1,2,1,2,2,1,2,-1,-1,1,2,2,0,2,2,1,2,1,1,2,0,0,1,2,2,0,2,2,0,1,0,0,0,0,0,1,1,2,1,2,2,0,1,0,1,0,2,2,0,0,2,2,0,1,1,2,0,2,2,1,1,2,2,0,1,2,12,1,2,1,1,1,0,1,1,2,0,0,1,2,2,1,2,2,2,1,2,2,1,2,2,2,0,2,2,1,0,2,2,1,0,2,2,0,1,2,1,0,1,2,1,0,1,2,1,0,1,1,0,1,2,1,0,1,2,0,0,2,1,0,2,2,1,0,0,21,1,2,1,2,2,2,0,0,0,1,0,2,0,2,0,1,2,2,0,0,2,1,0,1,2,2,2,0,2,0,1,2,1,2,1,2,1,2,2,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	public int[] FluB={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,0,1,0,1,1,1,1,1,0,1,1,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,0,1,0,0,0,0,1,1,1,1,0,1,0,1,0,0,0,0,0,1,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,1,1,0,1,0,1,0,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0,1,0,0,1,0,1,1,1,1,1,1,1,0,1,0,0,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,0,0,1,0,0,1,1,0,1,1,1,1,1,1,1,0,0,1,1,1,0,1,1,1,1,1,1,1,0,0,1,1,1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,1,1,0,1,0,1,0,1,1,0,0,1,1,0,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,0,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,1,1,0,0,1,1,1,1,1,1,1,1,0,0,0,1,0,1,0,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	public Hemagglutinin(){
		int[] tmp ={2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,0,1,0,0,0,0,0,1,1,1,1,1,1,0,1,0,0,1,1,1,1,0,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,1,1,1,1,0,1,0,0,0,0,1,1,1,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,1,0,0,1,1,1,1,1,0,1,1,0,1,1,1,2,1,1,0,1,0,0,1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,0,1,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,1,0,1,1,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1,1,0,1,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,0,1,1,0,1,0,1,1,0,1,0,0,1,0,0,1,1,1,1,1,2,2,2,1,1,0,0,1,0,0,1,1,1,0,1,1,0,1,1,1,0,1,1,1,0,0,0,0,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,0,0,1,1,0,0,1,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,1,0,1,1,0,0,1,1,0,1,1,1,0,1,1,1,0,1,1,1,1,1,0,1,0,1,0,1,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
		this.H1=tmp;
		int[] tmp2 ={2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,0,1,0,1,0,0,0,1,0,1,1,1,1,1,1,0,1,0,0,0,1,1,1,0,1,0,1,1,1,1,1,0,0,0,1,1,1,1,1,1,0,0,1,1,1,0,1,1,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,1,1,0,1,1,0,1,1,1,1,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0,1,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,0,1,0,1,1,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,1,0,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,0,1,1,1,1,0,1,0,0,1,0,0,1,0,1,1,1,2,2,2,2,1,1,1,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,0,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,0,1,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,1,0,1,0,1,1,0,1,0,1,1,1,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
		this.H3=tmp2;
		// for ha1 ha2 h1

		for(int i=0;i<HAH1.length;i++){
			if(i>=13 && i<=340){
				HAH1[i]=1;
			}else if(i>=341 && i<=503){
				HAH1[i]=0;
			} else {
				HAH1[i]=2;
			}
		}
		
		// for ha1 ha2 h3
		
		for(int i=0;i<HAH3.length;i++){
			if(i>=24 && i<=344){
				HAH3[i]=1;
			}else if(i>=345 && i<=516){
				HAH3[i]=0;
			} else {
				HAH3[i]=2;
			}
		}

	}
	
	


	
	
}
