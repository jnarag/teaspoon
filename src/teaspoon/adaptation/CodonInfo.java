package teaspoon.adaptation;

public class CodonInfo {
	int[] codonbases;
	int codonnumber;
	double probability;
	
	
	
	public CodonInfo(){
		codonbases=new int[3];
		codonnumber=0;
		probability=0.0;
	}
	
	
	public int getcodonnumber(int pos1, int pos2, int pos3){ //same method as in teaspoon.adaptation.SiteEstMulti
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
	
	public void setpos1(String x){
		codonbases[0]=Integer.parseInt(x);
	}
	public void setpos2(String x){
		codonbases[1]=Integer.parseInt(x);
	}
	public void setpos3(String x){
		codonbases[2]=Integer.parseInt(x);
	}
	public void setprob(String x){
		probability=Double.parseDouble(x);
	}
	
}
