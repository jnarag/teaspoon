package teaspoon.adaptation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Codonmodel {

	String loc;
	public Codonmodel(String location){
				loc=location;
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
	
	public CodonInfo[] makeModel(){
		CodonInfo[] CI = new CodonInfo[67];
		for(int i=0;i<CI.length;i++){
			CI[i] = new CodonInfo();
		}
	//	System.out.println("Tester" + CI[1].codonbases[0]);
		try {
		BufferedReader fh =
	        new BufferedReader(new FileReader(loc));
	       String s;
	       while ((s=fh.readLine())!=null){
	         String f[] = s.split("\t");
	         int pos1 = Integer.parseInt(f[0]);
	         int pos2 = Integer.parseInt(f[1]);
	         int pos3 = Integer.parseInt(f[2]);
	         int cn = getcodonnumber(pos1,pos2,pos3);
	         CI[cn].setpos1(f[0]);
	         CI[cn].setpos2(f[1]);
	         CI[cn].setpos3(f[2]);
	         CI[cn].setprob(f[3]);
	         
	 
	       }
	       fh.close();
	     }
		catch (IOException e) {
			System.err.println("Caught IOException: " +  e.getMessage());
		}
		return CI;

	}
}