package teaspoon.adaptation;

import java.io.*;
import java.util.*;

public class Split {

	File input;
	String lookfor;
	DataOutputStream dos;
	String fileName;
	String Outputnamea;
	String Outputnames;

	public Split(String Name, String taxon,  String splitNamess,String splitNamesa){
		input = new File(Name);  // The file object
		fileName = Name;
		lookfor = taxon;
		Outputnamea = splitNamesa;
		Outputnames = splitNamess;
	}

	public boolean writeToFile(String fileName, String dataLine, boolean isAppendMode, boolean isNewLine) {
		if (isNewLine) {
			dataLine = "\n" + dataLine;
		}

		try {
			File outFile = new File(fileName);
			if (isAppendMode) {
				dos = new DataOutputStream(new FileOutputStream(fileName, true));
			} else {
				dos = new DataOutputStream(new FileOutputStream(outFile));
			}

			dos.writeBytes(dataLine);
			dos.close();
		} catch (FileNotFoundException ex) {
			return (false);
		} catch (IOException ex) {
			return (false);
		}
		return (true);

	}

	public Vector fileToVector() {
		Vector<String> v = new Vector<String>();
		String inputLine;
		try {
			File inFile = new File(fileName);
			BufferedReader br = new BufferedReader(new InputStreamReader(
					new FileInputStream(inFile)));

			while ((inputLine = br.readLine()) != null) {
				v.addElement(inputLine.trim());
			}
			br.close();
		} // Try
		catch (FileNotFoundException ex) {
			//
		} catch (IOException ex) {
			//
		}
		return (v);
	}

	public int findAns(){
		Vector v = fileToVector();
		int location = 0;
		for(int i=0;i<v.size();i++){
			if(v.elementAt(i).toString().regionMatches(false, 0, lookfor, 0, lookfor.length())){
				location = i;
				break;
			}
		}
		return location;
	}

	public void SeqtoFile() {
		Vector v = fileToVector();
		int index = findAns();
		writeToFile(Outputnames+"s", (String) v.elementAt(0), true, true);
		for (int i = 1; i < v.size(); i++) {
			if(i != index){	
				writeToFile(Outputnames+"s", (String) v.elementAt(i), true, true);
			}
		}
	}
	
	public void AnstoFile() {
		Vector v = fileToVector();
	//	int index = findAns();
		int index = 50; /////cange back

		for (int i = 0; i < v.size(); i++) {
			if(i == index){	
				writeToFile(Outputnamea+"a", (String) v.elementAt(i), true, true);
			}
		}
	}
	





}
