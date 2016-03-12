package teaspoon.adaptation;

import java.io.File;

import java.util.ArrayList;
import java.util.Collections;


import org.jdom.Document;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;
//remember to change paths from /Users/sam/Desktop/Pybus bhatt web program/xmls/ to new path when mounting image

public class OutputXmlOther {
	public ArrayList<String> data = new ArrayList<String>();
	public ArrayList<String> ansdata = new ArrayList<String>();
	public ArrayList<String> names = new ArrayList<String>();
	private String[] tajimas_headers = new String[3];		// headers for tajimas D
	private String[] fuli_headers = new String[2];
	private String[] estimator_headers = new String[3];
	private String[] counts_headers = new String[3];
	private String[] faywu_headers = new String[2];
	private String[] mk_headers = new String[4];

	public OutputXmlOther(){
		// read directory for content and get all the sequences extracted by php
		File input = new File("/Users/sam/Desktop/Pybus bhatt web program/seq/");												// Input 
		File[] list = input.listFiles();	
		// print absolute path

		for(int i=0;i<list.length;i++){
			if(list[i].isHidden() == false && list[i].isDirectory() == false)	{
				data.add(list[i].getAbsolutePath());
			}
		}	
		Collections.sort(data);
		// read directory for content and get all the sequences extracted by php
		File ans = new File("/Users/sam/Desktop/Pybus bhatt web program/ans/");												// Input 
		File[] anslist = ans.listFiles();	

		for(int i=0;i<anslist.length;i++){
			if(anslist[i].isHidden() == false && anslist[i].isDirectory() == false)	{
				ansdata.add(anslist[i].getAbsolutePath());
			}
		}	
		Collections.sort(ansdata);


		//****************
		//  normal file names
		for(int i=0;i<list.length;i++){
			if(list[i].isHidden() == false && list[i].isDirectory() == false)	{
				names.add(list[i].getName());
			}
		}
		Collections.sort(names);

		// headers for Tajimas D, ful and li, estimators, Mk test and fay and wu
		tajimas_headers[0] = "D";tajimas_headers[1] = "highConfidence";tajimas_headers[2] = "lowConfidence";
		fuli_headers[0] = "D*";fuli_headers[1] = "F*";
		counts_headers[0] = "Number of Segregating Sites"; counts_headers[1] = "Number of Pairwise Difference"; counts_headers[2] = "Number of Singletons";
		estimator_headers[0] = "Theta-S";estimator_headers[1] = "Theta-K";estimator_headers[2] = "Theta-n";
		faywu_headers[0] = "D"; faywu_headers[1] = "H"; 
		mk_headers[0] = "Chi-squared probability";mk_headers[1] = "Chi-squared value";mk_headers[2] = "Cramer's V";
	}

	public void OutputXmlSeqNames(){
		try {
			Element root = new Element("Names");
			for(int m=0;m<names.size();m++){
				Element seqname = new Element("Values");
				seqname.addContent(names.get(m).toString());
				root.addContent(seqname);
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/Names.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	public void OutputToXmlCounts(){
		try{
			Element root = new Element("Counts");	
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<counts_headers.length;i++){
				sbuf.append(counts_headers[i]);							// use string buffer
				sbuf.append("  ");
			}
			headerNames.addContent(sbuf.toString());
			Element Alignments = new Element("Alignments");
			Alignments.addContent(String.valueOf(data.size()));
			root.addContent(Alignments);
			root.addContent(headerNames);	

			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Number",String.valueOf(m));
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
				Read_main read = new Read_main(data.get(m).toString());
				int[][] seq = read.read();
				FuAndLi FL = new FuAndLi(seq);
				double[] values_array = FL.counts();	// watterson estimates
				float[] printed_values = new float[values_array.length];
				for(int i=0;i<values_array.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<values_array.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}
				sequence.addContent(values);
				root.addContent(sequence);			
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/Counts.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}		
	public void OutputToXmlEstimators(){
		try{
			Element root = new Element("Estimators");	
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<estimator_headers.length;i++){
				sbuf.append(estimator_headers[i]);							// use string buffer
				sbuf.append("  ");
			}
			headerNames.addContent(sbuf.toString());
			Element Alignments = new Element("Alignments");
			Alignments.addContent(String.valueOf(data.size()));
			root.addContent(Alignments);
			root.addContent(headerNames);	
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
				Read_main read = new Read_main(data.get(m).toString());
				int[][] seq = read.read();
				FuAndLi FL = new FuAndLi(seq);
				double[] values_array = FL.wattersonEstimates();	// watterson estimates
				float[] printed_values = new float[values_array.length];
				for(int i=0;i<values_array.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<values_array.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}

				sequence.addContent(values);
				root.addContent(sequence);			
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/Estimators.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public void OutputToXmlTajima(){
		try{
			Element root = new Element("Tajimas-D-teaspoon.adaptation.Test-1989-Results");
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<tajimas_headers.length;i++){
				sbuf.append(tajimas_headers[i]);							// use string buffer
				sbuf.append("  ");
			}
			headerNames.addContent(sbuf.toString());
			Element Alignments = new Element("Alignments");
			Alignments.addContent(String.valueOf(data.size()));
			root.addContent(Alignments);
			root.addContent(headerNames);	
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
				Read_main read = new Read_main(data.get(m).toString());
				int[][] seq = read.read();
				DiversityStats TD = new DiversityStats(seq);
				double[] values_array = TD.TajimasDValues();
				float[] printed_values = new float[values_array.length];
				for(int i=0;i<values_array.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<values_array.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}
				sequence.addContent(values);
				root.addContent(sequence);			
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/TajimasD.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public void OutputToXmlFuLi(){
		try{
			Element root = new Element("Fu-Li-teaspoon.adaptation.Test-1992-Results");
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<fuli_headers.length;i++){
				sbuf.append(fuli_headers[i]);							// use string buffer
				sbuf.append("  ");
			}
			headerNames.addContent(sbuf.toString());
			Element Alignments = new Element("Alignments");
			Alignments.addContent(String.valueOf(data.size()));
			root.addContent(Alignments);
			root.addContent(headerNames);		
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
				Read_main read = new Read_main(data.get(m).toString());
				int[][] seq = read.read();
				FuAndLi FL = new FuAndLi(seq);
				double[] values_array = FL.FuLiTest();
				float[] printed_values = new float[values_array.length];
				for(int i=0;i<values_array.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<values_array.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}
				sequence.addContent(values);
				root.addContent(sequence);			
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/teaspoon.adaptation.FuAndLi.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public void OutputToXmlFayWu(String needconsensus){
		try{
			Element root = new Element("Fay-Wu-teaspoon.adaptation.Test-2000-Results");
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<faywu_headers.length;i++){
				sbuf.append(faywu_headers[i]);							// use string buffer
				sbuf.append("  ");
			}
			headerNames.addContent(sbuf.toString());
			Element Alignments = new Element("Alignments");
			Alignments.addContent(String.valueOf(data.size()));
			root.addContent(Alignments);
			root.addContent(headerNames);	
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
				Read_main read = new Read_main(data.get(m).toString());
				// this code is to check weather the ancestral sequence is an alignment or a single sequence
				int[] ans;
				if(needconsensus.equals("y")){	// if is consensus
					Read_ancestral_consensus readans = new Read_ancestral_consensus(ansdata.get(m).toString());
					ans = readans.read();
				} else {	//else 
					Read_ancestral readans = new Read_ancestral(ansdata.get(m).toString());
					ans = readans.read();
				}
				int[][] seq = read.read();
				FayAndWu FW = new FayAndWu(seq,ans);
				double[] values_array = FW.FayWu();
				float[] printed_values = new float[values_array.length];
				for(int i=0;i<values_array.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<values_array.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}
				sequence.addContent(values);
				root.addContent(sequence);			
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/teaspoon.adaptation.FayAndWu.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public void OutputToXmlMK(String needconsensus){
		try{
			Element root = new Element("McDonald-Kreitman-teaspoon.adaptation.Test-1991-Results");
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<mk_headers.length;i++){
				sbuf.append(mk_headers[i]);							// use string buffer
				sbuf.append("  ");
			}
			headerNames.addContent(sbuf.toString());
			Element Alignments = new Element("Alignments");
			Alignments.addContent(String.valueOf(data.size()));
			root.addContent(Alignments);
			root.addContent(headerNames);	
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
				Read_main read = new Read_main(data.get(m).toString());
				// this code is to check weather the ancestral sequence is an alignment or a single sequence
				int[] ans;
				if(needconsensus.equals("y")){	// if is consensus
					Read_ancestral_consensus readans = new Read_ancestral_consensus(ansdata.get(m).toString());
					ans = readans.read();
				} else {	//else 
					Read_ancestral readans = new Read_ancestral(ansdata.get(m).toString());
					ans = readans.read();
				}
				int[][] seq = read.read();
				McDonaldKreitman MK = new McDonaldKreitman(seq,ans);
				double[] values_array = MK.significanceTest();
				float[] printed_values = new float[mk_headers.length];
				for(int i=0;i<mk_headers.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<mk_headers.length;i++){		
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}
				sequence.addContent(values);
				root.addContent(sequence);			
			}

			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/MK.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}



}
