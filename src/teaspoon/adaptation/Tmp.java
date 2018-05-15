package teaspoon.adaptation;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

public class Tmp{

	public static void main(String[] args){

		ArrayList<String> data = new ArrayList<String>();
		ArrayList<String> names = new ArrayList<String>();
		String[] tajimas_headers = new String[3];
		String[] tajima_headers1 = new String[2];
		String[] estimator_headers = new String[3];
		String[] headers = new String[7];
		String[] mk_headers = new String[5];  // remember to change back
		File input = new File(args[4]);												// Input 
		File[] list = input.listFiles();
		ArrayList<String> ansdata = new ArrayList<String>();

		for(int i=0;i<list.length;i++){
			if(list[i].isHidden() == false && list[i].isDirectory() == false)	{
				data.add(list[i].getAbsolutePath());
			}
		}	
		Collections.sort(data);


		File ans = new File(args[5]);												// Input 
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
		tajima_headers1[0] = "Low Confidence - 5%";tajima_headers1[1] = "High Confidence - 95%";
		estimator_headers[0] = "Theta-S";estimator_headers[1] = "Theta-K";estimator_headers[2] = "Theta-n";
		mk_headers[0] = "Chi-squared probability";mk_headers[1] = "Chi-squared value";mk_headers[2] = "Cramer's V";
		headers[0] = "Number Silent";headers[1] = "Number Replacement";headers[2] = "Total";headers[3] = "Silent/Replacement Ratio";
		headers[4] = "Neutral Ratio"; headers[5] = "Number Non Neutral";headers[6] = "Proportion Non-Neutral";
		//**********************************************************************************************************************************************
		
		
		
		/*	try {
			Element root = new Element("Names");
			for(int m=0;m<names.size();m++){
				Element seqname = new Element("Values");
				seqname.addContent(names.get(m).toString());
				root.addContent(seqname);
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter(args[0]);	// write to results.xml  ***** input location
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		} */
		//**********************************************************************************************************************************************
	/*			try{
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
				teaspoon.adaptation.Read_main read = new teaspoon.adaptation.Read_main(data.get(m).toString());
				int[][] seq = read.read();
				teaspoon.adaptation.FuAndLi FL = new teaspoon.adaptation.FuAndLi(seq);
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
			java.io.FileWriter writer = new java.io.FileWriter(args[1]);	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}	*/
		//**********************************************************************************************************************************************
	/*				try{
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
				teaspoon.adaptation.Read_main read = new teaspoon.adaptation.Read_main(data.get(m).toString());
				int[][] seq = read.read();
				teaspoon.adaptation.DiversityStats TD = new teaspoon.adaptation.DiversityStats(seq);
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
			java.io.FileWriter writer = new java.io.FileWriter(args[2]);	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		} 
		//**********************************************************************************************************************************************		
		try{
				Element root = new Element("Tajimas-D-Confidence-Interval");	
				Element headerNames = new Element("Headers");			// print header names for columns
				// add headers
				StringBuffer sbuf = new StringBuffer();
				for(int i=0;i<tajima_headers1.length;i++){
					sbuf.append(tajima_headers1[i]);							// use string buffer
					sbuf.append("     ");
				}
				headerNames.addContent(sbuf.toString());					
				root.addContent(headerNames);
				for(int m=0;m<data.size();m++){		// Loop through all tiles ****
					teaspoon.adaptation.Read_main read = new teaspoon.adaptation.Read_main(data.get(m).toString());
					int[][] seq = read.read();
					int sample_size = seq.length;
					teaspoon.adaptation.DiversityStats TD = new teaspoon.adaptation.DiversityStats(seq);
					double theta  = TD.theta();	
					Element sequence = new Element("Simulation");	
					sequence.setAttribute("Sequence", String.valueOf(m));
					sequence.setAttribute("Sample-Size", String.valueOf(sample_size));
					sequence.setAttribute("Theta", String.valueOf(theta));
		//			theta=Double.valueOf(args[5]); //*******
					teaspoon.adaptation.CoalescentSim CS = new teaspoon.adaptation.CoalescentSim(sample_size,theta); //******
					double[] values_array = CS.TajimaCS();	
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
				java.io.FileWriter writer = new java.io.FileWriter(args[3]);	// write to results.xml
				out.output(document, writer);
				writer.flush();
				writer.close();

			}catch (Exception ex) {
				ex.printStackTrace();
			} */


	/*		try{
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
				teaspoon.adaptation.Read_main read = new teaspoon.adaptation.Read_main(data.get(m).toString());
				// this code is to check weather the ancestral sequence is an alignment or a single sequence
				int[] ansmat;

				teaspoon.adaptation.Read_ancestral readans = new teaspoon.adaptation.Read_ancestral(ansdata.get(m).toString());
				ansmat = readans.read();			
				int[][] seq = read.read();
				teaspoon.adaptation.McDonaldKreitman MK = new teaspoon.adaptation.McDonaldKreitman(seq,ansmat);
				double[] values_array = MK.significanceTest();
				double distance = MK.Distance();
				//double distance = MK.K2P();
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
				Element node = new Element("PWD");
				node.addContent(String.valueOf(distance));
				values.addContent(node);
				
				sequence.addContent(values);
				root.addContent(sequence);			
			}

			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter(args[6]);	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		} */


	/*	try{
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
			teaspoon.adaptation.Read_main read = new teaspoon.adaptation.Read_main(data.get(m).toString());
			// this code is to check weather the ancestral sequence is an alignment or a single sequence
			int[] ansmat;

			teaspoon.adaptation.Read_ancestral readans = new teaspoon.adaptation.Read_ancestral(ansdata.get(m).toString());
			ansmat = readans.read();			
			int[][] seq = read.read();
			teaspoon.adaptation.SiteEstMulti s_mid = new teaspoon.adaptation.SiteEstMulti(seq,ansmat);
	//		double[] L = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
	//		double[] H = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
			double[] L = {0.0,0.5};
			double[] H = {0.5,1.0};
			boolean[] Nvec = {false,false,false,true,true,true,true,false,false,false};
			double[][] bins = new double[2][L.length];
			for(int i=0;i<L.length;i++){
				bins[0][i]=L[i];
				bins[1][i]=H[i];
			}
			double values_array = s_mid.totalNoAdapt(Nvec, bins, 0.0);
			double distance = s_mid.Distance();
			//double distance = MK.K2P();
			float printed_values = (float) values_array;

			Element values = new Element("Values");	
			Element node = new Element("V"+String.valueOf(1));
			node.addContent(String.valueOf(printed_values));
			values.addContent(node);
			
			Element node2 = new Element("PWD");
			node.addContent(String.valueOf(distance));
			values.addContent(node2);
			
			sequence.addContent(values);
			root.addContent(sequence);			
		}

		Document document = new Document(root);				// create document on root
		XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
		java.io.FileWriter writer = new java.io.FileWriter(args[6]);	// write to results.xml
		out.output(document, writer);
		writer.flush();
		writer.close();

	}catch (Exception ex) {
		ex.printStackTrace();
	} */

/*			try{
			Element root = new Element("Eyre-Walker-teaspoon.adaptation.Test-1991-Results");
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<headers.length;i++){
				sbuf.append(headers[i]);							// use string buffer
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
				teaspoon.adaptation.Read_main read = new teaspoon.adaptation.Read_main(data.get(m).toString());
				// this code is to check weather the ancestral sequence is an alignment or a single sequence
				int[] ansmat;

				teaspoon.adaptation.Read_ancestral readans = new teaspoon.adaptation.Read_ancestral(ansdata.get(m).toString());
				ansmat = readans.read();

				int[][] seq = read.read();
				teaspoon.adaptation.Williamson W = new teaspoon.adaptation.Williamson(seq,ansmat);
				double[][] value_matrix = W.eyrewalker_method();
				float[][] printed_values = new float[value_matrix.length][value_matrix[0].length];
				for(int i=0;i<value_matrix.length;i++){
					for(int j=0;j<value_matrix[0].length;j++){
						printed_values[i][j] = (float) value_matrix[i][j];
					}
				}
				Element values = new Element("Values");	
	//			Element bin1 = new Element("B1");
	//			for(int i=0;i<value_matrix.length;i++){
	//				Element node = new Element("V"+String.valueOf(i));
	//				node.addContent(String.valueOf(printed_values[i][0]));
	//				bin1.addContent(node);
	//			}
	//			values.addContent(bin1);
				Element bin2 = new Element("B2");
				for(int i=0;i<value_matrix.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i][1]));
					bin2.addContent(node);
				}
				values.addContent(bin2);
				sequence.addContent(values);
				root.addContent(sequence);			
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter(args[6]);	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}  catch (Exception ex) {
			ex.printStackTrace();
		} */
/*		try{
			Element root = new Element("teaspoon.adaptation.Williamson-teaspoon.adaptation.Test-1991-Results");
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<headers.length;i++){
				sbuf.append(headers[i]);							// use string buffer
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
				teaspoon.adaptation.Read_main read = new teaspoon.adaptation.Read_main(data.get(m).toString());
				// this code is to check weather the ancestral sequence is an alignment or a single sequence
				int[] ansmat;

				teaspoon.adaptation.Read_ancestral readans = new teaspoon.adaptation.Read_ancestral(ansdata.get(m).toString());
				ansmat = readans.read();

				int[][] seq = read.read();
				teaspoon.adaptation.Williamson W = new teaspoon.adaptation.Williamson(seq,ansmat);
				double[][] value_matrix = W.williamson_method();
				float[][] printed_values = new float[value_matrix.length][value_matrix[0].length];
				for(int i=0;i<value_matrix.length;i++){
					for(int j=0;j<value_matrix[0].length;j++){
						printed_values[i][j] = (float) value_matrix[i][j];
					}
				}
				Element values = new Element("Values");	
				Element bin1 = new Element("sigma");
				for(int i=0;i<value_matrix[0].length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[0][i]));
					bin1.addContent(node);
				}
				values.addContent(bin1);
				Element bin2 = new Element("rho");
				for(int i=0;i<value_matrix[0].length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[1][i]));
					bin2.addContent(node);
				}
				values.addContent(bin2);
				Element bin3 = new Element("alpha");
				for(int i=0;i<value_matrix[0].length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[3][i]));
					bin3.addContent(node);
				}
				values.addContent(bin3);	
				sequence.addContent(values);
				root.addContent(sequence);
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter(args[6]);	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}  catch (Exception ex) {
			ex.printStackTrace();
		}*/
		
		
		
		
		

		// Start Formatting
/*		try{
			Element root = new Element("Pybus-Bhatt-results");		// root element
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<headers.length;i++){
				sbuf.append(headers[i]);							// use string buffer
				sbuf.append("  ");
			}
			headerNames.addContent(sbuf.toString());					
			Element Alignments = new Element("Alignments");
			Alignments.addContent(String.valueOf(data.size()));
			root.addContent(Alignments);
			root.addContent(headerNames);							// add header names to root
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
				Element datamat = new Element("Values");		// data matrix
				teaspoon.adaptation.Read_main read = new teaspoon.adaptation.Read_main(data.get(m).toString());
				// this code is to check weather the ancestral sequence is an alignment or a single sequence
				int[] ansmat;

				teaspoon.adaptation.Read_ancestral readans = new teaspoon.adaptation.Read_ancestral(ansdata.get(m).toString());
				ansmat = readans.read();

				int[][] seq = read.read();
	//			teaspoon.adaptation.SiteEstMulti Pb = new teaspoon.adaptation.SiteEstMulti(seq,ansmat,new Double(10.0));
	//			boolean[] vec = {false,false,false,false,true,true,false,false,false,false};
				teaspoon.adaptation.SiteEstMulti Pb = new teaspoon.adaptation.SiteEstMulti(seq,ansmat);
				Pb.SetNumBins(3.0);
				boolean[] vec = {true,false,false};
				double [][] vec2 = new double[2][3];
				vec2[0][0]=0.0;
				vec2[1][0]=0.5;
				vec2[0][1]=0.5;
				vec2[1][1]=0.9;
				vec2[0][2]=0.9;
				vec2[1][2]=1.0;

				
				double[][] matrix = Pb.value_matrix(vec,vec2);	// values
			//	double[][] intervals = Pb.intervals();
				double[][] intervals = vec2;
				float[][] printed_intervals = new float[intervals.length][intervals[0].length];
				float[][] printed_matrix = new float[matrix.length][matrix[0].length];
				//  create float values
				for(int i=0;i<matrix.length;i++){
					for(int j=0;j<matrix[0].length;j++){
						Double temp = new Double(matrix[i][j]);
						printed_matrix[i][j] = temp.floatValue();	// printed matrix with float values
					}
				}
				for(int i=0;i<printed_intervals.length;i++){
					for(int j=0;j<printed_intervals[0].length;j++){
						Double temp = new Double(intervals[i][j]);
						printed_intervals[i][j] = temp.floatValue();
					}
				}
				Element values = new Element("Values");	
				for(int j=0;j<printed_matrix[0].length;j++){
					Element node = new Element("V"+String.valueOf(j));
					// add matrix content
					node.addContent(String.valueOf(printed_matrix[6][j]));
					values.addContent(node);
					// add matrix
				}
				sequence.addContent(values);
				root.addContent(sequence);						// add to root
			}			
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter(args[6]);	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}*/
		
		

	}




}

