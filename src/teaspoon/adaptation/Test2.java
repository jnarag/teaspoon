package teaspoon.adaptation;

import java.util.*;
import java.lang.*;
//import umontreal.iro.lecuyer.probdistmulti.*;
//import umontreal.iro.lecuyer.probdist.GammaDist;
//import umontreal.iro.lecuyer.rng.*;
//import umontreal.iro.lecuyer.randvarmulti.*;
//import weka.core.ContingencyTables;
//import umontreal.iro.lecuyer.randvarmulti.DirichletGen;
//import umontreal.iro.lecuyer.rng.RandMrg;
//import cc.mallet.types.Dirichlet;
//import dr.app.tools.NexusToPhylip;
//import java.lang.System.*;
public class  Test2 {


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int count=0;
		DataSet D = new DataSet("/Users/sam/Desktop/shankarappa/patient2.fa");
		String[] s = {"s005","s012","s020","s030","s040","s051","s061","s073","s080","s085","s091","s091","s103","s126","s138","s144","s158","s164"};
		D.CreateDataFrame();
		for(int run=0;run<s.length;run++){
			Iterator<SequenceInfo> It =  D.DataMainFrame.iterator();
			ArrayList<SequenceInfo> ss = new ArrayList<SequenceInfo>();
			while(It.hasNext()){
				SequenceInfo element = It.next();
				if(element.Taxon.matches(".*"+s[run]+".*")){

					ss.add(element);
					count++;
				}
			}
			D.exportFASTA("/Users/sam/Desktop/tmp/seqHIV2/"+s[run],ss);
		}
		System.out.println("A total of  "+count);
		
	}

}
