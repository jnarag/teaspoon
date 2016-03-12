package teaspoon.adaptation;

import jebl.evolution.sequences.Sequence;
import jebl.evolution.sequences.SequenceType;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: jayna
 * Date: 12/08/2013
 * Time: 10:59
 * To change this template use File | Settings | File Templates.
 */
public class ProcessFiles {

    BufferedReader reader = null;
    String filelist = null;

    public ProcessFiles() {

        filelist = "filelist.txt";


    }

    public void readFiles() {

        try {
            reader = new BufferedReader(new FileReader(filelist));

            while(reader.ready()) {

                String fileName = reader.readLine().trim();

                BufferedReader reader2 = new BufferedReader(new FileReader(fileName));

                List<String> sequenceText = new ArrayList<String>();
                while(reader2.ready()) {

                    sequenceText.add(reader2.readLine().trim());


                }

                //write out Nexus file

                String newFileName = fileName.replaceAll("H1N1/","H1N1/processed_");
                BufferedWriter writer = new BufferedWriter(new FileWriter(new File(newFileName)));

                String [] parts = sequenceText.get(0).split("\\s");
                writer.write("#NEXUS\nBegin Data;\nDimensions ntax="+sequenceText.size()+" nchar="+parts[1].trim().length()+";\nFormat datatype=DNA;\nMatrix\n");

                for(String s: sequenceText) {

                    String [] strings = s.split("\\s");
                    String name = strings[0].trim();
                    String seq = strings[1].trim();

                    writer.write(name+"\t"+seq+"\n");
                    writer.flush();
                }

                writer.write(";\nEnd;\n");
                writer.close();



            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    public void concatenateNexusFiles(String concatenateList) {
        List<Sequence> allSeqs = new ArrayList<Sequence>();
        try {
            reader = new BufferedReader(new FileReader(concatenateList));
            BufferedWriter writer = new BufferedWriter(new FileWriter("concatenated_NS1.H1N1_2.nex"));


            while(reader.ready())  {

                String fileName = reader.readLine().trim().replaceAll("M2","NS1");

                BufferedReader nexusReader = new BufferedReader(new FileReader(fileName));
                jebl.evolution.io.NexusImporter importer = new jebl.evolution.io.NexusImporter(nexusReader);


                List<Sequence> seqs = importer.importSequences();
                allSeqs.addAll(seqs);


            }
            jebl.evolution.io.NexusExporter exporter = new jebl.evolution.io.NexusExporter(writer);
            exporter.exportSequences(allSeqs);

            reader.close();
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (jebl.evolution.io.ImportException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }

    public void phylipToNexus(String filelist) throws IOException {
        try {

            reader = new BufferedReader(new FileReader(filelist));


            while(reader.ready()) {

                String fileName = reader.readLine().trim();
                BufferedReader reader2 = new BufferedReader(new FileReader(fileName));

                BufferedWriter writer = new BufferedWriter(new FileWriter(fileName+".nex"));

                jebl.evolution.io.PhylipSequentialImporter phylipSequentialImporter = new jebl.evolution.io.PhylipSequentialImporter(reader2, SequenceType.NUCLEOTIDE, 360);


                List<Sequence> seqs = phylipSequentialImporter.importSequences();

                jebl.evolution.io.NexusExporter exporter = new jebl.evolution.io.NexusExporter(writer);
                exporter.exportSequences(seqs);


                writer.close();
                reader2.close();
            }


        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (jebl.evolution.io.ImportException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    public void fastaToNexus(String filelist) {

        try {

            reader = new BufferedReader(new FileReader(filelist));


            while(reader.ready()) {

                String fileName = reader.readLine().trim();
                BufferedReader reader2 = new BufferedReader(new FileReader(fileName));

                BufferedWriter writer = new BufferedWriter(new FileWriter(fileName+".nex"));

               jebl.evolution.io.FastaImporter fastaImporter = new jebl.evolution.io.FastaImporter(new File(fileName), SequenceType.NUCLEOTIDE); //(reader2, SequenceType.NUCLEOTIDE, 360);


                List<Sequence> seqs = fastaImporter.importSequences();

                jebl.evolution.io.NexusExporter exporter = new jebl.evolution.io.NexusExporter(writer);
                exporter.exportSequences(seqs);


                writer.close();
                reader2.close();
            }


        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (jebl.evolution.io.ImportException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }



    public static void main(String [] args) {

        ProcessFiles p = new ProcessFiles();
        try {
            p.phylipToNexus("filelist.txt");

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
}
