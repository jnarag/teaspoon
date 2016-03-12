package teaspoon.adaptation;

public class RunFile2 {

	public static void main(String[] args) { 
		Split s = new Split(args[0],args[1],args[2],args[3]);
		s.SeqtoFile();
		s.AnstoFile();
	}
}
