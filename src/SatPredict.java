import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

public class SatPredict {
	
	public static void main(String [] args) throws IOException{
		
		 if(args[0].equals("DataProcess")){
			 String dataFile = args[1];
			 String satFile = args[2];
			 String outputFile = args[3];
			 int slidingWindow = Integer.parseInt(args[4]);
			 DataPre.simplifyData(dataFile, satFile, outputFile, slidingWindow);
		 }
		 if(args[0].equals("MotifExtraction")){
			 String type = args[1];
			 int R = Integer.parseInt(args[2]);
			 String train = args[3];
			 String test = args[4];
			 int num = Integer.parseInt(args[5]);
			 ArrayList<Candidate> motifs = new ArrayList<Candidate>();
			 if(type.equals("Frequency")){
				 motifs = GenMotifs.genMotifsByFrequency(R,train,2*num);
			 }
			 if(type.equals("Distance")){
				 motifs = GenMotifs.genMotifsByDistance(R,train,2*num);
			 }
			 if(type.equals("Distribution")){
				 motifs = GenMotifs.genMotifsByDistribution(R,train,2*num);
			 }
			ArrayList<Candidate> pMotifs = GenMotifs.getTopMotifs(motifs, num);
			GenMotifs.genData(pMotifs, test,  args[6]);
		 }

		
		
	
	}
}
