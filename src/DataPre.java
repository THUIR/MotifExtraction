import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.crypto.Data;


/**
 * Description: Data pre-processing
 */

public class DataPre {
	
	/**
	 * get motif candidates
	 * @param inputDataFile: raw data file
	 * @param satFile: filename + Satisfaction score, 1 for DSAT ones and 2 for SAT ones
	 * @param outputFile: the output motif candidates, in the form of {x,y;}
	 * @param windowSize: length of sliding window, in millisecond
	 */
	public static void simplifyData(String inputDataFile,String satFile,String outputFile,int windowSize){
		try {
			ArrayList<String> SatData = new ArrayList<String>();
			ArrayList<String> DsatData = new ArrayList<String>();
			BufferedReader bReader = new BufferedReader(new FileReader(satFile));
			String line = bReader.readLine();
			while(line != null){
				String []t = line.split("\t");
				if(Double.parseDouble(t[1]) == 1)
					DsatData.add(t[0]);
				else {
					SatData.add(t[0]);
				}
				line = bReader.readLine();
			}
			bReader.close();
			
			File f = new File(inputDataFile);
			File []files = f.listFiles();
			for(int i = 0; i < files.length;i ++){
				bReader = new BufferedReader(new FileReader(files[i]));
				BufferedWriter bWriter;
				
				String filename = files[i].getName();
				ArrayList<DataNode> list = new ArrayList<DataNode>();
				line = bReader.readLine();
				while(line != null){
					String []tStrings = line.split("\t");	
					Pattern p = Pattern.compile(".*ACTION=([a-zA-Z_]*)");
					Matcher m = p.matcher(line);
					if(m.find()){
						String type = m.group(1);
						if(type.equals("OVER")){
							break;
						}
						else if(type.equals("MOUSE_MOVE") || type.equals("SCROLL")){
							DataNode node = new DataNode();
							node.time = Integer.parseInt(tStrings[2].substring(11));
							node.x = Integer.parseInt(tStrings[8].substring(2));
							node.y = Integer.parseInt(tStrings[9].substring(2));
							list.add(node);
						}
					}
					line = bReader.readLine();
				}
				bReader.close();
				
				ArrayList<Candidate> canList = new ArrayList<Candidate>();
				
				int start = 0;
				int end = 0;
				while(end < list.size()){
					for(;end <= list.size();end ++){
						if(end == list.size()){
							Candidate tempCandidate =calCandidate(list, start, end - 1,filename);	
							if (tempCandidate.time.length >= 3){//motifs with a length no longer than 2 are excluded 
								canList.add(tempCandidate);
							}		
							break;
						}
						int t = list.get(end).time;
						if((t - list.get(start).time) > windowSize){
							if(end > (start + 2)){
								Candidate tempCandidate =calCandidate(list, start, end - 1,filename);
								if (tempCandidate.time.length >= 3){
									canList.add(tempCandidate);
								}
							}
							while((t - list.get(start).time) > windowSize){
								start ++;
							}
							break;
						}
					}
				}
				
				if(SatData.contains(filename) || DsatData.contains(filename)){
					bWriter = new BufferedWriter(new FileWriter(outputFile + "//" + filename));
					for(int j = 0;j < canList.size();j ++){
						Candidate t = canList.get(j);
						for(int k = 0;k < t.length; k ++){
							bWriter.write(t.x[k] + "," + t.y[k] + ";");
						}
						if(SatData.contains(filename))
							bWriter.write("\t2\n");
						else {
							bWriter.write("\t1\n");
						}
					}
					bWriter.close();
				}		
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Description: calculate the information of a motif
	 */
	public static Candidate calCandidate(ArrayList<DataNode>list, int start,int end,String ind){
		Candidate c = new Candidate();
		c.index = ind;
		c.length = end - start + 1;
		c.x = new double[c.length];
		c.y = new double[c.length];
		c.time = new int[c.length];
		for (int i = 0;i < c.length; i ++){
			c.x[i] = list.get(i + start).x;
			c.y[i] = list.get(i + start).y;
			c.time[i] = list.get(i + start).time;
		}
		double tx = 0,ty = 0;
		for(int i = 0;i < c.length; i ++){
			tx += c.x[i];
			ty += c.y[i];
		}
		double meanX = tx / c.length;
		double meanY = ty / c.length;
		for(int i = 0;i < c.length;i ++){
			c.x[i] -= meanX;
			c.y[i] -= meanY; 
		}
		return c;
	}
}
