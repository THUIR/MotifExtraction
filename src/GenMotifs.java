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

import javax.swing.text.StyledEditorKit.BoldAction;

/**
 * Description: extract motifs from data file
 */
public class GenMotifs {
	
	/**
	 * get the minimum of a, b and c
	 */
	public static double get_min(double a,double b,double c){
		if((a <= b) && (a <= c))
			return a;
		if((b <= a) && (b <= c))
			return b;
		return c;
	}
	
	/**
	 * calculate the DTW distance between two motif candidates, based on Euclidean distance
	 */
	public static double DTW_Euclid(Candidate c1,Candidate c2){
		double res = 0;
		double x1[],x2[],y1[],y2[];
		x1 = c1.x;
		x2 = c2.x;
		y1 = c1.y;
		y2 = c2.y;
		double dis[][] = new double[x1.length][x2.length];
		for(int i = 0;i < x1.length;i ++)
			for(int j = 0;j < x2.length;j ++)
				dis[i][j] = Math.sqrt((x1[i] - x2[j]) * (x1[i] - x2[j]) + (y1[i] - y2[j]) * (y1[i] - y2[j]));
		
		double g[][] = new double[x1.length][x2.length];
		g[0][0] = 2 * dis[0][0];
		for(int i = 1;i < x1.length;i ++)
			g[i][0] = g[i - 1][0] + dis[i][0];
		for(int i = 1;i < x2.length;i ++)
			g[0][i] = g[0][i - 1] + dis[0][i];
		for(int i = 1;i < x1.length; i ++)
			for(int j = 1;j < x2.length;j ++){
				g[i][j] = get_min(g[i-1][j-1] + dis[i][j], g[i-1][j], g[i][j-1]) + dis[i][j];
			}
		return g[x1.length - 1][x2.length - 1];
	}

	/**
	 * calculate the DTW distance between motif candidates
	 */
	public static double[][] calDtwDis(ArrayList<Candidate> list){
		double d[][] = new double[list.size()][list.size()];
		for(int i = 0;i < list.size();i ++){
			for(int j = i + 1;j < list.size();j ++){
				d[i][j] = DTW_Euclid(list.get(i),list.get(j));
				d[j][i] = d[i][j];
			}
		}
		return d;
	}
	
	/**
	 * calculate the frequency of a motif
	 * @return
	 */
	public static ArrayList<Candidate> findMotifs(double dtw[][],ArrayList<Candidate> canList,double R){
		ArrayList<Candidate> motifs = new ArrayList<Candidate>();
		for(int i = 0;i < canList.size();i ++){ 
			ArrayList<Integer> similar = new ArrayList<Integer>();
			for(int j = 0;j < canList.size();j ++){
				if(i == j)
					continue;
				if(dtw[i][j] <= R)
					similar.add(j);
			}
			canList.get(i).score = similar.size();
			motifs.add(canList.get(i));
		}
		return motifs;
	}
	
	/**
	 * sort candidates by scores
	 */
	public static void sortMotifsByScore(ArrayList<Candidate> t){
		for(int i = 0;i < t.size();i ++)
			for(int j = i + 1;j < t.size();j ++){
				if(t.get(i).score < t.get(j).score){
					Candidate temp = t.get(i);
					t.set(i, t.get(j));
					t.set(j,temp);
				}
			}
	}
	
	/**
	 * extract motifs based on frequency
	 * @param candidateFile: raw data file
	 * @param num: the num of SAT / DSAT motifs 
	 */
	public static ArrayList<Candidate> genMotifsByFrequency(double R,String candidateFile,int num){
		ArrayList<Candidate> featureList = new ArrayList<Candidate>();
		try {
			File f = new File(candidateFile);
			File []files = f.listFiles();
			ArrayList<Candidate> satList = new ArrayList<Candidate>();
			ArrayList<Candidate> dsatList = new ArrayList<Candidate>();
			for(int i = 0;i < files.length;i ++){
					
				BufferedReader bReader = new BufferedReader(new FileReader(files[i]));
				String line = bReader.readLine();
				while(line != null){
					Candidate t = new Candidate();
					t.index = files[i].getName();
					String []tStrings = line.split("\t");
					String []pos = tStrings[0].split(";");
					t.length = pos.length;
					t.x = new double[t.length];
					t.y = new double[t.length];
					for(int j = 0;j < pos.length;j ++){
						String []xy = pos[j].split(",");
						t.x[j] = Double.parseDouble(xy[0]);
						t.y[j] = Double.parseDouble(xy[1]);
					}
					if(tStrings[1].equals("2")) 
						satList.add(t);
					else if(tStrings[1].equals("1"))
						dsatList.add(t);
					line = bReader.readLine();
				}
				bReader.close();
			}
			
			double satDTW[][] = calDtwDis(satList);
			double dsatDTW[][] = calDtwDis(dsatList);
						
			ArrayList<Candidate> sortedSatMotifs = findMotifs(satDTW, satList, R);
			ArrayList<Candidate> sortedDsatMotifs = findMotifs(dsatDTW, dsatList, R);
			sortMotifsByScore(sortedSatMotifs);
			sortMotifsByScore(sortedDsatMotifs);
			
			featureList = new ArrayList<Candidate>();
			
			for(int i = 0;i < num;i ++){
				featureList.add(sortedSatMotifs.get(i));
			}
			for(int i = 0;i < num;i ++){
				featureList.add(sortedDsatMotifs.get(i));
			}			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return featureList;
	}
	
	/**
	 * extract motifs based on distance
	 * @param candidateFile: raw data file
	 * @param num: the num of SAT / DSAT motifs 
	 */
	public static ArrayList<Candidate> genMotifsByDistance(double R,String candidateFile,int num){
		ArrayList<Candidate> featureList = new ArrayList<Candidate>();
		try {
			
			ArrayList<Candidate> raw_motifs = genMotifsByFrequency(R,candidateFile,num);
			
			ArrayList<Candidate> sortedSatMotifs = new ArrayList<Candidate>();
			ArrayList<Candidate> sortedDsatMotifs = new ArrayList<Candidate>();
			
			for(int k = 0;k < num;k ++){
				sortedSatMotifs.add(raw_motifs.get(k));
				sortedDsatMotifs.add(raw_motifs.get(k + num));
			}
		
			double dtw[][] = new double[sortedSatMotifs.size()][sortedDsatMotifs.size()];
			
			for(int i = 0;i < sortedSatMotifs.size();i ++)
				for(int j = 0;j < sortedDsatMotifs.size();j ++){
					dtw[i][j] = DTW_Euclid(sortedSatMotifs.get(i), sortedDsatMotifs.get(j));
				}
			
			for(int i = 0;i < sortedSatMotifs.size();i ++){
				double total = 0;
				for(int j = 0;j < sortedDsatMotifs.size();j ++){
					total += dtw[i][j];
				}
				sortedSatMotifs.get(i).score = (total / sortedDsatMotifs.size());
			}
			
			for(int i = 0;i < sortedDsatMotifs.size();i ++){
				double total = 0;
				for(int j = 0;j < sortedSatMotifs.size();j ++){
					total += dtw[j][i];
				}
				sortedDsatMotifs.get(i).score = (total / sortedSatMotifs.size());
			}
			
			sortMotifsByScore(sortedSatMotifs);
			sortMotifsByScore(sortedDsatMotifs);
			
			featureList = new ArrayList<Candidate>();
			
			for(int i = 0;i < num;i ++){
				featureList.add(sortedSatMotifs.get(i));
			}
			for(int i = 0;i < num;i ++){
				featureList.add(sortedDsatMotifs.get(i));
			}			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return featureList;
	}
	
	/**
	 * extract motifs based on distribution
	 * @param candidateFile: raw data file
	 * @param num: the num of SAT / DSAT motifs 
	 */
	public static ArrayList<Candidate> genMotifsByDistribution(double R,String candidateFile,int num){
		ArrayList<Candidate> featureList = new ArrayList<Candidate>();
		try {
			
			ArrayList<Candidate> raw_motifs = genMotifsByDistance(R,candidateFile,num);
			ArrayList<Candidate> sortedSatMotifs = new ArrayList<Candidate>();
			ArrayList<Candidate> sortedDsatMotifs = new ArrayList<Candidate>();
			
			for(int k = 0;k < num;k ++){
				sortedSatMotifs.add(raw_motifs.get(k));
				sortedDsatMotifs.add(raw_motifs.get(k + num));
			}
			
			// 2 SAT 1 DSAT
			//(d)satM_(d)sat are files to store middle results for motif selection
			calData(candidateFile, sortedSatMotifs, "satM_sat", num, "2");
			calData(candidateFile, sortedSatMotifs, "satM_dsat", num, "1");
			calData(candidateFile, sortedDsatMotifs, "dsatM_sat", num, "2");
			calData(candidateFile, sortedDsatMotifs, "dsatM_dsat", num, "1");
			
			ArrayList<Candidate> satMotifs = select_Motifs_Distr("satM_sat", "satM_dsat", num,sortedSatMotifs);
			ArrayList<Candidate> dsatMotifs = select_Motifs_Distr("dsatM_sat", "dsatM_dsat", num,sortedDsatMotifs);

			featureList = new ArrayList<Candidate>();
			
			for(int i = 0;i < satMotifs.size();i ++){
				featureList.add(satMotifs.get(i));
			}
			for(int i = 0;i < dsatMotifs.size();i ++){
				featureList.add(dsatMotifs.get(i));
			}		
		} catch (Exception e) {
			e.printStackTrace();
		}
		return featureList;
	}
	
	/**
	 * select motifs based on distribution difference
	 */
	public static ArrayList<Candidate> select_Motifs_Distr(String satDataFile,String dsatDataFile,int num,ArrayList<Candidate> raw_motifs){
		ArrayList<Candidate> motifs = new ArrayList<Candidate>();
		double DSAT[][];
		double SAT[][];
		
		try {
			double mean = 0;
			BufferedReader bReader = new BufferedReader(new FileReader(dsatDataFile));
			String line = bReader.readLine();
			String []tStrings = line.split("\t");
			DSAT = new double[num][tStrings.length];
			int count = 0;
			while(line != null){
				tStrings = line.split("\t");
				for(int i = 0;i < tStrings.length;i ++){
					DSAT[count][i] = Double.parseDouble(tStrings[i]);
					mean += DSAT[count][i];
				}
				count ++;
				line = bReader.readLine();
			}
			bReader.close();
			
			bReader = new BufferedReader(new FileReader(satDataFile));
			line = bReader.readLine();
			tStrings = line.split("\t");
			SAT = new double[num][tStrings.length];
			count = 0;
			while(line != null){
				tStrings = line.split("\t");
				for(int i = 0;i < tStrings.length;i ++){
					SAT[count][i] = Double.parseDouble(tStrings[i]);	
					mean += SAT[count][i];
				}
				count ++;
				line = bReader.readLine();
			}
			bReader.close();
			mean /= (num * (DSAT[0].length + SAT[0].length));
			
			count = 0;
			int sta[][] = new int[num][2];
			for(int i = 0;i < num;i ++)
				for(int j = 0;j < 2;j ++)
					sta[i][j] = 0;
			for(int i = 0;i < num;i ++){
				for(int j = 0;j < SAT[0].length;j++){
					double sim = SAT[i][j] / (3 * mean);
					if(sim < 0.15)
						sta[i][0] ++;
				}
			}
			for(int i = 0;i < num;i ++){
				for(int j = 0;j < DSAT[0].length;j++){
					double sim = DSAT[i][j] / (3 * mean);
					if(sim < 0.15)
						sta[i][1] ++;
				}
			}
			
			HashMap<Integer, Double> goodMotifs = new HashMap<Integer, Double>();
			for(int i = 0;i < num;i ++){
				double rate1 = (double)sta[i][0] / SAT[0].length;
				double rate2 = (double)sta[i][1] / DSAT[0].length;
				if((rate1 > 0.1) && (rate2 == 0)){
						goodMotifs.put(i, 2.0);
						continue;
				}
				if((rate2 > 0.1) && (rate1 == 0)){
					goodMotifs.put(i, 2.0);
					continue;
				}
				if((rate2 > 0.05) || (rate1 >0.05)){
					if(rate1 / rate2 > 1.0)
						goodMotifs.put(i, rate1/rate2);
					else
						goodMotifs.put(i, rate2/rate1);
				}
			}
			count = 0;
			for(int i = 0;i < raw_motifs.size();i ++){
				if(goodMotifs.containsKey(i)){
					Candidate m = raw_motifs.get(i);
					m.score = goodMotifs.get(i);
					motifs.add(m);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}		 
			
		return motifs;
	}
	
	/**
	 * calculate distance information
	 */
	public static void calData(String candidateFile,ArrayList<Candidate> motifs,String outputfile,int motif_num, String fileType) {
		try {
			BufferedWriter bWriter = new BufferedWriter(new FileWriter(outputfile));
			for(int i = 0;i < motif_num;i ++){
				double res = 0;
				File []files = new File(candidateFile).listFiles();
				for(int j = 0;j < files.length;j ++){

					BufferedReader bReader = new BufferedReader(new FileReader(files[j]));
					String line = bReader.readLine();
					String t = line.split("\t")[1];
					if(! t.equals(fileType))
						continue;
					double tdis = 30000000;
					while(line != null){
						Candidate temp = new Candidate();			
						String []ts = line.split("\t")[0].split(";");
						temp.length = ts.length;
						temp.x = new double[temp.length];
						temp.y = new double[temp.length];
						for(int k = 0;k < ts.length;k ++){
							String []xy = ts[k].split(",");
							temp.x[k] = Double.parseDouble(xy[0]);
							temp.y[k] = Double.parseDouble(xy[1]);
						}
						double ttdis = DTW_Euclid(temp, motifs.get(i));
						if( ttdis < tdis){
							tdis = ttdis;
						}
						line = bReader.readLine();
					}
					bReader.close();
					bWriter.write(tdis + "\t");
				}
				bWriter.write("\n");
			}
			bWriter.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}	
	}

	/**
	 * generate files for training and testing
	 * file format: index f1 f2 f3... SatScore;
	 */
	public static void genData(ArrayList<Candidate> motifList,String inputfile,String outputfile){
		try {			
			File f = new File(inputfile);
			File []files = f.listFiles();
			BufferedWriter bWriter = new BufferedWriter(new FileWriter(outputfile));
			for(int i = 0;i < files.length;i ++){
				double []fea = calFeatures(files[i].getAbsolutePath(),motifList);
				bWriter.write(files[i].getName() + "\t");
				for(int j = 0;j < fea.length;j ++){
					bWriter.write(fea[j] + "\t");
				}
				bWriter.write("\n");
			}
			bWriter.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static double[] calFeatures(String dataFile,ArrayList<Candidate> featureList){
		double []f = new double[featureList.size() + 1];
		double s = 0;
		for(int i = 0;i < f.length - 1;i ++){
			double res = 30000000;
			try {		
				Candidate fCandidate = featureList.get(i);
				BufferedReader bReader = new BufferedReader(new FileReader(dataFile));
				String line = bReader.readLine();
				s = Double.parseDouble(line.charAt(line.length() - 1) + "");
				while(line != null){
					Candidate temp = new Candidate();
					String []t = line.split("\t")[0].split(";");
					temp.length = t.length;
					temp.x = new double[temp.length];
					temp.y = new double[temp.length];
					for(int j = 0;j < t.length;j ++){
						String []xy = t[j].split(",");
						temp.x[j] = Double.parseDouble(xy[0]);
						temp.y[j] = Double.parseDouble(xy[1]);
					}
					Double tres;
					tres = DTW_Euclid(temp, fCandidate);
					if(tres < res)
						res = tres;
					line = bReader.readLine();
				}
				bReader.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			f[i] = res;
		}
		f[featureList.size()] = s;
		return f;
	}
	
	/**
	 * get the best K motifs
	 */
	public static ArrayList<Candidate> getTopMotifs(ArrayList<Candidate> raw_motifs, int K){
		sortMotifsByScore(raw_motifs);
		if(raw_motifs.size() > K){
			for(int j = raw_motifs.size() - 1;j >= K ;j --)
				raw_motifs.remove(j);
		}
		return raw_motifs;
	}
}
