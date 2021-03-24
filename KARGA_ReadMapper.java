import java.io.*;
import java.io.File;
import java.lang.*;
import java.math.*;
import java.util.*;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class KARGA_ReadMapper
{	
	public static String checkAndAmendRead(String s)
	{
		StringBuffer k = new StringBuffer();
		for (int i=0; i<s.length(); i++)
		{
			char c=s.charAt(i);
			if (c=='A' || c=='a') {k.append('A');}
			else
				if (c=='C' || c=='c') {k.append('C');}
				else
					if (c=='G' || c=='g') {k.append('G');}
					else
						if (c=='T' || c=='t') {k.append('T');}
							else 
								{k.append('N');}
		}
		return k.toString();		
	}
	
	public static String reverseComplement(String s)
	{
		char[] reverse = new char[s.length()];
		for (int i=0; i<s.length(); i++) 
		{
			char c = s.charAt(i);
			if (c=='A') {reverse[(reverse.length-1)-i]='T';}
			else
				if (c=='C') {reverse[(reverse.length-1)-i]='G';}
				else
					if (c=='G') {reverse[(reverse.length-1)-i]='C';}
					else
						if (c=='T') {reverse[(reverse.length-1)-i]='A';}
							else 
								if (c=='N') {reverse[(reverse.length-1)-i]='N';}
		}
		return String.valueOf(reverse);
	}
	
	public static String randomString(int n)
	{
		StringBuffer k = new StringBuffer();
		for (int i=0; i<n; i++)
		{
			double d = Math.random();
			if (d<0.000001d) k.append('N');
				else 
				{
					d = Math.random();
					if (d<0.25d) k.append('A');
						else if (d<0.5d) k.append('C');
							else if (d<0.75d) k.append('G');
								else k.append('T');
				}
		}
		return k.toString();
	}
	
	public static void main(String[] args) throws Exception
	{
		final int DEFAULT_BUFFER_SIZE=16384;
		long time0 = System.currentTimeMillis();
		long startTime = System.currentTimeMillis();
		long endTime = System.currentTimeMillis();
		long elapsedTime = endTime - startTime;
		float allram = (float)(Runtime.getRuntime().maxMemory());
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		
		int k = 17;
		String dbfile="megares_full_database_v2.00.fasta";
		String readfile="";
		
		for (int t=0; t<args.length; t++)
		{
			if (args[t].charAt(0)=='d') dbfile=args[t].split(":")[1];
			if (args[t].endsWith(".fastq") || args[t].endsWith(".gz")) readfile=args[t];
			if (args[t].charAt(0)=='f') readfile=args[t].split(":")[1];
			if (args[t].charAt(0)=='k') k=Integer.parseInt(args[t].split(":")[1]);
		}
		
		if (k%2==0) k=k+1; if (k<11) {System.out.println("Minimum value of k must be 11"); k=11;}
		if (readfile.equals("")) {System.out.println("Please specify a read file"); System.exit(0);}
		
		System.out.println("Reading AMR gene database, creating k-mer mapping (k="+k+")");
		System.out.print("\t");
		startTime = System.currentTimeMillis();
		HashMap<String,ArrayList<String>> kmerGeneMapping = new HashMap<String,ArrayList<String>>();
		BufferedReader r = new BufferedReader(new FileReader(dbfile));
		String line=r.readLine();
		long i=0;
		while(line!=null)
		{
			if (!line.startsWith(">")) {System.out.println("Wrong fasta format"); System.exit(0);}
			String header = line;
			String sequence = "";
			do {line=r.readLine(); sequence=sequence+line; if (line==null) break; if (line.equals("")) break;} while (!line.startsWith(">"));
			if (sequence==null) break;
			if (sequence.equals("")) break;
			if (sequence.length()>=k && header.indexOf("RequiresSNPConfirmation")==-1)
			{
				sequence = checkAndAmendRead(sequence);
				String rwd = reverseComplement(sequence);
				for (int g=0; g<sequence.length()-k+1; g++)
				{
					String fk = sequence.substring(g,g+k);
					ArrayList<String> al = kmerGeneMapping.get(fk);
					if (al==null)
					{
						al = new ArrayList<String>();
						al.add(header);
						kmerGeneMapping.put(fk,al);
					}
					else
					{
						al.add(header);
						kmerGeneMapping.put(fk,al);
					}
					String rk = rwd.substring(sequence.length()-(g+k),sequence.length()-g);
					al = kmerGeneMapping.get(rk);
					if (al==null)
					{
						al = new ArrayList<String>();
						al.add(header);
						kmerGeneMapping.put(rk,al);
					}
					else
					{
						al.add(header);
						kmerGeneMapping.put(rk,al);
					}	 
				}
			}
			i++;
			if (i%1000==0) System.out.print(i+"..");
		}
		r.close();
		
		/* this prints the k-mer gene mapping
		FileWriter fw = new FileWriter("kmerGeneMapping.txt");
		BufferedWriter w = new BufferedWriter(fw);
		Set<String> keys = kmerGeneMapping.keySet();
		for (String key : keys)
		{
			w.write(key+"; ");
			ArrayList<String> hl = kmerGeneMapping.get(key);
			while (hl.size()>0)
			{
				w.write(hl.remove(0)+", ");
			}
			w.write("\r\n");
		}
		w.close();
		*/
		
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.println("\r\n"+i+" genes read and k-mers mapped in "+elapsedTime/1000+" seconds");
		
		System.out.print("Estimating background/random k-mer match distribution");
		startTime = System.currentTimeMillis();
		if(readfile.endsWith(".gz"))
		{
			r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readfile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);
		}
		else
		{
			r=new BufferedReader(new FileReader(readfile),DEFAULT_BUFFER_SIZE);
		}
		i=0;
		double avg=0f;
		while((line=r.readLine())!=null || i<50000)
		{
			line=r.readLine();
			String fwd = line;
			if (fwd==null) break;
			avg=avg+(double)(fwd.length());
			r.readLine();
			r.readLine();
			i++;
		}
		avg=avg/(double)(i);
		System.out.println(" (average read length is "+Math.round(avg)+" bases)");
		int numT = 125000;
		int [] matchDist = new int [numT];
		System.out.print("\t");
		for (int y=0; y<numT; y++)
		{
			String fwd = randomString((int)(avg));
			for (int g=0; g<fwd.length()-k+1; g++)
			{
				String fk = fwd.substring(g,g+k);
				if (kmerGeneMapping.get(fk)!=null) {matchDist[y]=matchDist[y]+1;}
			}
			if (y%25000==0) System.out.print(y+"..");
		}
		System.out.println();
		Arrays.sort(matchDist);
		/*
		System.out.print(matchDist[0]+" "); 
		System.out.print(matchDist[1*numT/1000]+" "); 
		System.out.print(matchDist[10*numT/1000]+" "); 
		System.out.print(matchDist[25*numT/1000]+" "); 
		System.out.print(matchDist[5*numT/100]+" "); System.out.print(matchDist[10*numT/100]+" "); 
		System.out.print(matchDist[20*numT/100]+" "); System.out.print(matchDist[25*numT/100]+" "); 
		System.out.print(matchDist[50*numT/100]+" "); 
		System.out.print(matchDist[75*numT/100]+" "); System.out.print(matchDist[80*numT/100]+" "); 
		System.out.print(matchDist[90*numT/100]+" "); System.out.print(matchDist[95*numT/100]+" ");
		System.out.print(matchDist[975*numT/1000]+" "); 
		System.out.print(matchDist[99*numT/100]+" "); 
		System.out.print(matchDist[999*numT/1000]+" "); 
		System.out.print(matchDist[numT-1]+" "); System.out.println();
		*/
		int pvalthres=matchDist[99*numT/100];
		System.out.println("99th percentile of random k-mers match distribution is "+pvalthres+" (max is "+matchDist[numT-1]+")");
		r.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.println("Empirical distribution for "+numT+" random reads estimated in "+elapsedTime/1000+" seconds");
		
		System.out.println("Reading file and classifying reads");
		startTime = System.currentTimeMillis();
		if(readfile.endsWith(".gz"))
		{
			r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readfile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);
		}
		else
		{
			r=new BufferedReader(new FileReader(readfile),DEFAULT_BUFFER_SIZE);
		}
		FileWriter filewriter = new FileWriter(readfile.substring(0,readfile.indexOf("."))+"_ReadMapperKARGA.csv");
		BufferedWriter writer = new BufferedWriter(filewriter);
		writer.write("idx,type,class,mechanism,group,type_hits/mapped/total,class_hits/mapped/total,mecha_hits/mapped/total,group_hits/mapped/total\r\n");
		i=0;
		while((line=r.readLine())!=null)
		{
			String header = line;
			line=r.readLine();
			String fwd = line;
			i++;
			if (line==null) break;
			r.readLine();
			r.readLine();
			if (fwd.length()>k)
			{
				fwd = checkAndAmendRead(fwd);
				
				HashMap<String,Float> type = new HashMap<String,Float>();
				HashMap<String,Float> clas = new HashMap<String,Float>();
				HashMap<String,Float> mech = new HashMap<String,Float>();
				HashMap<String,Float> grou = new HashMap<String,Float>();
				
				float totKmersMapped = 0f;
				for (int g=0; g<fwd.length()-k+1; g++)
				{
					String fk = fwd.substring(g,g+k);
					ArrayList<String> kmerGenes = kmerGeneMapping.get(fk);
					if (kmerGenes!=null)
					{
						totKmersMapped = totKmersMapped+1f;
						for (int y=0; y<kmerGenes.size(); y++)
						{
							String key = kmerGenes.get(y);
							float frac = 1f/(float)(kmerGenes.size());
							String [] ont = key.split("\\|");
							String ty = ont[1];
							String cl = ont[2];
							String me = ont[3];
							String gr = ont[4];
							if (type.get(ty)==null) {type.put(ty,frac);} else {type.put(ty,type.get(ty)+frac);}
							if (clas.get(cl)==null) {clas.put(cl,frac);} else {clas.put(cl,clas.get(cl)+frac);}
							if (mech.get(me)==null) {mech.put(me,frac);} else {mech.put(me,mech.get(me)+frac);}
							if (grou.get(gr)==null) {grou.put(gr,frac);} else {grou.put(gr,grou.get(gr)+frac);}
						}
		 			
					}
				}
			
				writer.write(header+",");
				if (totKmersMapped>pvalthres)
				{
					Set<String> keys = type.keySet();
					float typeMaxI=0;
					String typeMaxN="";
					for (String key : keys)
					{
						float curr = type.get(key);
						if (curr>typeMaxI) {typeMaxI=curr;typeMaxN=key;}
					}
					keys = clas.keySet();
					float clasMaxI=0;
					String clasMaxN="";
					for (String key : keys)
					{
						float curr = clas.get(key);
						if (curr>clasMaxI) {clasMaxI=curr;clasMaxN=key;}
					}
					keys = mech.keySet();
					float mechMaxI=0;
					String mechMaxN="";
					for (String key : keys)
					{
						float curr = mech.get(key);
						if (curr>mechMaxI) {mechMaxI=curr;mechMaxN=key;}
					}
					keys = grou.keySet();
					float grouMaxI=0;
					String grouMaxN="";
					for (String key : keys)
					{
						float curr = grou.get(key);
						if (curr>grouMaxI) {grouMaxI=curr;grouMaxN=key;}
					}
					writer.write(typeMaxN+",");
					writer.write(clasMaxN+",");
					writer.write(mechMaxN+",");
					writer.write(grouMaxN+",");
					writer.write(Math.round(typeMaxI)+"/"+Math.round(totKmersMapped)+"/"+(fwd.length()-k+1)+",");
					writer.write(Math.round(clasMaxI)+"/"+Math.round(totKmersMapped)+"/"+(fwd.length()-k+1)+",");
					writer.write(Math.round(mechMaxI)+"/"+Math.round(totKmersMapped)+"/"+(fwd.length()-k+1)+",");
					writer.write(Math.round(grouMaxI)+"/"+Math.round(totKmersMapped)+"/"+(fwd.length()-k+1)+"\r\n");
				}
				else {writer.write("?,?,?,?,?,?,?,?\r\n");}
			}
			if (i%100000==0) 
			{
				System.gc();
				allram = (float)(Runtime.getRuntime().maxMemory());
				usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());	
				endTime = System.currentTimeMillis();
				elapsedTime = endTime - startTime;
				System.out.print("\t"+i+" reads processed; used RAM = "+100*usedram/allram+"%; elapsed time = "+elapsedTime/1000+" s \r\n");
			}
		}
		r.close();
		writer.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.print("Read file classified in = "+elapsedTime/1000+" s\r\n");
		
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - time0;
		System.out.print("Total time employed = "+elapsedTime/1000+" s\r\n");
	}
}