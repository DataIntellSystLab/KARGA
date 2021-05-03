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

class AMRGene
{
	public HashMap<String,Integer> kmerFreque;
	public HashMap<String,Float> kmerMapped;
	public AMRGene()
	{
		kmerFreque = new HashMap<String,Integer>();
		kmerMapped = new HashMap<String,Float>();
	}
}

public class KARGA
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
		//String separator = "\\|";
		String readfile="";
		boolean classifyReads = true;
		
		for (int t=0; t<args.length; t++)
		{
			if (args[t].startsWith("d:")) dbfile=args[t].split(":")[1];
			if (args[t].endsWith(".fastq") || args[t].endsWith(".gz")) readfile=args[t];
			if (args[t].startsWith("f:")) readfile=args[t].split(":")[1];
			if (args[t].startsWith("k:")) k=Integer.parseInt(args[t].split(":")[1]);
			//if (args[t].startsWith("s:")) separator=args[t].split(":")[1];
			if (args[t].equals("r:n") || args[t].equals("r:no")) classifyReads = false;
			if (args[t].equals("r:y") || args[t].equals("r:yes")) classifyReads = true;
		}
		
		if (k%2==0) k=k+1; if (k<11) {System.out.println("Minimum value of k must be 11"); k=11;}
		if (readfile.equals("")) {System.out.println("Please specify a read file"); System.exit(0);}
		
		System.out.println("Reading AMR gene database, creating k-mer mapping (k="+k+")");
		startTime = System.currentTimeMillis();
		HashMap<String,ArrayList<String>> kmerGeneMapping = new HashMap<String,ArrayList<String>>();
		HashMap<String,AMRGene> geneKmerMapping = new HashMap<String,AMRGene>();
		BufferedReader r = new BufferedReader(new FileReader(dbfile));
		String line=r.readLine();
		//int annot_num = line.split(separator).length;
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
				AMRGene amrgene = new AMRGene();
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
					if (amrgene.kmerFreque.get(fk)==null) {amrgene.kmerFreque.put(fk,1);} else {amrgene.kmerFreque.put(fk,amrgene.kmerFreque.get(fk)+1);}	 
					if (amrgene.kmerFreque.get(rk)==null) {amrgene.kmerFreque.put(rk,1);} else {amrgene.kmerFreque.put(rk,amrgene.kmerFreque.get(fk)+1);}	 
				}
				geneKmerMapping.put(header,amrgene);
			}
			i++;
			if (i%1000==0) System.out.print(i+"..");
		}
		
		r.close();
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
		int pvalthres=matchDist[99*numT/100];
		System.out.println("99th percentile of random k-mers match distribution is "+pvalthres+" (max is "+matchDist[numT-1]+")");
		r.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.println("Empirical distribution for "+numT+" random reads estimated in "+elapsedTime/1000+" seconds");
		
		System.out.println("Reading file and mapping resistome");
		startTime = System.currentTimeMillis();
		
		String readOutFile = readfile.substring(0,readfile.indexOf("."))+"_KARGA_mappedReads.csv";
		FileWriter rfilewriter = new FileWriter(readOutFile);
		BufferedWriter rwriter = new BufferedWriter(rfilewriter);
		rwriter.write("idx,");
		//for (int d=0; d<annot_num; d++) rwriter.write("annot"+d+",");
		rwriter.write("hits/mapped/total,");
		rwriter.write("annotation\r\n");
		
		if(readfile.endsWith(".gz"))
		{
			r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readfile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);
		}
		else
		{
			r=new BufferedReader(new FileReader(readfile),DEFAULT_BUFFER_SIZE);
		}
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
				
				ArrayList<String> kmerhits = new ArrayList<String>();
				HashMap<String,Float> genehits = new HashMap<String,Float>();
				
				for (int g=0; g<fwd.length()-k+1; g++)
				{
					String fk = fwd.substring(g,g+k);
					ArrayList<String> kmerGenes = kmerGeneMapping.get(fk);
					if (kmerGenes!=null)
					{
						kmerhits.add(fk);
						for (int y=0; y<kmerGenes.size(); y++)
						{
							String key = kmerGenes.get(y);
							float frac = 1f/(float)(kmerGenes.size());
							if (genehits.get(key)==null) {genehits.put(key,frac);} else {genehits.put(key,genehits.get(key)+frac);}
						}
		 			
					}
				}
			
				if (kmerhits.size()>pvalthres)
				{
					//Set<String> keys = genehits.keySet();
					List<String> keys = new ArrayList<>(genehits.keySet());
					Collections.shuffle(keys);
					float maxGeneFreq = 0;
					String maxGene="";
					for (String key : keys)
					{
						float curr = genehits.get(key);
						if (curr>maxGeneFreq) {maxGeneFreq=curr;maxGene=key;}
					}
					if (classifyReads)
					{
						rwriter.write(header+",");
						//String [] out = maxGene.split(separator);
						//for (int d=0; d<out.length; d++)
						//	rwriter.write(out[d]+",");
						rwriter.write(Math.round(maxGeneFreq)+"/"+kmerhits.size()+"/"+(fwd.length()-k+1)+",");
						rwriter.write(maxGene);
						rwriter.write("\r\n");
					}
					AMRGene genehit = geneKmerMapping.get(maxGene);
					for (int y=0; y<kmerhits.size(); y++)
					{
						String kh = kmerhits.get(y);
						if (genehit.kmerFreque.get(kh)!=null)
						{
							if (genehit.kmerMapped.get(kh)==null) {genehit.kmerMapped.put(kh,1f);}
							else {genehit.kmerMapped.put(kh,genehit.kmerMapped.get(kh)+1);}
						}
					}
				}
				else 
					if (classifyReads)
					{
						rwriter.write(header+",");
						//for (int d=0; d<annot_num; d++)
						//	rwriter.write("?,");
						rwriter.write("?/?/?,");
						rwriter.write("?");
						rwriter.write("\r\n");
					}
					
			}
			if (i%100000==0)
			{
				System.gc();
				allram = (float)(Runtime.getRuntime().maxMemory());
				usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());	
				endTime = System.currentTimeMillis();
				elapsedTime = endTime - startTime;
				System.out.print(i+" reads processed; used RAM = "+100*usedram/allram+"%; time = "+elapsedTime/1000+" s \r\n");
			}
		}
		r.close();
		rwriter.close();
		if (!classifyReads) {File f = new File(readOutFile); f.delete();}

		FileWriter filewriter = new FileWriter(readfile.substring(0,readfile.indexOf("."))+"_KARGA_mappedResistome.csv");
		BufferedWriter writer = new BufferedWriter(filewriter);
		writer.write("Gene,PercentGeneCovered,MedianKMerCoverage\r\n");
		Collection<String> keysc = geneKmerMapping.keySet();
		ArrayList<String> keys = new ArrayList<String>(keysc);
		Collections.sort(keys);
		for (String key : keys)
		{
			AMRGene ag = geneKmerMapping.get(key);
			Set<String> rk = ag.kmerFreque.keySet();
			Set<String> mk = ag.kmerMapped.keySet();
			float percCovered = (float)(mk.size())/(float)(rk.size());
			if (percCovered>0.001f)
			{
				writer.write(key+",");
				writer.write(100*percCovered+"%,");
				Collection<Float> kmc1 = ag.kmerMapped.values();
				ArrayList<Float> kmc = new ArrayList<Float>(kmc1);
				Collections.sort(kmc);
				Collection<Integer> krc1 = ag.kmerFreque.values();
				ArrayList<Integer> krc = new ArrayList<Integer>(krc1);
				Collections.sort(krc);
				float mkmc = kmc.get(kmc.size()/2)/(float)(krc.get(krc.size()/2));
				writer.write(mkmc+"\r\n");
			}
		}
		writer.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.print("Reads and resistome mapped in = "+elapsedTime/1000+" s\r\n");
		
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - time0;
		System.out.print("Total time employed  = "+elapsedTime/1000+" s\r\n");
	}
}