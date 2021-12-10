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
	public HashMap<String,String> kmerStrand;
	public HashMap<String,Integer> kmerFreque;
	public HashMap<String,Float> kmerMapped;
	public AMRGene()
	{
		kmerStrand = new HashMap<String,String>();
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
						if (c=='T' || c=='t' || c=='U' || c=='u') {k.append('T');}
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
	
	public static Comparator<HashMap.Entry<String,Float>> sortHashMapByValueFloat = new Comparator<HashMap.Entry<String,Float>>()
	{
		@Override
		public int compare(Map.Entry<String,Float> e1, Map.Entry<String,Float> e2)
		{
			Float f1 = e1.getValue();
			Float f2 = e2.getValue();
			return f2.compareTo(f1);
		}
	};
	
	public static void main(String[] args) throws Exception
	{
		long time0 = System.currentTimeMillis();
		long startTime = System.currentTimeMillis();
		long endTime = System.currentTimeMillis();
		long elapsedTime = endTime - startTime;
		final int DEFAULT_BUFFER_SIZE=16384;
		float allram = (float)(Runtime.getRuntime().maxMemory());
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		
		int k = 17;
		int numT = 125000;
		String dbfile="megares_full_database_v2.00.fasta";
		String readfile="";
		boolean classifyReads = true;
		boolean reportMultipleHits = false;
		
		for (int t=0; t<args.length; t++)
		{
			if (args[t].startsWith("d:")) dbfile=args[t].split(":")[1];
			if (args[t].endsWith(".fastq") || args[t].endsWith(".gz")) readfile=args[t];
			if (args[t].startsWith("f:")) readfile=args[t].split(":")[1];
			if (args[t].startsWith("k:")) k=Integer.parseInt(args[t].split(":")[1]);
			if (args[t].startsWith("i:")) numT=Integer.parseInt(args[t].split(":")[1]);
			if (args[t].equals("r:n") || args[t].equals("r:no")) classifyReads = false;
			if (args[t].equals("r:y") || args[t].equals("r:yes")) classifyReads = true;
			if (args[t].equals("m:n") || args[t].equals("m:no")) reportMultipleHits = false;
			if (args[t].equals("m:y") || args[t].equals("m:yes")) reportMultipleHits = true;
		}
		if (k%2==0) k=k+1; if (k<11) {System.out.println("Minimum value of k must be 11"); k=11;}
		if (readfile.equals("")) {System.out.println("Please specify a read file"); System.exit(0);}
		
		System.out.println("Reading AMR gene database, creating k-mer mapping (k="+k+")");
		startTime = System.currentTimeMillis();
		HashMap<String,ArrayList<String>> kmerGeneMapping = new HashMap<String,ArrayList<String>>();
		HashMap<String,AMRGene> geneKmerMapping = new HashMap<String,AMRGene>();
		BufferedReader r = new BufferedReader(new FileReader(dbfile));
		String header = r.readLine();
		long i=0;
		while(true)
		{
			if (!header.startsWith(">")) {System.out.println("Wrong fasta format"); System.exit(0);}
			if (header==null) break;
			String sequence = r.readLine();
			if (sequence==null) break;
			String nextl = r.readLine();
			if (nextl==null) break;
			while(nextl!=null && !nextl.startsWith(">")) {sequence=sequence+nextl; nextl=r.readLine();}
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
						if (!al.contains(header)) al.add(header);
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
						if (!al.contains(header)) al.add(header);
						kmerGeneMapping.put(rk,al);
					}
					if (amrgene.kmerFreque.get(fk)==null) {amrgene.kmerFreque.put(fk,1);} else {amrgene.kmerFreque.put(fk,amrgene.kmerFreque.get(fk)+1);}
					if (amrgene.kmerFreque.get(rk)==null) {amrgene.kmerFreque.put(rk,1);} else {amrgene.kmerFreque.put(rk,amrgene.kmerFreque.get(rk)+1);}
					if (amrgene.kmerStrand.get(fk)==null) {amrgene.kmerStrand.put(fk,rk);}
				}
				geneKmerMapping.put(header,amrgene);
			}
			header=nextl;
			if (nextl==null) break;
			i++;
			if (i%1000==0)
			{
				System.gc();
				allram = (float)(Runtime.getRuntime().maxMemory());
				usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());	
				System.out.println("\t"+i+" genes processed; used RAM = "+100*usedram/allram+"%");
			}
		}
		
		r.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.println(i+" genes read and k-mers mapped in "+elapsedTime/1000+" seconds");
		
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
		String line;
		while((line=r.readLine())!=null || i<numT)
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
		if ( avg<k ) {System.out.println("Avergae read length too short for the chosen k"); System.exit(0);}
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
			if (y%(numT/5)==0) System.out.print(y+"..");
		}
		System.out.println();
		Arrays.sort(matchDist);
		int pvalthres=matchDist[99*numT/100];
		System.out.println("99th percentile of random k-mers match distribution is "+pvalthres+" (max is "+matchDist[numT-1]+")");
		r.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.println("Empirical distribution for "+numT+" random reads estimated in "+elapsedTime/1000+" seconds");
		
		System.out.println("Reading file and mapping genes");
		startTime = System.currentTimeMillis();
		
		String readOutFile = readfile.substring(0,readfile.indexOf("."))+"_KARGA_mappedReads.csv";
		FileWriter rfilewriter = new FileWriter(readOutFile);
		BufferedWriter rwriter = new BufferedWriter(rfilewriter);
		rwriter.write("Idx,");
		rwriter.write("GeneProbability/KmersHitsOnGene/KmersHitsOnAllGenes/KmersTotal,");
		rwriter.write("GeneAnnotation");
		if (reportMultipleHits) rwriter.write("s,...");
		rwriter.write("\r\n");
		
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
			header = line;
			line=r.readLine();
			String fwd = line;
			i++;
			if (line==null) break;
			r.readLine();
			r.readLine();
			fwd = checkAndAmendRead(fwd);
			if (fwd.length()>k)
			{
				ArrayList<String> kmerhits = new ArrayList<String>();
				HashMap<String,Float> geneHitsWeighted = new HashMap<String,Float>();
				HashMap<String,Integer> geneHitsUnweighted = new HashMap<String,Integer>();
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
							if (geneHitsWeighted.get(key)==null) {geneHitsWeighted.put(key,frac);} else {geneHitsWeighted.put(key,geneHitsWeighted.get(key)+frac);}
							if (geneHitsUnweighted.get(key)==null) {geneHitsUnweighted.put(key,1);} else {geneHitsUnweighted.put(key,geneHitsUnweighted.get(key)+1);}
						}
		 			
					}
				}
			
				if (kmerhits.size()>pvalthres)
				{
					if (!reportMultipleHits)
					{
						List<String> keys = new ArrayList<>(geneHitsWeighted.keySet());
						Collections.shuffle(keys);
						float maxGeneFreq = 0;
						String maxGene="";
						for (String key : keys)
						{
							float curr = geneHitsWeighted.get(key);
							if (curr>maxGeneFreq) {maxGeneFreq=curr;maxGene=key;}
						}
						if (classifyReads)
						{
							rwriter.write(header+",");
							float fr = (float)Math.round(maxGeneFreq*100)/100;
							fr = fr/kmerhits.size();
							fr = (float)Math.round(fr*100)/100;
							rwriter.write(fr+"/"+geneHitsUnweighted.get(maxGene)+"/"+kmerhits.size()+"/"+(fwd.length()-k+1)+",");
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
								else {genehit.kmerMapped.put(kh,genehit.kmerMapped.get(kh)+1f);}
							}
						}
					}
					if (reportMultipleHits)
					{
						ArrayList<HashMap.Entry<String,Float>> genehitsarr = new ArrayList<HashMap.Entry<String,Float>>();
						for (HashMap.Entry<String,Float> e: geneHitsWeighted.entrySet()) {genehitsarr.add(e);}
						Collections.sort(genehitsarr,sortHashMapByValueFloat);
						if (classifyReads)
						{
							rwriter.write(header+",");
							float cumul = 0f;
							for (int y=0; y<genehitsarr.size(); y++)
							{
								float fr = genehitsarr.get(y).getValue();
								fr = (float)(fr)/(float)(kmerhits.size());
								float fp = (float)Math.round(fr*100)/100;
								rwriter.write(fp+"/"+geneHitsUnweighted.get(genehitsarr.get(y).getKey())+"/"+kmerhits.size()+"/"+(fwd.length()-k+1)+",");
								rwriter.write(genehitsarr.get(y).getKey());
								cumul = cumul+fr;
								if (y>19 || cumul>0.95f) break;
								rwriter.write(",");
							}
							rwriter.write("\r\n");
						}
						float cumul = 0f;
						for (int y=0; y<genehitsarr.size(); y++)
						{
							if (y<=19 && cumul<=0.95f)
							{
								AMRGene genehit = geneKmerMapping.get(genehitsarr.get(y).getKey());
								float fr = genehitsarr.get(y).getValue();
								fr = (float)(fr)/(float)(kmerhits.size());
								cumul = cumul+fr;
								fr = 1;
								for (int c=0; c<kmerhits.size(); c++)
								{
									String kh = kmerhits.get(c);
									if (genehit.kmerFreque.get(kh)!=null)
									{
										if (genehit.kmerMapped.get(kh)==null) {genehit.kmerMapped.put(kh,fr);}
										else {genehit.kmerMapped.put(kh,genehit.kmerMapped.get(kh)+fr);}
									}
								}
							}
						}
					}
					
				}
				else 
					if (classifyReads)
					{
						rwriter.write(header+",");
						rwriter.write("?/?/?/?,");
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

		FileWriter filewriter = new FileWriter(readfile.substring(0,readfile.indexOf("."))+"_KARGA_mappedGenes.csv");
		BufferedWriter writer = new BufferedWriter(filewriter);
		writer.write("GeneIdx,PercentGeneCovered,AverageKMerDepth\r\n");
		Collection<String> keysc = geneKmerMapping.keySet();
		ArrayList<String> keys = new ArrayList<String>(keysc);
		Collections.sort(keys);
		for (String key : keys)
		{
			AMRGene ag = geneKmerMapping.get(key);
			Set<String> actualKmers = ag.kmerStrand.keySet();
			int totKmers = actualKmers.size();
			double percCovered = 0;
			double kmerDepth = 0;
			for (String fk : actualKmers)
			{
				String rk = ag.kmerStrand.get(fk);
				if (ag.kmerMapped.get(fk)!=null || ag.kmerMapped.get(rk)!=null)
				{
					percCovered = percCovered + 1d;
					double dd = 0;
					if (ag.kmerMapped.get(fk)!=null) dd+=(double)(ag.kmerMapped.get(fk));
					if (ag.kmerMapped.get(rk)!=null) dd+=(double)(ag.kmerMapped.get(rk));
					if (ag.kmerStrand.get(rk)!=null) {dd=dd/(double)(ag.kmerFreque.get(fk)+ag.kmerFreque.get(rk));} else {dd=dd/(double)(ag.kmerFreque.get(fk));}
					kmerDepth = kmerDepth + dd;
				}
			}
			percCovered = percCovered/(double)(totKmers);
			kmerDepth = kmerDepth/(double)(totKmers);
			if (percCovered>0.001f)
			{
				writer.write(key+",");
				writer.write(100*percCovered+"%,");
				writer.write(kmerDepth+"\r\n");
			}
		}
		writer.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.print("Reads and genes mapped in = "+elapsedTime/1000+" s\r\n");
		
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - time0;
		System.out.print("Total time employed  = "+elapsedTime/1000+" s\r\n");
	}
}