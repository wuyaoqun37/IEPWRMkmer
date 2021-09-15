import java.io.BufferedReader;
import java.io.File;  
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.*;
import java.math.*;


public class position {
	

	public static void main ( String[] args) throws FileNotFoundException, IOException
	{
		int n=7;
		for (int kermas=6;kermas<n;kermas++)
		{

			String filepath="F:/";
			String filedataOutpath="F:/"+kermas+".txt";//
			String fileOutpath="F:/";
	                     readfile(filepath,fileOutpath,filedataOutpath,kermas);
		}
	}
	
	  
         public static  boolean readfile(String filepath,String fileOutpath,String filedataOutpath, int kermas) throws FileNotFoundException, IOException {
                File file = new File(filepath);
                Map<Integer, Double> temp =new HashMap();
                Map<Integer, Double> temp1 =new HashMap();
                  if (file.isDirectory()) {
				
				        File[] tempLilstFile =file.listFiles();
				        String s =new String();
				        int n= tempLilstFile.length;
				        double [][]K=new double[n][(int) (Math.pow(4, kermas)+1)];
				        double [][]K2=new double[n][(int) (Math.pow(4, kermas)+1)];
				        int[]K1=new int[(int)Math.pow(4, kermas)];
				        for(int i=0;i<n;i++)
				        {
				        	File f=tempLilstFile[i];
				        	s=readString(f);
				        	 temp= ComputWeightKermasbegin(s,kermas);
				        	 for (Map.Entry<Integer, Double> entry : temp.entrySet()) {
				                 K1[entry.getKey()]=K1[entry.getKey()]+1;
				                 K[i][entry.getKey()]=entry.getValue();
				             }

				        }
				        double [][]P=new double[n][(int) (Math.pow(4, kermas)+1)];
				        P=ComputWeightProbability(K,K1,kermas,n);
				        String dis[][]=new String[n-1][n-1];
				        String dis2[][]=new String[n-1][n-1];
				       dis= DisCalculationMD(P,n,(int) (Math.pow(4, kermas)));
				        SelectK(dis,kermas,n);
				    PrintMegOut(fileOutpath,dis,tempLilstFile,n);   
				}
                return true;
        }
         
         public static void SelectK( String dis [][],int kermas,int n) {
			
        	   double s2=0.0;
		        for(int i = 0;i < dis.length;i++){   
		        	 Double s1 =0.00; 
	     			for(int j = 0;j < dis[i].length;j++){
	     				if(dis[i][j]!=null)
	     				{
	     					s2=Double.parseDouble(dis[i][j])*Math.log(Double.parseDouble(dis[i][j]))+(1.0-Double.parseDouble(dis[i][j]))*Math.log(1.00-  Double.parseDouble(dis[i][j]))+s2;
	     				}
	     			}
		        }
		        DecimalFormat df = new DecimalFormat("#.########");
				  String s3;
				  s3=df.format((-1)*(s2)/n);
				
		
		}
         public static String readString(File file){
             StringBuilder result = new StringBuilder();
             String sequence=new String();
             try{
                 BufferedReader br = new BufferedReader(new FileReader(file));
                 String s = null;
                 while((s = br.readLine())!=null){
                if (!s.startsWith(">")) {
                	 result.append(s.trim());
                }
                 }
                 sequence =result.toString().toUpperCase();
                 br.close();  
             }catch(Exception e){
                 e.printStackTrace();
             }
			return sequence;
         } 

         public static int[] seqTonum(String sequence) {
        	 
        	 int[] seq =new int[sequence.length()];
        	  for(int i=0;i<sequence.length();i++){
        		 char st = sequence.charAt(i);
        		 switch (Character.toUpperCase(st)){
        		 case 'A':
        			 seq[i] = 0;
        			 break;
        		 case 'C':
        			 seq[i] = 1;
        			 break;
        		 case 'G':
        			 seq[i] = 2;
        			 break;
        		 case 'T':
        			 seq[i] = 3;
        			 break;
        		 default:
        			 seq[i]=0;
        			 break;
        		 }
        	
        	  }
        	
			return seq;
		}
        
         public static Map<Integer,Double> ComputWeightKermasbegin(String sequence, int kermas){
     		int[] strnum=seqTonum(sequence);
     		int countKmerT=strnum.length-kermas+1;
     		int L=strnum.length;
     		Map<Integer,Double> map=new HashMap();
     		Map<Integer,Integer> map1=new HashMap();
     		Map<Integer,Double> map2=new HashMap();

     		for(int i=0;i<countKmerT;i++){
     	    	int m=0;
     	    	for(int j=0;j<kermas;j++){
     	    		m=4*m+strnum[i+j];
     	          }

     	    	if(map.containsKey(m)){
     	    			Double newCount1=(double)(i+1)/L;
    					Double newCount2=map.get(m);
    					int Count1=map1.get(m);
    					Count1=Count1+1;
  
    					map1.put(m, Count1);
    					Double NC=newCount2+newCount1;
  
    					map.put(m,NC);
     	    		 }
     	    	 else {
 
     	    			Double NC=(double)(i+1)/L;
        				map.put(m,NC);
        				int Count=1;
        				map1.put(m, Count);	
     	    		 }
				}

     		map2=ComputFreq(map,map1,countKmerT,L);		
     	    return  map2;     		
     	}


         public static Map<Integer,Double> ComputCumulativeFreq( Map<Integer, Double> map,Map<Integer,Integer> map1, int countKmerT) {
        	    double P=0.0;
        	 for (Integer key : map.keySet()) {
     			  int Count2=map1.get(key);
     			  Double nc= map.get(key)/Count2;
     			  P =P+(double)Count2/countKmerT; 
     			 nc=nc*P;
     			 double nc1=new BigDecimal(nc).setScale(8,BigDecimal.ROUND_HALF_UP).doubleValue();
     			 map.put(key, nc1);
        	 }	
        	 return map;
		}
     
         public static Map<Integer,Double> ComputFreq( Map<Integer, Double> map,Map<Integer,Integer> map1, int countKmerT,int L) {
        	 for (Integer key : map.keySet()) {
     			  Double nc= map.get(key)/countKmerT;
     			 double nc1=new BigDecimal(nc).setScale(8,BigDecimal.ROUND_HALF_UP).doubleValue();
     			 map.put(key, nc1);
        	 }	
        	 return map;
		}
         

         public static double[][] ComputWeightProbability(double K[][],int K1[],int kermas,int n)
         {
		        double H1=0;
		        double H2=0;
		        for(int i=0;i<n;i++)
		        { 
		        	for(int j=0;j<Math.pow(4, kermas);j++)
			        {
                      double F= (double)K1[j]/n;
                      if(F==0 || F==1 )
                      {
                      	H2=0;
                      }
                      else {
                      	 H1= -(F*(Math.log(F)/Math.log(2))+(1-F)*(Math.log(1-F)/Math.log(2)));
                      	 H2=new BigDecimal(H1).setScale(8,BigDecimal.ROUND_HALF_UP).doubleValue();
						}
		        		K[i][j] = (K[i][j])*H2;
			        }
		        }
		        

        	 return K;
         }

         public static  String[][] DisCalculationMD (double K[][],int n,int m){   //N为文件的长度,m为k-mer的个数

        	 String dis[][]=new String[n-1][m-1];
        	 for(int i=0;i<n;i++)
        	 {
				  for(int k=1;k<n-i;k++)
				  { 
        	
					  double d=0;
					  for(int j=0;j<m;j++)
					   {
						    d= Math.abs(K[i][j]-K[k+i][j])+d;
						
					   }
					
					  DecimalFormat df = new DecimalFormat("#.########");
					  String s1;
					  s1=df.format(d);

					  dis[i+k-1][i]=s1;
				  }
			   }
        	
        	 return  dis;
		}
         
    
         private static String big2(double d) {
             BigDecimal d1 = new BigDecimal(Double.toString(d));
             BigDecimal d2 = new BigDecimal(Integer.toString(1));
             return d1.divide(d2,8,BigDecimal.ROUND_HALF_UP).toString();
         }
         

         public static boolean PrintOut(String filedataOutpath, double dis [][], File[]tempLilstFile,int n)
        	 {
				try {
					 PrintWriter fw =new  PrintWriter(new File(filedataOutpath ));
					 for(int k=0;k<n;k++)
	        		 {
	        			 String s =tempLilstFile[k].getName();
	        			
	        		 }
					 for(int i = 0;i < dis.length;i++){   
						
						 for(int j = 0;j < dis[i].length-1;j++){
			     				if(dis[i][j]!=0)
			     				{
			     					 DecimalFormat df = new DecimalFormat("#.########");
			   					     String s1;
			   					     s1=df.format(dis[i][j]);
			   					
			     					fw.print(s1);
			     					fw.print(" ");
			     				}
			     				else {
			     					fw.print("0.00000000");
			     					fw.print(" ");
								}
			     			}
			     			fw.println("");
			     			fw.flush();
					 }
					 fw.close();
				} catch (Exception e) {
					// TODO: handle exception
				}
        	 
        	 return true;
        	 
         }
    
         public static boolean PrintMegOut(String fileOutpath, String dis [][], File[]tempLilstFile,int n){
        	 try {
        		 PrintWriter  fw = new PrintWriter(new File(fileOutpath));
        		 fw.println("#mega");
        		 fw.println("!TITLE;");
        		 for(int k=0;k<n;k++)
        		 {
        			 String s =tempLilstFile[k].getName();
        			 fw.println("#"+s.substring(0, s.indexOf(".")));
        		 }
				for(int i = 0;i < dis.length;i++){   
	     			for(int j = 0;j < dis[i].length;j++){
	     				if(dis[i][j]!=null)
	     				{
	     					fw.print(dis[i][j]);
	     					fw.print(" ");
	     				}
	     			}
	     			fw.println("");
	     			fw.flush();
	     			
				}
				fw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        	 return true; 
         }
     
}
