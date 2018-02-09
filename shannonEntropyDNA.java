/*
 * Author:  Christopher McFall
 * shannonEntropyDNA allows DNA sequence data in the form of a .fasta file to be analyzed for character frequencies, zero order Markov chain probabilities, and single
 * symbol predictability using the Shannon entropy information theoretic.
 */

package shannonEntropyDNA;

import java.io.*;
import java.util.*;
import java.io.File;

public class shannonEntropyDNA {
	public static void main(String[] args) throws java.io.IOException {

		File a = new File("");  //PATH to DNA .fasta file excluding header information
   		FileReader b = new FileReader(a);
   		BufferedReader c = new BufferedReader(b);
   	
   			ArrayList<String> aar = new ArrayList<String>();
		
   				String d;
        	
   				while ((d = c.readLine()) != null) {
                		aar.add(d);
   				}
   					String e = aar.toString();
   					//System.out.println("\n" + e + "\n");
   					String regex = e.replaceAll("[\\s+\\-\\+\\.\\[\\]\\^:,]","");
   					c.close();
		
   						StringBuilder sb = new StringBuilder();
				
   							sb.append(regex);
	      
   								String pukeItUp = sb.toString();
	    
	        			char [] dnascanner = pukeItUp.toCharArray();

	        		String f = String.valueOf(dnascanner);
	        		System.out.println("\n"+"DNA Scanner Output: \n\n" + f);
	        		
	        		Map < Character, Integer > baseMap = new HashMap < Character, Integer > ();

	        			for (int i = 0; i <= 0; i++) {
	        				
							for(char nt:dnascanner) {
					    	    if (baseMap.containsKey(nt)) {
					    	        int ntCount = baseMap.get(nt);
					    	        baseMap.put(nt, ntCount + 1);
					    	    }
					    	    else {
					    	        baseMap.put(nt, 1);
					    	    }
					    	}

					     float charCount = baseMap.get('A') + baseMap.get('T') + baseMap.get('G') + baseMap.get('C');
					     System.out.println("\n" + "Sequence Length: " + charCount);

						    		System.out.print("\nCharacter Frequencies: \n");
							    	for (char nt : baseMap.keySet()) {
							    	    System.out.print("\n" + nt + " = " + baseMap.get(nt));
							    	}

							    	System.out.print("\n\nZero Order Markov Chain Probabilities: \n");
							    	for (char nt : baseMap.keySet()) {
							    		float zeroOrder = baseMap.get(nt)/charCount;
							    	    System.out.print("\n" + nt + " = " + zeroOrder);
							    	}
								
								double shannonEntropy = 0.0;

								System.out.print("\n\nShannon Entropy of Sequence: \n");
								for (char nt : baseMap.keySet()) {
									double Api = -(baseMap.get('A')/charCount * (Math.log (baseMap.get('A')/charCount) / Math.log(2)));
									double Tpi = -(baseMap.get('T')/charCount * (Math.log (baseMap.get('T')/charCount) / Math.log(2)));
									double Gpi = -(baseMap.get('G')/charCount * (Math.log (baseMap.get('G')/charCount) / Math.log(2)));
									double Cpi = -(baseMap.get('C')/charCount * (Math.log (baseMap.get('C')/charCount) / Math.log(2)));
									shannonEntropy = Api + Tpi + Gpi + Cpi;
									
								}
							System.out.print(shannonEntropy + "\n");
					
	        			}
				System.out.print("\n");
		}

}