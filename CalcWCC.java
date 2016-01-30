package CCOSCD;

import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class CalcWCC {
	
	public static void main(String[] args) throws IOException {		
		
		/*String pathToGraph = args[0];
		String pathToComms = args[1;
		String betas = ParseDoubleArray(args[2]);
		String outputDirPath = args[3];*/
	
		String pathToGraph = "C:/Users/t-amirub/Desktop/amazon/y/com-amazon.ungraph.txt";
		String pathToCommsDir = "C:/Users/t-amirub/Desktop/amazon/CCOSCD/";
		double[] betas = {1.5,2.0,3.0,5.0,10.0,15.0,20.0};
		String outputDirPath = "C:/Users/t-amirub/Desktop/amazon/CCOSCD/wcc-";
		
		System.out.println("reading graph");
		UndirectedUnweightedGraph g = new UndirectedUnweightedGraph(Paths.get(pathToGraph));
		String pathToComms;
		for(double beta : betas){
			pathToComms = pathToCommsDir+beta+ ".txt";
			System.out.println("reading comms beta: " + beta  );
			Map<Integer, Set<Integer>> firstPart = GetPartitionFromFile(pathToComms);	
			System.out.println("building metaData");
			SCDGraphMetaData metaData = new SCDGraphMetaData(g,firstPart,true);
			System.out.println("calc WCC");
			double wcc = WCC(metaData);
			System.out.println(wcc);
			
			// Write to file
			PrintWriter writer = new PrintWriter(outputDirPath + beta + ".txt", "UTF-8");
			writer.println("WCC " + wcc);
			writer.close();
		}
	}
	
	private static Map<Integer,Set<Integer>> GetPartitionFromFile(String partFile) throws IOException{		
		Map<Integer,Set<Integer>> comm2Nodes= new HashMap<Integer,Set<Integer>>();
		List<String> lines= Files.readAllLines(Paths.get(partFile));		
	    int commID=0;
	    for (String line : lines){
	        String[] nodes = line.split(" |\t");	        		 
	    	if (nodes.length >2){
	    		Set<Integer> comm = new HashSet<>();
	    		for (String node : nodes){
	    			comm.add(Integer.parseInt(node.trim()));
	    		}
	            comm2Nodes.put(commID, comm);
	            commID ++;
	    	}
	    }
	    return comm2Nodes;
	}

	public static double WCC(SCDGraphMetaData metaData){
		double dividor = 0;
		for (Set<Integer> comms : metaData.node2coms.values()){
			dividor = dividor + comms.size();					
		}
		double sum = 0;
		int counter = 0;
		int size = metaData.com2nodes.keySet().size();
		for (Integer commId : metaData.com2nodes.keySet()){
			sum = 	sum + WCCPerComm(commId, metaData);
		}	
		
			
		return sum/(double)dividor;
	}

	private static double WCCPerComm(Integer commId, SCDGraphMetaData metaData) {
		double wcc = 0;
		for (Integer node : metaData.com2nodes.get((Integer)commId)){
						
			wcc = wcc + WCCPerNode(node,commId, metaData);
		}		
		return wcc;
	}
	
	private static double WCCPerNode(Integer x, Integer comm, SCDGraphMetaData metaData) {
		Set<Integer> commMembers = metaData.com2nodes.get((Integer)comm);
		long TxV = metaData.T.get((Integer)x);	    
	    if (TxV==0){
	        return 0;
	    }
	    
		long TxC = calcT(commMembers, x, metaData);	    
		if(TxC == 0){
			return 0;
		}
		BigDecimal partA = new BigDecimal(TxC).divide(new BigDecimal(TxV),10, BigDecimal.ROUND_DOWN); 
	    
	    int VTxV = metaData.VT.get((Integer)x).size();
	    if(VTxV == 0){
			return 0;
		}
	    int VTxVnoC = calcVTWithoutComm(commMembers, x, metaData);	    
	    double divesor = (double)(commMembers.size() +(VTxVnoC));	    
	    if (divesor==0){
	        return 0;
	    }	    
	    BigDecimal partB = new BigDecimal(VTxV).divide(new BigDecimal(divesor),10, BigDecimal.ROUND_DOWN);	
	    double ans = (partA.multiply(partB)).doubleValue();
	    
	    return ans;
	    
	}
	
	private static double[] ParseDoubleArray(String string) {
		String[] parts = string.split(",");
		double[] ans= new double[parts.length];
	    int i=0;
	    for(String str:parts){
	    	ans[i]=Double.parseDouble(str);
	        i++;
	    }
	    return ans;
	}
	
	private static int calcVTWithoutComm(Set<Integer> commMembers, int node, SCDGraphMetaData metaData) {		
		Set<Integer> nodesWithTriangle = metaData.VT.get(node);
		return nodesWithTriangle.size() - Utills.IntersectionSize(nodesWithTriangle, commMembers);
	}
	
	private static long calcT(Set<Integer> commMembers, int node, SCDGraphMetaData metaData) {
		UndirectedUnweightedGraph g = metaData.g;
		long t=0;
	    Set<Integer> neighbours = g.neighbors(node);
	    Set<Integer> neighInComm = Utills.Intersection(commMembers, neighbours);
	    for (int v : neighInComm){
	        for (int u : neighInComm){
	            if (u > v && g.get_edge_weight(u,v)>0){
	                t++;
	            }
	        }
	    }
	    return t;
	}

	
}
