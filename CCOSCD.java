package CCOSCD;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.management.RuntimeErrorException;

import java.util.Map.Entry;

public class CCOSCD {
	UndirectedUnweightedGraph g;
	public double[] betas;
	public double alpha;
	public String outputPath;
	public int iteratioNumToStartMerge;
	public int maxIterationsToRun;
	public int percentageOfStableNodes;
	public String pathToGraph;
	public SCDGraphMetaData ORIGINALmetaData;
	public SCDGraphMetaData metaData;
	
	public CCOSCD(String pathToGraph, double[]betas, double alpha, String outputPath, int iteratioNumToStartMerge, int maxIterationsToRun, int percentageOfStableNodes, int firstPartMode) throws IOException{
		this.percentageOfStableNodes= percentageOfStableNodes;
		this.betas= betas;
		this.alpha = alpha;
		this.outputPath =outputPath;
		this.iteratioNumToStartMerge = iteratioNumToStartMerge;
		this.maxIterationsToRun = maxIterationsToRun;
		this.pathToGraph = pathToGraph;
		this.g = new UndirectedUnweightedGraph(Paths.get(pathToGraph));	
		Map<Integer, Set<Integer>> firstPart;
		if (firstPartMode == 0){
			firstPart = GetFirstPartition(g);
		}
		else if (firstPartMode == 3){
			firstPart = GetFirstPartitionCliques3(g);
		}
		else if (firstPartMode == 4){
			firstPart = GetFirstPartitionCliques4(g);
		}
		else{
			throw new RuntimeException("param firstPartMode must be on of 0=CC, 3=clique 3, 4=clique 4");
		}
		
		this.ORIGINALmetaData = new SCDGraphMetaData(g,firstPart);
		this.metaData = this.ORIGINALmetaData; 
	}
	
	public CCOSCD(String pathToGraph, String pathToPartition, double[]betas, double alpha, String outputPath, int iteratioNumToStartMerge, int maxIterationsToRun, int percentageOfStableNodes) throws IOException{
		this.percentageOfStableNodes= percentageOfStableNodes;
		this.betas= betas;
		this.alpha = alpha;
		this.outputPath =outputPath;
		this.iteratioNumToStartMerge = iteratioNumToStartMerge;
		this.maxIterationsToRun = maxIterationsToRun;
		this.pathToGraph = pathToGraph;
		this.g = new UndirectedUnweightedGraph(Paths.get(pathToGraph));
		Map<Integer, Set<Integer>> firstPart = GetPartitionFromFile(pathToPartition);	
		System.out.println(firstPart.entrySet());
		this.ORIGINALmetaData = new SCDGraphMetaData(g,firstPart,true);
		this.metaData = this.ORIGINALmetaData; 
	}
	
	public void FindCommunities() throws IOException{
		for (double betta : betas){
			System.out.println("");
			System.out.println("                       Input: " + pathToGraph);
			System.out.println("                       betta: " + betta);
			// Create a copy of the original meta data
			metaData = new SCDGraphMetaData(ORIGINALmetaData);			
			Map<Integer,Set<Integer>> comms = FindCommunities(betta);
			
			WriteToFile(comms, betta);
		}
	}
	
	private Map<Integer,Set<Integer>> FindCommunities(double betta) throws FileNotFoundException, UnsupportedEncodingException {
		int numOfIterations = 0;
    	int amountOfDone = 0;
    	boolean haveMergedComms = false; 	
    	    	
		int nonEmptyComms = metaData.com2nodes.keySet().size();
    	while(numOfIterations<2 ||( amountOfDone < nonEmptyComms && numOfIterations < maxIterationsToRun)){
    		
    		System.out.println(""); 
    		System.out.println("numOfIterations: " + numOfIterations);
    		System.out.println("amount of communities: "+ nonEmptyComms);
    		System.out.print("    progress in current iteration: ");
    		nonEmptyComms = 0;
    		amountOfDone = 0;
    		numOfIterations++;
    		// Go over all comms
    		List<Integer> comms = new ArrayList(metaData.com2nodes.keySet());
    		System.out.println("amount of communities: "+ comms.size());
    		int n = comms.size();
    	    int tenPercent = n/10+1;
    	    int commCounter=0;
    		for(Integer comm : comms){
    			commCounter++;
    		
    			if(metaData.com2nodes.get(comm)== null || metaData.com2nodes.get(comm).size()<1){
    				continue;
    			}
    			nonEmptyComms++;
    			if ((commCounter%tenPercent) == 0){
	        		System.out.print(commCounter/tenPercent*10 + "%  ");
	        	}
    			
        		Map<Integer, Double> nodes_inc = FindNewNodes(metaData.g.neighbors, metaData.g.size(), comm, betta);
        		
        		if(nodes_inc.size()==0)
        			continue;
        		double maxInc = Collections.max(nodes_inc.values());
        		List<Integer> newNodes =FilterBestNodes(nodes_inc, betta);
        		boolean nodesHaveChanged = CheckIfOtherNodesAreToBeRemovedFromComm(newNodes, comm, metaData.g.neighbors, metaData.g.size(), betta, maxInc);
        		
    			Map<Integer,Map<Integer, Double>> commsCouplesIntersectionRatio = AddNodesToComm(comm, newNodes, alpha);    			
    			
    			if(numOfIterations>iteratioNumToStartMerge){
    				haveMergedComms = FindAndMergeComms(commsCouplesIntersectionRatio);   			
    			}    			
	            
    			if (!haveMergedComms && !nodesHaveChanged ){
    				amountOfDone++;
    			}
    			    			
    		}
    	}
    	WriteNode2Comm(betta);
    	System.out.println("numOfIterations: " + numOfIterations);
		return metaData.com2nodes;
	}
	
	private Map<Integer,Map<Integer, Double>> AddNodesToComm( 
			Integer comm,
			List<Integer> newNodes, double alpha) {
		Map<Integer,Map<Integer, Double>> ans = new HashMap<>();
		Set<Integer> commNodes = metaData.com2nodes.get(comm);
		for (Integer node : newNodes){
			//add to community
			commNodes.add(node);
			// add to nodes list
			metaData.node2coms.get(node).add(comm);
						
			for(Integer otherComm : metaData.node2coms.get(node)){
				if(otherComm!=comm){
					Integer highComm = comm;
					Integer lowComm = otherComm;
					if(highComm<lowComm){
						highComm = otherComm;
						lowComm = comm;
					}
					Set<Integer> comm1 = commNodes;
				    Set<Integer> comm2 = metaData.com2nodes.get(otherComm);				
		            int intersectionSize = Utills.IntersectionSize(comm1,comm2);
		            if (intersectionSize < 3){
		            	intersectionSize= 0;
		            }	            
		            
		            double intersectionRatio = (double)intersectionSize/(double)Math.min(comm1.size(), comm2.size());
		            if (intersectionRatio > alpha){
			            Map<Integer,Double> lowMap = ans.getOrDefault(lowComm, new HashMap<Integer,Double>());
			            lowMap.put(highComm,intersectionRatio);
			            ans.put(lowComm, lowMap);
		            }
				}
			}
		}	
		
		return ans;
	}
		
	private boolean CheckIfOtherNodesAreToBeRemovedFromComm(List<Integer> newNodes, Integer comm,			
			Map<Integer, Set<Integer>> a, double e, double betta, double maxInc) {
		boolean ans = false;
		List<Integer> nodesToRemove = new ArrayList<Integer>();
		if(metaData.com2nodes.get(comm).size()<3){
			return false;
		}
		//go over all nodes in comm
		for(Integer node :metaData.com2nodes.get(comm)){
			//only check for other nodes
			if(!newNodes.contains(node)){
				double wcc= Calc_WCC(comm, node);
				if (wcc*betta<maxInc){
					ans = true;
					nodesToRemove.add(node);					
				}
			}
		}
		for(Integer node : nodesToRemove){
			metaData.RemoveCommForNode(node, comm);
		}
		return ans;
	}	
		
	private Map<Integer,Double> FindNewNodes( Map<Integer, Set<Integer>> a, double e,
			int comm, double betta) {
		Map<Integer,Double> nodes_inc = new HashMap<Integer,Double>(); 
		// Go over all border nodes
		List<Integer> borderNodes = FindBorderNodes(comm, a);
        if (borderNodes.size() == 0){
        	return new HashMap<Integer,Double>();
        }
        for (int borderNode : borderNodes){
        	double inc= Calc_WCC(comm, borderNode);            
        	nodes_inc.put(borderNode, inc);
        }		
        return nodes_inc;
	}
	
	private List<Integer> FindBorderNodes(
			Integer comm, 
			Map<Integer, Set<Integer>> a
			) {
		Set<Integer> ans = new HashSet<>();
		Set<Integer> nodes = metaData.com2nodes.get(comm);
		for(Integer node:nodes){
			Set<Integer> neis = a.get(node);
			for (int nei : neis){
				if(!nodes.contains(nei)){
					ans.add(nei);
				}
			}
		}
		List<Integer> ret = new ArrayList<>();
		ret.addAll(ans);
		return ret;
	}
	
	private List<Integer> FilterBestNodes(
			Map<Integer, Double> nodes_inc, double betta) {
		List<Integer> ans = new ArrayList<Integer>(); 
		double max = Collections.max(nodes_inc.values());
		if (max <=0){
			return ans;
		}
		
		for (Integer k : nodes_inc.keySet()){
			if(nodes_inc.get(k)*betta >= max){
				ans.add(k);
			}
		}
		return ans;
	}
		
	
	private void WriteNode2Comm(double betta) throws FileNotFoundException, UnsupportedEncodingException {
		PrintWriter writer = new PrintWriter(outputPath + betta + "node2Comm.txt", "UTF-8");
		PrintWriter writer1 = new PrintWriter(outputPath + betta + "node2Comm1.clu", "UTF-8");
		int numOfNodes= metaData.node2coms.size();
		//BuildNode2CommsFromComm2Nodes();
		int maxNumOfComms = numOfNodes;
		writer1.println("*Vertices " + numOfNodes); 
		
		for ( Integer node: metaData.node2coms.keySet()){
			if (metaData.node2coms.get((Integer)node).size()>0){
				writer.print(node + " " + metaData.node2coms.get(node).toArray()[0]);
				writer.println("");	
				writer1.println(metaData.node2coms.get(node).toArray()[0]);
			}
			else{				
				writer.print(node + " " + maxNumOfComms);
				writer.println("");	
				writer1.println(maxNumOfComms);
				maxNumOfComms++;
			}
		}		
		writer.close();
		writer1.close();	
		
	}

	private void BuildNode2CommsFromComm2Nodes() {
		for(Integer node : metaData.node2coms.keySet()){
			metaData.node2coms.put(node, new HashSet<>());			
		}
		
		for(Integer commId : metaData.com2nodes.keySet()){
			for(Integer node : metaData.com2nodes.get((Integer)commId)){
				metaData.node2coms.get((Integer)node).add(commId);
			}			
		}		
	}

	private boolean FindAndMergeComms (Map<Integer,Map<Integer, Double>> commsCouplesIntersectionRatio){
	    boolean haveMergedComms = false;
	    //Set<Integer> commsToClean = new HashSet<Integer>();
	    for (Entry<Integer,Map<Integer, Double> > c1c2intersectionRate : commsCouplesIntersectionRatio.entrySet()){	    	
	        Integer c1 = c1c2intersectionRate.getKey();
	        Map<Integer, Double> values = c1c2intersectionRate.getValue();
	        for(Entry<Integer, Double > c2intersection: values.entrySet()){	        	
		    	if(c2intersection.getValue()>alpha){
		        	Integer[] c1c2 = new Integer[]{c1,c2intersection.getKey()};
		        	//commsToClean.add(c1c2[0]);
		        	MergeComms(c1c2);
		        	haveMergedComms = true;
		    	}
	        }
	    }
	    //ClearNodesFromComms(commsToClean);
	    return haveMergedComms;
	}

	private void MergeComms(Integer[] commsToMerge){
		Integer c1 = commsToMerge[0];
		Integer c2 = commsToMerge[1];
		List<Integer> copyOfC1= new ArrayList<>(metaData.com2nodes.get(c1));
		List<Integer> copyOfC2= new ArrayList<>(metaData.com2nodes.get(c2));
	    for (Integer node : copyOfC1){
	    	metaData.RemoveCommForNode(node,c1);
	        if(!copyOfC2.contains(node)){
	        	metaData.AddCommForNode(node,c2);
	        }	        
	    }
	}
	
	private Set<Integer> Keep_Best_Communities(Map<Integer, Double>comms_imps,double betta){
	    double bestImp = 0;
	    for( double imp : comms_imps.values()){
	    	bestImp = Math.max(bestImp, imp);
	    }
	    Set<Integer> bestComs = new HashSet<Integer>();
	    for(Entry<Integer, Double> entry: comms_imps.entrySet()){
	    		 if (entry.getValue() == bestImp){
	    				 bestComs.add(entry.getKey());
	    		 }
	    }
	    return bestComs;
	}	

	private Set<Integer> Find_Neighbor_Comms(Integer node){
	    Set<Integer>neighborComms = new HashSet<Integer>();
	    for (Integer neighbor : g.neighbors(node)){
	        neighborComms.addAll(metaData.node2coms.get(neighbor));
	    }
    return neighborComms;
    }
	
	private void WriteToFile(Map<Integer, Set<Integer>> comms, double betta) throws FileNotFoundException, UnsupportedEncodingException {
		PrintWriter writer = new PrintWriter(outputPath + betta + ".txt", "UTF-8");
		for ( Set<Integer> listOfNodes : comms.values()){
			if(listOfNodes.size()>0){
				for(int node : listOfNodes){
					writer.print(node + " ");
				}
				writer.println("");
			}
		}		
		writer.close();	
	}

	public static Map<Integer,Set<Integer>> GetFirstPartition(UndirectedUnweightedGraph G){
		Map<Integer,Set<Integer>> result = new HashMap<>();
		Map<Integer, Double> CC = G.Clustring();		
	    Map<Integer, Double> sorted_CC = MapUtil.sortByValue(CC);
	    double maxSeenSoFar=1.0;    
	    boolean[] isVisited = new boolean[G.maxNodeId()+1];	    
	    int commID=0;	    
	    for (int v : sorted_CC.keySet()){
	    	if(maxSeenSoFar<CC.get(v)){
	    		throw(new RuntimeException(String.format("sortedCC was not sorted. node: %1$d.", v)));
	    	}	    	
	        if (!isVisited[v]){
	            isVisited[v]= true;
	            Set<Integer> vSet = new HashSet<>();
	            vSet.add(v);
	            result.put(commID, vSet);
	            for(int  neigh : G.neighbors(v)){
	            	if (!isVisited[neigh]){
	            		isVisited[neigh]= true;
	                    result.get(commID).add(neigh);
	                }
	            }
	            commID+=1;
	        }
	    }
	    
	    return result;
	}
	
	public double OLD_WCC(int comm, int  node){
		int n =(int) g.number_of_nodes(); 
		Set<Integer> commMembers = metaData.com2nodes.get(comm);
		 long t = calcT(commMembers, node);
		 double divesor = commMembers.size() + n*((metaData.T.get(node) - t));
		 if (divesor==0)
		        return 0;
		 return (double)n*t/(double)divesor;
	}
				   
				   
	public double Calc_WCC(int comm, int  x){	    
		Set<Integer> commMembers = metaData.com2nodes.get(comm);
		long TxV = metaData.T.get(x);	    
	    if (TxV==0){
	        return 0;
	    }
	    
		long TxC = calcT(commMembers, x);	    
		if(TxC == 0){
			return 0;
		}
		BigDecimal partA = new BigDecimal(TxC).divide(new BigDecimal(TxV),10, BigDecimal.ROUND_DOWN); 
	    
	    int VTxV = metaData.VT.get(x).size();
	    if(VTxV == 0){
			return 0;
		}
	    int VTxVnoC = calcVTWithoutComm(commMembers, x);	    
	    double divesor = (double)(commMembers.size() +(VTxVnoC));	    
	    if (divesor==0){
	        return 0;
	    }	    
	    BigDecimal partB = new BigDecimal(VTxV).divide(new BigDecimal(divesor),10, BigDecimal.ROUND_DOWN);	
	    double ans = (partA.multiply(partB)).doubleValue();
	    
	    return ans;
	    
	}

	private int calcVTWithoutComm(Set<Integer> commMembers, int node) {		
		Set<Integer> nodesWithTriangle = metaData.VT.get(node);
		return nodesWithTriangle.size() - Utills.IntersectionSize(nodesWithTriangle, commMembers);
	}

	private long calcT(Set<Integer> commMembers, int node) {
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
	
	public static Map<Integer,Set<Integer>> GetFirstPartitionCliques4(UndirectedUnweightedGraph G){
		Set<Integer> hasComm = new HashSet<>();
		boolean vHasComm=false;
		Map<Integer,Set<Integer>> result = new HashMap<>();
			    
	    int commID=0;
	    Set<Integer> nodes = G.nodes();
	    for (int v :nodes){
	    	if(hasComm.contains(v)){
	    		continue;
	    	}
	    	vHasComm = false;
	    	Set<Integer> vNeigh = G.neighbors(v);
	    	for(int u:vNeigh){
	    		if(vHasComm){
	    			break;
	    		}
	    		if(!hasComm.contains(u)){
	    			Set<Integer> UVNeigh= Utills.Intersection(vNeigh, G.neighbors(u));
	    			for(int w:UVNeigh){
	    				if(vHasComm){
	    	    			break;
	    	    		}
	    				if(!hasComm.contains(w)){	
	    					for(int z : Utills.Intersection(UVNeigh, G.neighbors(w))){
	    						if(vHasComm){
	    			    			break;
	    			    		}
	    						if(!hasComm.contains(z)){	    					
			    					Set<Integer> comm = new HashSet<>();
			    					comm.add(v);
			    					comm.add(u);
			    					comm.add(w);
			    					comm.add(z);
			    					result.put(commID, comm);
			    					commID+=1;
			    					hasComm.add(v);
			    					hasComm.add(u);
			    					hasComm.add(w);
			    					hasComm.add(z);
			    					vHasComm = true;
			    					break;
			    				}	    						
	    					}
	    				}
	    			}
	    		}
	    	}
	    	/*if(!vHasComm){
	    		Set<Integer> comm = new HashSet<>();
				comm.add(v);
	    		result.put(commID, comm);
	    		commID+=1;
	    		hasComm.add(v);
	    	}*/
	    }
	    return result;
	}
	
	
	public static Map<Integer,Set<Integer>> GetFirstPartitionCliques3(UndirectedUnweightedGraph G){
		Set<Integer> hasComm = new HashSet<>();
		boolean vHasComm=false;
		Map<Integer,Set<Integer>> result = new HashMap<>();
			    
	    int commID=0;
	    Set<Integer> nodes = G.nodes();
	    for (int v :nodes){
	    	if(hasComm.contains(v)){
	    		continue;
	    	}
	    	vHasComm = false;
	    	Set<Integer> vNeigh = G.neighbors(v);
	    	for(int u:vNeigh){
	    		if(vHasComm){
	    			break;
	    		}
	    		if(!hasComm.contains(u)){
	    			Set<Integer> UVNeigh= Utills.Intersection(vNeigh, G.neighbors(u));
	    			for(int w:UVNeigh){
	    				if(vHasComm){
	    	    			break;
	    	    		}
	    				if(!hasComm.contains(w)){		    					   					
	    					Set<Integer> comm = new HashSet<>();
	    					comm.add(v);
	    					comm.add(u);
	    					comm.add(w);
	    					result.put(commID, comm);
	    					commID+=1;
	    					hasComm.add(v);
	    					hasComm.add(u);
	    					hasComm.add(w);
	    					vHasComm = true;
	    					break;
	    				}
	    			}
	    		}
	    	}
	    	/*if(!vHasComm){
	    		Set<Integer> comm = new HashSet<>();
				comm.add(v);
	    		result.put(commID, comm);
	    		commID+=1;
	    		hasComm.add(v);
	    	}*/
	    }
	    return result;
	}
	
	public Map<Integer,Set<Integer>> GetPartitionFromFile(String partFile) throws IOException{		
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
}





