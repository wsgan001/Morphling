/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dataStructures;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.File;
import java.io.FileWriter;
import java.util.*;
import utils.MemoryLogger;
/**
 *
 * @author jiadonglin
 */
public class SequenceDatabase {
    
    private long startTime;
    private long endTime;
    
    /** saves all sequences */
    private final List<Sequence> sequences = new ArrayList<Sequence>();
    List<List<SuperItem>> superitemInput = new ArrayList<>();
    Map<Integer, String> seqChromMap = new HashMap<>();
    
    private double weightRatioThresh = 0.2;     
    
    public void loadSample(String sampleDir){
        File folder = new File(sampleDir);
        File[] superitemFiles = folder.listFiles();
        for(File file : superitemFiles){
            if (file.isFile()){               
                String[] tokens = file.getName().split("\\.");
                String chromName = tokens[4];
                System.out.println("Loading chrom: " + chromName);
                String filePath = file.getAbsolutePath();
                loadSequencesFromFile(filePath);
            }
        }
        databaseFormatter(superitemInput);
        endTime = System.currentTimeMillis();
        
        printDatabaseStats();
    }
    
    
    /**
     * load each superitem from file
     * @param superitemFilePath 
     */
    
    public void loadSequencesFromFile(String superitemFilePath){
//        System.out.println("Start loading mutational database.....");
        
        String thisLine;
        startTime = System.currentTimeMillis();
        MemoryLogger.getInstance().reset();
        
        
        
        BufferedReader myInput = null;
        try {
            
            FileInputStream fin = new FileInputStream(new File(superitemFilePath));
            myInput = new BufferedReader(new InputStreamReader(fin));
            // skip header
            myInput.readLine();
            int curSeq = -1;
            while((thisLine = myInput.readLine()) != null){
                String[] tokens = thisLine.split("\t");   
                String chrom = tokens[1];
                int ChromIdx = Integer.parseInt(chrom);
                if (!seqChromMap.containsKey(ChromIdx)){
                    curSeq = superitemInput.size();
                    superitemInput.add(new ArrayList<>());
                    seqChromMap.put(superitemInput.size() - 1, chrom);                   
                }else{
                    curSeq = ChromIdx;
                }
                double ratio = Double.parseDouble(tokens[7]);
                int splitAlignPos = Integer.parseInt(tokens[4]);
                // if the item is not discordant read-pair, consider the weight ratio threshold
                if (!tokens[0].contains("ARP") && (ratio > 0.1)){
                    if (ratio > 0.2){
                        SuperItem superitem = new SuperItem(tokens);
                        superitem.setChromName(seqChromMap.get(curSeq));
                        superitemInput.get(curSeq).add(superitem);
                    }else if (splitAlignPos != -1){
                        SuperItem superitem = new SuperItem(tokens);
                        superitem.setChromName(seqChromMap.get(curSeq));
                        superitemInput.get(curSeq).add(superitem);
                    }
                                                            
                }else if (tokens[0].contains("ARP") && (ratio > 0.15)){   
//                    if(ratio > 0.15){
//                        SuperItem superitem = new SuperItem(tokens);
//                        superitem.setChromName(seqChromMap.get(curSeq));
//                        superitemInput.get(curSeq).add(superitem);
//                    }else if (splitAlignPos != -1){
//                        SuperItem superitem = new SuperItem(tokens);
//                        superitem.setChromName(seqChromMap.get(curSeq));
//                        superitemInput.get(curSeq).add(superitem);
//                    }
                    SuperItem superitem = new SuperItem(tokens);
                    superitem.setChromName(seqChromMap.get(curSeq));
                    superitemInput.get(curSeq).add(superitem);
                }               
            }
        } catch (IOException e) {
            System.err.println(e);
        }
        databaseFormatter(superitemInput);
        endTime = System.currentTimeMillis();
        
        printDatabaseStats();
       
    }
    /**
     * This is used to re-format the database, in case of <(AB)(C)> happens.
     * @param superitemSeqs 
     */
    private void databaseFormatter(List<List<SuperItem>> superitemSeqs) {       
        System.out.println("\nStart formatting the database...");
//        superitemTypeCount(SIsequences);
        int prePos = 0;
        for (int i = 0; i < superitemSeqs.size(); i ++){
            Sequence sequence = new Sequence(i);
            List<SuperItem> smallIndels = new ArrayList<>();
            List<SuperItem> SIlist = superitemSeqs.get(i);
            Collections.sort(SIlist);
            int sequenceSize = 0;
            List<SuperItem> itemset = new ArrayList<>();
            for (SuperItem si : SIlist){
                
                if (si.isSmallIndel()){
                    smallIndels.add(si);
                    continue;
                } 
                if (itemset.isEmpty()){
                    itemset.add(si);
                }else{
                    sequenceSize += 1;
                    sequence.addItemset(itemset);
                    itemset = new ArrayList<>();  
                    itemset.add(si);
                }
                                              
            }
            sequence.addItemset(itemset);
            sequence.setSmallIndelSuperitem(smallIndels);
            sequences.add(sequence);
//            System.out.println("Sequence " + sequence.getId() + " processed, " + " size: " + sequenceSize);           
        }
    }
    
    public int size(){
        return sequences.size();
    }
    
    public List<Sequence> getSequences(){
        return sequences;
    }
    
    public Sequence getSequenceByID(int seqID){
        return sequences.get(seqID);
    }

   
    
    private void printDatabaseStats(){
        System.out.println("\n============  Sequence Database STATS ==========");
	System.out.println("Number of sequences : " + sequences.size());
        long size = 0;
        for(Sequence sequence : sequences){
            size += sequence.size();
            
        }
        double meansize = ((float)size) / ((float)sequences.size());
        System.out.println("Average sequence size : " + meansize);
        System.out.println("Time: " + (endTime - startTime) + "ms");
        System.out.println("Memory usage: " + MemoryLogger.getInstance().getMaxMemory());
        System.out.println("================================================\n");
        
    }
}
