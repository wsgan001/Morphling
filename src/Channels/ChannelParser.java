/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Channels;
import htsjdk.samtools.*;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import dataStructures.SuperItem;

/**
 *
 * @author jiadonglin
 */
public class ChannelParser {
    private List<SuperItem> superitems = new ArrayList<>();
    // number of normal reads is not set for some break superitems in previous window, will be set in the following window.
    private List<SuperItem> unSetSuperitems = new ArrayList<>();
    // save abnormal mapped reads - mutation signals 
    private List<MutSignal> mutSignals = new ArrayList<>();
    // Read-pair signals are separated into forward and reverse sub-channels.
    private List<MutSignal> forwardMutSignals = new ArrayList<>();
    private List<MutSignal> reverseMutSignals = new ArrayList<>();
    // some break reads may indicate small indels are saved to a separate channel while dealing with break reads
    private List<MutSignal> smallIndelMutSignals = new ArrayList<>();
    

    private final int maximumDeletionSize = 1000000;
    // number of superitems generated in this channel
    private static int SuperItemCount = 0;


    // Abnormal signal writer
    BufferedWriter writer = null; 
    private boolean isARPchannel = false;
    private final int BufferSize = 2000;
    private String channelName;
    private int rpClusterMaxDist;
    
    public ChannelParser(int maxDist, String name, String signalOutPath){
        rpClusterMaxDist = maxDist;
        this.channelName = name;   
        if (!signalOutPath.equals("")){
            String fileName = signalOutPath + channelName + ".abnormal.signals.txt";
            try {
                writer = new BufferedWriter(new FileWriter(fileName));
                writer.write("qname\tmut_type\tmut_pos\tse_pos\tref\tinsert_size\tcigar\n");
            } catch (IOException e) {
                System.out.println(e);
            }            
        }
        
    }
       
    public void addSignals(MutSignal signal, int fragMean, int readLen){
        String mutSignalType = signal.getMutSignalType();     
        if (writer != null){
            try {
                writer.write(signal.toString());
                writer.newLine();
            } catch (Exception e) {
                e.printStackTrace();
            }            
        }
        if (signal.isARPSignal()){
            isARPchannel = true;            

            if (signal.insertSize <= maximumDeletionSize){
                int forwardDistToLast = 0;
                int reverseDistToLast = 0;

                if (forwardMutSignals.isEmpty() && signal.getMutSignalOri().equals("+")){                
                    forwardMutSignals.add(signal);
                }
                if (reverseMutSignals.isEmpty() && signal.getMutSignalOri().equals("-")){
                    reverseMutSignals.add(signal);
                }
                else{
                    if (signal.getMutSignalOri().equals("+")){
                        MutSignal forwardLastSignal = forwardMutSignals.get(forwardMutSignals.size() - 1);

                        forwardMutSignals.add(signal);
                        forwardDistToLast = signal.getMutPos() - forwardLastSignal.getMutPos();
                    }else{
                        MutSignal reverseLastSignal = reverseMutSignals.get(reverseMutSignals.size() - 1);
                        reverseMutSignals.add(signal);
                        reverseDistToLast = signal.getMutPos() - reverseLastSignal.getMutPos();
                    }

                }
                if ((forwardMutSignals.size() > BufferSize && forwardDistToLast > rpClusterMaxDist) ){                
                    List<SuperItem> createdSuperItems = signalLinearClustering(forwardMutSignals, rpClusterMaxDist);
                    if (! createdSuperItems.isEmpty()){
                        for (SuperItem si : createdSuperItems){
                            addSuperItem(si);
    //                        writeSuperItems(si, superitemWriter);
                        }
                    }                
//                    forwardMutSignals.subList(0, forwardMutSignals.size() - 1).clear();
                    MutSignal lastSignal = forwardMutSignals.get(forwardMutSignals.size() - 1);
                    mutSignalListClear(forwardMutSignals, lastSignal);
                }

                if ((reverseMutSignals.size() > BufferSize && reverseDistToLast > rpClusterMaxDist) ){
                    List<SuperItem> createdSuperItems = signalLinearClustering(reverseMutSignals, rpClusterMaxDist);
                    if (! createdSuperItems.isEmpty()){
                        for (SuperItem si : createdSuperItems){
                            addSuperItem(si);
    //                        writeSuperItems(si, superitemWriter);
                        }
                    }
//                    reverseMutSignals.subList(0, reverseMutSignals.size() - 1).clear();
                    MutSignal lastSignal = reverseMutSignals.get(reverseMutSignals.size() - 1);
                    mutSignalListClear(reverseMutSignals, lastSignal);
                }
            }                        
        } else {
            
            if ((mutSignalType.contains("I") || mutSignalType.contains("D")) && signal.isIsizeNormal()){
                MutSignal lastSignal = signal;
                if (!smallIndelMutSignals.isEmpty()){
                    lastSignal = smallIndelMutSignals.get(smallIndelMutSignals.size() - 1);
                }
                smallIndelMutSignals.add(signal);
                int distToLastSignal = signal.getMutPos() - lastSignal.getMutPos();
                
                if (smallIndelMutSignals.size() > BufferSize && distToLastSignal > 0){                    
                    List<SuperItem> createdSuperItems = signalLinearClustering(smallIndelMutSignals, 0);
                    if (!createdSuperItems.isEmpty()){
                        for(SuperItem si : createdSuperItems){
                            addSuperItem(si);
//                            writeSuperItems(si, superitemWriter);
                        }
                    }
                    // clear the list only keep the last element.
//                    smallIndelMutSignals.subList(0, smallIndelMutSignals.size() - 1).clear();
                    MutSignal restSignal = smallIndelMutSignals.get(smallIndelMutSignals.size() - 1);
                    mutSignalListClear(smallIndelMutSignals, restSignal);
                    
                }
            }else if (!mutSignalType.contains("I") && !mutSignalType.contains("D")){

                MutSignal lastSignal = signal;
                if (! mutSignals.isEmpty()){
                    lastSignal = mutSignals.get(mutSignals.size() - 1);
                }
                
                int distToLastSignal = signal.getMutPos() - lastSignal.getMutPos();
                mutSignals.add(signal);
                if (mutSignals.size() > BufferSize && distToLastSignal > 0 ){
                    List<SuperItem> createdSuperItems = signalLinearClustering(mutSignals, 0);
                    if (!createdSuperItems.isEmpty()){
                        for(SuperItem si : createdSuperItems){
                            addSuperItem(si);
//                            writeSuperItems(si, superitemWriter);
                        }
                    }
//                    mutSignals.subList(0, mutSignals.size() - 1).clear();
                    MutSignal restSignal = mutSignals.get(mutSignals.size() - 1);
                    mutSignalListClear(mutSignals, restSignal);
                }
            }
        }
        
    }
    
    public void mutSignalListClear(List<MutSignal> listToClear, MutSignal signalToAdd){
        listToClear.clear();
        listToClear.add(signalToAdd);
    }
    
    public int getSuperitemCount(){
        return SuperItemCount;
    }
    
    private void addSuperItem(SuperItem superitem){
        SuperItemCount += 1;                                    
        superitems.add(superitem);                        
    }
    
    
    
    /**
     * get number of normal read per base.
     * @param readDepthContainer
     * @param windowStart
     * @param windowSize
     * @param preReadDepthBuffer 
     */
    public void setSuperitemWeightRatio(int[] readDepthContainer, int windowStart, int windowSize, int[] preReadDepthBuffer){
        List<SuperItem> setSuperItem = new ArrayList<>();
        for (SuperItem si : superitems){
            int pos = si.getPos();
            int indexOfArray = pos - windowStart;
            if (indexOfArray < windowSize && indexOfArray >= 0){
                int readDepthAtPos = readDepthContainer[indexOfArray];
                si.setSuperitemReadDepth(readDepthAtPos);                   
                setSuperItem.add(si);   

            }else if (indexOfArray < 0){
                int indexInPreBuffer = preReadDepthBuffer.length + indexOfArray;
                int readDepthAtPos = preReadDepthBuffer[indexInPreBuffer];
                si.setSuperitemReadDepth(readDepthAtPos);                   
                setSuperItem.add(si); 
            }
            else{
                unSetSuperitems.add(si);
            }                                                
        }
        superitems.clear();
        superitems = setSuperItem;        
    }
    
    public void writeSuperItemsInChannel(BufferedWriter superitemWriter) throws IOException{
        
        if (isARPchannel){          
            for (SuperItem si : superitems){

                String strOut = si.toString();
                superitemWriter.write(strOut);
                superitemWriter.newLine();   
                
            }       
            superitems.clear();
            SuperItemCount = 0;
        }else{
           for (SuperItem si : superitems){
                String strOut = si.toString();
                superitemWriter.write(strOut);
                superitemWriter.newLine();            
            }       
            superitems.clear();
            for (SuperItem item : unSetSuperitems){
                superitems.add(item);
            }
            unSetSuperitems.clear();
            SuperItemCount = 0;
        }
        
    }        
    
    public void processFinalSignals(int fragMean, int readLen){

        if (isARPchannel){
            if (!forwardMutSignals.isEmpty()){
                List<SuperItem> createdSuperItems = signalLinearClustering(forwardMutSignals, rpClusterMaxDist);
                if (! createdSuperItems.isEmpty()){
                    for (SuperItem si : createdSuperItems){
                        addSuperItem(si);
//                        writeSuperItems(si, superitemWriter);
                    }
                }
                forwardMutSignals.clear();
            }
            if (!reverseMutSignals.isEmpty()){
                List<SuperItem> createdSuperItems = signalLinearClustering(reverseMutSignals, rpClusterMaxDist);
                if (! createdSuperItems.isEmpty()){
                    for (SuperItem si : createdSuperItems){
                        addSuperItem(si);
//                        writeSuperItems(si, superitemWriter);
                    }
                }
                reverseMutSignals.clear();
            }
        }else{
            if (!mutSignals.isEmpty()){
                List<SuperItem> createdSuperItems = signalLinearClustering(mutSignals, 0);
                if (!createdSuperItems.isEmpty()){
                    for(SuperItem si : createdSuperItems){
                        addSuperItem(si);
//                        writeSuperItems(si, superitemWriter);
                    }
                }
                mutSignals.clear();
            }
            if (!smallIndelMutSignals.isEmpty()){
                List<SuperItem> createdSuperItems = signalLinearClustering(smallIndelMutSignals, 0);
                if (!createdSuperItems.isEmpty()){
                    for(SuperItem si : createdSuperItems){
                        addSuperItem(si);
//                        writeSuperItems(si, superitemWriter);
                    }
                }
                smallIndelMutSignals.clear();
            }
        }
        
    }

    private List<SuperItem> signalLinearClustering(List<MutSignal> signals, int maxDist){
                           
        List<SuperItem> superitemList = new ArrayList<>();
        if (signals.size() < 3){
            return superitemList;
        }
        // sort mutation signals in ascending order
        Collections.sort(signals);
        Iterator<MutSignal> iter = signals.iterator();
        List<List<MutSignal>> clusters = new ArrayList<>();
        List<MutSignal> cluster = new ArrayList<>();
        cluster.add(iter.next());
        clusters.add(cluster);
        
        while (iter.hasNext()){
            MutSignal mutSignal = iter.next();

            List<MutSignal> lastCluster = clusters.get(clusters.size() - 1);
            if (mutSignal.withinDistance(lastCluster.get(lastCluster.size() - 1), maxDist)){
                lastCluster.add(mutSignal);
            }else{
                if (lastCluster.size() < 3){
                    clusters.remove(clusters.size() - 1);
                }else{
                    
                    SuperItem superItem = new SuperItem(lastCluster);
                    
                    superitemList.add(superItem);
                    clusters.remove(clusters.size() - 1);
                }
                List<MutSignal> newCluster = new ArrayList<>();
                newCluster.add(mutSignal);
                clusters.add(newCluster);
            }            
            
        }
        
        if (clusters.get(clusters.size() - 1).size() >= 3) {
            SuperItem superitem = new SuperItem(clusters.get(clusters.size() - 1));
            superitemList.add(superitem);
        }
        
        return superitemList;
    }
    
//    public void channelClean(){
//        forwardMutSignals.clear();
//        reverseMutSignals.clear();
//        smallIndelMutSignals.clear();
//        mutSignals.clear();
//    }
    public int getNumOfMut(){
        return forwardMutSignals.size() + reverseMutSignals.size() + mutSignals.size() + smallIndelMutSignals.size();
    }
    public void setARPSuperItemRatio(int[] readDepthContainer, int windowStart, int windowSize, int[] preReadDepthBuffer){
        for (SuperItem superItem : superitems){
            int superItemIntervalStart = superItem.getSuperitemRegion().start - windowStart - 1;
            int superItemIntervalEnd = superItem.getSuperitemRegion().end - windowStart;
            int nreadSum = 0;
            double superItemCov = 0;
            if (superItemIntervalEnd < windowSize && superItemIntervalStart >= 0){
                 for (int i = superItemIntervalStart; i< superItemIntervalEnd; i++){
                    nreadSum += readDepthContainer[i];
                }
                superItemCov = (double) nreadSum / (superItemIntervalEnd - superItemIntervalStart);
            }else if (superItemIntervalEnd < 0 && superItemIntervalEnd > -1000000){                
                int endIndexPreBuffer = preReadDepthBuffer.length + superItemIntervalEnd;
                int startIndexPreBuffer = preReadDepthBuffer.length + superItemIntervalStart;
                for (int i = startIndexPreBuffer; i < endIndexPreBuffer ; i ++){
                    nreadSum += preReadDepthBuffer[i];
                }
                superItemCov = (double) nreadSum / (superItemIntervalEnd - superItemIntervalStart);
            }else{
                superItemCov = 0;
            }                                        
            superItem.setARPsuperitemRatio(superItemCov);
        }        
    }
}
