/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dataStructures;

import contiguousfspm.pseudoSuperItem;
import java.util.Map.Entry;
import java.util.*;
/**
 *
 * @author jiadonglin
 */
public class SequentialPattern {
    
    private final List<Itemset> itemsets = new ArrayList<Itemset> ();
    // sequence id
    private int id;
    private int numTypes = 0;
    private List<ItemSeqIdentifier> itemAppear;
    
    /** 
     * Test: alternative way to extract sequential pattern
     * Idea is that, each sequential pattern is a sub-sequence, thus only need to keep
     * the initial prefix start and from which sequence.While adding more items, just 
     * use a variable to keep the length.Having these three we can easily extract the 
     * sub-sequence.
     */
    
    List<Integer[]> initialPrefixTracker;
//    private int initialPrefixStartIdx;
//    private int initialPrefixStartSeqID;
    
    
    public SequentialPattern(int id) {
        this.id = id;
    }
    public int getNumOfTypes(){
        return numTypes;
    }
    public double patternEntropy(){
        Map<String ,Integer> itemTypeCount = new HashMap<>();
        for (Itemset itemset : itemsets){
            for (String itemType : itemset.getItems()){
                if (!itemTypeCount.containsKey(itemType)){
                    numTypes += 1;
                    itemTypeCount.put(itemType, 1);
                }else{
                    int count = itemTypeCount.get(itemType);
                    count += 1;
                    itemTypeCount.put(itemType, count);
                }
            }
        }
        List<Integer> countList = new ArrayList<>();
        for (Entry<String, Integer> entry : itemTypeCount.entrySet()){
            countList.add(entry.getValue());
        }
        double entropy = patternEntropy(countList, length());
        return entropy;
    }
    
    private double patternEntropy(List<Integer> countList, int contentLength){
        double entropy = 0;
        for(int i = 0; i < countList.size(); i++){
            double prob = (double) countList.get(i) / (double) contentLength;
            double val = prob * (Math.log(prob)/Math.log(2));
            entropy -= val;
        }
        if (numTypes == 1){
            entropy = 1;
        }
        return entropy;
    }
    
    public int length(){
        return itemsets.size();
    }
    
    public void addItemset(Itemset itemset){
        itemsets.add(itemset);
    }
    
    public boolean isARPCandidatePattern(){

        boolean hasARPitem = false;
        int arpItemCount = 0;
        int arpOEMcount = 0;        
        for (Itemset itemset: itemsets){
            List<String> items = itemset.getItems();
            for (String item : items){                
                if (item.contains("ARP")){
                    hasARPitem = true;
                    arpItemCount += 1;
                    if (item.contains("OEM")){
                        arpOEMcount += 1;
                    }                    
                }              
            }
        }
        // if a pattern only contains ARP_OEM, discard here
        if (arpOEMcount == arpItemCount){
            hasARPitem = false;
        }
        return hasARPitem;
    }
    
    
    
    public void setItemAppear(List<ItemSeqIdentifier> itemAppearIdx){
        this.itemAppear = itemAppearIdx;
    }
    public List<ItemSeqIdentifier> getItemAppear(){
        return this.itemAppear;
    }
    public int size(){
        return itemsets.size();
    }
    
    public List<Itemset> getItemsets(){
        return itemsets;
    }
    
    public int getSupport(){
        return itemAppear.size();
    }
    public Itemset get(int index) {
        return itemsets.get(index);
    }
    public int getId(){
        return id;
    }
    public SequentialPattern cloneSequence(){
        SequentialPattern newSequence = new SequentialPattern(getId());
        for (Itemset itemset : itemsets){
            newSequence.addItemset(itemset);
        }
        return newSequence;
    }
    
    public String toString(){
        StringBuilder r = new StringBuilder("");
        for(Itemset itemset:itemsets){
            r.append('(');
            for(String item : itemset.getItems()){
                r.append(item);
                r.append(' ');
            }
            r.append(')');
        }
        r.append(" Sup: ");
        r.append(getItemAppear().size());
        return r.toString();
    }
    
    public String patternDetailsString(SequenceDatabase database){
        StringBuilder r = new StringBuilder("");
        
        for(Itemset itemset:itemsets){
            r.append('(');
            for(String item : itemset.getItems()){
                r.append(item);
                r.append(' ');
            }
            r.append(')');
        }
        r.append(" Sup: ");
        r.append(getItemAppear().size());
        r.append("\n");
        
//        for(ItemSeqIdentifier itemIdentity : getItemAppear()){
//            
//            // write the sequence id
//            int seqId = itemIdentity.getSeqID();
//            r.append("seqID: " + seqId);
//            // write the start position of this pattern in the sequence.
//            r.append(' ');
//            int superitemSetStartIdx = itemIdentity.getSubSeqID() - length() + 1;
//            
//            StringBuilder patternInfo = new StringBuilder();
//            for (int i = 0; i < length(); i ++){
//                int superitemSetIdx = superitemSetStartIdx + i;
//                for(int j = 0; j < database.getSequenceByID(seqId).getItemsets().get(superitemSetIdx).size();j ++){
//                    SuperItem curSuperItem = database.getSequenceByID(seqId).superItemAtPos(superitemSetIdx, j);
//                    patternInfo.append('(');
//                    patternInfo.append(curSuperItem.toConciseString());
//                    patternInfo.append(')');
//                    patternInfo.append(' ');
//                }
//                
//            }
//            r.append(patternInfo.toString());                    
//            r.append('\n');
//        }
        return r.append("  ").toString();
    }
    
}
