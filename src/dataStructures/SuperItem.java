/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dataStructures;

import java.util.*;
import java.util.Map.Entry;
import Channels.MutSignal;
import htsjdk.samtools.QueryInterval;
import java.text.DecimalFormat;
import htsjdk.samtools.util.*;
/**
 *
 * @author jiadonglin
 */
public class SuperItem implements Comparable<SuperItem>{
    
    protected int genomePos;
    protected int weight;
    protected String type;
    protected String ori;
    protected QueryInterval superitemInterval;
    protected QueryInterval superitemMateInterval;

    protected int chromIdx;
    protected String chromName;
    protected int nread = 0;
    protected double weightRatio = 0;
    protected List<byte[]> qnames = new ArrayList<>();
    protected int splitAlignedPos = -1;
    
//    DecimalFormat df = new DecimalFormat("0.000");
    
    public SuperItem(){
        
    }
    
    public SuperItem(String[] tokens){
        type = tokens[0];
        chromIdx = Integer.parseInt(tokens[1]);
        
        nread = Integer.parseInt(tokens[2]);
        genomePos = Integer.parseInt(tokens[3]);
        splitAlignedPos = Integer.parseInt(tokens[4]);
        ori = tokens[5];
        weight = Integer.parseInt(tokens[6]);        
        weightRatio = Double.parseDouble(tokens[7]);
        
        String superitemRegion = tokens[8];
        superitemInterval = decodeIntervalFromString(superitemRegion);
        
        String superitemMateRegion = tokens[9];
        superitemMateInterval = decodeIntervalFromString(superitemMateRegion);
        
        String qNameColumn = tokens[10];
        if (!qNameColumn.equals("*")){
            String[] qNameTokens = qNameColumn.split(",");
            for(String qName : qNameTokens){
                qnames.add(qName.getBytes());
            }
        }
    }
    public SuperItem(List<MutSignal> mutSignals){
        createSuperItem(mutSignals);
    }
    
    @Override
    public int compareTo(SuperItem otherSuperItem){
        return this.genomePos - otherSuperItem.getPos();
    }
    /**
     * Consider the equality of superitem. If all attribute considered, too much identical superitems.
     * @param obj
     * @return 
     */
    public boolean equals(Object obj){
        if (obj instanceof SuperItem){
            SuperItem si = (SuperItem) obj;
            return (si.type.equals(this.type));
        }else{
            return false;
        }
    }
    
    
    public boolean isEqual(SuperItem other){
        return ((type.equals(other.type)) && (genomePos == other.genomePos));
    }
    public boolean isSmallIndel(){
        if (!type.contains("ARP") && (type.contains("I")||type.contains("D"))){
            return true;
        }
        else return false;
    }
    // just for debug useage
    public String toConciseString(){
        StringBuilder sb = new StringBuilder();
        sb.append("id: " + chromIdx);
        sb.append(" t: " + type);
        sb.append(" p: " + genomePos);
        sb.append(" w: " + weight);
//        sb.append( " I: ");
//        sb.append(superitemInterval);
//        sb.append(" mateI: ");
//        sb.append(superitemMateInterval.toString());
        return sb.toString();
    }    
    private String qNamesToString(){
        String str = "*";
        if (!qnames.isEmpty()){
            StringBuilder sb = new StringBuilder();
            for (byte[] array : qnames){
                sb.append(new String(array));
                sb.append(",");
            }
            str = sb.substring(0, sb.length() - 1);
        }        
        return str;
    }
    
    public String toString(){
        StringBuilder sb = new StringBuilder();
        sb.append(type);
        sb.append("\t");
        sb.append(chromIdx);
        sb.append("\t");
        sb.append(nread);
        sb.append("\t");
        sb.append(genomePos);
        sb.append("\t");        
        sb.append(splitAlignedPos);
        sb.append("\t");              
        sb.append(ori);
        sb.append("\t");
        sb.append(weight);
        sb.append("\t");
        sb.append(weightRatio);
        sb.append("\t");
        sb.append(superitemInterval.toString());
        sb.append("\t");
        sb.append(superitemMateInterval.toString());
        sb.append("\t");
        sb.append(qNamesToString());
        return sb.toString();
    }
    @Override
    public int hashCode(){
        return type.hashCode();
    }
    public List<byte[]> getQNames(){
        return qnames;
    }
    public boolean isARPsuperitem(){
        return type.contains("ARP");
    }
    public QueryInterval getSuperitemRegion(){
        return superitemInterval;
    }
    public QueryInterval getSuperitemMateRegion(){
        return superitemMateInterval;
    }
    public double getWeightRatio(){
        return this.weightRatio;
    }
    public String getType(){
        return this.type;
    }
    public int getPos(){
        return this.genomePos;
    }
    public int getWeight(){
        return this.weight;
    }
    public int getChromIdx(){
        return chromIdx;
    }
    public void setChromName(String chrom) {
        chromName = chrom;
    }
    public String getChromName(){
        return chromName;
    }
    public String getOri(){
        return ori;
    }
    public int getSplitAlignPos(){
        return splitAlignedPos;
    }
    
    /**
     * Use the dominant signals
     * @param mutList sorted mutational signals by mutPos
     */
    private void createSuperItem(List<MutSignal> mutList){
        
        Map<MutSignal, Integer> typeCountMap = new HashMap<>();
        
        List<Integer> superitemRecordPos = new ArrayList<>();
        List<Integer> superitemMateRecordPos = new ArrayList<>();
        
        Set<String> queryName = new HashSet<>();
        List<byte[]> qNameByteList = new ArrayList<>();
        
//        List<String> splitAligns = new ArrayList<>();
//        List<Integer> splitAlignPos = new ArrayList<>();
        
        for (MutSignal signal : mutList){              
            if(signal.isSplitAlign()){
//                splitAligns.add(signal.getSplitAlignInfo());
//                splitAlignPos.add(signal.getSplitAlignPos());
                splitAlignedPos = signal.getSplitAlignPos();
            }
            
            String signalQname = signal.getqName();
            superitemRecordPos.add(signal.getRecordPos());
            superitemMateRecordPos.add(signal.getMateRecordPos());
            
            if (!queryName.contains(signalQname)){
                qNameByteList.add(signalQname.getBytes());
                queryName.add(signalQname);
            }
            
            if (!typeCountMap.containsKey(signal)){
                typeCountMap.put(signal, 1);
            }else{
                int count = typeCountMap.get(signal);
                count += 1;
                typeCountMap.put(signal, count);
            }
        }
        // decide the superitem interval and its mate interval.
        Collections.sort(superitemRecordPos);
        Collections.sort(superitemMateRecordPos);
                     
        
        MutSignal majoritySignal = new MutSignal();
        int maxCount = 0;
        for (Entry<MutSignal, Integer> entry : typeCountMap.entrySet()){
            if (entry.getValue() > maxCount){
                maxCount = entry.getValue();
                majoritySignal = entry.getKey();
            }
        }
           
        type = majoritySignal.getMutSignalType();
        
        ori = majoritySignal.getMutSignalOri();
        chromIdx = majoritySignal.getSignalChromIdx();
        chromName = majoritySignal.getSignalRef();
        genomePos = majoritySignal.getMutPos();
        weight = mutList.size();
        
        // If this is a discordant read-pair based superitem, adjust the position and assign query names for further process.
        if (type.contains("ARP")){
            qnames = qNameByteList;
            superitemInterval = new QueryInterval(chromIdx, superitemRecordPos.get(0), superitemRecordPos.get(superitemRecordPos.size() - 1));
            superitemMateInterval = new QueryInterval(chromIdx, superitemMateRecordPos.get(0), superitemMateRecordPos.get(superitemMateRecordPos.size() - 1));
            
            
            if (ori.equals("+")){
                genomePos = mutList.get(mutList.size() - 1).getMutPos() - 1;
            }else{
                genomePos = mutList.get(0).getMutPos() + 1;
            }
        }else{
            qnames = qNameByteList;
            // make sure there is a potential second majority superitem
            if (type.equals("SMS") && typeCountMap.size() > 1){
                refineSuperItemType(typeCountMap);
            }    
            superitemInterval = new QueryInterval(chromIdx, genomePos, genomePos);
            superitemMateInterval = new QueryInterval(chromIdx, superitemMateRecordPos.get(0), superitemMateRecordPos.get(superitemMateRecordPos.size() - 1));                        
        }        
    }   
    
    
    public void setSuperitemReadDepth(int rd){
        nread = rd;
        weightRatio = (double) weight / (weight + nread);
    }
    public void setARPsuperitemRatio(double avgCov){
        nread = (int) avgCov;
        weightRatio = weight / (weight + avgCov);
    }
    /**
     * Refine superitem of type SMS
     * @param signalCountMap 
     */
    private void refineSuperItemType(Map<MutSignal, Integer> signalCountMap){
        MutSignal secondMajoritySignal = new MutSignal();
        int maxCount = 0;
       for(Entry<MutSignal, Integer> entry : signalCountMap.entrySet()){
           String signalType = entry.getKey().getMutSignalType();
           if (!signalType.equals("SMS")){
               if (entry.getValue() > maxCount){
                   maxCount = entry.getValue();
                   secondMajoritySignal = entry.getKey();
               }
           }
       }
       if (secondMajoritySignal != null){
           type = secondMajoritySignal.getMutSignalType();
       }
       
    }
    private QueryInterval decodeIntervalFromString(String str){
        int colonIndex = str.lastIndexOf(':');
        int start = 0;
        int end = 0;
        if (colonIndex == -1){
            start = 1;
            end = Integer.MAX_VALUE;
        }else{
            int dashIndex = str.indexOf('-', colonIndex);
            if (dashIndex == - 1){
                System.err.println("Incorrect interval format !");
            }else{
                start = Integer.parseInt(str.substring(colonIndex + 1, dashIndex).replaceAll(",", ""));
                end = Integer.parseInt(str.substring(dashIndex + 1).replaceAll(",", ""));               
            }
        } 
        return new QueryInterval(chromIdx, start, end);
    }

    
}
