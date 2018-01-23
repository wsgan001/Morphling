/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package contiguousfspm;

import htsjdk.samtools.QueryInterval;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.io.FileWriter;
import java.util.Map.Entry;

import dataStructures.*;
import superitemGenerator.svOutInfo;
import utils.MemoryLogger;
import utils.Linkage;


/**
 *
 * @author jiadonglin
 */
public class ContiguousFSPM {
    
    private long startTime;
    private long endTime;
    
    private int minsuppAbsolute;

    BufferedWriter writer = null;
    BufferedWriter idWriter;
    BufferedWriter intMeaWriter;
    
    boolean showPatternIdentifiers = false;
    
    
    /** original sequence count **/
    int sequenceCount = 0;

    private int arpPatternCount = 0;
//    private int nonARPpatternCount = 0;
    private SequentialPatterns patterns = null;
    
    private List<List<pseudoSequentialPattern>> ARPPatternCandidates = new ArrayList<>();
    private List<List<pseudoSequentialPattern>> clippedPatternCandidates = new ArrayList<>();
    
    private int patternSpanMaxRegion;
    private SequenceDatabase database;
    private Map<String, List<ItemSeqIdentifier>> itemAppearMap;
    
    /**
     * Default constructor
     */
    public ContiguousFSPM(int minSup, int maxRegionSpan){
        this.minsuppAbsolute = minSup;
        this.patternSpanMaxRegion = maxRegionSpan;
    }
    
    public SequentialPatterns runAlgorithm(SequenceDatabase database, String outputPath, String superItemPatternOut, String svRegionOut) throws IOException{
        MemoryLogger.getInstance().reset();
        
        this.database = database;
        startTime = System.currentTimeMillis();
        
        prefixSpan(outputPath);

        // doing post-processing of FSPMs
        wgsPatternMerge(superItemPatternOut, svRegionOut);
                        
        if (writer != null){
            writer.close();
        }
        endTime = System.currentTimeMillis();
        return patterns;
    }
    private void prefixSpan(String outputPath) throws IOException{
        // do some post-processing of generated sequential patterns
//        patterns = new SequentialPatterns();
        if (outputPath == null){
            writer = null;            
        }
        // save results in file
        else{
            patterns = null;
            writer = new BufferedWriter(new FileWriter(outputPath));
        }
        // Infomation of single superitem
        System.out.println("Collect information of single superitem .....");
        itemAppearMap = findSequencesContainItems(database);
        
        /**
         * Start creating initial pseudo-projected database.
         */
        System.out.println("Start creating initial pseudo sequence database ....");
        List<PseudoSequence> initialContext = new ArrayList<>();
        for (Sequence sequence : database.getSequences()){
            
            if (sequence.size() != 0){
                initialContext.add(new PseudoSequence(sequence, 0, 0, 0));
            }
        }
        for (Entry<String, List<ItemSeqIdentifier>> entry : itemAppearMap.entrySet()){
            String item = entry.getKey();
            if (entry.getValue().size() >= minsuppAbsolute){
                List<ItemSeqIdentifier> itemAppearIdx = entry.getValue();
                List<PseudoSequence> projectedDatabase = buildProjectedContext(new SequentialPattern(0), item, initialContext, false);
                
                // Create a prefix with initial sequence ID 0
                SequentialPattern prefix = new SequentialPattern(0);
                prefix.addItemset(new Itemset(item));
                prefix.setItemAppear(itemAppearIdx);
                
                depthFirstRecursion(prefix, projectedDatabase);
            }
        }
        
    }
    
    
    private void depthFirstRecursion(SequentialPattern prefix, List<PseudoSequence> database) throws IOException{
        Set<Pair> pairs = itemCountsInProjectedDB(prefix, database);
        for (Pair pair : pairs){
            if (pair.getCount() >= minsuppAbsolute){
                SequentialPattern newPrefix;
                // If the frequent item is of form (_A), append it to the last itemset of the current prefix. 
                if(pair.isPostfix()){
                    newPrefix = appendItemToPrefixLastItemset(prefix, pair.getItem());
                }else{
//                    pair.printItemAppears();
                    newPrefix = appendItemToPrefixSequence(prefix, pair.getItem());
                }
                // Build pseudo-projected database of appended item.
                List<PseudoSequence> projectedDB = buildProjectedContext(newPrefix, pair.getItem(), database, pair.isPostfix());
                newPrefix.setItemAppear(pair.getItemAppear());
                
                double entropy = newPrefix.patternEntropy();
                
                if (entropy >= 1){
                    
                    saveARPPatternCandidate(newPrefix);                      
                    // save patterns in file, not superitems
                    if (newPrefix.isARPCandidatePattern()){
                        savePattern(newPrefix);
                    }
//                    else{
//                        saveClippedPatternCandidate(newPrefix);
//                    }

                }
                depthFirstRecursion(newPrefix, projectedDB);
            }
        }
        MemoryLogger.getInstance().checkMemory();
    }
    
     
   
    
    /**
     * Pair is used to record and separate two conditions:
     * 1) (_A) -> (true, A)
     * 2) (A) -> (false, A)
     * Only count the item at start of the sequence, since the pattern growth in a continuous fashion.
     * @param prefixPattern
     * @param sequences
     * @return a set of pairs
     */
    private Set<Pair> itemCountsInProjectedDB(SequentialPattern prefixPattern, List<PseudoSequence> sequences){       
        Map<Pair, Pair> mapPairs = new HashMap<Pair, Pair>();
        for (PseudoSequence sequence : sequences){
            Sequence oriSequence = this.database.getSequenceByID(sequence.getId());
            for (int j = 0; j < sequence.getSizeOfItemsetAt(0, oriSequence); j++){
                String item = sequence.getItemAtItemsetAt(j, 0, oriSequence).getType();
                Pair paire = new Pair(sequence.isPostfix(0), item);
                Pair oldPaire = mapPairs.get(paire);
                if (oldPaire == null){
                    mapPairs.put(paire, paire);
                }
                /** 
                * same item found, use the previous one. 
                * previous pair object record the item appearance index.
                */
                else{
                    paire = oldPaire;
                }
                // Update item index for each item.
//                paire.addItemAppearIdx(sequence.getId(), 0);
                paire.addItemAppearIdx(sequence.getId(), sequence.getFirstItemsetIdx(), 0, j);
            }
        }
        return mapPairs.keySet();
    }
    
    
    private Map<String, List<ItemSeqIdentifier>> findSequencesContainItems(SequenceDatabase database) {
        Map<String, List<ItemSeqIdentifier>> itemAppearMap = new HashMap<>();
        for(Sequence sequence : database.getSequences()){
            int sequenceID = sequence.getId();
            for(int i = 0; i < sequence.getItemsets().size();i++){
                List<SuperItem> itemset = sequence.getItemsets().get(i);
                for(int j = 0; j < itemset.size(); j ++){
                    SuperItem superitem = itemset.get(j);
                    ItemSeqIdentifier itemIdentity = new ItemSeqIdentifier(sequenceID, sequenceID, i, j);
                    String sitype = superitem.getType();
                    List<ItemSeqIdentifier> itemAppearIdx = itemAppearMap.get(sitype);
                    if(itemAppearIdx == null){
                        itemAppearIdx = new ArrayList<>();
                        itemAppearIdx.add(itemIdentity);
                        itemAppearMap.put(sitype, itemAppearIdx);
                    }else{
                        itemAppearMap.get(sitype).add(itemIdentity);
                    }
                }
            }
        }
        return itemAppearMap;
    }
    

    /**
     * Create a pseudo projected database. 
     * For the tree, root is empty. Level one is the initial projection, where you have to keep all possible suffix string of a give prefix.
     * 
     * @param item superitem type, but has to locate exact object while doing projection
     * @param database a pseudo database
     * @param isSuffix if the item is a suffix or not
     * @return 
     */
 
    private List<PseudoSequence> buildProjectedContext(SequentialPattern prefix, String item, List<PseudoSequence> database, boolean inSuffix){
        List<PseudoSequence> newPseudoProjectedDatabase = new ArrayList<>();
        
        for (PseudoSequence psSequence : database) {
            // simple check if this is level one projection. Single item projection.
            if (psSequence.getOriSeqSize() == psSequence.getSize()){
                for (int i = 0; i < psSequence.getSize();i++){
                    int seqId = psSequence.getId();
                    int itemsetIdx = i;
                    Sequence oriSequence = this.database.getSequenceByID(seqId);
                    
                    int idxOfItemInItemset = psSequence.indexOf(itemsetIdx, item, oriSequence);
                    
                    if (idxOfItemInItemset != -1 && psSequence.isPostfix(itemsetIdx) == inSuffix){
                        SuperItem curSuperItem = oriSequence.superItemAtPos(itemsetIdx, idxOfItemInItemset);                                               
                                        
                        if (idxOfItemInItemset != psSequence.getSizeOfItemsetAt(itemsetIdx, oriSequence) - 1){
                            PseudoSequence newSequence = new PseudoSequence(psSequence, itemsetIdx, idxOfItemInItemset + 1);
                            newSequence.setGenomeStartPos(curSuperItem.getPos());
                            SuperItem nextSuperItem = oriSequence.superItemAtPos(i + 1, 0);
                            boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, newSequence);
                            
//                            if (curSuperItem.isARPsuperitem() && nextSuperItem.isARPsuperitem()){
//                                boolean isLinkable = patternGrowthLinkage(curSuperItem, nextSuperItem);
//                                nextSuperItemInRange = isLinkable;
//                            }
                            if (curSuperItem.isARPsuperitem()){
                                if (i + 3 < psSequence.getSize()){
                                    int deeperSearchRange = i + 3;
                                    for (int k = i + 1; k < deeperSearchRange; k++){
                                        SuperItem si = oriSequence.superItemAtPos(k, 0);
                                        boolean isLinkable = patternGrowthLinkage(curSuperItem, si);
                                        if (isLinkable) {
                                            nextSuperItemInRange = isLinkable;
                                            break;
                                        }
                                    }                                                                        
                                }    
                            }
                                                        
                            if (newSequence.getSize() > 0 && nextSuperItemInRange){
                                newPseudoProjectedDatabase.add(newSequence);
                            }
                            
                        }
                        else if (itemsetIdx != psSequence.getSize() - 1){
                            PseudoSequence newSequence = new PseudoSequence(psSequence, itemsetIdx + 1, 0);
                            newSequence.setGenomeStartPos(curSuperItem.getPos());
                            SuperItem nextSuperItem = oriSequence.superItemAtPos(i + 1, 0);
                            
                            boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, newSequence);
                            
//                            if (curSuperItem.isARPsuperitem() && nextSuperItem.isARPsuperitem()){
//                                boolean isLinkable = patternGrowthLinkage(curSuperItem, nextSuperItem);
//                                nextSuperItemInRange = isLinkable;
//                            }
                            if (curSuperItem.isARPsuperitem()){
                                if (i + 3 < psSequence.getSize()){
                                    int deeperSearchRange = i + 3;
                                    for (int k = i + 1; k < deeperSearchRange; k++){
                                        SuperItem si = oriSequence.superItemAtPos(k, 0);
                                        boolean isLinkable = patternGrowthLinkage(curSuperItem, si);
                                        if (isLinkable){
                                            nextSuperItemInRange = isLinkable;
                                            break;
                                        }
                                    }
                                    
                                }                                                                
                            }
                            if (newSequence.getSize() > 0 && nextSuperItemInRange){
                                newPseudoProjectedDatabase.add(newSequence);
                            }
                            
                        }
                        
                    }
                }         
            }
            // In the deeper level of the tree, only build pseudo-sequence of the start item. Otherwise, the pattern is not consecutive.
            else{
                int seqId = psSequence.getId();
                int itemsetIdx = 0;
                Sequence oriSequence = this.database.getSequenceByID(seqId);

                int idxOfItemInItemset = psSequence.indexOf(itemsetIdx, item, oriSequence);
                
                if (idxOfItemInItemset != -1 && psSequence.isPostfix(itemsetIdx) == inSuffix){
                    SuperItem curSuperItem = oriSequence.superItemAtPos(itemsetIdx, idxOfItemInItemset); 
                    
                    if (idxOfItemInItemset != psSequence.getSizeOfItemsetAt(itemsetIdx, oriSequence) - 1){
                        int nextSuperItemIdx = psSequence.getFirstItemsetIdx() + 1;
                        SuperItem nextSuperItem = oriSequence.superItemAtPos(nextSuperItemIdx, 0);

                        boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, psSequence);
                        
//                        if (curSuperItem.isARPsuperitem() && nextSuperItem.isARPsuperitem()){
//                            boolean isLinkable = patternGrowthLinkage(curSuperItem, nextSuperItem);
//                            nextSuperItemInRange = isLinkable;
//                        }
                        if (curSuperItem.isARPsuperitem()){
                            if (nextSuperItemIdx + 3 < psSequence.getSize()){
                                int deeperSearchRange = nextSuperItemIdx + 3;
                                for (int k = nextSuperItemIdx + 1; k < deeperSearchRange; k++){
                                    SuperItem si = oriSequence.superItemAtPos(k, 0);
                                    boolean isLinkable = patternGrowthLinkage(curSuperItem, si);
                                    if (isLinkable){
                                        nextSuperItemInRange = isLinkable;
                                        break;
                                    }
                                }

                            }                                                                
                        }
                        if (nextSuperItemInRange){
                            PseudoSequence newSequence = new PseudoSequence(psSequence, itemsetIdx, idxOfItemInItemset + 1);
                            if (newSequence.getSize() > 0){
                                newPseudoProjectedDatabase.add(newSequence);
                            }
                        }                 
                    }
                    else if (itemsetIdx != psSequence.getSize() - 1){
                        int nextSuperItemIdx = psSequence.getFirstItemsetIdx() + 1;
                        SuperItem nextSuperItem = oriSequence.superItemAtPos(nextSuperItemIdx, 0);
                        boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, psSequence);

//                        if (curSuperItem.isARPsuperitem() && nextSuperItem.isARPsuperitem()){
//                            boolean isLinkable = patternGrowthLinkage(curSuperItem, nextSuperItem);
//                            nextSuperItemInRange = isLinkable;
//                        }                                                
                        if (curSuperItem.isARPsuperitem()){
                            if (nextSuperItemIdx + 3 < psSequence.getSize()){
                                int deeperSearchRange = nextSuperItemIdx + 3;
                                for (int k = nextSuperItemIdx + 1; k < deeperSearchRange; k++){
                                    SuperItem si = oriSequence.superItemAtPos(k, 0);
                                    boolean isLinkable = patternGrowthLinkage(curSuperItem, si);
                                    if (isLinkable){
                                        nextSuperItemInRange = isLinkable;
                                        break;
                                    }
                                }

                            }                                                                
                        }
                        if (nextSuperItemInRange){
                            PseudoSequence newSequence = new PseudoSequence(psSequence, itemsetIdx + 1, 0);
                            if (newSequence.getSize() > 0){
                                newPseudoProjectedDatabase.add(newSequence);
                            }
                        }                 
                    }                   
                }                     
            }                        
        }
        return newPseudoProjectedDatabase;
    }
    /**
     * Genome start position of current sequence, limit the pattern span region.
     * The idea is that, I will check it while building the pseudo-projection.
     * 1) Once the genome cord of the left most superitem of the new projection minus the genome start cord of 
     * its corresponding initial projection is beyond the max region threshold, then this new pseudo-projection will be discarded.
     * 2) Though the distance calculate in 1) is beyond the threshold, if these two superitems are connected by read-pairs, I wont
     * discard the new appended superitem.
     * @return 
     */
    private boolean ableToBuildProjection(SuperItem nextSuperItem, PseudoSequence psSequence){
        boolean isNextItemInRange = false;
        int genomeStartPos = psSequence.getGenomeStartPos();
        int spannedRegion = nextSuperItem.getPos() - genomeStartPos;
        
        if (spannedRegion <= patternSpanMaxRegion){
            isNextItemInRange = true;
        }
        
        return isNextItemInRange;
    }
    
    private boolean patternGrowthLinkage(SuperItem curSuperItem, SuperItem nextSuperItem){
        Linkage linker = new Linkage();
        return linker.supportARPs(curSuperItem, nextSuperItem);
    }
    
    private SequentialPattern appendItemToPrefixSequence(SequentialPattern prefix, String item){
        SequentialPattern newPrefix = prefix.cloneSequence(); 
        newPrefix.addItemset(new Itemset(item)); 
        return newPrefix;
    }
    
    private SequentialPattern appendItemToPrefixLastItemset(SequentialPattern prefix, String item){
        SequentialPattern newPrefix = prefix.cloneSequence();
        Itemset itemset = newPrefix.get(newPrefix.size() - 1);
        itemset.addItem(item);
        return newPrefix;
    }
    
    private void savePattern(SequentialPattern prefix) throws IOException{
        // increase the pattern count
//        patternCount ++;
        if (writer != null){
            String outStr = prefix.patternDetailsString(database);
            writer.write(outStr);
            writer.newLine();
        }
    }
    
    private void saveClippedPatternCandidate(SequentialPattern prefix){
        
        List<ItemSeqIdentifier> itemSeqIdentifiers = prefix.getItemAppear();
        int patternLength = prefix.length();
        
        for (ItemSeqIdentifier identifier : itemSeqIdentifiers){
            boolean patternHasSplitInfo = false;
            boolean patternHasARP = false;
            List<pseudoSuperItem> curPattern = new ArrayList<>();
            int seqId = identifier.getSeqID();
            
            int superitemSetStartIdx = identifier.getSubSeqID() - patternLength + 1;  
            for (int i = 0; i < patternLength; i++){
                int superitemSetIdx = superitemSetStartIdx + i;
                int length = database.getSequenceByID(seqId).getItemsets().get(superitemSetIdx).size();
                for(int j = 0; j < length;j ++){   
                    SuperItem superitem = database.getSequenceByID(seqId).superItemAtPos(superitemSetIdx, j);
                    if (superitem.getSplitAlignPos() != -1){
                        patternHasSplitInfo = true;
                    }
                    if (superitem.getType().contains("ARP")){
                        patternHasARP = true;                        
                    }
                    pseudoSuperItem psItem = new pseudoSuperItem(seqId, superitemSetIdx, j);                        
                    psItem.setPsSuperitemLeftPos(database);
                    curPattern.add(psItem);
                }                
            }
            if (patternHasSplitInfo && !patternHasARP){
                pseudoSequentialPattern pattern = new pseudoSequentialPattern(curPattern, database, "clipped");
                System.out.println(pattern.toString(database));
                int chromIdx = pattern.ChromId;
                clippedPatternCandidates.get(chromIdx).add(pattern);
            }
            
        }
        
      
        
    }
    
    /**
     * patterns contain ARP patterns can be used to build connections.
     *
     * @param prefix     
     */
    private void saveARPPatternCandidate(SequentialPattern prefix){
        int numOfTypes = prefix.getNumOfTypes();
        int patternLength = prefix.length();
        boolean isARPPattern = prefix.isARPCandidatePattern();
                
        if (isARPPattern){           
            List<ItemSeqIdentifier> itemSeqIdentifiers = prefix.getItemAppear();
            for (ItemSeqIdentifier itemIdentity : itemSeqIdentifiers){
                List<pseudoSuperItem> curPattern = new ArrayList<>();
                int seqId = itemIdentity.getSeqID();
                int superitemSetStartIdx = itemIdentity.getSubSeqID() - patternLength + 1;                
                for (int i = 0; i < patternLength; i ++){
                    int superitemSetIdx = superitemSetStartIdx + i;
                    int length = database.getSequenceByID(seqId).getItemsets().get(superitemSetIdx).size();
                    for(int j = 0; j < length;j ++){                    
                        pseudoSuperItem psItem = new pseudoSuperItem(seqId, superitemSetIdx, j);                        
                        psItem.setPsSuperitemLeftPos(database);
                        curPattern.add(psItem);
                    }                
                }
                pseudoSequentialPattern pattern = new pseudoSequentialPattern(curPattern, database);

                if (patternLength >= 2 && numOfTypes == 1){
                    int supportArps = numOfSupportARPs(pattern);
                    if (supportArps > 1){
//                        int patternLeftMostPos = pattern.patternLeftMostPos;                
                        int chromIdx = pattern.ChromId;
                        while (ARPPatternCandidates.size() < chromIdx + 1){
                            ARPPatternCandidates.add(new ArrayList<>());
                        }
                        ARPPatternCandidates.get(chromIdx).add(pattern);
                        arpPatternCount ++;
                    }
                }else{
                    
                    int chromIdx = pattern.ChromId;
                    while (ARPPatternCandidates.size() < chromIdx + 1){
                        ARPPatternCandidates.add(new ArrayList<>());
                    }
                    ARPPatternCandidates.get(chromIdx).add(pattern);
                    arpPatternCount ++;
                }                                
                
            }                       
        }         
    }
    
    private void wgsPatternMerge(String linkerOut, String svRegionOut) throws IOException{
        System.out.println("\nStart pattern post-processing, total arp patterns: " + arpPatternCount);
        
        BufferedWriter regionWriter = new BufferedWriter(new FileWriter(svRegionOut));
        BufferedWriter linkWriter = new BufferedWriter(new FileWriter(linkerOut));
        
        int numChrs = ARPPatternCandidates.size();
        for (int i = 0; i < numChrs; i ++){
            List<pseudoSequentialPattern> arpPatterns = ARPPatternCandidates.get(i);
            Map<Integer, List<Integer>> indexMap = getPatternStartIndexMap(arpPatterns);
            List<pseudoSequentialPattern> mergedPatterns = oneChromMerge(i, arpPatterns, indexMap);
            oneChrARPPatternLinkageAnalysis(linkWriter, regionWriter, mergedPatterns);
        }
        regionWriter.close();
        linkWriter.close();
    }
    
    private Map<Integer, List<Integer>> getPatternStartIndexMap(List<pseudoSequentialPattern> arpPatterns){
        Map<Integer, List<Integer>> indexMap = new HashMap<>();
        int numOfPatterns = arpPatterns.size();
        for (int i = 0; i < numOfPatterns ; i++){
            int patternLeftMostPos = arpPatterns.get(i).patternLeftMostPos;
            List<Integer> indexList = indexMap.get(patternLeftMostPos);
            if (indexList == null) {
                indexList = new ArrayList<>();
                indexMap.put(patternLeftMostPos, indexList);
            }
            indexList.add(i);
        }
        return indexMap;
    }
    
    
    private List<pseudoSequentialPattern> oneChromMerge(int chrom, List<pseudoSequentialPattern> arpPatterns, Map<Integer, List<Integer>> patternStartIndexMap){
        int chr = chrom + 1;
        System.out.println("\nProcess Chr:"+chr+" ARP pattern before merge: " + ARPPatternCandidates.get(chrom).size());
        List<pseudoSequentialPattern> mergedARPPatternCandidates = new ArrayList<>();
//        List<Entry<Integer, List<Integer>>> patternIndexEntrys = new ArrayList<>(patternStartAndIndexMap.entrySet());
        List<Entry<Integer, List<Integer>>> patternIndexEntrys = new ArrayList<>(patternStartIndexMap.entrySet());
        Collections.sort(patternIndexEntrys, new Comparator<Entry<Integer, List<Integer>>>(){
            @Override
            public int compare(Entry<Integer, List<Integer>> o1, Entry<Integer, List<Integer>> o2){
                return o1.getKey().compareTo(o2.getKey());
            }
        
        });
        
        int entrysSize = patternIndexEntrys.size();
        Set<Integer> tracker = new HashSet<>();
        for (int i = 0; i < entrysSize - 1; i ++){
            int candidateSize = mergedARPPatternCandidates.size();
            Entry<Integer, List<Integer>> entry = patternIndexEntrys.get(i);
            int pos = entry.getKey();
            
            List<Integer> patternIndex = entry.getValue();

            Entry<Integer, List<Integer>> nextEntry = patternIndexEntrys.get(i + 1);                
            List<Integer> nextPatternIndex = nextEntry.getValue();
            int nextPos = nextEntry.getKey();
            
            pseudoSequentialPattern mergedPattern = mergePatternList(arpPatterns, patternIndex);
            pseudoSequentialPattern nextMergedPattern = mergePatternList(arpPatterns, nextPatternIndex);
    
//            System.out.println(entry.getKey() + ": " + patternIndex.toString() + "\t" + mergedPattern.toString(database));
//            System.out.println(nextEntry.getKey() + ": " + nextPatternIndex.toString() + "\t" + nextMergedPattern.toString(database));
            List<pseudoSuperItem> mergedSuperItems = mergedPattern.mergeTwoPattern(nextMergedPattern, database);
            if (!mergedSuperItems.isEmpty()){
                tracker.add(pos);
                tracker.add(nextPos);
                
                pseudoSequentialPattern newMergedPattern = new pseudoSequentialPattern(mergedSuperItems, database);
//                System.out.println("merged: " + newMergedPattern.toString(database));

                // the new pattern might be merged with the last pattern in the candidate list.
                if (!mergedARPPatternCandidates.isEmpty()){
                    
                    List<pseudoSuperItem> superitems = secondaryMerge(mergedARPPatternCandidates, newMergedPattern);
                    if (!superitems.isEmpty()){
                        pseudoSequentialPattern secondaryMergedPattern = new pseudoSequentialPattern(superitems, database);
//                        System.out.println("Added pattern: " + secondaryMergedPattern.toString(database));
                        mergedARPPatternCandidates.remove(candidateSize - 1);
                        mergedARPPatternCandidates.add(secondaryMergedPattern);
                    }else{
                        
//                        System.out.println("Added pattern: " + newMergedPattern.toString(database));
                        mergedARPPatternCandidates.add(newMergedPattern);
                    }
                }else{
//                    System.out.println("Added pattern: " + newMergedPattern.toString(database));

                    mergedARPPatternCandidates.add(newMergedPattern);
                }
                               
            }else{
                if (!tracker.contains(pos)){
                    tracker.add(pos);
                    if (! mergedARPPatternCandidates.isEmpty()){
                        List<pseudoSuperItem> superitems = secondaryMerge(mergedARPPatternCandidates, mergedPattern);
                        if (!superitems.isEmpty()){
                            pseudoSequentialPattern secondaryMergedPattern = new pseudoSequentialPattern(superitems, database);
                            mergedARPPatternCandidates.remove(candidateSize - 1);
                            mergedARPPatternCandidates.add(secondaryMergedPattern);
                        }else{
                            mergedARPPatternCandidates.add(mergedPattern);
                        }
                    }else{
                        mergedARPPatternCandidates.add(mergedPattern);
                    }                                        
                }                                                               
            }
                                
        }

        System.out.println("ARP pattern after merge: " + mergedARPPatternCandidates.size());
        return mergedARPPatternCandidates;
    }
    
    private List<pseudoSuperItem> secondaryMerge(List<pseudoSequentialPattern> ARPCandidates, pseudoSequentialPattern aPattern){
        int candidateSize = ARPCandidates.size();
        pseudoSequentialPattern lastSPInCandidates = ARPCandidates.get(candidateSize - 1);
        List<pseudoSuperItem> superitems = lastSPInCandidates.mergeTwoPattern(aPattern, database);
        return superitems;
    }
    
    private void oneChrARPPatternLinkageAnalysis(BufferedWriter linkWriter, BufferedWriter regionWriter, List<pseudoSequentialPattern> mergedPatterns) throws IOException{
        Linkage linkageAnalyzer = new Linkage();
        
        Collections.sort(mergedPatterns);
     
       
        int ARPPatternCount = mergedPatterns.size();
        
        Set<Integer> linkedPatternIdx = new HashSet<>();
        Set<Integer> selfUnLinkablePatternIdx = new HashSet<>();
        for (int i = 0; i < ARPPatternCount; i ++){
            StringBuilder sb = new StringBuilder();
            pseudoSequentialPattern ARPpattern = mergedPatterns.get(i);
//            if (ARPpattern.patternLeftMostPos == 80928556){
//                System.out.println(ARPpattern.toString(database));
//            }
            
            boolean isSelfLinked = ARPpattern.isSelfLinked(database);
            
            if (isSelfLinked){
                
                linkedPatternIdx.add(i);
                int[] bpPos = ARPpattern.selfLinkedPatternBP(database);                                                                
                // not a multi-breakpoint case
                if (bpPos[0] != 0){ 
                    int linkSupport = ARPpattern.getLinkSupport();
                    svOutInfo svInoInfo = new svOutInfo(bpPos[0], bpPos[1], ARPpattern.toString(database), -1, linkSupport, ARPpattern.getWeights(), ARPpattern.getPos());
                    sb.append(ARPpattern.ChromId);
                    sb.append("\t");
//                    sb.append(bpPos[0]);
//                    sb.append("\t");
//                    sb.append(bpPos[1]);    
                    sb.append((svInoInfo.toString()));
                    regionWriter.write(sb.toString());
                    regionWriter.newLine();
                    
                }
                sb = new StringBuilder();
                sb.append(ARPpattern.toString(database));
                linkWriter.write(sb.toString());
                linkWriter.newLine();

            }else{
                Map<Integer, Integer> indexMap = new HashMap<>();
                List<pseudoSequentialPattern> removedPatternCandidates = linkageAnalyzer.minusItemAndCopyWithIndexMap(mergedPatterns, i, indexMap);
                int[] vals = searchMatePattern(removedPatternCandidates, ARPpattern, linkageAnalyzer);
                int mateARPindex = vals[0];
                int linkSup = vals[1];
                
                if (mateARPindex != -1 && !linkedPatternIdx.contains(i)){
                    int orignialIndex = indexMap.get(mateARPindex);
                    if (selfUnLinkablePatternIdx.contains(orignialIndex)){
                        selfUnLinkablePatternIdx.remove(orignialIndex);
                    }
                    if (linkedPatternIdx.contains(orignialIndex)){
                        continue;
                    }
                    pseudoSequentialPattern matePattern = removedPatternCandidates.get(mateARPindex);
//                    System.out.println(matePattern.toString(database));
                    boolean selfLinked = matePattern.isSelfLinked(database);
                    if (selfLinked){
                        linkedPatternIdx.add(i);
//                        int[] bpPos = ARPpattern.selfLinkedPatternBreakpoint(database);
                        int[] bpPos = matePattern.selfLinkedPatternBP(database);
                        // not a multi-breakpoint case
                        if (bpPos[0] != 0){
                            int linkSupport = ARPpattern.getLinkSupport();
                            svOutInfo svInoInfo = new svOutInfo(bpPos[0], bpPos[1], ARPpattern.toString(database), -1, linkSupport, ARPpattern.getWeights(), ARPpattern.getPos());
                            sb.append(ARPpattern.ChromId);
//                            sb.append("\t");
//                            sb.append(bpPos[0]);
//                            sb.append("\t");
//                            sb.append(bpPos[1]);  
                            sb.append(svInoInfo.toString());
                            regionWriter.write(sb.toString());
                            regionWriter.newLine();                            
                        }
                        sb = new StringBuilder();
                        sb.append(ARPpattern.toString(database));
                        linkWriter.write(sb.toString());
                        linkWriter.newLine();
                    }else{                        
//                        if(selfUnLinkablePatternIdx.contains(orignialIndex)){
//                            continue;
//                        }
                        linkedPatternIdx.add(i);                        
                        linkedPatternIdx.add(orignialIndex);
//                        linkedPatternIdx.add(mateARPindex);

                        boolean mateIsAfter = matePattern.patternLeftMostPos > ARPpattern.patternLeftMostPos ? true : false;

                        if (mateIsAfter){
                            int[] coords = ARPpattern.estimateBreakpointPos(matePattern, database);
                            
                            if (coords[0] != -1){
                                String patternStr = ARPpattern.toString(database) + "<>" + matePattern.toString(database);
                                List<Integer> weights = ARPpattern.getWeights();
                                List<Integer> weightsTwo = matePattern.getWeights();
                                weights.addAll(weightsTwo);
                                List<Integer> pos = ARPpattern.getPos();
                                pos.addAll(matePattern.getPos());
                                
                                svOutInfo svInoInfo = new svOutInfo(coords[0], coords[1], patternStr, 1, linkSup, weights, pos);
                                sb.append(ARPpattern.ChromId);
                                sb.append("\t");
//                                sb.append(coords[0]);
//                                sb.append("\t");
//                                sb.append(coords[1]);
                                sb.append(svInoInfo.toString());
                                regionWriter.write(sb.toString());
                                regionWriter.newLine();
//                                
                            }
                            sb = new StringBuilder();
                            sb.append(ARPpattern.toString(database));
                            sb.append(" -- ");
                            sb.append(matePattern.toString(database));
                            linkWriter.write(sb.toString());
                            linkWriter.newLine(); 
                        }else{
                            int[] coords = matePattern.estimateBreakpointPos(ARPpattern, database);
                            if (coords[0] != -1){
                                String patternStr = matePattern.toString(database) + "<>" + ARPpattern.toString(database);
                                List<Integer> weights = matePattern.getWeights();
                                List<Integer> weightsTwo = ARPpattern.getWeights();
                                weights.addAll(weightsTwo);
                                
                                List<Integer> pos = matePattern.getPos();
                                pos.addAll(ARPpattern.getPos());
                                
                                svOutInfo svInoInfo = new svOutInfo(coords[0], coords[1], patternStr, 1, linkSup, weights, pos);
                                sb.append(ARPpattern.ChromId);
//                                sb.append("\t");
//                                sb.append(coords[0]);
//                                sb.append("\t");
//                                sb.append(coords[1]); 
                                sb.append(svInoInfo.toString());
                                regionWriter.write(sb.toString());
                                regionWriter.newLine();
                                
                            }
                            sb = new StringBuilder();
                            sb.append(matePattern.toString(database));
                            sb.append(" -- ");
                            sb.append(ARPpattern.toString(database));
                            linkWriter.write(sb.toString());
                            linkWriter.newLine(); 

                        }
                    }                          
                }else{      
                    if (!linkedPatternIdx.contains(i)){
                        selfUnLinkablePatternIdx.add(i);
                    }
                }       
            }
           
        }
        for (Integer idx : selfUnLinkablePatternIdx){       
            StringBuilder sb = new StringBuilder();
            pseudoSequentialPattern ARPpattern = mergedPatterns.get(idx);

            int[] coords = ARPpattern.noneLinkablePatternBreakpoint(database);
            if (coords[0] != 0 && coords[1] != 0){
                int linkSupport = ARPpattern.getLinkSupport();
                svOutInfo svInoInfo = new svOutInfo(coords[0], coords[1], ARPpattern.toString(database), 0, linkSupport, ARPpattern.getWeights(), ARPpattern.getPos());
                sb.append(ARPpattern.ChromId);
                sb.append("\t");
//                sb.append(coords[0]);
//                sb.append("\t");
//                sb.append(coords[1]);
                sb.append(svInoInfo.toString());
                regionWriter.write(sb.toString());
                regionWriter.newLine();
            }
            sb = new StringBuilder();
            sb.append(ARPpattern.toString(database));
            linkWriter.write(sb.toString());
            linkWriter.newLine();
        }
        
    }
   /**
     * For arp pattern of length larger than 2 with same items, check connections
     * @param pattern
     * @return 
     */
    public int numOfSupportARPs(pseudoSequentialPattern pattern){
        List<SuperItem> superitemList = pattern.getSuperItemsOfPattern(database);
        int length = superitemList.size();
        int maximuSup = 0;
        for (int i = 0; i < length; i++){
            int supportedARPs = 0;
            for (int j = 0; j < length; j ++){
                if (i != j){
                    SuperItem superitemOne = pattern.getSuperItemFromOriginal(database, i);
                    SuperItem superitemTwo = pattern.getSuperItemFromOriginal(database, j);

                    List<byte[]> qnameOneByteList = superitemOne.getQNames();
                    List<byte[]> qnameTwoByteList = superitemTwo.getQNames();
                    Set<String> uniqueQName = new HashSet<>();

                    
                    for (byte[] array : qnameOneByteList){
                        String qname = new String(array);              
                        uniqueQName.add(qname);
                    }       
                    for (byte[] array : qnameTwoByteList){
                        String qname = new String(array);
                        if (uniqueQName.contains(qname)){
                            supportedARPs += 1;                
                        }
                        uniqueQName.add(qname);
                    }
                }
                
            }
            if (supportedARPs > maximuSup){
                maximuSup = supportedARPs;
            }
            
        }
        
        return maximuSup;
    }
  
      
    private int[] searchMatePattern(List<pseudoSequentialPattern> sortedARPPatterns, pseudoSequentialPattern targetPattern, Linkage linker){        
        List<QueryInterval> targetPatternMateInterval = targetPattern.superitemMateInterval;
        int length = sortedARPPatterns.size();
        int startIdx = 0;
        int endIdx = length - 1;
        int mateIndex = -1;
        int noQueryInterval = targetPatternMateInterval.size();
        
        int linkSup = -1;
        int targetPatternMatchSuperItemIdx = -1;
        int matchedMatePatternSuperItemIdx = -1;
        for (int i = 0; i < noQueryInterval ; i++){
            QueryInterval interval = targetPatternMateInterval.get(i);
            while (startIdx <= endIdx){
                int midIdx = startIdx + (endIdx - startIdx) / 2;

                pseudoSequentialPattern pattern = sortedARPPatterns.get(midIdx);
                List<QueryInterval> sortedIntervals = pattern.superitemInterval;

                int overlapAt = hasOverlap(interval, sortedIntervals);
                if (overlapAt != -1){
                    mateIndex = midIdx;                    
                    targetPatternMatchSuperItemIdx = i;
                    matchedMatePatternSuperItemIdx = overlapAt;
                    break;
                }
                if (isAfterInterval(interval, sortedIntervals)){
                    startIdx = midIdx + 1;                    
                }
                if (isAheadInterval(interval, sortedIntervals)){
                    endIdx = midIdx - 1;
                }
                else if(overlapAt == -1 && !isAfterInterval(interval, sortedIntervals) && !isAheadInterval(interval, sortedIntervals)){
                    startIdx = midIdx + 1;
                }
            }
            if (mateIndex != -1){
                break;
            }
            // reset start and end for next interval match
            startIdx = 0;
            endIdx = length - 1;
            
        }
        if (mateIndex != -1){
            SuperItem superitemOne = targetPattern.getSuperItemOfPatternAtPos(database, targetPatternMatchSuperItemIdx);
            pseudoSequentialPattern matchedSequentialPattern = sortedARPPatterns.get(mateIndex);
            SuperItem superitemTwo = matchedSequentialPattern.getSuperItemOfPatternAtPos(database, matchedMatePatternSuperItemIdx);
            boolean isEnoughARPs = linker.supportARPs(superitemOne, superitemTwo);
            linkSup = linker.getSupLink();
            if (!isEnoughARPs){
                mateIndex = -1;
            }
        }
        int[] returnVals = new int[2];
        returnVals[0] = mateIndex;
        returnVals[1] = linkSup;
        return returnVals;
    }
    
    
    private pseudoSequentialPattern mergePatternList(List<pseudoSequentialPattern> arpPatterns, List<Integer> patternIndex){
        List<pseudoSequentialPattern> patterns = new ArrayList<>();
        for (Integer idx : patternIndex){
            patterns.add(arpPatterns.get(idx));
        }
        int maxLength = 0;
        int maxLengthPatternIndex = 0;
        int patternsSize = patterns.size();
        for (int i = 0; i < patternsSize; i ++){                
            if (patterns.get(i).patternLength > maxLength){
                maxLength = patterns.get(i).patternLength;
                maxLengthPatternIndex = i;
            }
        }
        return patterns.get(maxLengthPatternIndex);
    }
    private int hasOverlap(QueryInterval targetInterval, List<QueryInterval> intervals){
        int overlapAtIdx = -1;
        for (int i = 0; i < intervals.size(); i++){
            QueryInterval interval = intervals.get(i);
            if (interval != null){
//                if(targetInterval.overlaps(interval)){
//                    overlapAtIdx = i;
//                    break;
//                }
                if (reciprocalOverlap(targetInterval, interval)){
                    overlapAtIdx = i;
                    break;
                }
            }
        }
        return overlapAtIdx;
    }
    private boolean reciprocalOverlap(QueryInterval a, QueryInterval b){
        int aSize = a.end - a.start;
        int bSize = b.end - b.start;
        boolean isOverlapped = false;
        
        if (b.start < a.end && b.start >= a.start){
            int overlapSize = a.end - b.start;
            double aOverlapRatio = (double) overlapSize / aSize;
            double bOverlapRatio = (double) overlapSize / bSize;
            if (aOverlapRatio >= 0.1 && bOverlapRatio >= 0.1){
                isOverlapped = true;
            }
        }else if (a.start < b.end && a.start >= b.start){
            int overlapSize = a.end - b.start;
            double aOverlapRatio = (double) overlapSize / aSize;
            double bOverlapRatio = (double) overlapSize / bSize;
            if (aOverlapRatio >= 0.1 && bOverlapRatio >= 0.1){
                isOverlapped = true;
            }
        }
        
        return isOverlapped;
    }
    private boolean isAheadInterval(QueryInterval targetInterval, List<QueryInterval> intervals){
        boolean isAhead = false;
        QueryInterval leftMostInterval = intervals.get(0);
        if (targetInterval.end < leftMostInterval.start){
            isAhead = true;
        }
        return isAhead;        
    }
    
    private boolean isAfterInterval(QueryInterval targetInterval, List<QueryInterval> intervals){
        boolean isAfter = false;
        QueryInterval lastInterval = intervals.get(intervals.size() - 1);
        if (targetInterval.start > lastInterval.end){
            isAfter = true;
        }
        return isAfter;
    }
    
    public void printAlgoStatistics(){
        StringBuilder r = new StringBuilder(200);
        r.append("\n=============  Contiguous-FSPM =============\n Total time ~ ");
        r.append(endTime - startTime);
        r.append(" ms\n");     
        r.append(" Max memory (mb) : " );
        r.append(MemoryLogger.getInstance().getMaxMemory());        
        r.append('\n');
        r.append("===================================================\n");
        System.out.println(r.toString());
    }
    
        
}
