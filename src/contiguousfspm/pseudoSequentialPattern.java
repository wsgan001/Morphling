/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package contiguousfspm;

import dataStructures.SequenceDatabase;
import dataStructures.SuperItem;
import utils.Linkage;
import htsjdk.samtools.QueryInterval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;

/**
 *
 * @author jiadonglin
 */
public class pseudoSequentialPattern implements Comparable<pseudoSequentialPattern>{
    
        List<pseudoSuperItem> superitems = new ArrayList<>();
        List<pseudoSuperItem> ARPSuperItems = new ArrayList<>();
        List<QueryInterval> superitemInterval;
        List<QueryInterval> superitemMateInterval;
//        int leftMostNotNull = -1;
        int ChromId = -1;
        String chromName;
        int patternLength;
        int numOfBPs;
        int patternLeftMostPos; // position of the first superitem in the pattern
        int patternLeftMostIntervalStart; // read position of the left most interval of the pattern
        int[] linkedPatterns = new int[2]; // for a self linked pattern, save its linked coord for boundary estimation
        int numOfLinkedEvidence;
        List<Integer> weights = new ArrayList<>();
        List<Integer> postions = new ArrayList<>();
        
        Map<QueryInterval, List<Integer>> indexMap = new HashMap<>();
        
        @Override
        public int compareTo(pseudoSequentialPattern other){                                          
            return patternLeftMostIntervalStart - other.patternLeftMostIntervalStart;
        }
        
        public List<QueryInterval> getSortedIntervals(){
            List<QueryInterval> sortedInterval = new ArrayList<>();
            for (QueryInterval val : superitemInterval){
                sortedInterval.add(val);
            }
            return sortedInterval;
        }
        
        public pseudoSequentialPattern(List<pseudoSuperItem> itemset, SequenceDatabase database, String clipped){
            superitems = itemset;
            patternLength = itemset.size();
            for (int i = 0; i < patternLength ; i ++){
                SuperItem superitem = superitems.get(i).getSuperItem(database);    
                if (ChromId == -1){
                    ChromId = superitem.getChromIdx();
                    break;
                }
            }  
        }
        
        public pseudoSequentialPattern (List<pseudoSuperItem> itemset, SequenceDatabase database){

            superitems = itemset;
            patternLength = itemset.size();
            superitemInterval = new ArrayList<>();
            superitemMateInterval = new ArrayList<>();
            
            patternLeftMostPos = itemset.get(0).getSuperItem(database).getPos();
            
            for (int i = 0; i < patternLength ; i ++){
                SuperItem superitem = superitems.get(i).getSuperItem(database);    
                weights.add(superitem.getWeight());
                postions.add(superitem.getPos());
                if (ChromId == -1){
                    ChromId = superitem.getChromIdx();
                }
                if (superitem.isARPsuperitem() && ! superitem.getType().contains("OEM")){                      
                    ARPSuperItems.add(superitems.get(i));                                         
                    superitemInterval.add(superitem.getSuperitemRegion());
                    superitemMateInterval.add(superitem.getSuperitemMateRegion());
                }else{
                    numOfBPs += 1;
                }
            }   
            Collections.sort(ARPSuperItems);
            Collections.sort(superitemInterval);
            patternLeftMostIntervalStart = superitemInterval.get(0).start;
        } 
        
        public List<Integer> getIndex(QueryInterval aInterval){
            return indexMap.get(aInterval);
        }
        public List<Integer> getWeights(){
            return weights;
        }
        public List<Integer> getPos(){
            return postions;
        }
        public int getLinkSupport(){
            return numOfLinkedEvidence;
        }
        
        public int getNumOfBPs(){
            return numOfBPs;
        }
        
        public String toString(SequenceDatabase database){
            StringBuilder sb = new StringBuilder();
//            sb.append("IS: " + superitemInterval.get(0).start);
//            sb.append(" - ");
            for (pseudoSuperItem item : superitems){
                SuperItem superitem = item.getSuperItem(database);
                sb.append('(');
                sb.append(superitem.getType());
                sb.append(')');
//                sb.append(';');
            }
//            String str = sb.substring(0, sb.length() - 1);
            return sb.toString();
        }
        public List<SuperItem> getSuperItemsOfPattern(SequenceDatabase database){
            List<SuperItem> superItemsList = new ArrayList<>();
            for (pseudoSuperItem item : superitems){
                SuperItem superitem = item.getSuperItem(database);
                superItemsList.add(superitem);
            }
            return superItemsList;
        }
        
        public boolean isSelfLinked(SequenceDatabase database){
            
            boolean linked = false;
            boolean hasEnoughARPs = false;
            Linkage linker = new Linkage();
            
            int Arps = ARPSuperItems.size();
            int curSuperItemIdx = -1;
            List<pseudoSuperItem> searchSpace;
            int machtedSuperItemIdx = -1;
            
            Set<Integer> matchedItem = new HashSet<>();
            
            for (int i = 0; i < Arps; i++){
                if (matchedItem.contains(i)){
                    continue;
                }
                pseudoSuperItem target = ARPSuperItems.get(i);
                Map<Integer, Integer> idxMap = new HashMap<>();
                searchSpace = minusSelf(ARPSuperItems, i, idxMap);
                
                int mateIndex = linker.mateSuperItemSearch(searchSpace, target, database);
                if (mateIndex != -1){    
                    
                    curSuperItemIdx = i;
                    machtedSuperItemIdx = mateIndex;
                    int originalIdx = idxMap.get(machtedSuperItemIdx);
                    matchedItem.add(originalIdx);
                    SuperItem superItemOne = ARPSuperItems.get(curSuperItemIdx).getSuperItem(database);
//                    SuperItem superItemTwo = searchSpace.get(machtedSuperItemIdx).getSuperItem(database);
                    SuperItem superItemTwo = ARPSuperItems.get(originalIdx).getSuperItem(database);
                    
                    boolean isEnoughARPs = linker.supportARPs(superItemOne, superItemTwo);
                    numOfLinkedEvidence = linker.getSupLink();
                    if (isEnoughARPs){
                        hasEnoughARPs = true;
                        int superItemOnePos = superItemOne.getPos();
                        int superItemTwoPos = superItemTwo.getPos();
                        if (superItemTwoPos > superItemOnePos){
                            linkedPatterns[0] = curSuperItemIdx;
                            linkedPatterns[1] = idxMap.get(machtedSuperItemIdx);
                        }else{
                            linkedPatterns[0] = idxMap.get(machtedSuperItemIdx);
                            linkedPatterns[1] = curSuperItemIdx;
                        }
//                        break;
                    }
                    
                }
            }
            if (hasEnoughARPs){ 
                linked = true;              
            }
            return linked;
        }
//        private int depthSearch(List<pseudoSuperItem> psSuperItems, int start){
//            
//        }
        private List<pseudoSuperItem> minusSelf(List<pseudoSuperItem> psSuperItems, int idx, Map<Integer, Integer> idxMap){
            List<pseudoSuperItem> newSuperItems= new ArrayList<>();
            int length = psSuperItems.size();
            for (int i = 0; i < length; i ++){
                if (i != idx){
                    newSuperItems.add(psSuperItems.get(i));
                    idxMap.put(newSuperItems.size() - 1, i);
                }
                
            }
            return newSuperItems;
        }
        
        public SuperItem getSuperItemFromOriginal(SequenceDatabase database, int idx){
            pseudoSuperItem item = superitems.get(idx);
            return item.getSuperItem(database);
        }
        public SuperItem getSuperItemOfPatternAtPos(SequenceDatabase database, int idx){
//            Collections.sort(ARPSuperItems);
            pseudoSuperItem item = ARPSuperItems.get(idx);
            return item.getSuperItem(database);
        }
        public List<pseudoSuperItem> mergeTwoPattern(pseudoSequentialPattern aPattern, SequenceDatabase database){          
            
            List<pseudoSuperItem> mergedSuperitems = new ArrayList<>();
            List<SuperItem> patternOneSuperItems = getSuperItemsOfPattern(database);
            List<SuperItem> patternTwoSuperItems = aPattern.getSuperItemsOfPattern(database);

            int lengthOne = patternOneSuperItems.size();
            int lengthTwo = patternTwoSuperItems.size();

            int matchedIndexAtPatternOne = -1;           
            int lastMatchedIndexAtPatternOne = lengthOne; 
            
            SuperItem patternTwoStartSuperItem = patternTwoSuperItems.get(0);
            SuperItem patternTwoLastSuperItem = patternTwoSuperItems.get(lengthTwo - 1);
            
            for (int i = 0; i < lengthOne ;i++){
                SuperItem patternOneSuperItem = patternOneSuperItems.get(i);
                if (patternTwoStartSuperItem.isEqual(patternOneSuperItem)){
                    matchedIndexAtPatternOne = i;                    
                }
                if (patternTwoLastSuperItem.isEqual(patternOneSuperItem)){
                    lastMatchedIndexAtPatternOne = i;
                }
            }
                                                      
            if (matchedIndexAtPatternOne != - 1){
                if (lastMatchedIndexAtPatternOne == lengthOne){
                    List<pseudoSuperItem> subListOfPatternOne = superitems.subList(0, matchedIndexAtPatternOne);
                    mergedSuperitems.addAll(subListOfPatternOne);
                    mergedSuperitems.addAll(aPattern.superitems);
                }else{
                    mergedSuperitems = superitems;
                }                
            }
            
            return mergedSuperitems;            
        }

    public int[] estimateBreakpointPos(pseudoSequentialPattern matePattern, SequenceDatabase database){
        int[] coords = new int[2];
        
        // Process the patten ahead.
        int rightMostARPsuperitemIdx = -1;
        String ARPsuperItemOri = "";
        int ARPsuperItemPos = 0;
        for (int k = patternLength - 1; k >=0 ;k--){
            SuperItem superitem = getSuperItemFromOriginal(database, k);
            if (superitem.isARPsuperitem()){
                coords[0] = superitem.getPos();
                rightMostARPsuperitemIdx = k;
                ARPsuperItemOri = superitem.getOri();
                ARPsuperItemPos = superitem.getPos();
                break;
            }
        }        
        if (ARPsuperItemOri.equals("+")){
            int closestBP = 500;
            if (rightMostARPsuperitemIdx < patternLength){
                for (int i = rightMostARPsuperitemIdx; i < patternLength; i ++){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("MS")){
                        coords[0] = superitem.getPos();
                        int dist = Math.abs(superitem.getPos() - ARPsuperItemPos);
                        if (dist < closestBP){
                            closestBP = dist;
                        }
                        break;
                    }
                }
            }if(rightMostARPsuperitemIdx > 0){
                for (int i = rightMostARPsuperitemIdx; i>=0 ; i--){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("MS")){
                        int dist = Math.abs(superitem.getPos() - ARPsuperItemPos);
                        if (dist < closestBP){                            
                            coords[0] = superitem.getPos();
                            break;
                        }
                    }
                }
            }
            
        }

        

        // Process the mate pattern
        int leftMostARPsuperitemIdx = -1;
        int matePatternLength = matePattern.patternLength;
        String mateARPsuperItemOri = "";
        int mateARPsuperItemPos = 0;
        for (int k = 0; k < matePatternLength; k ++){
            SuperItem superitem = matePattern.getSuperItemFromOriginal(database, k);
            if (superitem.isARPsuperitem()){
                leftMostARPsuperitemIdx = k;
                mateARPsuperItemOri = superitem.getOri();
                coords[1] = superitem.getPos();
                mateARPsuperItemPos = superitem.getPos();
                break;
            }
        }
        if (mateARPsuperItemOri.equals("-")){
            int closestBP = 500;
            if (leftMostARPsuperitemIdx > 0){
                for (int i = leftMostARPsuperitemIdx; i >= 0; i--){
                    SuperItem superitem = matePattern.getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("SM")){
                        coords[1] = superitem.getPos();
                        int dist = Math.abs(mateARPsuperItemPos - superitem.getPos());
                        if (dist < closestBP){
                            closestBP = dist;
                        }
                        break;
                    }
                }
            }
            if (leftMostARPsuperitemIdx < matePatternLength){
                for (int i = leftMostARPsuperitemIdx; i < matePatternLength; i++){
                    SuperItem superitem = matePattern.getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("SM")){
                        int dist = Math.abs(mateARPsuperItemPos - superitem.getPos());
                        if (dist < closestBP){
                            coords[1] = superitem.getPos();
                            break;
                        }
                        
                    }
                }
            }
        }
        return coords;
    }
    private int getARPpatternIdx(SuperItem superitem, SequenceDatabase database){
        int idx = -1;
        for (int i = 0; i < patternLength; i++){
            SuperItem curSuperItem = getSuperItemFromOriginal(database, i);
            if (curSuperItem.isEqual(superitem)){
                idx = i;
                break;
            }
        }
        return idx;
    }
    
    public int[] noneLinkablePatternBreakpoint(SequenceDatabase database){
        // use split-alignment to get the BP pos. 
        // If there is no split-alignment, use traditional way to estimate, this is not accurate. 
        int[] pos = splitAlignForBP(database);                        
        if (pos[0] == 0 || pos[1] == 0){
            int leftIdx = -1; // left breakpoint SuperItem idx
            int leftBp = -1; // estimated left breakpoint position
            int leftARPpos = -1; // left ARP SuperItem position
            int leftARPidx = -1; // left ARP SuperItem idx
            for (int i = 0; i < patternLength; i++){
                SuperItem superitem = getSuperItemFromOriginal(database, i);
                if (superitem.isARPsuperitem() && superitem.getOri().equals("+") && leftARPpos == -1){
                    leftARPpos = superitem.getPos();
                    leftARPidx = i;
                    break;
                }
            }

            int distToLeftARPpos = 500;
            int searchLeftBp = 0;
            if (leftARPidx != -1){
                searchLeftBp = leftARPidx;
            }
            for (int i = searchLeftBp; i < patternLength ;i ++){
                SuperItem superitem = getSuperItemFromOriginal(database, i);
                if (!superitem.isARPsuperitem() && superitem.getType().equals("MS")){

                    if (leftARPpos != -1){
                        int dist = Math.abs(superitem.getPos() - leftARPpos);
                        if (dist < distToLeftARPpos){
                            leftBp = superitem.getPos();
                            distToLeftARPpos = dist;
                            leftIdx = i;
                            break;
                        }
                    }else{
                        leftBp = superitem.getPos();
                        leftIdx = i;
                        break;
                    }               
                }
            }
            if(leftIdx == -1 && leftARPidx != -1){
                for(int i = 0;i < leftARPidx;i++){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                        if (!superitem.isARPsuperitem() && superitem.getType().equals("MS")){

                        if (leftARPpos != -1){
                            int dist = Math.abs(superitem.getPos() - leftARPpos);
                            if (dist < distToLeftARPpos){
                                leftBp = superitem.getPos();
                                distToLeftARPpos = dist;
                                leftIdx = i;
                                break;
                            }
                        }else{
                            leftBp = superitem.getPos();
                            leftIdx = i;
                            break;
                        }
                    }
                }
            }


            int rightIdx = -1;
            int rightBp = -1;
            int rightARPpos = -1;
            int rightARPidx = -1;

            for (int i = patternLength - 1;i>=0;i--){
                SuperItem superitem = getSuperItemFromOriginal(database, i);
                if (superitem.isARPsuperitem() && superitem.getOri().equals("-") && rightARPpos == -1){
                    rightARPpos = superitem.getPos();
                    rightARPidx = i;
                    break;
                }
            }

            int distToRightARPpos = 500;
            int searchRightBp = patternLength - 1;
            if (rightARPidx != -1){
                searchRightBp = rightARPidx;
            }

            for (int i = searchRightBp; i >=0; i--){
                SuperItem superitem = getSuperItemFromOriginal(database, i);

                if (!superitem.isARPsuperitem() && superitem.getType().equals("SM")){

                    if (rightARPpos != -1){
                        int dist = Math.abs(superitem.getPos() - rightARPpos);
                        if (dist < distToRightARPpos){
                            distToRightARPpos = dist;
                            rightBp = superitem.getPos();
                            rightIdx = i;
                            break;
                        }
                    }else{
                        rightBp = superitem.getPos();
                        rightIdx = i;
                        break;
                    }                

                }
            }
            if (rightIdx == -1 && rightARPidx != -1 && rightARPidx < patternLength){
                for (int i=rightARPidx;i<patternLength;i++){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);

                    if (!superitem.isARPsuperitem() && superitem.getType().equals("SM")){

                        if (rightARPpos != -1){
                            int dist = Math.abs(superitem.getPos() - rightARPpos);
                            if (dist < distToRightARPpos){
                                distToRightARPpos = dist;
                                rightBp = superitem.getPos();
                                rightIdx = i;
                                break;
                            }
                        }else{
                            rightBp = superitem.getPos();
                            rightIdx = i;
                            break;
                        }                

                    }
                }

            }

            if (leftIdx != -1 && rightIdx != -1 && (rightIdx > leftIdx)){
                pos[0] = leftBp;
                pos[1] = rightBp;
            }
        }
        
        
        return pos;
    }   
    private int[] splitAlignForBP(SequenceDatabase database){
        int[] pos = new int[2];
        
        for (int i = 0;i < patternLength;i++){
            SuperItem superitem = getSuperItemFromOriginal(database, i);
            if (superitem.isARPsuperitem() || superitem.getWeight() < 5){
                continue;
            }
            int splitAlignPos = superitem.getSplitAlignPos();
            if (splitAlignPos != -1){
                if (superitem.getPos() < splitAlignPos){
                    pos[0] = superitem.getPos();
                    pos[1] = splitAlignPos;
                }else{
                    pos[0] = splitAlignPos;
                    pos[1] = superitem.getPos();
                }
//                pos[0] = superitem.getPos();
//                pos[1] = splitAlignPos;
                break;
            }
            
        }
        return pos;
    }
    
    
    public int[] selfLinkedPatternBP(SequenceDatabase database){
        int[] pos = new int[2];
        
        SuperItem leftARPpattern = ARPSuperItems.get(linkedPatterns[0]).getSuperItem(database);
        int leftARPpatternSuperItemPos = leftARPpattern.getPos();
        SuperItem rigthARPpattern = ARPSuperItems.get(linkedPatterns[1]).getSuperItem(database);
        int rightARPpatternSuperItemPos = rigthARPpattern.getPos();

        
        int leftClosestToARP = -1;
        int rightClosestToARP = -1;
        int maxDist = 500;
        for (int i = 0; i < patternLength; i++){
            SuperItem superitem = getSuperItemFromOriginal(database, i);
            if (!superitem.isARPsuperitem() && superitem.getType().equals("MS")){
                int dist = Math.abs(superitem.getPos() - leftARPpatternSuperItemPos);
                if (dist < maxDist){
                    maxDist = dist;
                    leftClosestToARP = i;
                }
            }            
        }
        if (leftClosestToARP != -1){
            SuperItem superitem = getSuperItemFromOriginal(database, leftClosestToARP);
            pos[0] = superitem.getPos();
            pos[1] = superitem.getSplitAlignPos();
            if (pos[1] == -1){
                pos[1] = rightARPpatternSuperItemPos;
            }
        }else{
            for(int i = patternLength - 1; i>=0; i--){
                SuperItem superitem = getSuperItemFromOriginal(database, i);
                if (!superitem.isARPsuperitem() && superitem.getType().equals("SM")){
                    int dist = Math.abs(superitem.getPos() - rightARPpatternSuperItemPos);
                    if (dist < maxDist){
                        maxDist = dist;
                        rightClosestToARP = i;
                    }
                } 
            }
            if(rightClosestToARP != -1){
                SuperItem superitem = getSuperItemFromOriginal(database, rightClosestToARP);
                pos[0] = superitem.getSplitAlignPos();
                if (pos[0] == -1){
                    pos[0] = leftARPpatternSuperItemPos;
                }
                pos[1] = superitem.getPos();
            }else{                
                pos[0] = leftARPpatternSuperItemPos;
                pos[1] = rightARPpatternSuperItemPos;
            }
        }
        
        
        return pos;
        
    }
    
    
    public int[] selfLinkedPatternBreakpoint(SequenceDatabase database){
        int[] pos = splitAlignForBP(database);

        SuperItem leftARPpattern = ARPSuperItems.get(linkedPatterns[0]).getSuperItem(database);
        int leftARPpatternSuperItemPos = leftARPpattern.getPos();
        SuperItem rigthARPpattern = ARPSuperItems.get(linkedPatterns[1]).getSuperItem(database);
        int rightARPpatternSuperItemPos = rigthARPpattern.getPos();
        int leftArpPatternIdx = getARPpatternIdx(leftARPpattern, database);
        int rightArpPatternIdx = getARPpatternIdx(rigthARPpattern, database);
        
        if (pos[0] == 0 && pos[1] == 0){
            pos[0] = leftARPpatternSuperItemPos;
            pos[1] = rightARPpatternSuperItemPos;
        }
        
        boolean rightBPnotRefined = true;

        if (rightArpPatternIdx - leftArpPatternIdx > 1){
            int leftBPidx = -1;

            if (leftARPpattern.getOri().equals("+")){
                for (int i = leftArpPatternIdx; i < rightArpPatternIdx; i++){
                    SuperItem superItem = getSuperItemFromOriginal(database, i);
                    int dist = Math.abs(superItem.getPos() - leftARPpatternSuperItemPos);
                    if (superItem.getType().equals("MS") && dist < 500){
                        pos[0] = superItem.getPos();
                        if (superItem.getSplitAlignPos() != -1){
                            pos[1] = superItem.getSplitAlignPos();
                        }
                        leftBPidx = i;
                        break;
                    }
                }
                if (leftBPidx == -1){
                    for (int i=leftArpPatternIdx;i>=0;i--){
                        SuperItem superItem = getSuperItemFromOriginal(database, i);
                        int dist = Math.abs(leftARPpatternSuperItemPos - superItem.getPos());
                        if (superItem.getType().equals("MS") && dist < 500){
                            pos[0] = superItem.getPos();
                            if (superItem.getSplitAlignPos() != -1){
                                pos[1] = superItem.getSplitAlignPos();
                                rightBPnotRefined = false;
                            }
                            leftBPidx = i;
                            break;
                        }
                    }
                }
            }
            
            if (rightBPnotRefined){
                int rightBPidx = -1;
                if (rigthARPpattern.getOri().equals("-")){
                    for (int k = rightArpPatternIdx; k > leftArpPatternIdx ;k--){
                        SuperItem superItem = getSuperItemFromOriginal(database, k);
                        int dist = Math.abs(rightARPpatternSuperItemPos - superItem.getPos());
                        if (superItem.getType().equals("SM") && dist < 500){   
                            rightBPidx = k;                        
                            pos[1] = superItem.getPos();
                            break;

                        }
                    }
                    if (rightBPidx == -1){
                        for (int k = rightArpPatternIdx; k < patternLength ;k++){
                            SuperItem superItem = getSuperItemFromOriginal(database, k);
                            int dist = superItem.getPos() - rightARPpatternSuperItemPos;
                            if (superItem.getType().equals("SM") && dist < 500){   
                                rightBPidx = k;
                                pos[1] = superItem.getPos();
                                break;
                            }
                        }
                    }
                }
            }
            
//            // BP superitem cross-over
//            if (rightBPidx <= leftBPidx) {
//                pos[0] = leftARPpattern.getPos();
//                pos[1] = rigthARPpattern.getPos();
//            }

        }else{
            if (leftARPpattern.getOri().equals("+")){
                for (int i=leftArpPatternIdx;i>=0;i--){
                    SuperItem superItem = getSuperItemFromOriginal(database, i);
                    int dist = Math.abs(superItem.getPos() - leftARPpatternSuperItemPos);

                    if (superItem.getType().equals("MS") && dist < 500){
                        pos[0] = superItem.getPos();
                        if (superItem.getSplitAlignPos() != -1){
                            pos[1] = superItem.getSplitAlignPos();
                        }
                        break;
                    }
                }
            }
//            if (rigthARPpattern.getOri().equals("-")){
//                for (int k = leftArpPatternIdx; k < patternLength ;k++){
//                    SuperItem superItem = getSuperItemFromOriginal(database, k);
//                    int dist = Math.abs(superItem.getPos() - rightARPpatternSuperItemPos);
//                    if (superItem.getType().equals("SM") && dist < 500){   
//                        pos[1] = superItem.getPos();
//                        break;
//                    }
//                }
//            }
        }
//        if (ARPbasedBP){
//            pos = splitAlignForBP(database);
//        }
                      
        return pos;
    }
}
