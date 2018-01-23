/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package superitemGenerator;
import htsjdk.samtools.*;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;

import java.io.BufferedWriter;
import java.io.FileWriter;

import Channels.*;
import dataStructures.SuperItem;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
/**
 *
 * @author jiadonglin
 */
public class SignalReader {
    
    private long startTime;
    private long endTime;
    
    private BufferedWriter writer = null;
    
    private String CHROM;
    private int chromStart = -1;
    private int chromEnd = -1;
    private Map<Integer, String> chromInfoMap;
    private List<Integer> chromLength;
    
        
    private int isizeUpper;
    private int isizeLower;
    private int readLen;
    private int fragMean;
    
    private int minMapQ;
    private int clusteringDist; 
    
//    private int superitemMaxSpanRange;
    
    
    private Map<String, List<SAMRecord>> rpTracker;
    private ChannelParser breakChannel;
    private ChannelParser isizeLargeChannel;
    private ChannelParser isizeSmallChannel;
    private ChannelParser oemChannel;
    private ChannelParser oriChannel;
    
    private int superitemCount;
    // keep normal read per base every 1Mb length.
    private int readDepthContainerBuffer = 1000000;
    private int[] readDepthContainer = new int[readDepthContainerBuffer];
    
    private int numOfARPs;
    private int numOfRPs;
    
    
    
    public SignalReader(int fragMean, int fragStd, int cutStd, int readLen, int maxDist, int minMapQ){
        isizeUpper = fragMean + cutStd * fragStd;
        isizeLower = fragMean - cutStd * fragStd;
        
        this.readLen = readLen;
        this.fragMean = fragMean;
        this.minMapQ = minMapQ;
        this.clusteringDist = maxDist;
        if (maxDist == -1){
            this.clusteringDist = fragMean - 2 * readLen;
        }
//        superitemMaxSpanRange = 2 * fragMean;
        
        
    }
    /**
     * 
     * @param bamFile the alignment file
     * @param fastaFile the reference index file
     * @param chrom a user specified chromosome (optional)
     * @param chromStart start of a genome region (optional)
     * @param chromEnd end of a genome region (optional)
     * @param superitemOutPath output path of the created superitems
     * @param abnormalSigOut output path of the abnormal signals (optional)
     * @throws IOException 
     */
    public void doWork(String bamFile, String fastaFile, String chrom, int chromStart, int chromEnd, String superitemOutPath, String abnormalSigOut) throws IOException{
        startTime = System.currentTimeMillis();
        
        // the channel constructor need a name, it can be whatever you like, just used for naming some output file.
        breakChannel = new ChannelParser(clusteringDist, "break", abnormalSigOut);
        isizeLargeChannel = new ChannelParser(clusteringDist, "isize_large", abnormalSigOut);
        isizeSmallChannel = new ChannelParser(clusteringDist, "isize_small", abnormalSigOut);
        oemChannel = new ChannelParser(clusteringDist, "oem", abnormalSigOut);
        oriChannel = new ChannelParser(clusteringDist, "ori", abnormalSigOut);
        
        
        if (!superitemOutPath.isEmpty()){
            createSuperItemWriter(superitemOutPath);            
        }
        extractSignalsFromBAM(bamFile, fastaFile, chrom, chromStart, chromEnd);
                
        endTime = System.currentTimeMillis();
        writer.close();
        printSuperItemGeneratorStats();
            
    }
    
    private void extractSignalsFromBAM(String bamFile, String fastaFile, String chrom, int regionStart, int regionEnd) throws IOException{        
               
        rpTracker = new HashMap<>();
        CHROM = chrom;
        chromStart = regionStart;
        chromEnd = regionEnd;
        
        boolean singleChrom = false;
        if (!chrom.isEmpty()){
            singleChrom = true;
        }
        final SamReader samReader = openBAMReader(bamFile, ValidationStringency.SILENT, false);
               
        // access user specified region
        if (singleChrom){  
            readFastaFile(fastaFile);
            SAMFileHeader samFileHeader = samReader.getFileHeader();
            
            SAMSequenceDictionary sequenceDictionary = samFileHeader.getSequenceDictionary();
            SAMSequenceRecord refSequenceRecord = sequenceDictionary.getSequence(chrom);
            int refSequenceLength = refSequenceRecord.getSequenceLength();           
            int nWindows = refSequenceLength / readDepthContainerBuffer;
            int windowStart = 0;
            int windowEnd = 0;
            
            if (chromStart != 0 && chromEnd != 0){
                refSequenceLength = chromEnd - chromStart + 1;
                if (refSequenceLength <= readDepthContainerBuffer){
                    nWindows = 1;
                }
                windowStart = chromStart;
                windowEnd = chromEnd;
            }
            
            int[] readDepthPreStepBuffer = new int[readDepthContainerBuffer];
            
            for (int i = 0; i < nWindows; i++){  
                windowEnd = windowStart + readDepthContainerBuffer;                 
                SAMRecordIterator iterator = samReader.query(chrom, windowStart, windowEnd, false);               
                
                analysisAlignment(iterator, windowStart);                
                processRemainingSignals();
                int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
                System.out.println("processed region: [" + windowStart + ", " + windowEnd + "] " + "#superitems: " + curBinSuperitemCount);
                windowStart = windowEnd;
                
//                readDepthPreStepBuffer = copyFromReadDepthBuffer();
                readDepthPreStepBuffer = readDepthContainer;
                
                readDepthContainer = new int[readDepthContainerBuffer];
                if (writer != null){
                    writeAllSuperItems();
                }    
            }
            // process remaining alignment in BAM
            SAMRecordIterator iterator = samReader.query(chrom, windowStart, refSequenceLength, false);
            analysisAlignment(iterator, windowStart);
            processRemainingSignals();
            int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
            System.out.println("processed region: [" + windowStart + ", " + refSequenceLength + "] " + "#superitems: " + curBinSuperitemCount);
            if (writer != null){
                writeAllSuperItems();
            }
           
        } 
        // read whole genome
        else{
            readFastaFile(fastaFile);
            int length = chromLength.size();
            SAMRecordIterator iterator;
            for (int i = 0;i < length; i ++){
                int refSequenceLength = chromLength.get(i);
                String curChrom = chromInfoMap.get(i);
                System.out.println("Start processing chrom: " + curChrom + ", chrom length: " + refSequenceLength);
                
                int nWindows = refSequenceLength / readDepthContainerBuffer;
                int windowStart = 0;
                int windowEnd = 0;
                int[] readDepthPreStepBuffer = new int[readDepthContainerBuffer];
                for (int k = 0; k < nWindows; k++){
                    windowEnd = windowStart + readDepthContainerBuffer;                 
                    iterator = samReader.query(curChrom, windowStart, windowEnd, false);
                    
                    analysisAlignment(iterator, windowStart);                
                    processRemainingSignals();
                    int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
                    System.out.println("processed region: [" + windowStart + ", " + windowEnd + "] " + "#superitems: " + curBinSuperitemCount);
                    windowStart = windowEnd;

                    readDepthPreStepBuffer = copyFromReadDepthBuffer();
                    readDepthContainer = new int[readDepthContainerBuffer];
                    if (writer != null){
                        writeAllSuperItems();
                        
                    }
                    
                    
                }
                iterator = samReader.query(curChrom, windowStart, refSequenceLength, false);
                analysisAlignment(iterator, windowStart);
                processRemainingSignals();
                int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
                System.out.println("processed region: [" + windowStart + ", " + refSequenceLength + "] " + "#superitems: " + curBinSuperitemCount);
                if (writer != null){
                    writeAllSuperItems();
                }  
                
                int mutSignalInChannels = breakChannel.getNumOfMut();
                mutSignalInChannels += isizeLargeChannel.getNumOfMut();
                mutSignalInChannels += isizeSmallChannel.getNumOfMut();
                mutSignalInChannels += oriChannel.getNumOfMut();
                mutSignalInChannels += oemChannel.getNumOfMut();
                System.out.println("Chrom: " + curChrom + ", #mutSignals remain:" + mutSignalInChannels);
                System.out.println("\n");
            }
        }                                               
    }        
      
    private void analysisAlignment(SAMRecordIterator iterator, int windowStart){
//        CigarOps corasenCigar = new CigarOps();
        while(iterator.hasNext()){
            SAMRecord record = iterator.next();
            int mapq = record.getMappingQuality();
            // Discard some alignment, including seconday, supplementary, low mapping quality.
                      
            
//            if (record.isSecondaryOrSupplementary() || mapq <= minMapQ){
//                continue;
//            }
//            if (record.getReadName().equals("H06HDADXX130110:2:2203:7358:54872")){
//                System.out.println("ssss");
//            }
            if (mapq <= minMapQ){
                continue;
            }
            List<CigarElement> cigarElements = record.getCigar().getCigarElements();
            // discard reads with clipped sequence longer than 70% of read length.
            if (badReads(cigarElements)){
                continue;
            }
            // count the number of normal read per base
            int isGoodAlign = goodAlignment(cigarElements);
            if (isGoodAlign != -1){
                updateReadDepthArray(record.getAlignmentStart(), isGoodAlign, windowStart);
            }                        
            if (!isInterChrom(record)){                
                
                SEClippedParser(record, cigarElements);   
                RPUnmappedParser(record);
                int incorrectSize = RPisizeParser(record, cigarElements);

                if (!rpTracker.containsKey(record.getReadName())){
                    List<SAMRecord> records = new ArrayList<>();
                    records.add(record);
                    rpTracker.put(record.getReadName(), records);
                }else{
                    List<SAMRecord> records = rpTracker.get(record.getReadName());
                    records.add(record);
                    numOfRPs += 1;
                    int incorrectOri = RPoriParser(records, cigarElements);
                    if (incorrectOri == -1){
//                        oriChannel.updateNumOfARPs();
                        numOfARPs += 1;
                    }
                    if (incorrectSize == -1){
//                        isizeLargeChannel.updateNumOfARPs();
                        numOfARPs += 1;
                    }
                    rpTracker.remove(record.getReadName());
                }                
            }
            
        } 
        iterator.close();
    }
    private int[] copyFromReadDepthBuffer(){
        int[] newBuffer = new int[readDepthContainerBuffer];
//        int startPosToCopy = readDepthContainerBuffer - readLen;
        for (int i = 0; i < readDepthContainerBuffer; i ++){
            int val = readDepthContainer[i];
//            newBuffer[i - startPosToCopy] = val;
            newBuffer[i] = val;
        }
        return newBuffer;
    }
    
    private boolean badReads(List<CigarElement> cigarElements){
        int totalClipped = 0;
        boolean isBad = false;
        for (CigarElement element : cigarElements){
            String operation = element.getOperator().toString();
            int optLength = element.getLength();
            if (operation.equals("S") || operation.equals("H")){
                totalClipped += optLength;
            }
        }
        if (totalClipped > 0.7 * readLen){
            isBad = true;
        }
        return isBad;
    }
    
    private int goodAlignment(List<CigarElement> cigarElements){
        if (cigarElements.size() == 1){
            String cigarOperation = cigarElements.get(0).getOperator().toString();
            int opLength = cigarElements.get(0).getLength();
            return cigarOperation.equals("M") ? opLength : -1;
        }
        else return -1;
    }
    
    private void updateReadDepthArray(int pos, int length, int windowStart){
        
        for (int i = 0; i < length ; i ++){
            if ( (pos + i - windowStart) >= readDepthContainerBuffer){
                continue;
            }
            else if (pos + i < windowStart) {
                continue;
            }
            else{
                readDepthContainer[pos + i - windowStart] += 1;            
            }
            
        }
        
    }
    private int isOverlapRP(List<SAMRecord> records){
        int overlap = -1;
        SAMRecord leftMostAlign = records.get(0);
        int leftMostAlignStart = leftMostAlign.getAlignmentStart();
        SAMRecord mateAlign = records.get(1);
        int mateAlignStart = mateAlign.getAlignmentStart();
        int tmp = mateAlignStart - leftMostAlignStart;
        if (tmp < readLen){
            overlap = tmp;
        }
        return overlap;
    }   
    
    private void SEClippedParser(SAMRecord record, List<CigarElement> cigarElements){
        // For a mapped read
        if (!record.getReadUnmappedFlag()){
            
            String firstOperation = cigarElements.get(0).getOperator().toString();
                        
            CigarOps coarsenCigar = new CigarOps();
            coarsenCigar.calQueryPosFromCigar(cigarElements, 1, record.getReadNegativeStrandFlag(),readLen);
            int qsPos = coarsenCigar.getqsPos();
            
            String cigarStr = coarsenCigar.getCigarStr();
            int mutCoord = record.getAlignmentStart();
            
            if (!cigarStr.equals("M") && !cigarStr.isEmpty()){
                if (firstOperation.equals("M")){                    
                    mutCoord += qsPos;                
                }
                if (coarsenCigar.isCoIDread() && firstOperation.equals("S")){
                    mutCoord += qsPos;
                }
                
                String ori = record.getReadNegativeStrandFlag() ? "-" : "+";                               
                MutSignal mutSignal = new MutSignal(record, cigarStr, mutCoord, ori);
//                MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                        record.getInferredInsertSize(), cigarStr, mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
                mutSignal.setIsizeNormal(isizeUpper, isizeLower);

                breakChannel.addSignals(mutSignal, fragMean, readLen);
            }
        }
        
    }
     /**
     * One end unmapped read
     * @param record
     * @return 
     */
    private int RPUnmappedParser(SAMRecord record){
        int isUnmappd = 0;
        // read unmapped
        if (record.getReadUnmappedFlag()){
            int mutCoord = record.getMateAlignmentStart();
            String ori = record.getMateNegativeStrandFlag() ? "-": "+";
            
            MutSignal mutSignal = new MutSignal(record, "ARP_OEM", mutCoord, ori);
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
           
            oemChannel.addSignals(mutSignal, fragMean, readLen);
            isUnmappd = -1;
//            oemChannel.addSignalsTest(mutSignal);
        }else if (record.getMateUnmappedFlag()){
            int mutCoord = record.getMateAlignmentStart();
            String ori = record.getMateNegativeStrandFlag() ? "-": "+";
            
            MutSignal mutSignal = new MutSignal(record, "ARP_OEM", mutCoord, ori);
//                MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                        record.getInferredInsertSize(), "ARP_OEM", mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
            
            
            oemChannel.addSignals(mutSignal, fragMean, readLen);
            isUnmappd = -1;
            
//            oemChannel.addSignalsTest(mutSignal);
        }
        return isUnmappd;
    }
    
    private int RPisizeParser(SAMRecord record, List<CigarElement> cigarElements){
        // only process read-pair mapped on the same chrom.
        
        int incorrectSize = 0;
        CigarElement leftMostCigarElement = cigarElements.get(0);
        String leftMostCigarOperator = leftMostCigarElement.getOperator().toString();

        int mutCoord = record.getAlignmentStart();
        if (leftMostCigarOperator.equals("M") && !record.getReadNegativeStrandFlag()){
            mutCoord += leftMostCigarElement.getLength();
        }

        int insertSize = record.getInferredInsertSize();

        String ori = record.getReadNegativeStrandFlag() ? "-" : "+";
        if (Math.abs(insertSize) >= isizeUpper){    
//                System.out.println(record.getReadName() + " isize: " +Math.abs(insertSize));

            MutSignal mutSignal = new MutSignal(record, "ARP_LARGE_INSERT", mutCoord, ori);               
//            MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                        record.getInferredInsertSize(), "ARP_LARGE_INSERT", mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
            isizeLargeChannel.addSignals(mutSignal, fragMean, readLen);
            incorrectSize = -1;
//                isizeLargeChannel.addSignalsTest(mutSignal);
        }
        else if (Math.abs(insertSize) <= isizeLower && insertSize != 0){

            MutSignal mutSignal = new MutSignal(record, "ARP_SMALL_INSERT", mutCoord, ori);       
//            MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                record.getInferredInsertSize(), "ARP_SMALL_INSERT", mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);

            isizeSmallChannel.addSignals(mutSignal, fragMean, readLen); 
//                isizeSmallChannel.addSignalsTest(mutSignal);
            incorrectSize = -1;
        }
        return incorrectSize;
        
    }
    
    private int RPoriParser(List<SAMRecord> records, List<CigarElement> cigarElements){
        int incorrectOri = 0;
        SAMRecord leftMostRecord = records.get(0);
        SAMRecord rightMostRecord = records.get(records.size() - 1);
        // For read-pair, it should be proper paired. Its read and mate are all mapped.
        if (leftMostRecord.getReadPairedFlag() && !leftMostRecord.getReadUnmappedFlag() && !leftMostRecord.getMateUnmappedFlag()){
            int mutCoord = leftMostRecord.getAlignmentStart();
            if (leftMostRecord.getReadNegativeStrandFlag()== leftMostRecord.getMateNegativeStrandFlag()){
                String mutType;
                String ori = leftMostRecord.getReadNegativeStrandFlag() ? "-":"+";
                if (leftMostRecord.getReadNegativeStrandFlag()){
                    mutType = "ARP_RR";
                }else{
                    CigarElement leftMostCigarElement = cigarElements.get(0);
                    String leftMostCigarOperation = leftMostCigarElement.getOperator().toString();

                    if (leftMostCigarOperation.equals("M")){
                        mutCoord += leftMostCigarElement.getLength();
                    }
                    mutType = "ARP_FF";
                }
                MutSignal readMutSignal = new MutSignal(leftMostRecord, mutType, mutCoord, ori);
//                MutSignal readMutSignal = new MutSignal(leftMostRecord.getReadName(), leftMostRecord.getReferenceIndex(), leftMostRecord.getReferenceName(), 
//                        leftMostRecord.getInferredInsertSize(), mutType, mutCoord, ori, leftMostRecord.getAlignmentStart(), leftMostRecord.getMateAlignmentStart());
                readMutSignal.setIsizeNormal(isizeUpper, isizeLower);

                MutSignal mateMutSignal = new MutSignal(leftMostRecord, mutType, rightMostRecord.getMateAlignmentStart(), ori);
//                MutSignal mateMutSignal = new MutSignal(leftMostRecord.getReadName(), leftMostRecord.getReferenceIndex(), leftMostRecord.getReferenceName(), 
//                        leftMostRecord.getInferredInsertSize(), mutType, mutCoord, ori, leftMostRecord.getAlignmentStart(), leftMostRecord.getMateAlignmentStart());
                mateMutSignal.setIsizeNormal(isizeUpper, isizeLower);   

                oriChannel.addSignals(readMutSignal, fragMean, readLen);
                oriChannel.addSignals(mateMutSignal, fragMean, readLen);
                incorrectOri = -1;
            }
            else if (leftMostRecord.getReadNegativeStrandFlag() && !leftMostRecord.getMateNegativeStrandFlag() && isOverlapRP(records) == -1){
                String mutType = "ARP_RF";

                MutSignal readMutSignal = new MutSignal(leftMostRecord, mutType, mutCoord, "-");
//                MutSignal readMutSignal = new MutSignal(leftMostRecord.getReadName(), leftMostRecord.getReferenceIndex(), leftMostRecord.getReferenceName(), 
//                        leftMostRecord.getInferredInsertSize(), mutType, mutCoord, "-", leftMostRecord.getAlignmentStart(), leftMostRecord.getMateAlignmentStart());
                readMutSignal.setIsizeNormal(isizeUpper, isizeLower);

//                oriChannel.addSignalsTest(readMutSignal);
                MutSignal mateMutSignal = new MutSignal(leftMostRecord, mutType, rightMostRecord.getMateAlignmentStart(), "+");
//                MutSignal mateMutSignal = new MutSignal(leftMostRecord.getReadName(), leftMostRecord.getReferenceIndex(), leftMostRecord.getReferenceName(), 
//                leftMostRecord.getInferredInsertSize(), mutType, mutCoord, "-", leftMostRecord.getAlignmentStart(), leftMostRecord.getMateAlignmentStart());
                mateMutSignal.setIsizeNormal(isizeUpper, isizeLower);

                oriChannel.addSignals(readMutSignal, fragMean, readLen);
                oriChannel.addSignals(mateMutSignal, fragMean, readLen);
//                oriChannel.addSignalsTest(mateMutSignal);
                incorrectOri = -1;
            }
        }
        return incorrectOri;
    }

    
    
    private boolean isInterChrom(SAMRecord record){
        
        if (!record.getReferenceName().equals(record.getMateReferenceName())){
            return true;
        }
        else return false;
    }

    private SamReader openBAMReader(String bamFile, ValidationStringency stringency, boolean includeFileSource) throws IOException{
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(stringency).enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
        if(includeFileSource){
            samReaderFactory.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS);            
        }
//        samReaderFactory.samRecordFactory(DefaultSAMRecordFactory.getInstance());
        final SamReader samReader = samReaderFactory.open(new File(bamFile));
        return samReader;
    }
    
    private void readFastaFile(String fastaIndexFile) throws IOException{
        chromInfoMap = new TreeMap<>();
        chromLength = new ArrayList<>();
        FileInputStream fin = new FileInputStream(new File(fastaIndexFile));
        BufferedReader myInput = new BufferedReader(new InputStreamReader(fin));
        String thisLine;
        while ((thisLine = myInput.readLine()) != null){
            String[] tokens = thisLine.split("\t");
            String refSeq = tokens[0];
            // escape "M", "MT", "chrM"
            if (refSeq.contains("M") || refSeq.contains("_")){
                continue;
            }
                    
            int refLength = Integer.parseInt(tokens[1]);
            chromLength.add(refLength);
            chromInfoMap.put(chromLength.size() - 1, refSeq);
            // Read until Y chromosome, ignore other contigs.
            if (refSeq.equals("Y") || refSeq.equals("chrY")){
                break;
            }
        }
        
    }
       
    private void processRemainingSignals() {
        breakChannel.processFinalSignals(fragMean, readLen);
        isizeLargeChannel.processFinalSignals(fragMean, readLen);
        isizeSmallChannel.processFinalSignals(fragMean, readLen);
        oemChannel.processFinalSignals(fragMean, readLen);
        oriChannel.processFinalSignals(fragMean, readLen);
                             
    }
    
    private void writeAllSuperItems() throws IOException{
        breakChannel.writeSuperItemsInChannel(writer);
        isizeLargeChannel.writeSuperItemsInChannel(writer);
        isizeSmallChannel.writeSuperItemsInChannel(writer);
        oemChannel.writeSuperItemsInChannel(writer);
        oriChannel.writeSuperItemsInChannel(writer);
    }
    /**
     * calculate normal read aligned at a specific position and the number of superitems that generated within this window.
     * @param windowStart
     * @param windowSize
     * @param preReadDepthBuffer
     * @return 
     */
    private int assignReadDepthAndCountSuperItem(int windowStart, int windowSize, int[] preReadDepthBuffer){
        breakChannel.setSuperitemWeightRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        isizeLargeChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        isizeSmallChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        oriChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        oemChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        
        int superitemCount = 0;
        superitemCount += breakChannel.getSuperitemCount();
        superitemCount += isizeLargeChannel.getSuperitemCount();
        superitemCount += isizeSmallChannel.getSuperitemCount();
        superitemCount += oemChannel.getSuperitemCount();
        superitemCount += oriChannel.getSuperitemCount();
        
        return superitemCount;
    }
    private void createSuperItemWriter(String superitemOutPath) throws IOException{
        writer = new BufferedWriter(new FileWriter(superitemOutPath));
        writer.write("type\tchromIdx\tnread\tpos\tsplitAlignPos\tori\tweight\tratio\tregion\tmateRegion\tqnames\n");
    }
    
    public int getWGARPNum(){
        return numOfARPs;
    }    
    public void printSuperItemGeneratorStats(){
        StringBuilder sb = new StringBuilder();
        sb.append("\n==============  SuperItem Generation =============\n");
        sb.append("Time: " + (endTime - startTime) + "ms");  
        sb.append("\nTotal superitems: " + superitemCount);   
        sb.append("\nTotal number of read-pairs: " + numOfRPs);
        sb.append("\nTotal number of doiscordant read-pairs: " + numOfARPs);
        sb.append("\n======================================================");
        
        System.out.println(sb.toString());
    }
    
}
