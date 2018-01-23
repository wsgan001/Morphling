/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package contiguousfspm;

import java.io.IOException;

import superitemGenerator.SignalReader;
import dataStructures.*;
import htsjdk.samtools.QueryInterval;
import java.io.OutputStream;
import java.util.*;


/**
 *
 * @author jiadonglin
 */
public class TestContiguousFSPM {
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException{
        // TODO code application logic here
                                                    
        String workingDir = "/Users/jiadonglin/SV_data/GIAB/";
//        String workingDir = "/Users/jiadonglin/SV_data/NA12878/";

        
        String bamFile = workingDir + "10X/HG001.hs37d5.10X.bam";
//        String bamFile = workingDir + "30X/NA12878.30X.bam";
        String fastaIndexFile = "/Users/jiadonglin/SV_data/1K_Project_Ref/hs37d5.fa.fai";
        String superitemOut = workingDir + "10X/FSP/all.superitems.unsorted.mapq10.txt";
//        String superitemOut = workingDir + "50X/NA12878/all.superitems.unsorted.mapq10.chr1.txt";
        
        String patternOut = workingDir + "10X/FSP/frequent.patterns.mapq10.sup30.newBP.txt";
        String superitemPatternOut = workingDir + "10X/FSP/frequent.patterns.linked.sup30.newBP.txt";
        String svRegionOut = workingDir + "10X/FSP/svRegion.sup30.out";
       
        
        String syntheticDir = "/Users/jiadonglin/SV_data/synthetics/";
        
        String syntheticBamFile = syntheticDir + "rd50_align/venter.rd50.bwamem.aln.sorted.bam";
        String syntheticSuperItemOut = syntheticDir + "rd10_align/all.superitems.unsorted.mapq10.txt";
        String syntheticPatternOut = syntheticDir + "rd10_align/frequent.patterns.mapq10.sup30.txt";
        String syntheticSuperItemPatternOut = syntheticDir + "rd10_align/frequent.patterns.linked.sup30.txt";
        String syntheticRegionOut = syntheticDir + "rd10_align/svRegion.sup30.bed";
        
//        String bamFile = "";
//        String fastaIndexFile = "";
//        String superitemOut = "";
        String signalOut = "";
//        
        int fragMean = 556;
//        int fragMean = 425; // 1000GP
//        int fragMean = 500; // synthetic

        int readLen = 148;
//        int readLen = 250; // 1000GP
//        int readLen = 150; // synthetic
    
        int fragStd = 160; 
//        int fragStd = 147; // 1000GP
//        int fragStd = 50; // synthetic

        int cutStd = 3;
        int clusteringDist = readLen;

        int minMapQ = 10;
        
        String chr = "";      
        int regionS = 0;
        int regionE = 0;
        
        for (int i = 0; i < args.length; i++){            
            String[] argTokens = args[i].split("=");
            if (argTokens[0].equals("readLen")){
                readLen = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("fragMean")){
                fragMean = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("fragStd")){
                fragStd = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("cutStd")){
                cutStd = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("maxD")){
                clusteringDist = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("minQ")){
                minMapQ = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("chrom")){
                chr = argTokens[1];
            }
            if (argTokens[0].equals("bamFile")){
                bamFile = argTokens[1];
            }
            if (argTokens[0].equals("faiFile")){
                fastaIndexFile = argTokens[1];
            }
            if (argTokens[0].equals("siOut")){
                superitemOut = argTokens[1];
            }
            if (argTokens[0].equals("sigOut")){
                signalOut = argTokens[1];
            }
        }
        System.out.println("fragMean: " + fragMean + " fragStd: " + fragStd + " readLen: " + readLen + " clustD: " + clusteringDist);

//        SignalReader myReader = new SignalReader(fragMean, fragStd, cutStd, readLen, clusteringDist, minMapQ);
//        myReader.doWork(bamFile, fastaIndexFile, chr, regionS, regionE, superitemOut, signalOut); 
       


        SequenceDatabase sequenceDatabase = new SequenceDatabase(); 
        int minSup = 30;
        int patternMaxRegionSpan = 3 * fragMean;
           
        sequenceDatabase.loadSequencesFromFile(superitemOut);
//        sequenceDatabase.loadSample(workingDir + "50X/NA12878");
        
        ContiguousFSPM algoContiguousFSPM = new ContiguousFSPM(minSup, patternMaxRegionSpan);
        algoContiguousFSPM.runAlgorithm(sequenceDatabase, patternOut, superitemPatternOut, svRegionOut);
        algoContiguousFSPM.printAlgoStatistics();

    }
    
}
