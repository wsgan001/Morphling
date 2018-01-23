/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dataStructures;

/**
 * A class to keep identify the sequence.
 * @author jiadonglin
 */
public class ItemSeqIdentifier {
    
    private int seqID;
    private int subSeqStart;
    // Item index in the sequence in terms of itemset index.
    
    private int itemSetIdx;
    private int itemIdx;
    // Allow sequence identifier of strings.
//    private String sampleName;

    public ItemSeqIdentifier(int seqID, int subseqStart, int itemsetIdx, int itemIdx) {
        this.seqID = seqID;
        this.subSeqStart = subseqStart;
        this.itemSetIdx = itemsetIdx;
        this.itemIdx = itemIdx;
    }
    public int getSeqID(){
        return this.seqID;
    }
    public int getSubSeqID(){
        return this.subSeqStart;
    }
    public int getItemSetIdx(){
        return this.itemSetIdx;
    }
    public int getItemIdx(){
        return this.itemIdx;
    }
    
}
