/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package superitemGenerator;

import htsjdk.samtools.CigarElement;
import java.util.List;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author jiadonglin
 */
public class CigarOps {
    
    private String cigarString;
    private int qsPos;
    // indicate co-occurence of I and D on a single read
    private boolean coID = false;
    public CigarOps(){
        
    }
    
    public String getCigarStr(){
        return this.cigarString;
    }
    public int getqsPos(){
        return this.qsPos;
    }
    // not a good way, but it is ok to reverse cigar operations which is usually short.
    private List reverseList(List myList) {
        List invertedList = new ArrayList<>();
        for (int i = myList.size() - 1; i >= 0; i--) {
            invertedList.add(myList.get(i));
        }
        return invertedList;
    }
    
    public void calQueryPosFromCigar(List<CigarElement> cigarElements, int nIgnore, boolean isReverse, int readLen){
//        List<CigarElement> cigarElements = cigar.getCigarElements();
//        if (isReverse){
//            cigarElements = reverseList(cigarElements);
//        }
        int qsPos = 0;        
        String cigarStr = "";
        int IDnum = 0;  // number of 'I' 'D' operation in cigar
        int lastIDindex = -1;
        int firstIDindex = -1;
        int longestMatch = 0;
        int operators = cigarElements.size();
        for (int i = 0; i < operators; i++){
            CigarElement cigarElement = cigarElements.get(i);
            String operation = cigarElement.getOperator().toString();
            
            operation = operation.equals("H") ? "S":operation;
            int opLength = cigarElement.getLength();             
            longestMatch = (operation.equals("M") && opLength > longestMatch) ? opLength : longestMatch; 
            if (i == 0){                
                cigarStr += operation;
                if (operation.equals("M")){
                    qsPos += opLength;
                }                
            }
            else if (i > 0 && (operation.equals("I")||operation.equals("D"))){
                
                lastIDindex = i;
                if (opLength <= nIgnore){
                    qsPos += opLength;
                }else{
                    cigarStr += operation;
                    IDnum += 1;                    
                }                    
                
                if (firstIDindex == -1) {
                    firstIDindex = i;
                }
            }else if (!cigarStr.isEmpty()){
                Character lastChar = cigarStr.charAt(cigarStr.length() - 1);
                Character curChar = operation.charAt(0);
                // = 0, two char same
                if (curChar.compareTo(lastChar) == 0){
                    continue;
                }else if (opLength > nIgnore){
                    cigarStr += operation;
                }
            }
        }
 
        
        if (IDnum > 1 && longestMatch > 0.5 * readLen){
            coID = true;
            if (isReverse){
                int distToLastID = 0;
                for (int i = 0; i < lastIDindex ; i++){
                    CigarElement cigarElement = cigarElements.get(i);
                    String operation = cigarElement.getOperator().toString();
                    operation = operation.equals("H") ? "S":operation;
                    int opLength = cigarElement.getLength(); 
                    if (i == 0 && operation.equals("M")){
                        distToLastID += opLength;
                    }else{
                        distToLastID += opLength;
                    }                    
                }
                cigarStr = "SM";
                qsPos = distToLastID;
            }  
            else{
                cigarStr = "MS";
            }
        }
        

        this.cigarString = cigarStr; 
        this.qsPos = qsPos;
    }
    /**
     * 'I' 'D' operation occur together in Cigar, which is not good.
     * @return 
     */
    public boolean isCoIDread(){
        return coID;
    }
    
    
    
}
