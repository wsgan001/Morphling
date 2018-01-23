/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package superitemGenerator;
import java.util.List;
import java.util.ArrayList;

/**
 *
 * @author jiadonglin
 */
public class svOutInfo {
     
    int start;
    int end;
    String pattern;
    int linkType; // 1 mate linked. -1 self linked. 0 non-linkable
    int linkSup = -1;
    List<Integer> weights;
    List<Integer> postions;
    
    public svOutInfo(int s, int e, String patternStr, int linkFlag, int supLink, List<Integer> weights, List<Integer> pos){
        start = s;
        end = e;
        pattern = patternStr;
        linkType = linkFlag;
        linkSup = supLink;
        this.weights = weights;
        postions = pos;
    }
    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        sb.append(start);
        sb.append("\t");
        sb.append(end);
        sb.append("\t");
        
        sb.append(infoToString());
        return sb.toString();
    }
    private String infoToString(){
        StringBuilder sb = new StringBuilder();
        sb.append("LinkType=");
        if (linkType == 1){
            sb.append("Mate");        
        }
        if(linkType == 0){
            sb.append("None");
        }
        if(linkType == -1){
            sb.append("Self");
        }
        sb.append(";");
        sb.append("Pattern=");
        sb.append(pattern);
        sb.append(";");
        sb.append("LinkSup=");
        if (linkType == 0){
            sb.append("Supp");
        }else{
            sb.append(linkSup);
        }        
        sb.append(";");
        sb.append("Weights=");
        sb.append(ListToString(weights));
        sb.append(";");
        sb.append("Pos=");
        sb.append(ListToString(postions));
        return sb.toString();
    }
    private String ListToString(List<Integer> alist){
        StringBuilder sb = new StringBuilder();
        for (Integer ele : alist){
            sb.append(ele);
            sb.append(",");
        }
        String str = sb.toString();
        String outStr = str.substring(0, str.length() - 1);
        return outStr;
    }
    
//    private String ratioToString(){
//        StringBuilder sb = new StringBuilder();
//    }
}
