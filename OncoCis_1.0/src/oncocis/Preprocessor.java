/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package OncoCis;

/**
 *
 * @author Dilmi
 */
public class Preprocessor {
    private String mutations;
    private String geneExpression;
    private String outputFolder;
    private String mutationFileError;
    private String geneExpFileError;
    private int totalNumofLines = 0;

    public Preprocessor(String cMutations, String cGeneExpession, String cOutputFolder) {
        mutations = cMutations;
        geneExpression = cGeneExpession;
        outputFolder = cOutputFolder;
        System.out.println("Files being checked for errors");
        if(isTextFile(mutations, outputFolder)){
            if(isTextFile(geneExpression, outputFolder)){
                if(isBed()){
                    mutationFileError = "No Errors";  
                }
                if(checkGeneExpFile()){
                    geneExpFileError = "No Errors";  
                }
            }
            else{
            geneExpFileError = "Gene Expression file: Incorrect file format";
                System.out.println(""+geneExpFileError);
            }
        }
        else{
            mutationFileError = "Mutations file: Incorrect file format";
            System.out.println(""+mutationFileError);
        }
        
    }
    
    public Preprocessor(String cMutations, String cOutputFolder) {
        mutations = cMutations;
        outputFolder = cOutputFolder;
        if(isTextFile(mutations, outputFolder)){
          if(isBed()){
          mutationFileError = "No Errors";  
            }
          else{
              System.out.println(""+mutationFileError);
          }
        }
        else{
            mutationFileError = "Mutations file is not a text file";
            System.out.println(""+mutationFileError);
        }
        
        
    }
    
    public static void main(String[] args) {
        
//        System.out.println(Preprocessor.isTextFile("/home/bioinfo/Darkmatter-data/DarkMatter/UserInput/test/Mutation-input-test.bed", 
//                "/home/bioinfo/Darkmatter-data/DarkMatter/UserInput/test"));
//        
//        Preprocessor eh = new Preprocessor(
//        "/home/bioinfo/Darkmatter-data/DarkMatter/UserInput/test/Mutation-input-test.bed",
//        "home/bioinfo/Darkmatter-data/DarkMatter/UserInput/test/GeneExpression-test.bed",
//        "/home/bioinfo/Darkmatter-data/DarkMatter/UserInput/test");
//        if(eh.isBed()){
//            System.out.println("Is bed");
//        }
//        else{
//            System.out.println(""+eh.getMutationFileError());
//        }
//        
        
    }
    
    public static boolean isTextFile(String file, String outputFolder){
        boolean isText = false;
        UtilitiesLinux.shellCommandExecuter("file " +file+ " > " +outputFolder+"temp-fileType.txt");
        String fileFormat[] = DataSetReader.readLines(outputFolder+"temp-fileType.txt");
        System.out.println(""+fileFormat[0]);
        if(fileFormat[0].contains("ASCII text")
        ||fileFormat[0].contains("UTF-8 text")){
            isText=true;
        }
        
       
        return isText;
    }
    
    
    
    public boolean isBed(){
        boolean isBed = false;
        int numOfCols = getNumOfCols(mutations, outputFolder);
        System.out.println(""+numOfCols);
        String mutationsFile[][] = null;
        if(numOfCols==7){
            if(checkNumOfCols(mutations, outputFolder)){
                mutationsFile=DataSetReader.readDataSet(this.mutations, (numOfCols+1), "\t");        
                totalNumofLines = mutationsFile.length;
                for (int i = 0; i < mutationsFile.length; i++) {

                    //Chromosome
                    if(mutationsFile[i][0].regionMatches(0, "chr", 0, 3) &&(
                            (isNumericInt(mutationsFile[i][0].substring(3)) && 
                            Integer.parseInt(mutationsFile[i][0].substring(3))>0 &&
                            Integer.parseInt(mutationsFile[i][0].substring(3))<23)||
                            mutationsFile[i][0].substring(3).equals("X")||
                            mutationsFile[i][0].substring(3).equals("Y"))){                    

                        //Start and Stop
                        if(isNumericInt(mutationsFile[i][1])&&isNumericInt(mutationsFile[i][2])){

                            //Type
                            if(mutationsFile[i][4].equalsIgnoreCase("Substitution")||
                                mutationsFile[i][4].equalsIgnoreCase("Insertion")||
                                mutationsFile[i][4].equalsIgnoreCase("Deletion")){                            

                                //Reference Sequence
                                
                                if((mutationsFile[i][6].equals("."))&& mutationsFile[i][4].equalsIgnoreCase("Insertion")){
                                        isBed=true;
                                }                                
                                else{
                                    char refseq[] = mutationsFile[i][5].toCharArray();
                                    for (int j = 0; j < refseq.length; j++) {
                                       if(refseq[j]=='A'||refseq[j]=='C'||refseq[j]=='G'||refseq[j]=='T'){
                                            isBed=true;
                                        }
                                       else{
                                           isBed=false;
                                           mutationFileError="Error in line " + (i+1) + " of Mutation file";
                                           break;
                                       }
                                    }
                                }
                                
                                if(isBed){

                                    //Mutated Sequence
                                    if((mutationsFile[i][6].equals("."))&& mutationsFile[i][4].equalsIgnoreCase("Deletion")){
                                        isBed=true;
                                    }
                                    else{
                                        char mutseq[] = mutationsFile[i][6].toCharArray();
                                    
                                        for (int j = 0; j < mutseq.length; j++) {
                                           if(mutseq[j]=='A'||mutseq[j]=='C'||mutseq[j]=='G'||mutseq[j]=='T'){
                                                isBed=true;
                                            }
                                           else{
                                               isBed=false;
                                               mutationFileError="Error in line " + (i+1) + " of Mutation file";
                                               break;
                                           }

                                        }                                    
                                    }
                                    
                                    if(isBed){
                                        String temp = mutationsFile[i][6];
                                        mutationsFile[i][7]=temp;
                                        temp = mutationsFile[i][5];
                                        mutationsFile[i][6]=temp;
                                        temp = mutationsFile[i][4];
                                        mutationsFile[i][5]=temp;
                                        temp = mutationsFile[i][3];
                                        mutationsFile[i][4]=temp;
                                        mutationsFile[i][3] = String.valueOf(i+1);
                                        
                                    }
                                    else{
                                        mutationFileError="Error in line " + (i+1) + " of Mutation file";
                                        break;
                                    }

                                }
                                else{
                                    isBed=false;
                                    mutationFileError="Error in line " + (i+1) + " of Mutation file";
                                    break;
                                }
                            }
                            //type
                            else{
                                isBed=false;
                                mutationFileError="Error in line " + (i+1) + " of Mutation file";
                                break;
                            }

                        }
                        //Start and Stop
                        else{
                            isBed=false;
                            mutationFileError="Error in line " + (i+1) + " of Mutation file";
                            break;
                        }
                    }
                    //chr
                    else{
                        isBed=false;
                        mutationFileError="Error in line " + (i+1) + " of Mutation file";
                        break;

                    }

                }

            }
            
        }
        else{
          isBed = false;
          mutationFileError="Incorrect number of columns in Mutation file";  
        }
        if(isBed){
            DataSetWriter.writeToFile2D(mutations,mutationsFile);
        }
        return isBed;
    }
    
    public boolean checkGeneExpFile(){
        boolean geneExpFileCheck = false;
        int numofCols = getNumOfCols(geneExpression, outputFolder);
        String geneExpFile[][];
        if(checkNumOfCols(geneExpression, outputFolder, numofCols)){            
            geneExpFile = DataSetReader.readDataSet(this.geneExpression, numofCols, "\t");            
            for (int i = 1; i < geneExpFile.length; i++) {
                for (int j = 1; j < numofCols; j++) {
                    if(!isNumericDouble(geneExpFile[i][j])){
                       geneExpFileCheck = false;
                       break; 
                    }
                    else{
                        geneExpFileCheck = true;
                    }
                }
                if(!geneExpFileCheck){
                    break;
                }
            }
        }
        
        
        
        return geneExpFileCheck;
        
    }
    
    public static boolean isNumericInt(String str){  
          try  
          {  
            int i = Integer.parseInt(str);  
          }  
          catch(NumberFormatException nfe)  
          {  
            return false;  
          }  
          return true;  
        }
    
    public static boolean isNumericDouble(String str){  
          try  
          {  
            double d = Double.parseDouble(str);  
          }  
          catch(NumberFormatException nfe)  
          {  
            return false;  
          }  
          return true;  
        }

    public boolean checkNumOfCols(String file, String outputFolder, int numOfCols){
        boolean checkNumOfCols = false;
        UtilitiesLinux.shellCommandExecuter("awk '{print NF}' < "+file+ " > " + outputFolder+"numofCols");
        String colNum[][] = DataSetReader.readDataSet(outputFolder+"numofCols", 1, "\t");
        numOfCols=Integer.parseInt(colNum[0][0]);
        for (int i = 0; i < colNum.length; i++) {
            if(numOfCols==Integer.parseInt(colNum[i][0])){
                checkNumOfCols = true;
                continue;
            }
            else{
                checkNumOfCols=false;
                geneExpFileError = "Incorrect number of columns in line "+i+" of gene expression file";
                break;
            }
        }
        return checkNumOfCols;
    }
    public boolean checkNumOfCols(String file, String outputFolder){
        boolean checkNumOfCols = false;
        UtilitiesLinux.shellCommandExecuter("awk '{print NF}' < "+file+ " > " + outputFolder+"numofCols");
        String colNum[][] = DataSetReader.readDataSet(outputFolder+"numofCols", 1, "\t");
        int numOfCols=7;
        for (int i = 0; i < colNum.length; i++) {
            if(numOfCols==Integer.parseInt(colNum[i][0])){
                checkNumOfCols = true;
                continue;
            }
            else{
                checkNumOfCols=false;
                mutationFileError = "Incorrect number of columns in line "+i+" of Mutations file";
                break;
            }
        }
        return checkNumOfCols;
    }
    
    public static int getNumOfCols(String file, String outputFolder){
        int numOfCols = 0;
        UtilitiesLinux.shellCommandExecuter("awk '{print NF}' < "+file+ " > " + outputFolder+"numofCols");
        String colNum[][] = DataSetReader.readDataSet(outputFolder+"numofCols", 1, "\t");
        numOfCols=Integer.parseInt(colNum[0][0]);
//        for (int i = 0; i < colNum.length; i++) {
//            if(numOfCols==Integer.parseInt(colNum[i][0])){
//                continue;
//            }
//            else{
//                numOfCols=-1;
//                break;
//            }
//        }
        return numOfCols;
    }
    
    public static void removeRows(){
        UtilitiesLinux.shellCommandExecuter("");
    }

    public String getMutationFileError() {
        return mutationFileError;
    }

    public void setMutationFileError(String mutationFileError) {
        this.mutationFileError = mutationFileError;
    }

    public String getGeneExpFileError() {
        return geneExpFileError;
    }

    public void setGeneExpFileError(String geneExpFileError) {
        this.geneExpFileError = geneExpFileError;
    }

    public int getTotalNumofLines() {
        return totalNumofLines;
    }

    public void setTotalNumofLines(int totalNumofLines) {
        this.totalNumofLines = totalNumofLines;
    }
    
    
}
