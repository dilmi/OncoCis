/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package OncoCis;
import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

/**
 *
 * @author dilmi
 */
public class OncoCis {
    public static String MAIN_FOLDER;//e.g "/home/usr/OncoCis/";
    private String OUTPUT_FOLDER;
    private String MUTATIONS;      
    private String DNASE1HS;
    private String H3K4ME1;
    private String H3K4ME3;    
    private String H3K27AC;
    private String CONSERVATION = MAIN_FOLDER + "default/PhastCons/BigWig/";//default
    private static String HG19fa = MAIN_FOLDER + "default/hg19.fa";//default
    private static String GENEMAP = MAIN_FOLDER + "default/regDoms1000kb.bed";//default
    private boolean HASGENEEXP = false;
    private String GENEEXP;
    private String CELLTYPE;
//    private String requestID;
       
    public static void main(String[] args){
        dbConnect.connectDriver();
         //Automatically get current file location i.e. links.txt file is in the same directory as the .jar file
        File directory = new File (".");
        String link = "";
        try {
            //Defaults files
            link = MAIN_FOLDER+"/links.txt";             
            //Get user inputs
            String links[][] = DataSetReader.readDataSet(link, 1, "\t");
            
            
            String request[]=new String[10];
            request[0]=System.currentTimeMillis() + "";
            request[1]=System.currentTimeMillis() + "";
            request[2] = links[0][0];
            request[3] = "";
            if(links[1][0].equals("-")){
               request[4] = "0";
               request[5] = "null";
            }
            else{
               request[4] = "1";
               request[5]=links[1][0]; 
            }
            request[6] = System.currentTimeMillis()+"";
            System.out.print(""+request[0]);//id-queue
        System.out.print("\t"+request[1]);//request id
        System.out.print("\t"+request[2]);//mutation file
        System.out.print("\t"+request[3]);//cell type
        System.out.print("\t"+request[4]);//gene exp flag
        System.out.print("\t"+request[5]);//gene exp file
        System.out.println("\t"+request[6]);//timestamp
        
        init(request, links);
        }
        catch(Exception e) {
            System.out.println("The following error occured : "+e.getMessage());
        }
        
//       String request[] = getRequest();
//        writeOutput("dark_matter_temp.Final_results", "myTest", "/var/www/darkmatter/Output/");
       
//        String flag=null;
//        while(true){
//            String request[] = getRequest();
//            flag=request[8];         
//            if(flag!=null && flag.equals("0")){
//                System.out.println("Executing Request ID: "+request[1]);
//                boolean check1 = new File(MAIN_FOLDER + "Input/", request[2]).exists();
//                boolean check2 = new File(MAIN_FOLDER + "Input/", request[5]).exists();
////                if(!check1||!(request[4].equals("1")&&check2)){
////                   UtilitiesLinux.shellCommandExecuter("cp /home/bioinfo/Darkmatter-data/DarkMatter/default/nosuchfile.bed /var/www/OncoCis/Summary/"+request[1]);
////                    
////                }
////                else{
//                init(request);
////            }
//                execute("UPDATE dark_matter.request_system SET flag='1' WHERE id=" +request[0]+";");
//                System.out.println("Execution of request ID "+request[1]+ " complete. Flag set to 1");
//                
//                
//            }
//            else{
//                try{
//                    Thread.sleep(10000);
//                }
//                catch(Exception e){
//                    Thread.currentThread().interrupt();
//                }
//            }
//            if(request[1]==null){
//                System.out.println("Waiting for new requests");
//            }
//           
//        }
//        
    }
    
    public static void init(String request[],String links[][]){
//        System.out.println("cleardb");
        dbConnect.connectDriver();
        
        clearDB();

        OncoCis oc = new OncoCis();
        System.out.println("true");
        oc.OUTPUT_FOLDER = createOPFolder(request[1]); 
        oc.DNASE1HS = links[2][1];
        oc.H3K4ME3 = links[3][1];
        oc.H3K4ME1 = links[4][1];
        oc.H3K27AC = links[5][1];
        
//        oc.requestID = request[1];
//        oc.OUTPUT_FOLDER = createOPFolder(request[1]); 
        UtilitiesLinux.shellCommandExecuter("mv " + MAIN_FOLDER + "Input/" + request[2]
        +" "+oc.OUTPUT_FOLDER+"MutationsInput.bed");
        oc.MUTATIONS = oc.OUTPUT_FOLDER+"MutationsInput.bed";
        oc.CELLTYPE = request[3];
        if(request[4].equals("1")){
           oc.HASGENEEXP = true;
           UtilitiesLinux.shellCommandExecuter("mv " + MAIN_FOLDER + "Input/" + request[5]
        +" "+oc.OUTPUT_FOLDER+"GeneExpUserInput.bed");
           oc.GENEEXP = oc.OUTPUT_FOLDER+"GeneExpUserInput.bed"; 
        }
        
        
        
        
      
//        String data[] = getdata(oc.CELLTYPE);
//        System.out.println(data.length);
//        oc.DNASE1HS = MAIN_FOLDER + "Cell_type_specific_data/"+data[1];
//        oc.H3K4ME1 = MAIN_FOLDER + "Cell_type_specific_data/"+data[2];
//        oc.H3K4ME3 = MAIN_FOLDER + "Cell_type_specific_data/"+data[3];
//        oc.H3K27AC = MAIN_FOLDER + "Cell_type_specific_data/"+data[4];
//        System.out.println("Mutations "+dm.MUTATIONS);
//        System.out.println("Cell Type "+dm.CELLTYPE);
//        System.out.println("Gene Expression data "+dm.GENEEXP);
//        System.out.println("Output Folder "+dm.OUTPUT_FOLDER);
//        System.out.println("DNase1 "+dm.DNASE1HS);
//        System.out.println("H3K4me1 "+dm.H3K4ME1);
//        System.out.println("H3K4me3 "+dm.H3K4ME3);
//        System.out.println("H3K27ac "+dm.H3K27AC);    
        UtilitiesLinux.fromDosCheck(oc.MUTATIONS);
        
        Preprocessor preP =null;
        String error = "No Errors";
        if(oc.HASGENEEXP){UtilitiesLinux.fromDosCheck(oc.GENEEXP);
            preP = new Preprocessor(oc.MUTATIONS, oc.GENEEXP, oc.OUTPUT_FOLDER);}
        else{
            preP = new Preprocessor(oc.MUTATIONS, oc.OUTPUT_FOLDER); 
        }
        
        if(!preP.getMutationFileError().equals("No Errors")){
            error = preP.getMutationFileError();
            
        }
        else if(oc.HASGENEEXP&&(!preP.getGeneExpFileError().equals("No Errors"))){
            error = preP.getGeneExpFileError();
        }
        
        if(error.equals("No Errors")){
            System.out.println("True");
        
            Mutations mutations = new Mutations(oc.MUTATIONS, oc.OUTPUT_FOLDER);

            DNase1 dnase1 = new DNase1(oc.DNASE1HS, oc.OUTPUT_FOLDER, mutations);

            Histones histones = new Histones(oc.H3K4ME1, oc.H3K4ME3, oc.H3K27AC, mutations);

            GeneMapping genemap = new GeneMapping(oc.GENEMAP, oc.OUTPUT_FOLDER, mutations);                     

//            FimoMotif motifs = new FimoMotif(dm.HG19fa, dm.MOTIFDB, dm.OUTPUT_FOLDER, mutations);
            PossumMotif motifs = new PossumMotif(oc.HG19fa, oc.OUTPUT_FOLDER,mutations);

            Conservation  cons = new Conservation(oc.CONSERVATION, oc.OUTPUT_FOLDER, mutations);

            Fantom5 f5 = new Fantom5(oc.OUTPUT_FOLDER, mutations);
            
            updateFinalResTable();
            if(oc.HASGENEEXP){
                GeneExpression ge = new GeneExpression(oc.GENEEXP,oc.OUTPUT_FOLDER);
            }
            getFinalResults("FinalOutput/MutationAnnotations-"+System.currentTimeMillis(), preP.getTotalNumofLines(), oc.CELLTYPE,oc.HASGENEEXP);
        }
        else{
            String errorSummary[][] = new String[1][10];
            for (int i = 0; i < errorSummary[0].length; i++) {
                if(i==9){errorSummary[0][i]=error;}
                else {errorSummary[0][i]="0";}
            }
            DataSetWriter.writeToFile2D(MAIN_FOLDER + "FinalOutput/"+"Errors/"+System.currentTimeMillis(),errorSummary);
        }
    
                
    }
    
    public static void updateFinalResTable(){
        dbConnect.execute("INSERT INTO dark_matter_temp.Final_results ( \n" +
                            "Chromosome, Start, Stop, SampleID, MutationID, \n" +
                            "Gene, Druggable, Distance_to_TSS, Fold_Change, pValue, DHS, H3K4me1, H3K4me3, H3K27ac,\n" +
                            "Conservation_of_mutation, Background_conservation, " +
                            "Motifs_Created, Motifs_Removed, Fantom_Promoter,Fantom_Enhancer) \n" +
                            "SELECT Chromosome, Start, Stop, SampleID, MutationID," +
                            "'NONE', '0','1000000', '9999', '9999', '0', '0','0','0'," +
                            "'0','0', '-','-','0','0' FROM dark_matter_temp.mutations;");
        
        dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.mutations_to_dhs\n" +
                        "ON (Final_results.MutationID=mutations_to_dhs.MutationID)\n" +
                        "SET Final_results.DHS = mutations_to_dhs.DHS;");
        
        dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.mutations_to_histones\n" +
                        "ON (Final_results.MutationID=mutations_to_histones.MutationID)\n" +
                        "SET Final_results.H3K4me1 = mutations_to_histones.H3K4me1,\n" +
                        "Final_results.H3K4me3 = mutations_to_histones.H3K4me3,\n" +
                        "Final_results.H3K27ac = mutations_to_histones.H3K27ac;");
        
        dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.mutations_to_genes\n" +
                        "ON (Final_results.MutationID=dark_matter_temp.mutations_to_genes.MutationID)\n" +
                        "SET Final_results.Gene= dark_matter_temp.mutations_to_genes.Gene,\n" +
                        "Final_results.Distance_to_TSS = dark_matter_temp.mutations_to_genes.DistanceToTSS;");
        
        dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.mutations_to_conservation\n" +
                        "ON (Final_results.MutationID=dark_matter_temp.mutations_to_conservation.MutationID)\n" +
                        "SET Final_results.conservation_of_mutation= dark_matter_temp.mutations_to_conservation.Score;");
        
        dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.mutations_to_backgroundConservation\n" +
                        "ON (Final_results.MutationID=dark_matter_temp.mutations_to_backgroundConservation.MutationID)\n" +
                        "SET Final_results.background_conservation= dark_matter_temp.mutations_to_backgroundConservation.Score;");
        
        dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.mutations_to_motifs_created\n" +
                        "ON (Final_results.MutationID=dark_matter_temp.mutations_to_motifs_created.MutationID)\n" +
                        "SET Final_results.Motifs_Created= dark_matter_temp.mutations_to_motifs_created.Motifs;");
        
        dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.mutations_to_motifs_removed\n" +
                        "ON (Final_results.MutationID=dark_matter_temp.mutations_to_motifs_removed.MutationID)\n" +
                        "SET Final_results.Motifs_Removed= dark_matter_temp.mutations_to_motifs_removed.Motifs;");
        
        dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.mutations_to_Fantom\n" +
                        "ON (Final_results.MutationID=mutations_to_Fantom.MutationID)\n" +
                        "SET Final_results.Fantom_Promoter = mutations_to_Fantom.Promoter,\n" +
                        "Final_results.Fantom_Enhancer = mutations_to_Fantom.Enhancer;");
        
        dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.gene_hg19_1\n" +
                        "ON (Final_results.Gene=gene_hg19_1.GeneSymbol)\n" +
                        "SET Final_results.Druggable = '1' where gene_hg19_1.link_drug!='';");
    }
    
    public static void getFinalResults(String outputFile, int total, String cellType, boolean hasGeneExpData){
//        if(hasGeneExpData){
//           dbConnect.updateDB("ALTER TABLE `dark_matter_temp`.`Final_results` "
//                   + "CHANGE COLUMN `Fold_Change` `Fold_Change` DOUBLE NULL DEFAULT '0'  , "
//                   + "CHANGE COLUMN `pValue` `pValue` DOUBLE NULL DEFAULT '0'  ;");
//           
//           dbConnect.execute("INSERT INTO dark_matter_temp.Final_results ( \n" +
//                        "Chromosome, Start, Stop, SampleID, MutationID, \n" +
//                        "Gene, Distance_to_TSS, Fold_Change, pValue, DHS, H3K4me1, H3K4me3, H3K27ac,\n" +
//                        "Conservation_of_mutation, Background_conservation, \n" +
//                        "Motifs_Created, Motifs_Removed\n" +
//                        ") \n" +
//                        "SELECT Chromosome, Start, Stop, SampleID, MutationID,\n" +
//                        "'NONE', '1000000', 0, 0, '0', '0','0','0',\n" +
//                        "'0','0', '-','-'\n" +
//                        "FROM dark_matter_temp.mutations");
//        }
//        
//        else{
//            dbConnect.updateDB("ALTER TABLE `dark_matter_temp`.`Final_results` "
//                    + "CHANGE COLUMN `Fold_Change` `Fold_Change` VARCHAR(10) NULL DEFAULT 'NaN'  , "
//                    + "CHANGE COLUMN `pValue` `pValue` VARCHAR(10) NULL DEFAULT 'NaN'  ;");
//            
//            dbConnect.execute("INSERT INTO dark_matter_temp.Final_results ( \n" +
//                            "Chromosome, Start, Stop, SampleID, MutationID, \n" +
//                            "Gene, Distance_to_TSS, Fold_Change, pValue, DHS, H3K4me1, H3K4me3, H3K27ac,\n" +
//                            "Conservation_of_mutation, Background_conservation, \n" +
//                            "Motifs_Created, Motifs_Removed\n" +
//                            ") \n" +
//                            "SELECT Chromosome, Start, Stop, SampleID, MutationID,\n" +
//                            "'NONE', '1000000', 'NaN', 'NaN', '0', '0','0','0',\n" +
//                            "'0','0', '-','-'\n" +
//                            "FROM dark_matter_temp.mutations");
////        }
        
        
        if(hasGeneExpData){
            dbConnect.execute("UPDATE dark_matter_temp.Final_results " +
                        "SET \n" +
                        "Final_results.Fold_Change = '-1',\n" +
                        "Final_results.pValue = '-1';");
            
            dbConnect.execute("UPDATE dark_matter_temp.Final_results\n" +
                        "inner join dark_matter_temp.mutations_to_geneExp\n" +
                        "ON (Final_results.MutationID=mutations_to_geneExp.MutationID)\n" +
                        "SET \n" +
                        "Final_results.Fold_Change = mutations_to_geneExp.Foldchange,\n" +
                        "Final_results.pValue = mutations_to_geneExp.pValue;");
            
            dbConnect.execute("UPDATE dark_matter_temp.Final_results " +
                        "SET \n" +
                        "Final_results.Fold_Change = '-1',\n" +
                        "Final_results.pValue = '-1'"
                    + "where Final_results.Fold_Change = '0';");
        }      
        
        
        writeOutput("dark_matter_temp.Final_results", outputFile, MAIN_FOLDER+"FinalOutput/",hasGeneExpData);
        
        dbConnect.execute("Insert into dark_matter_temp.Summary\n" +
                        "(Total, Non_coding, DNase1, H3K4me1, H3K4me3, H3K27ac, Mapped_to_Gene, Possible_Promoter, Highly_cons, Cell_Type)\n" +
                        "Select "+total+", count(Final_results.MutationID), sum(Final_results.DHS), \n" +
                        "sum(Final_results.H3K4me1), sum(Final_results.H3K4me3), sum(Final_results.H3K27ac),\n" +
                        "sum(case when Gene='NONE' then 0 else 1 end),\n" +
                        "sum(case when Distance_to_TSS>-10000 and Distance_to_TSS<1000  then 1 else 0 end),\n" +
                        "sum(case when Conservation_of_mutation>0.8 then 1 else 0 end), \n" + "'"+cellType + "'"+
                        " from dark_matter_temp.Final_results;");
        
        dbConnect.exportTableToFile("dark_matter_temp.Summary", outputFile, MAIN_FOLDER+"FinalOutput/Summary/");
//        dbConnect.exportDataToFile("SELECT Chromosome, Start, Start, concat('Mutation-',MutationID) FROM dark_matter_temp.mutations_reduced", outputFile + "old", "/var/www/OncoCis/BED/");
//        UtilitiesLinux.shellCommandExecuter("sed '1i\\track name=coords description=\"Chromosome coordinates list\" visibility=2' " + "/var/www/OncoCis/BED/"+outputFile + "old > " + "/var/www/OncoCis/BED/"+outputFile);
        //sed '1i\ track name=coords description="Chromosome coordinates list" visibility=2' file_name > new_filename

        
}
    
    public static void writeOutput(String tableName, String fileName, String outputLocation, boolean geneexp){
        UtilitiesLinux.shellCommandExecuter("mkdir /tmp/mysql");
        UtilitiesLinux.shellCommandExecuter("chmod 777 /tmp/mysql");
        String sqlCommand = "";
        if(geneexp){
            sqlCommand = "SELECT Chromosome, Start, Stop, SampleID,"
                + " Gene, Druggable, Distance_to_TSS, FORMAT(Fold_Change,3), pValue, DHS,"
                + " H3K4me1, H3K4me3, H3K27ac,FORMAT(Conservation_of_mutation,3),"
                + " FORMAT(Background_conservation,3), Motifs_Created, Motifs_Removed, Fantom_Promoter, Fantom_Enhancer"
                + " FROM " + tableName + " INTO OUTFILE '" + "/tmp/mysql/"+fileName + "';";
                }
        else{
            sqlCommand = "SELECT Chromosome, Start, Stop, SampleID,"
                + " Gene, Druggable, Distance_to_TSS, 9999, 9999, DHS,"
                + " H3K4me1, H3K4me3, H3K27ac,FORMAT(Conservation_of_mutation,3),"
                + " FORMAT(Background_conservation,3), Motifs_Created, Motifs_Removed, Fantom_Promoter, Fantom_Enhancer"
                + " FROM " + tableName + " INTO OUTFILE '" + "/tmp/mysql/"+fileName + "';";
        
        }
        dbConnect.connectDriver();
        dbConnect.queryDB(sqlCommand);
//        UtilitiesLinux.outputFolder = outputLocation;
        UtilitiesLinux.shellCommandExecuter("sed -i '1iChrom\\tStart\\tEnd\\tSampleID"
                + "\\tGene\\tDruggable\\tDistanceToTSS"
                + "\\tGeneExpression_FoldChange\\tGeneExpression_p-value"
                + "\\tDHS\\tH3K4me1\\tH3K4me3\\tH3K27ac"
                + "\\tConservation_at_mutation\\tBackgroundConservation"
                + "\\tMotifs_Created\\tMotifs_removed"
                + "\\tFantom5_Promoter\\tFantom5_Enhancer'"+" /tmp/mysql/"+fileName);
        UtilitiesLinux.shellCommandExecuter("mv /tmp/mysql/"+fileName+ " "+ outputLocation + "." + "\n" + "rm /tmp/mysql/"+fileName);
        
        
    }
    
    public static void execute(String sqlCommand){       
       Connection conn = null;
       Statement stmt = null;
        try { 
            conn =
               DriverManager.getConnection("jdbc:mysql://localhost:3306/dark_matter" ,"root","");
            stmt = conn.createStatement();
//            stmt.executeUpdate("SELECT * FROM mutations_mapped_to_dhs INTO OUTFILE '/tmp/MutationsMappedToDHS4.bed';");
            stmt.execute(sqlCommand);
           
           


        } catch (SQLException ex) {
            // handle any errors
            System.out.println("SQLException: " + ex.getMessage());
            System.out.println("SQLState: " + ex.getSQLState());
            System.out.println("VendorError: " + ex.getErrorCode());
        }
       try{
                
                stmt.close();
                conn.close();
            }
            catch(SQLException ex){
                System.out.println("SQLException: " + ex.getMessage());
            System.out.println("SQLState: " + ex.getSQLState());
            System.out.println("VendorError: " + ex.getErrorCode());
            }
    }
    
    public static String[] getRequest(){
        String request[]=new String[9];        
        String query = "SELECT * FROM dark_matter.request_system;";
        Connection conn = null;
        Statement stmt = null;
        ResultSet rs= null;
        try {
                // jdbc:mysql://localhost:3306/sample 
                conn = DriverManager.getConnection("jdbc:mysql://localhost:3306/dark_matter" ,"root","");
                stmt = conn.createStatement();
                rs = stmt.executeQuery(query);            
                if(rs.next()){                    
                    request[0] = rs.getString(1);
                    request[1] = rs.getString(2);
                    request[2] = rs.getString(3);
                    request[3] = rs.getString(4);
                    request[4] = rs.getString(5);
                    request[5] = rs.getString(6);
                    request[6] = rs.getString(7);
                    request[7] = rs.getString(8);
                    request[8] = rs.getString(9);
                }
            
        } 
        catch (SQLException ex) {
            // handle any errors
            System.out.println("SQLException: " + ex.getMessage());
            System.out.println("SQLState: " + ex.getSQLState());
            System.out.println("VendorError: " + ex.getErrorCode());
        }
        
            try{
                rs.close();
                stmt.close();
                conn.close();
            }
            catch(SQLException ex){
                System.out.println("SQLException: " + ex.getMessage());
            System.out.println("SQLState: " + ex.getSQLState());
            System.out.println("VendorError: " + ex.getErrorCode());
            }
        
        return request;
    }
    
    public static String createOPFolder(String timestamp){        
        String folder = MAIN_FOLDER + "Output/"+timestamp;
        UtilitiesLinux.shellCommandExecuter("mkdir "+folder);
        return folder + "/";
    }

    public static String[] getdata(String cellType){
           String request[]=new String[7];
           dbConnect.connectDriver();
           String query = "SELECT * FROM dark_matter_temp.file_locations where Cell_type = '"+ cellType + "';";
           Connection conn = null;
           Statement stmt = null;
           ResultSet rs= null;
           try {
               // jdbc:mysql://localhost:3306/sample 
               conn =
                  DriverManager.getConnection("jdbc:mysql://localhost:3306/dark_matter_temp" ,"root","");
               stmt = conn.createStatement();
               rs = stmt.executeQuery(query);

               if(rs.next()){
               request[0] = rs.getString(2);
               request[1] = rs.getString(3);
               request[2] = rs.getString(4);
               request[3] = rs.getString(5);
               request[4] = rs.getString(6);
               request[5] = rs.getString(7);
   //            System.out.println(""+request[0]);
               }

           } 
           catch (SQLException ex) {
               // handle any errors
               System.out.println("SQLException: " + ex.getMessage());
               System.out.println("SQLState: " + ex.getSQLState());
               System.out.println("VendorError: " + ex.getErrorCode());
           }
           try{
                rs.close();
                stmt.close();
                conn.close();
            }
            catch(SQLException ex){
                System.out.println("SQLException: " + ex.getMessage());
            System.out.println("SQLState: " + ex.getSQLState());
            System.out.println("VendorError: " + ex.getErrorCode());
            }
           return request;
       }
    
    public static void clearDB(){
        String request[]=new String[50];
           dbConnect.connectDriver();
           String query = "SHOW TABLES IN dark_matter_temp;";
           Connection conn = null;
           Statement stmt = null;
           ResultSet rs= null;
           int count=0;
           try {
               // jdbc:mysql://localhost:3306/sample 
               conn =
                  DriverManager.getConnection("jdbc:mysql://localhost:3306/dark_matter_temp" ,"root","");
               stmt = conn.createStatement();
               rs = stmt.executeQuery(query);
               

               while(rs.next()){
                   request[count] = rs.getString(1);
                   count++;
               }

           } 
           catch (SQLException ex) {
               // handle any errors
               System.out.println("SQLException: " + ex.getMessage());
               System.out.println("SQLState: " + ex.getSQLState());
               System.out.println("VendorError: " + ex.getErrorCode());
           }
           try{
                rs.close();
                stmt.close();
                conn.close();
            }
            catch(SQLException ex){
                System.out.println("SQLException: " + ex.getMessage());
            System.out.println("SQLState: " + ex.getSQLState());
            System.out.println("VendorError: " + ex.getErrorCode());
            }
           
           for (int i = 0; i < count; i++) {
               
               if(!request[i].equals("file_locations")){
//                   System.out.println("Cleared "+request[i]);
                   dbConnect.clearTable(request[i]);
               }
                
            }
           System.out.println("Number of tables cleared " + count);
      
        
    }   

     public String getOUTPUT_FOLDER() {
        return OUTPUT_FOLDER;
    }

    public void setOUTPUT_FOLDER(String OUTPUT_FOLDER) {
        this.OUTPUT_FOLDER = OUTPUT_FOLDER;
    }

    public String getMUTATIONS() {
        return MUTATIONS;
    }

    public void setMUTATIONS(String MUTATIONS) {
        this.MUTATIONS = MUTATIONS;
    }

    public String getDNASE1HS() {
        return DNASE1HS;
    }

    public void setDNASE1HS(String DNASE1HS) {
        this.DNASE1HS = DNASE1HS;
    }

    public String getH3K4ME1() {
        return H3K4ME1;
    }

    public void setH3K4ME1(String H3K4ME1) {
        this.H3K4ME1 = H3K4ME1;
    }

    public String getH3K4ME3() {
        return H3K4ME3;
    }

    public void setH3K4ME3(String H3K4ME3) {
        this.H3K4ME3 = H3K4ME3;
    }

    public String getH3K27AC() {
        return H3K27AC;
    }

    public void setH3K27AC(String H3K27AC) {
        this.H3K27AC = H3K27AC;
    }

    
    public String getCONSERVATION() {
        return CONSERVATION;
    }

    public void setCONSERVATION(String CONSERVATION) {
        this.CONSERVATION = CONSERVATION;
    }

       public String getHG19fa() {
        return HG19fa;
    }

    public void setHG19fa(String HG19fa) {
        this.HG19fa = HG19fa;
    }

    public String getGENEMAP() {
        return GENEMAP;
    }

    public void setGENEMAP(String GENEMAP) {
        this.GENEMAP = GENEMAP;
    }

    public String getGENEEXP() {
        return GENEEXP;
    }

    public void setGENEEXP(String GENEEXP) {
        this.GENEEXP = GENEEXP;
    }

}

