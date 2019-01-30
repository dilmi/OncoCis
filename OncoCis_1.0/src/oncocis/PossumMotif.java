/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package OncoCis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

/**
 *
 * @author Dilmi
 */
public class PossumMotif {
    public static String hg19fa = "";
    public static String motifDB = "";
    public static String normalSequence = "";
    public static String substitutionFile = "";
    public static String mutatedSequence = "";
    public static String normalMotifs = "";
    public static String mutatedMotifs = "";
    public static String possumOutputFolder1 = "";
    public static String possumOutputFolder2 = "";
    
    Mutations mutations;
//    SELECT * FROM dark_matter_temp.mutated_motifs
//WHERE (Motif, Location, Start, Stop, Strand, Score, pValue, qValue, Matched_sequence) NOT IN
//( SELECT *
//  FROM dark_matter_temp.normal_motifs
//) ;
//Insert INTO dark_matter_temp.motifs_temp
//SELECT 
//Motif, Location, Score, pValue FROM dark_matter_temp.mutated_motifs
//WHERE (Motif, Location, Start, Stop, Strand, Score, pValue, qValue, Matched_sequence) NOT IN
//( SELECT *
//  FROM dark_matter_temp.normal_motifs
//) ;
    public PossumMotif(String cHg19fa,  String cOutputFolder, Mutations m) {
        hg19fa = cHg19fa;
        normalSequence = cOutputFolder + "NormalSequence.fasta";
        substitutionFile = cOutputFolder + "substitutions.bed";
        mutatedSequence = cOutputFolder + "MutatedSequence.fasta";
        normalMotifs = cOutputFolder + "normalMotifs.bed";
        mutatedMotifs = cOutputFolder + "mutatedMotifs.bed";
        mutations = m;
        UtilitiesLinux.shellCommandExecuter("mkdir " + cOutputFolder + "/Motifs1");
        possumOutputFolder1 = cOutputFolder + "Motifs1";
        UtilitiesLinux.shellCommandExecuter("mkdir " + cOutputFolder + "/Motifs2");
        possumOutputFolder2 = cOutputFolder + "Motifs2";
        System.out.println("Get sequence");
        this.getSequence();
        System.out.println("Get mutated sequence");
        this.getMutatedSequence();
        System.out.println("get motifs");
        getMotifs();
        System.out.println("get altered motifs");
        getAlteredMotifsFinal(normalMotifs, mutatedMotifs, cOutputFolder + "motifsRemoved.bed", cOutputFolder + "motifsCreated.bed");
        updateMotifsInDB();
    }  
    
    public PossumMotif(String cHg19fa,  String cOutputFolder) {
        hg19fa = cHg19fa;
        normalSequence = cOutputFolder + "NormalSequence.fasta";
        substitutionFile = cOutputFolder + "substitutions.bed";
        mutatedSequence = cOutputFolder + "MutatedSequence.fasta";
        normalMotifs = cOutputFolder + "normalMotifs.bed";
        mutatedMotifs = cOutputFolder + "mutatedMotifs.bed";
        
        UtilitiesLinux.shellCommandExecuter("mkdir " + cOutputFolder + "/Motifs1");
        possumOutputFolder1 = cOutputFolder + "Motifs1";
        UtilitiesLinux.shellCommandExecuter("mkdir " + cOutputFolder + "/Motifs2");
        possumOutputFolder2 = cOutputFolder + "Motifs2";
        System.out.println("Get sequence");
        this.getSequence();
        System.out.println("Get mutated sequence");
        this.getMutatedSequence();
        System.out.println("get motifs");
        getMotifs();
        System.out.println("get altered motifs");
        getAlteredMotifsFinal(normalMotifs, mutatedMotifs, cOutputFolder + "motifsRemoved.bed", cOutputFolder + "motifsCreated.bed");
//        getMotifCreatedOrRemoved();
//        getMotifRemoved();
    }

    public static void main(String[] args) {
//        PossumMotif.getMotifCreatedOrRemoved();
//        testMutatedSequence();
        OncoCis.clearDB();
        Preprocessor p=new Preprocessor("/home/bioinfo/Darkmatter-data/DarkMatter/Output/95IZ0/Mutations-input.bed", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/95IZ0/");
        Mutations mut=new Mutations("/home/bioinfo/Darkmatter-data/DarkMatter/Output/95IZ0/Mutations-input.bed", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/95IZ0/");
        PossumMotif pm=new PossumMotif("/home/bioinfo/Darkmatter-data/DarkMatter/default/hg19.fa", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/95IZ0/",mut);
//////        
        
         
//        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==7 {if($2!='null'){ print $1,$2;}}' "+"/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/Motifs1/Possum-edited.bed > /home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/Motifs1/normalMotifs.bed");
//         
//        dbConnect.importDataFromFile("dark_matter_temp.possum_normal", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/Motifs1/normalMotifs-copy.bed");
//        
//        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==7 {if($2!='null'){ print $1,$2;}}' "+"/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/Motifs2/Possum-edited.bed > /home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/Motifs2/mutatedMotifs.bed");
//         dbConnect.importDataFromFile("dark_matter_temp.possum_mutated", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/Motifs2/mutatedMotifs-copy.bed");
         
//           dbConnect.importDataFromFile("dark_matter_temp.motifs_removed", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/MotifsRemoved.bed");
//           dbConnect.importDataFromFile("dark_matter_temp.motifs_created", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/MotifsCreated.bed");

//           dbConnect.execute(  "Insert into dark_matter_temp.mutations_to_motifs_created\n" +
//                            "(select Chromosome, mutations_extended.Start, mutations_extended.Stop, MutationID, count(Motif) as NumOfMotifs, group_concat(Distinct Motif)as Motifs  \n" +
//                            "from dark_matter_temp.motifs_created join dark_matter_temp.mutations_extended \n" +
//                            "Where Location=concat(Chromosome,'-',mutations_extended.Start,'-',mutations_extended.Stop)\n" +
//                            "Group by\n" +
//                            "MutationID);");
//        
//        dbConnect.execute(  "Insert into dark_matter_temp.mutations_to_motifs_removed\n" +
//                            "(select Chromosome, mutations_extended.Start, mutations_extended.Stop, MutationID, count(Motif) as NumOfMotifs, group_concat(Distinct Motif)as Motifs  \n" +
//                            "from dark_matter_temp.motifs_removed join dark_matter_temp.mutations_extended \n" +
//                            "Where Location=concat(Chromosome,'-',mutations_extended.Start,'-',mutations_extended.Stop)\n" +
//                            "Group by\n" +
//                            "MutationID);");
////        //dbConnect.connectDriver();
//////        dbConnect.execute(  "Select \n" +
//                            "Motif, Location, Score, pValue FROM dark_matter_temp.mutated_motifs\n" +
//                            "WHERE (Motif, Location, Start, Stop, Strand, Score, pValue, qValue, Matched_sequence) NOT IN\n" +
//                            "( SELECT *\n" +
//                            "  FROM dark_matter_temp.normal_motifs\n" +
//                            ")\n" +
//                            "into outfile 'motifs_temp.bed';"); 
//        UtilitiesLinux.shellCommandExecuter("mv /tmp/"+"motifs_temp.bed"+ "/home/bioinfo/Darkmatter-data/DarkMatter/Output/" + ".");
       
//        String jasperMotifDB[] = DataSetReader.readLines("/home/bioinfo/Darkmatter-data/DarkMatter/default/JASPAR_CORE_2009.meme");
//       
//       String temp[];
//        for (int i = 0; i < jasperMotifDB.length; i++) {
////            System.out.println("true");
//           if (jasperMotifDB[i].contains("MOTIF")) {
//                temp = jasperMotifDB[i].trim().split(" ");
//                System.out.println(""+temp[0]+" " + temp[1]+" " +temp[2]);
//                jasperMotifDB[i] = temp[0]+" "+temp[2];
//                System.out.println(""+jasperMotifDB[i]);
//                
//                
//            } 
//        }
        
//       String jasperMotifDB[] = DataSetReader.readLines("/home/bioinfo/Darkmatter-data/DarkMatter/default/JASPAR_CORE_2009_EDITED.meme");
//        ArrayList<Integer> lineNum = new ArrayList<Integer>();
//       String temp[];
//       int i=0;
//       while (i < jasperMotifDB.length) {
//           if (jasperMotifDB[i].contains("MOTIF")) {
//                temp = jasperMotifDB[i].trim().split(" ");
//                if(Character.isLowerCase(temp[1].charAt(0))||temp[1].length()==1){
//                    System.out.println(""+jasperMotifDB[i]);
//                    if(i+1<jasperMotifDB.length){
//                        i++;
//                        
//                        while(!jasperMotifDB[i].contains("MOTIF")){
//                            if(i+1<jasperMotifDB.length)i++;
//                            else{break;}
//                        }
//                    }
//                    else{
//                        break;
//                    }
//                }
//                else{
//                    lineNum.add(i); 
//                    i++;
//                }
//            }
//           else{
//               lineNum.add(i);
//               i++;
//           }
//        }
//        String jasperMotifDBnew[] = new String [lineNum.size()];
//        for (int j = 0; j < jasperMotifDBnew.length; j++) {
//            jasperMotifDBnew[j] = jasperMotifDB[lineNum.get(j)];
//        }
//        
//       DataSetWriter.writeToFile1D("/home/bioinfo/Darkmatter-data/DarkMatter/default/JASPAR_CORE_2009_EDITED_new1.meme", jasperMotifDBnew);
//       
    }
    
    public void getSequence(){         
        UtilitiesLinux.bedtoolsGetFasta(hg19fa, mutations.getMutations20bpExtented() , normalSequence);
//    UtilitiesLinux.bedtoolsGetFasta(hg19fa, "/home/bioinfo/Darkmatter-data/DarkMatter/attachments(1)/allmutations40kb.bed" , normalSequence);
    }
    
    public static void test(String[] args) {
        String seq = "ttgcaaaagctcaatgaatctcaaccaggattttctgttca";
        System.out.println(""+seq.substring(20,21));
        System.out.println(""+seq.substring(0,20));
        System.out.println(""+seq.substring(seq.length()-20));
        
    }
    
    public void getMutatedSequence(){
        String sequence[][] = DataSetReader.readDataSet(normalSequence, 1, "\t");
        
        //NOTE : Mutation.mutations has a header!!!!!!!!!!!
        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==8 { print $7, $8, $4, $6;}' " + mutations.getMutations() + " > "+substitutionFile);
//        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==8 { print $7, $8, $4, $6;}' " + "/home/bioinfo/Darkmatter-data/DarkMatter/attachments(1)/BreastCancerDataset-OncoCisinput.bed" + " > "+substitutionFile);
        
        String substitutions[][] = DataSetReader.readDataSet(substitutionFile,4,"\t");              
        String temp="";
        String temp1="";
        int seqLenth = 0;
        if((sequence.length/2)!=substitutions.length){
            System.out.println(""+sequence.length);
            System.out.println(""+substitutions.length);
                System.out.println("Error creating mutated sequence - substitution file and sequence file don't match"); 
            }
        else{
        for (int i = 0; i < sequence.length; i++) {            
                seqLenth = sequence[i][0].length();
                if(i%2!=0){
                    if(substitutions[(i-1)/2][3].equalsIgnoreCase("substitution")){
                        //if(substitutions[(i-1)/2][0].equals(sequence[i][0].substring(20,21))){
                        if(true){
                            temp=sequence[i][0].substring(0,20);
                            temp1=sequence[i][0].substring(seqLenth-20);                    
                            sequence[i][0]=temp+substitutions[(i-1)/2][1]+temp1;
                        }
                        else if(substitutions[(i-1)/2][0].equalsIgnoreCase(sequence[i][0].substring(20,21))){
                            temp=sequence[i][0].substring(0,20);
                            temp1=sequence[i][0].substring(seqLenth-20);                    
                            sequence[i][0]=temp+substitutions[(i-1)/2][1]+temp1;
                        }                       
                        else{
                            sequence[i][0]="ERROR";
                            System.out.println("Error"+i);
                        }

                    }
                    else if(substitutions[(i-1)/2][3].equalsIgnoreCase("deletion")){                        
                            temp=sequence[i][0].substring(0,20);
                            temp1=sequence[i][0].substring(seqLenth-20);                    
                            sequence[i][0]=temp+temp1;                       
                                              
                        
                                                                
                    }
                    else if(substitutions[(i-1)/2][3].equalsIgnoreCase("insertion")){                        
                            temp=sequence[i][0].substring(0,20);
                            temp1=sequence[i][0].substring(seqLenth-20);                    
                            sequence[i][0]=temp+substitutions[(i-1)/2][1]+temp1;
                        
                                                                
                    }
                }
            }
            
        }
        DataSetWriter.writeToFile2D(mutatedSequence, sequence);
        
    }
    
    public static void testMutatedSequence(){
        String sequence[][] = DataSetReader.readDataSet("/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/NormalSequence.fasta", 1, "\t");
        
        //NOTE : Mutation.mutations has a header!!!!!!!!!!!
        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==7 { print $6, $7, $3, $5;}' " + "/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/BreastCancerDataset-OncoCisinput.bed" + " > "+"/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/substitutionFile.bed");
//        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==8 { print $7, $8, $4, $6;}' " + "/home/bioinfo/Darkmatter-data/DarkMatter/attachments(1)/BreastCancerDataset-OncoCisinput.bed" + " > "+substitutionFile);
        
        String substitutions[][] = DataSetReader.readDataSet("/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/substitutionFile.bed",4,"\t");              
        String temp="";
        String temp1="";
        int seqLenth = 0;
        if((sequence.length/2)!=substitutions.length){
            System.out.println(""+sequence.length);
            System.out.println(""+substitutions.length);
                System.out.println("Error creating mutated sequence - substitution file and sequence file don't match"); 
            }
        else{
        for (int i = 0; i < sequence.length; i++) {            
                seqLenth = sequence[i][0].length();
                if(i%2!=0){
                    if(substitutions[(i-1)/2][3].equalsIgnoreCase("substitution")){
                        //if(substitutions[(i-1)/2][0].equals(sequence[i][0].substring(20,21))){
                        if(true){
                            temp=sequence[i][0].substring(0,20);
                            temp1=sequence[i][0].substring(seqLenth-20);                    
                            sequence[i][0]=temp+substitutions[(i-1)/2][1]+temp1;
                        }
                        else if(substitutions[(i-1)/2][0].equalsIgnoreCase(sequence[i][0].substring(20,21))){
                            temp=sequence[i][0].substring(0,20);
                            temp1=sequence[i][0].substring(seqLenth-20);                    
                            sequence[i][0]=temp+substitutions[(i-1)/2][1]+temp1;
                        }                       
                        else{
                            sequence[i][0]="ERROR";
                            System.out.println("Error"+i);
                        }

                    }
                    else if(substitutions[(i-1)/2][3].equalsIgnoreCase("deletion")){                        
                            temp=sequence[i][0].substring(0,20);
                            temp1=sequence[i][0].substring(seqLenth-20);                    
                            sequence[i][0]=temp+temp1;                       
                                              
                        
                                                                
                    }
                    else if(substitutions[(i-1)/2][3].equalsIgnoreCase("insertion")){                        
                            temp=sequence[i][0].substring(0,20);
                            temp1=sequence[i][0].substring(seqLenth-20);                    
                            sequence[i][0]=temp+substitutions[(i-1)/2][1]+temp1;
                        
                                                                
                    }
                }
            }
            
        }
        DataSetWriter.writeToFile2D("/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/MutatedSequence.fasta", sequence);
        
    }
    
    public void getMotifs(){
        UtilitiesLinux.shellCommandExecuter("sed -i 's/:/-/g' " + normalSequence);
        UtilitiesLinux.shellCommandExecuter("sed -i 's/:/-/g' " + mutatedSequence);
        UtilitiesLinux.possum( normalSequence, possumOutputFolder1);
        normalMotifs = possumOutputFolder1+"/normalMotifs.bed";
//        normalMotifs = fimoOutputFolder1+"/normalMotifs.bed";
//        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==9 { print $1,$2,$5,$6,$7,$8,$9;}' "+fimoOutputFolder1+"/fimo.txt"+" > "+normalMotifs);
//        System.out.println("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==9 { print $1,$2,$5,$6,$7,$8,$9;}' "+fimoOutputFolder1+"/fimo.txt"+" > "+normalMotifs);
        editPossumOutput(possumOutputFolder1+"/Possum.bed", possumOutputFolder1+"/Possum-edited.bed");
        motifFilter(possumOutputFolder1+"/Possum-edited.bed", possumOutputFolder1+"/Possum-edited-filtered.bed");
        UtilitiesLinux.shellCommandExecuter("bedtools groupby -i  "+possumOutputFolder1+"/Possum-edited-filtered.bed"+" -g 2 -c 1 -o collapse > " + normalMotifs);
        
//        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==7 {if($2!='null'){ print $1,$2;}}' "+possumOutputFolder1+"/Possum-edited.bed > "+normalMotifs);
         
//        dbConnect.importDataFromFile("dark_matter_temp.possum_normal", normalMotifs);
        
        UtilitiesLinux.possum( mutatedSequence, possumOutputFolder2);
        mutatedMotifs = possumOutputFolder2+"/mutatedMotifs.bed";
//        mutatedMotifs = fimoOutputFolder2+"/mutatedMotifs.bed";
//         UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==9 { print $1,$2,$5,$6,$7,$8,$9;}' "+fimoOutputFolder2+"/fimo.txt"+" > "+mutatedMotifs);
//        System.out.println("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==9 { print $1,$2,$5,$6,$7,$8,$9;}' "+fimoOutputFolder2+"/fimo.txt"+" > "+mutatedMotifs);
         editPossumOutput(possumOutputFolder2+"/Possum.bed", possumOutputFolder2+"/Possum-edited.bed"); 
//         UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==7 {if($2!='null'){ print $1,$2;}}' "+possumOutputFolder2+"/Possum-edited.bed > "+mutatedMotifs);
         motifFilter(possumOutputFolder2+"/Possum-edited.bed", possumOutputFolder2+"/Possum-edited-filtered.bed");
        UtilitiesLinux.shellCommandExecuter("bedtools groupby -i  "+possumOutputFolder2+"/Possum-edited-filtered.bed"+" -g 2 -c 1 -o collapse > " + mutatedMotifs);
        
         
//         dbConnect.importDataFromFile("dark_matter_temp.possum_mutated", mutatedMotifs);
        
        
        
    }
    
    public static void editPossumOutput(String possumOutput, String editedOutput){
        String data[][] = DataSetReader.readDataSet(possumOutput, 6, "\t");
        int count=0;
        String output[][] = new String[data.length][3];
         for (int i = 0; i <data.length; i++) {
             
             if(data[i][0].startsWith(">")){
//                 System.out.println("Sequence = "+data[i][0]);
//                 String temp[] = data[i][0].substring(1).split("-");
                 String location = data[i][0].substring(1).trim();
                 for (int j = i+1; j < data.length-1; j++) {                     
                     if(data[j][0].contains(">")){
//                         System.out.println("EOM");
                         break;
                     }
                     else{
                         if(!(data[j][0].equals("null")||data[j][0].equals(""))){
                            
//                             output[count][0]=temp[0].trim();
//                             output[count][1]=temp[1].trim();
//                             output[count][2]=temp[2].trim();
//                             output[count][3]= data[j][0];
//                             output[count][4]=data[j+1][1];
//                             output[count][5]=data[j+1][2];
//                             output[count][6]=data[j+1][3];
//                             output[count][7]=data[j+1][4];
//                             output[count][8]=data[j+1][5];
                             
                             output[count][0]= data[j][0];
//                             output[count][1]= location + "&"+data[j+1][1]+ "&" +data[j+1][2];
                             output[count][1]= location;
                             output[count][2]=data[j+1][4];
//                             output[count][4]=data[j+1][3];
//                             output[count][5]=data[j+1][5];
//                             output[count][6]=data[j+1][4];
//                             System.out.print(""+output[count][0]);
//                             System.out.print("\t"+output[count][1]);
//                             System.out.print("\t"+output[count][2]);
//                             System.out.print("\t"+output[count][3]);
//                             System.out.print("\t"+output[count][4]);
//                             System.out.print("\t"+output[count][5]);
//                             System.out.print("\t"+output[count][6]);
//                             System.out.print("\t"+output[count][7]);
                             count++;
                         }
                     }
                     
                 }
                 
             }
             else{
                 continue;
             }
         }
        DataSetWriter.writeToFile2D(editedOutput,output);
//        UtilitiesLinux.shellCommandExecuter("grep '.'"+editedOutput+"pre"+" > "+editedOutput);
    } 
                
    public static void motifFilter(String possumEditedOutput, String filteredMotifs) {
        String data[][] = DataSetReader.readDataSet(possumEditedOutput, 4, "\t");
        String motifFilter[][] = DataSetReader.readDataSet("/home/bioinfo/Darkmatter-data/DarkMatter/default/MotifFilter-JASPAR2014-with-RevComplement.txt", 3, "\t");
//        String motifList[][] = new String[367][7];
        HashMap hm = new HashMap<String, Vector>();
        Vector sequence;
        int max=0;
        for (int i = 0; i < motifFilter.length; i++) {            
            if(hm.containsKey(motifFilter[i][0])) {
                sequence=(Vector)hm.get(motifFilter[i][0]);
                sequence.add(motifFilter[i][1]);
                sequence.add(motifFilter[i][2]);
                if(sequence.size()>max){max=sequence.size();}
            }
            else{
//                motifList[hm.size()][0]=motifFilter[i][0];
                sequence=new Vector();
                sequence.add(motifFilter[i][1]);
                sequence.add(motifFilter[i][2]);
                hm.put(motifFilter[i][0], sequence);
                if(sequence.size()>max){max=sequence.size();}
            }
        }
//        System.out.println(""+hm.size());
//    
//        for (int i = 0; i < motifList.length; i++) {
//            sequence=(Vector)hm.get(motifList[i][0]);
//            for (int j = 0; j < sequence.size(); j++) {
//                motifList[i][j+1]=(String)sequence.get(j);
//            }
//        }
////        System.out.println(""+max);
//        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Oncocis\\motifFilter\\MotifList.bed", motifList);
        boolean mismatch;
        boolean found;
        for (int i = 0; i < data.length; i++) {
            if(data[i][0].equals("null")){break;}
            found=false;
            sequence=(Vector)hm.get(data[i][0]);            
            data[i][3]="0";
            for (int j = 0; j < sequence.size(); j++) {                
                if(data[i][2].length() == sequence.get(j).toString().length()){ 
                   char motif[]=data[i][2].toUpperCase().toCharArray();
                   char motifDB[]=sequence.get(j).toString().toCharArray();                   
                   for (int k = 0; k < motifDB.length; k++) {
                      if(motifDB[k]=='*'){
                          motif[k]='*';
                      }                     
                   }                   
                   if(String.valueOf(motif).equalsIgnoreCase(String.valueOf(motifDB))){
                       found=true;
                       data[i][3]="1";
                       break;
                   }
                }
                if(found){
                    break;
                }
            }
        }
        DataSetWriter.writeToFile2D(filteredMotifs+"pre", data);
        String cmd = "awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==4 { if($4==\"1\"){print $1, $2;} }' " + filteredMotifs+"pre" +" > " + filteredMotifs ;
        UtilitiesLinux.shellCommandExecuter(cmd);
    }
    
    public static void getAlteredMotifsFinal(
            String normalMotifsCollapsed,
            String mutatedMotifsCollapsed,
            String motifsRemoved,
            String motifsCreated
            ){ 
       HashMap mp=new HashMap<String, HashSet<String>>(); 
        HashSet<String> set;
       
        String normalMotifs[][] = DataSetReader.readDataSet(normalMotifsCollapsed, 3, "\t");
        String mutatedMotifs[][] = DataSetReader.readDataSet(mutatedMotifsCollapsed, 3, "\t");
       
        for (int i = 0; i < mutatedMotifs.length; i++) {
            set=new HashSet<String>();
            String temp[]=mutatedMotifs[i][1].replace("\"", "").split(",");
            for (int j = 0; j < temp.length; j++) {
               set.add(temp[j]);
            }
            mp.put(mutatedMotifs[i][0], set);
          
        }
       
        
            for (int j = 0; j < normalMotifs.length; j++) {
                normalMotifs[j][2]="";
                String temp[]=normalMotifs[j][1].replace("\"", "").split(",");
                set=(HashSet<String>)mp.get(normalMotifs[j][0]);
                if(j==0){
               
                    //System.out.print(set.toString() + "\t");
               
                }

                for (int i = 0; i < temp.length; i++) {
                    try{
                        if(set.contains(temp[i])){

                    }
                        else{
                           //System.out.println("true");
                            String temp2=normalMotifs[j][2];
                            if(temp2.equals("")){
                                normalMotifs[j][2]=temp[i];
                            }
                            else{
                                normalMotifs[j][2]=temp2+";"+temp[i];
                            }

                        }
                    }
                    catch(Exception e){
                        //System.out.println("true");
                        String temp2=normalMotifs[j][2];
                        if(temp2.equals("")){
                            normalMotifs[j][2]=temp[i];
                        }
                        else{
                            normalMotifs[j][2]=temp2+";"+temp[i];
                        }
                    }

                }           
            }
        
               
        DataSetWriter.writeToFile2D(motifsRemoved+ "pre",normalMotifs); 
        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==3 { print $1, $3;}' " + motifsRemoved+"pre > "+motifsRemoved);
        dbConnect.importDataFromFile("dark_matter_temp.possum_removed", motifsRemoved);
        
       
        mp=new HashMap<String, HashSet<String>>();
        for (int i = 0; i < normalMotifs.length; i++) {
            set=new HashSet<String>();
           
           
            String temp[]=normalMotifs[i][1].replace("\"", "").split(",");
            for (int j = 0; j < temp.length; j++) {
               set.add(temp[j]);
            }
            mp.put(normalMotifs[i][0], set);
             if(i==0){
               
                    //System.out.print(set.toString() + "\t");
               
                }
            //System.out.println(""+set.size());
        }
        //System.out.println("HashMap"+mp.size());
        
        for (int j = 0; j < mutatedMotifs.length; j++) {
            mutatedMotifs[j][2]="";
            String temp[]=mutatedMotifs[j][1].replace("\"", "").split(",");
            set=(HashSet<String>)mp.get(mutatedMotifs[j][0]);
            for (int i = 0; i < temp.length; i++) {
                try{
                    if(set.contains(temp[i])){
                   
                }
                else{
                   //System.out.println("true");
                    String temp2=mutatedMotifs[j][2];
                    if(temp2.equals("")){
                            mutatedMotifs[j][2]=temp[i];
                        }
                        else{
                            mutatedMotifs[j][2]=temp2+";"+temp[i];
                        } 
                }}
                catch(Exception e){
                    String temp2=mutatedMotifs[j][2];
                    if(temp2.equals("")){
                            mutatedMotifs[j][2]=temp[i];
                        }
                        else{
                            mutatedMotifs[j][2]=temp2+";"+temp[i];
                        } 
                }
               
            }           
        }
        DataSetWriter.writeToFile2D(motifsCreated+"pre",mutatedMotifs);
        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==3 { print $1, $3;}' " + motifsCreated+"pre > "+motifsCreated);
        dbConnect.importDataFromFile("dark_matter_temp.possum_created", motifsCreated);     
    
    }                    
//        HashMap mp=new HashMap<String, HashSet<String>>();  
//        HashSet<String> set; 
//        
//        String normalMotifs[][] = DataSetReader.readDataSet(normalMotifsCollapsed, 3, "\t"); 
//        String mutatedMotifs[][] = DataSetReader.readDataSet(mutatedMotifsCollapsed, 3, "\t");
//        
//        for (int i = 0; i < mutatedMotifs.length; i++) {
//            set=new HashSet<String>();
//            String temp[]=mutatedMotifs[i][1].replace("\"", "").split(",");
//            for (int j = 0; j < temp.length; j++) {
//               set.add(temp[j]);
//            }
//            mp.put(mutatedMotifs[i][0], set);
//           
//        }
//        
//        try{
//            for (int j = 0; j < normalMotifs.length; j++) {
//                normalMotifs[j][2]="";
//                String temp[]=normalMotifs[j][1].replace("\"", "").split(",");
//                set=(HashSet<String>)mp.get(normalMotifs[j][0]);
//                if(j==0){
//                
//                    System.out.print(set.toString() + "\t");
//                
//                }
//
//                for (int i = 0; i < temp.length; i++) {
//                    if(set.contains(temp[i])){
//
//                    }
//                    else{
//                       System.out.println("true");
//                        String temp2=normalMotifs[j][2];
//                        if(temp2.equals("")){
//                            normalMotifs[j][2]=temp[i];
//                        }
//                        else{
//                            normalMotifs[j][2]=temp2+";"+temp[i];
//                        }
//                          
//                    }
//
//                }            
//            }
//        }
//        
//        catch(Exception e){}
//                
//        DataSetWriter.writeToFile2D(motifsRemoved,normalMotifs);  
//        
//        mp=new HashMap<String, HashSet<String>>();
//        for (int i = 0; i < normalMotifs.length; i++) {
//            set=new HashSet<String>();
//            
//            
//            String temp[]=normalMotifs[i][1].replace("\"", "").split(",");
//            for (int j = 0; j < temp.length; j++) {
//               set.add(temp[j]);
//            }
//            mp.put(normalMotifs[i][0], set);
//             if(i==0){
//                
//                    System.out.print(set.toString() + "\t");
//                
//                }
//            System.out.println(""+set.size());
//        }
//        System.out.println("HashMap"+mp.size());
//        try{
//        for (int j = 0; j < mutatedMotifs.length; j++) {
//            mutatedMotifs[j][2]="";
//            String temp[]=mutatedMotifs[j][1].replace("\"", "").split(",");
//            set=(HashSet<String>)mp.get(mutatedMotifs[j][0]);
//            for (int i = 0; i < temp.length; i++) {
//                if(set.contains(temp[i])){
//                    
//                }
//                else{
//                   System.out.println("true");
//                    String temp2=mutatedMotifs[j][2];
//                    if(temp2.equals("")){
//                            mutatedMotifs[j][2]=temp[i];
//                        }
//                        else{
//                            mutatedMotifs[j][2]=temp2+";"+temp[i];
//                        }  
//                }
//                
//            }            
//        }}
//        catch(Exception e){}
//        DataSetWriter.writeToFile2D(motifsCreated,mutatedMotifs);
//              
//    }
    
    
    public static void updateMotifsInDB(){
        //Created temp
        //awk 'BEGIN { OFS = "\t" } NR>0 && NF==7 { if($2!='null'){print $1, $2;}}' Possum-edited-normal.bed > Possum-edited-normal-reduced.bed
//        dbConnect.execute(  "Insert into dark_matter_temp.motifs_temp_created\n" +
//                            "(Select distinctrow Motif, Location, Start, Stop FROM dark_matter_temp.mutated_motifs\n" +
//                            "WHERE (\n" +
//                            "Motif, Location, Start, Stop, Strand, Score,  Matched_sequence) NOT IN\n" +
//                            "(SELECT * FROM dark_matter_temp.normal_motifs));");
//
//        //Removed temp
//        dbConnect.execute(  "Insert into dark_matter_temp.motifs_temp_removed\n" +
//                            "(Select distinctrow Motif, Location, Start, Stop FROM dark_matter_temp.normal_motifs\n" +
//                            "WHERE (\n" +
//                            "Motif, Location, Start, Stop, Strand, Score, Matched_sequence) NOT IN\n" +
//                            "(SELECT * FROM dark_matter_temp.mutated_motifs));");
//        
//        //Created
//        dbConnect.execute(  "Insert into dark_matter_temp.motifs_created\n" +
//                            "Select distinctrow Motif, Location\n" +
//                            "FROM dark_matter_temp.possum_mutated\n" +
//                            "WHERE (\n" +
//                            "Motif, Location) NOT IN\n" +
//                            "(SELECT Motif, Location FROM dark_matter_temp.possum_normal)");
//        
//        //Removed
//        dbConnect.execute(  "Insert into dark_matter_temp.motifs_removed\n" +
//                            "Select distinctrow Motif, Location\n" +
//                            "FROM dark_matter_temp.possum_normal\n" +
//                            "WHERE (\n" +
//                            "Motif, Location) NOT IN\n" +
//                            "(SELECT Motif, Location FROM dark_matter_temp.possum_mutated)");
//        
        dbConnect.execute(  "Insert into dark_matter_temp.mutations_to_motifs_created\n" +
                            "(select Chromosome, mutations_extended.Start, \n" +
                            "mutations_extended.Stop, MutationID, Motif  \n" +
                            "from dark_matter_temp.possum_created \n" +
                            "join dark_matter_temp.mutations_extended\n" +
                            "Where Location=concat(Chromosome,'-',\n" +
                            "mutations_extended.Start,'-',mutations_extended.Stop));");
        
        dbConnect.execute(  "Insert into dark_matter_temp.mutations_to_motifs_removed\n" +
                            "(select Chromosome, mutations_extended.Start, \n" +
                            "mutations_extended.Stop, MutationID, Motif  \n" +
                            "from dark_matter_temp.possum_removed \n" +
                            "join dark_matter_temp.mutations_extended\n" +
                            "Where Location=concat(Chromosome,'-',\n" +
                            "mutations_extended.Start,'-',mutations_extended.Stop));");
        
       
//Insert into dark_matter_temp.motifs_created (Select Chromosome, Start, Stop, MutationID, count(Motif) as NumOfMotifs, group_concat(Distinct Motif)as Motifs FROM 
//dark_matter_temp.motifs_temp_created
//Join
//dark_matter_temp.mutations_extended
//Where Location=concat(Chromosome,'-',Start,'-',Stop)
//Group by
//MutationID);

//        dbConnect.execute("Truncate table dark_matter_temp.motifs_temp;");
    }
    
    public void getMotifRemoved(){
        
        
        dbConnect.execute("Truncate table dark_matter_temp.motifs_temp;");
    }
    
    public void compareMotifs(){
        //chr1:228535231-228535272
        String mutationList[][] = DataSetReader.readDataSet(mutations.getMutationsReduced(), 4, "\t");
        String normal[][] = DataSetReader.readDataSet(normalMotifs,9,"\t");
        String mutated[][] = DataSetReader.readDataSet(mutatedMotifs,9,"\t");
        String sequenceName = "";
        double normalScore = 0;
        double mutatedScore = 0;
        for (int i = 0; i < mutationList.length; i++) {
            sequenceName = mutationList[i][0]+":"+mutationList[i][1]+"-"+mutationList[i][2];
            for (int j = 0; j < normal.length; j++) {
                if(normal[j][1].equals(sequenceName)){
                    normalScore+=Double.parseDouble(normal[j][5]);
                }
            }
            for (int j = 0; j < mutated.length; j++) {
                if(mutated[j][1].equals(sequenceName)){
                    mutatedScore+=Double.parseDouble(mutated[j][5]);
                }
            }
            mutationList[i][3] = String.valueOf(mutatedScore-normalScore);
        }
        
    
    }
             
}
