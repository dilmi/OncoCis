/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package OncoCis;

import java.util.ArrayList;

/**
 *
 * @author Dilmi
 */
public class FimoMotif {
    public static String hg19fa = "";
    public static String motifDB = "";
    public static String normalSequence = "";
    public static String substitutionFile = "";
    public static String mutatedSequence = "";
    public static String normalMotifs = "";
    public static String mutatedMotifs = "";
    public static String fimoOutputFolder1 = "";
    public static String fimoOutputFolder2 = "";
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
    public FimoMotif(String cHg19fa, String cMotifDB, String cOutputFolder, Mutations m) {
        hg19fa = cHg19fa;
        motifDB = cMotifDB;
        normalSequence = cOutputFolder + "NormalSequence.fasta";
        substitutionFile = cOutputFolder + "substitutions.bed";
        mutatedSequence = cOutputFolder + "MutatedSequence.fasta";
        normalMotifs = cOutputFolder + "normalMotifs.bed";
        mutatedMotifs = cOutputFolder + "mutatedMotifs.bed";
        mutations = m;
        UtilitiesLinux.shellCommandExecuter("mkdir " + cOutputFolder + "/Motifs1");
        fimoOutputFolder1 = cOutputFolder + "Motifs1";
        UtilitiesLinux.shellCommandExecuter("mkdir " + cOutputFolder + "/Motifs2");
        fimoOutputFolder2 = cOutputFolder + "Motifs2";
        this.getSequence();
        this.getMutatedSequence();
        getMotifs();
        getMotifCreatedOrRemoved();
//        getMotifRemoved();
    }  

    public static void main(String[] args) {
//        dbConnect.connectDriver();
//        dbConnect.execute(  "Select \n" +
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
    
    public void getMotifs(){
        UtilitiesLinux.shellCommandExecuter("sed -i 's/:/-/g' " + normalSequence);
        UtilitiesLinux.shellCommandExecuter("sed -i 's/:/-/g' " + mutatedSequence);
        UtilitiesLinux.fimo(motifDB, normalSequence, fimoOutputFolder1);
        normalMotifs = fimoOutputFolder1+"/fimo.txt";
//        normalMotifs = fimoOutputFolder1+"/normalMotifs.bed";
//        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==9 { print $1,$2,$5,$6,$7,$8,$9;}' "+fimoOutputFolder1+"/fimo.txt"+" > "+normalMotifs);
//        System.out.println("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==9 { print $1,$2,$5,$6,$7,$8,$9;}' "+fimoOutputFolder1+"/fimo.txt"+" > "+normalMotifs);
        dbConnect.importDataFromFile("dark_matter_temp.normal_motifs", normalMotifs);
        
        UtilitiesLinux.fimo(motifDB, mutatedSequence, fimoOutputFolder2);
        mutatedMotifs = fimoOutputFolder2+"/fimo.txt";
//        mutatedMotifs = fimoOutputFolder2+"/mutatedMotifs.bed";
//         UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==9 { print $1,$2,$5,$6,$7,$8,$9;}' "+fimoOutputFolder2+"/fimo.txt"+" > "+mutatedMotifs);
//        System.out.println("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==9 { print $1,$2,$5,$6,$7,$8,$9;}' "+fimoOutputFolder2+"/fimo.txt"+" > "+mutatedMotifs);
         dbConnect.importDataFromFile("dark_matter_temp.mutated_motifs", mutatedMotifs);
        
        
        
    }
    
    public void getMotifCreatedOrRemoved(){
        //Created temp
        dbConnect.execute(  "Insert into dark_matter_temp.motifs_temp_created\n" +
                            "(Select distinctrow Motif, Location, Start, Stop FROM dark_matter_temp.mutated_motifs\n" +
                            "WHERE (\n" +
                            "Motif, Location, Start, Stop, Strand, Score, pValue, qValue, Matched_sequence) NOT IN\n" +
                            "(SELECT * FROM dark_matter_temp.normal_motifs));");

        //Removed temp
        dbConnect.execute(  "Insert into dark_matter_temp.motifs_temp_removed\n" +
                            "(Select distinctrow Motif, Location, Start, Stop FROM dark_matter_temp.normal_motifs\n" +
                            "WHERE (\n" +
                            "Motif, Location, Start, Stop, Strand, Score, pValue, qValue, Matched_sequence) NOT IN\n" +
                            "(SELECT * FROM dark_matter_temp.mutated_motifs));");
        
        //Created
        dbConnect.execute(  "Insert into dark_matter_temp.motifs_created\n" +
                            "Select distinctrow Motif, Location\n" +
                            "FROM dark_matter_temp.motifs_temp_created\n" +
                            "WHERE (\n" +
                            "Motif, Location) NOT IN\n" +
                            "(SELECT Motif, Location FROM dark_matter_temp.motifs_temp_removed)");
        
        //Removed
        dbConnect.execute(  "Insert into dark_matter_temp.motifs_removed\n" +
                            "Select distinctrow Motif, Location\n" +
                            "FROM dark_matter_temp.motifs_temp_removed\n" +
                            "WHERE (\n" +
                            "Motif, Location) NOT IN\n" +
                            "(SELECT Motif, Location FROM dark_matter_temp.motifs_temp_created)");
        
        dbConnect.execute(  "Insert into dark_matter_temp.mutations_to_motifs_created\n" +
                            "(select Chromosome, mutations_extended.Start, mutations_extended.Stop, MutationID, count(Motif) as NumOfMotifs, group_concat(Distinct Motif)as Motifs  \n" +
                            "from dark_matter_temp.motifs_created join dark_matter_temp.mutations_extended \n" +
                            "Where Location=concat(Chromosome,'-',mutations_extended.Start,'-',mutations_extended.Stop)\n" +
                            "Group by\n" +
                            "MutationID);");
        
        dbConnect.execute(  "Insert into dark_matter_temp.mutations_to_motifs_removed\n" +
                            "(select Chromosome, mutations_extended.Start, mutations_extended.Stop, MutationID, count(Motif) as NumOfMotifs, group_concat(Distinct Motif)as Motifs  \n" +
                            "from dark_matter_temp.motifs_removed join dark_matter_temp.mutations_extended \n" +
                            "Where Location=concat(Chromosome,'-',mutations_extended.Start,'-',mutations_extended.Stop)\n" +
                            "Group by\n" +
                            "MutationID);");
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
