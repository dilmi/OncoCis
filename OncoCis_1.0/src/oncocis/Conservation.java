/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package OncoCis;

/**
 *he
 * @author Dilmi
 */
public class Conservation {
    private String scoreLocation = "";
    private String mutationsMappedToScores = "";
    private String outputFolder;
    private Mutations mutations;
    
    public Conservation(String conservationScores, String cOutputFolder, Mutations m){
        scoreLocation = conservationScores;
        outputFolder = cOutputFolder;
        mutations = m;
        this.getScores();
    }
    public static void main(String[] args) {
        Mutations m = new Mutations("/home/bioinfo/Darkmatter-data/DarkMatter/UserInput/Mutations-input.bed", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/");
        Conservation c = new Conservation("/home/bioinfo/Darkmatter-data/DarkMatter/default/PhastCons/BigWig/", m.getOutputLocation(), m);
    }
    
    public void getScores(){
//       String regions[][] = DataSetReader.readDataSet(mutations.getMutationsReduced(), 4, "\t");
//       String output[][] = new String[regions.length][1]; 
//        for (int i = 0; i < regions.length; i++){
//            output[i][0] =  "./bigWigToBedGraph -chrom="+regions[i][0]+" -start="+regions[i][1]+" -end="+regions[i][2]+" "+ scoreLocation + regions[i][0] + ".phastCons46way.placental.bw "+outputFolder+"out.bedGraph"+"\n"+"cat " + outputFolder+"out.bedGraph";     
//        }
//        DataSetWriter.writeToFile2D(outputFolder+"convervation-temp.bed", output);
//        UtilitiesLinux.conservationScore(outputFolder+"convervation-temp.bed", outputFolder+"convervation-scores-temp.bed");
        
       String regions[][] = DataSetReader.readDataSet(mutations.getMutations20bpExtented(), 4, "\t");
       String output[][] = new String[regions.length][1]; 
       output = new String[regions.length][1]; 
        for (int i = 0; i < regions.length; i++){
            output[i][0] =  "./bigWigToBedGraph -chrom="+regions[i][0]+" -start="+regions[i][1]+" -end="+regions[i][2]+" "+ scoreLocation + regions[i][0] + ".phastCons46way.placental.bw "+outputFolder+"out-extended.bedGraph"+"\n"+"cat " + outputFolder+"out-extended.bedGraph";     
        }
        DataSetWriter.writeToFile2D(outputFolder+"convervation-temp-extended.bed", output);
        UtilitiesLinux.conservationScore(outputFolder+"convervation-temp-extended.bed", outputFolder+"convervation-scores-temp-extended.bed");
        dbConnect.importDataFromFile("dark_matter_temp.conservation_scores", outputFolder+"convervation-scores-temp-extended.bed");
//        dbConnect.execute(  "INSERT INTO dark_matter_temp.mutations_to_conservation\n" +
//                            "SELECT mutations_reduced.Chromosome, mutations_reduced.Start, mutations_reduced.Stop,\n" +
//                            "mutations_reduced.MutationID, conservation_scores.Score FROM dark_matter_temp.mutations_reduced JOIN dark_matter_temp.conservation_scores \n" +
//                            "WHERE (conservation_scores.Chromosome=mutations_reduced.Chromosome\n" +
//                            "AND conservation_scores.Start<=mutations_reduced.Start\n" +
//                            "AND conservation_scores.Stop>=mutations_reduced.Stop);");
//       
        dbConnect.execute(  "Insert into dark_matter_temp.Conservation_Score_temp Select distinctRow Chromosome, Start, Stop, Score*(Stop-Start) as newScore, (Stop-Start) as gap FROM dark_matter_temp.conservation_scores;" );
        
        dbConnect.execute(  "Insert into dark_matter_temp.mutations_to_conservation\n" +
                            "	SELECT mutations_reduced.Chromosome,\n" +
                            "    mutations_reduced.Start,\n" +
                            "    mutations_reduced.Stop,\n" +
                            "    mutations_reduced.MutationID,\n" +
                            "	sum(Conservation_Score_temp.Score)/sum(Conservation_Score_temp.Gap)\n" +
                            "FROM\n" +
                            "    dark_matter_temp.mutations_reduced\n" +
                            "        JOIN\n" +
                            "    dark_matter_temp.Conservation_Score_temp\n" +
                            "WHERE\n" +
                            "    (Conservation_Score_temp.Chromosome = mutations_reduced.Chromosome\n" +
                            "        AND Conservation_Score_temp.Start <= mutations_reduced.Start\n" +
                            "        AND Conservation_Score_temp.Stop >= mutations_reduced.Stop)\n" +
                            "GROUP BY\n" +
                            "    MutationID;");
        
        dbConnect.execute(  "UPDATE dark_matter_temp.mutations_to_conservation \n" +
                            "SET Score='1' WHERE Score>1;");
        
        dbConnect.execute(  "Insert into dark_matter_temp.mutations_to_backgroundConservation\n" +
                            "	SELECT mutations_extended.Chromosome,\n" +
                            "    mutations_extended.Start,\n" +
                            "    mutations_extended.Stop,\n" +
                            "    mutations_extended.MutationID,\n" +
                            "	sum(Conservation_Score_temp.Score)/sum(Conservation_Score_temp.Gap)\n" +
                            "FROM\n" +
                            "    dark_matter_temp.mutations_extended\n" +
                            "        JOIN\n" +
                            "    dark_matter_temp.Conservation_Score_temp\n" +
                            "WHERE\n" +
                            "    (Conservation_Score_temp.Chromosome = mutations_extended.Chromosome\n" +
                            "        AND Conservation_Score_temp.Start >= mutations_extended.Start\n" +
                            "        AND Conservation_Score_temp.Stop <= mutations_extended.Stop)\n" +
                            "GROUP BY\n" +
                            "    MutationID;");
        
        dbConnect.execute(  "UPDATE dark_matter_temp.mutations_to_backgroundConservation \n" +
                            "SET Score='1' WHERE Score>1;");
        

    }

    public String getScoreLocation() {
        return scoreLocation;
    }

    public void setScoreLocation(String scoreLocation) {
        this.scoreLocation = scoreLocation;
    }

    public String getMutationsMappedToScores() {
        return mutationsMappedToScores;
    }

    public void setMutationsMappedToScores(String mutationsMappedToScores) {
        this.mutationsMappedToScores = mutationsMappedToScores;
    }

    public String getOutputFolder() {
        return outputFolder;
    }

    public void setOutputFolder(String outputFolder) {
        this.outputFolder = outputFolder;
    }

    public Mutations getMutations() {
        return mutations;
    }

    public void setMutations(Mutations mutations) {
        this.mutations = mutations;
    }
    
    
}

//        "INSERT INTO dark_matter_temp.mutations_to_conservation
//SELECT mutations_reduced.Chromosome, mutations_reduced.Start, mutations_reduced.Stop,
//mutations_reduced.MutationID, conservation_scores.Score FROM dark_matter_temp.mutations_reduced JOIN dark_matter_temp.conservation_scores 
//WHERE (conservation_scores.Chromosome=mutations_reduced.Chromosome
//AND conservation_scores.Start<=mutations_reduced.Start
//AND conservation_scores.Stop>=mutations_reduced.Stop)"
//        
//        Insert into dark_matter_temp.Conservation_Score_temp Select distinctRow Chromosome, Start, Stop, Score*(Stop-Start) as newScore, (Stop-Start) as gap FROM dark_matter_temp.conservation_scores
//SELECT
//	mutations_extended.Chromosome,
//    mutations_extended.Start,
//    mutations_extended.Stop,
//    mutations_extended.MutationID,
//	sum(Conservation_Score_temp.Score)/sum(Conservation_Score_temp.Gap)
//FROM
//    dark_matter_temp.mutations_extended
//        JOIN
//    dark_matter_temp.Conservation_Score_temp
//WHERE
//    (Conservation_Score_temp.Chromosome = mutations_extended.Chromosome
//        AND Conservation_Score_temp.Start >= mutations_extended.Start
//        AND Conservation_Score_temp.Stop <= mutations_extended.Stop)
//GROUP BY
//    MutationID;