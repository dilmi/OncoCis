/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package OncoCis;

/**
 *
 * @author Dilmi
 */
public class Histones {
    private String mutations;
    private String outputLocation;
    private String H3K4me1;
    private String H3K4me3;
    private String H3K27ac;
    
    public Histones(String cH3K4me1, String cH3K4me3, String cH3K27ac, Mutations cMutations) {
        H3K4me1 = cH3K4me1;
        H3K4me3 = cH3K4me3;
        H3K27ac = cH3K27ac;
        UtilitiesLinux.fromDosCheck(H3K27ac);
//        mutations = cMutations.getMutationsReduced();
        outputLocation = cMutations.getOutputLocation();
        this.generateFlankingRegions();
        this.intersect();
    }

    public static void main(String[] args) {
        dbConnect.exportDataToFile("SELECT * FROM dark_matter_temp.Final_results", "Final_Results.bed", "~/Darkmatter-data/DarkMatter/Output/EIP6Y/");
    }

    public void generateFlankingRegions(){
        System.out.println("Flanking regions");
        dbConnect.exportDataToFile("SELECT `DHS_Chr`,`DHS_Start`,`DHS_Stop`,"
                + "`MutationID`,`DHS` FROM dark_matter_temp.mutations_mapped_to_dhs","DHS-for-Flanks.bed",outputLocation);
        String cmd1="awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==5 {if(($3-$2)%2==1){ print $1, (($2+$3)+1)/2, ((($2+$3)+1)/2)+1, $4, $5;} else{print $1, ($2+$3)/2, (($2+$3)/2)+1, $4, $5;}}' " + outputLocation+ "DHS-for-Flanks.bed" + " > "+ outputLocation + "mutations-in-DHS-forFlanks.bed";
//        System.out.println(""+cmd1);
        UtilitiesLinux.shellCommandExecuter(cmd1);
        dbConnect.exportDataToFile("SELECT * FROM dark_matter_temp.mutations_to_dhs WHERE DHS=0", "mutations-not-in-DHS.bed", outputLocation);
        mutations = outputLocation + "mutationsForFlankGen.bed";
        UtilitiesLinux.shellCommandExecuter("cat "+outputLocation+"mutations-in-DHS-forFlanks.bed "+outputLocation+"mutations-not-in-DHS.bed > "+ outputLocation+"mutationsForFlankGen.bed");
        String cmd = "awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==5 { print $1, ($2-500), ($3-150), $4;}' " + mutations + " > "+ outputLocation + "leftFlank.bed";
        UtilitiesLinux.shellCommandExecuter(cmd);
        cmd = "awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==5 { print $1, ($2+150), ($3+500), $4;}' " + mutations + " > "+ outputLocation + "rightFlank.bed";
        UtilitiesLinux.shellCommandExecuter(cmd);
        cmd = "cat "+ outputLocation + "leftFlank.bed" + " " + outputLocation + "rightFlank.bed > "+outputLocation+"flankingRegion.bed";
        UtilitiesLinux.shellCommandExecuter(cmd);
    }
    
    public void intersect(){
      UtilitiesLinux.bedtoolsIntersect(outputLocation+"flankingRegion.bed", H3K4me1,outputLocation+"flankingRegionOverlapH3K4me1.bed", "-c" );  
      UtilitiesLinux.bedtoolsIntersect(outputLocation+"flankingRegionOverlapH3K4me1.bed", H3K4me3,outputLocation+"flankingRegionOverlapH3K4me1H3K4me3.bed", "-c" );  
      UtilitiesLinux.bedtoolsIntersect(outputLocation+"flankingRegionOverlapH3K4me1H3K4me3.bed", H3K27ac,outputLocation+"flankingRegionOverlapH3K4me1H3K4me3H3K27ac.bed", "-c" );  
      
        System.out.println("intersect complete");
        dbConnect.importDataFromFile("dark_matter_temp.mutations_to_H3K4me1_H3K4me3_H3K27ac", outputLocation+"flankingRegionOverlapH3K4me1H3K4me3H3K27ac.bed");
        System.out.println("imported data");
        dbConnect.execute("INSERT INTO dark_matter_temp.temptable SELECT MutationID, SUM(H3K4me1), SUM(H3K4me3), SUM(H3K27ac) "
              + "FROM dark_matter_temp.mutations_to_H3K4me1_H3K4me3_H3K27ac GROUP BY MutationID;");
      dbConnect.execute("INSERT INTO dark_matter_temp.mutations_to_histones \n" +
                        "SELECT mutations_reduced.Chromosome, mutations_reduced.Start, mutations_reduced.Stop, \n" +
                        "mutations_reduced.MutationID, temptable.H3K4me1, temptable.H3K4me3, temptable.H3K27ac \n" +
                        "FROM dark_matter_temp.temptable INNER JOIN dark_matter_temp.mutations_reduced \n" +
                        "ON dark_matter_temp.temptable.MutationID=dark_matter_temp.mutations_reduced.MutationID;");
      
       
        dbConnect.execute("Update dark_matter_temp.mutations_to_histones Set H3K4me1=1 where H3K4me1>0;");
        dbConnect.execute("Update dark_matter_temp.mutations_to_histones Set H3K4me3=1 where H3K4me3>0;");
        dbConnect.execute("Update dark_matter_temp.mutations_to_histones Set H3K27ac=1 where H3K27ac>0;");

    }

    public String getMutations() {
        return mutations;
    }

    public void setMutations(String mutations) {
        this.mutations = mutations;
    }

    public String getOutputLocation() {
        return outputLocation;
    }

    public void setOutputLocation(String outputLocation) {
        this.outputLocation = outputLocation;
    }

    public String getH3K4me1() {
        return H3K4me1;
    }

    public void setH3K4me1(String H3K4me1) {
        this.H3K4me1 = H3K4me1;
    }

    public String getH3K4me3() {
        return H3K4me3;
    }

    public void setH3K4me3(String H3K4me3) {
        this.H3K4me3 = H3K4me3;
    }

    public String getH3K27ac() {
        return H3K27ac;
    }

    public void setH3K27ac(String H3K27ac) {
        this.H3K27ac = H3K27ac;
    }
    
    
    
}
