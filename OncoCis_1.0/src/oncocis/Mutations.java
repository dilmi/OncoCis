/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package OncoCis;

/**
 *
 * @author dilmi
 */
class Mutations {
    private String mutations="";
    private String mutationsReduced="";
    private String mutations20bpExtented="";
    private String mutations500bpExtented="";
    private String outputLocation;

    public Mutations(String cmutations, String cOutputLocation) {
        outputLocation = cOutputLocation;
        mutations=cmutations;
//        mutations = overlapWithExons(cmutations);
        mutationsReduced = cOutputLocation+"Mutation-reduced.bed";
        dbConnect.importDataFromFile("dark_matter_temp.mutations", mutations);
        dbConnect.execute("Delete from dark_matter_temp.mutations where MutationID=0;");
        this.generatereducedMutationFile();
        mutations20bpExtented = cOutputLocation+"Mutation-20bpExtended.bed";
        mutations500bpExtented = cOutputLocation+"Mutation-500bpExtended.bed";
        this.generateExtendedMutationFile();
        
    }
    
    public String overlapWithExons(String mutationFile){
        String newMutFile=mutationFile;
        UtilitiesLinux.shellCommandExecuter("mv "+mutationFile 
        +" "+outputLocation+"allMutations.bed");
        UtilitiesLinux.bedtoolsIntersect(outputLocation+"allMutations.bed", OncoCis.MAIN_FOLDER+"default/RefSeqGenes-Exons.bed", mutationFile, "-v");
        return newMutFile;
    }
    
    public void generatereducedMutationFile(){
        //" + MUTATIONS + " > "+ OUTPUT_FOLDER+"mutations-reduced.bed
        String cmd = "awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==8 { print $1, $2, $3, $4;}' " + mutations + " > "+ mutationsReduced;
        UtilitiesLinux.shellCommandExecuter(cmd);
        dbConnect.importDataFromFile("dark_matter_temp.mutations_reduced", mutationsReduced);
    } 
    
    public void generateExtendedMutationFile(){
        String cmd = "awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==4 { print $1, ($2-20), ($3+20), $4;}' " + mutationsReduced + " > "+ mutations20bpExtented;
        UtilitiesLinux.shellCommandExecuter(cmd);
        dbConnect.importDataFromFile("dark_matter_temp.mutations_extended", mutations20bpExtented);
        cmd = "awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==4 { print $1, ($2-500), ($3+500), $4;}' " + mutationsReduced + " > "+ mutations500bpExtented;
        UtilitiesLinux.shellCommandExecuter(cmd);
//        dbConnect.importDataFromFile("dark_matter_temp.mutations_extended", mutations20bpExtented);
    }
    
   

    public String getMutations() {
        return mutations;
    }

    public void setMutations(String mutations) {
        this.mutations = mutations;
    }

    public String getMutationsReduced() {
        return mutationsReduced;
    }

    public void setMutationsReduced(String mutationsReduced) {
        this.mutationsReduced = mutationsReduced;
    }

    public String getMutations20bpExtented() {
        return mutations20bpExtented;
    }

    public void setMutations20bpExtented(String mutations20bpExtented) {
        this.mutations20bpExtented = mutations20bpExtented;
    }

    public String getOutputLocation() {
        return outputLocation;
    }

    public void setOutputLocation(String outputLocation) {
        this.outputLocation = outputLocation;
    }

    public String getMutations500bpExtented() {
        return mutations500bpExtented;
    }

    public void setMutations500bpExtented(String mutations500bpExtented) {
        this.mutations500bpExtented = mutations500bpExtented;
    }

//    String getMutations150bpExtented() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
    
    
}
