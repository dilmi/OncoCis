/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package OncoCis;

/**
 *
 * @author Dilmi
 */
public class DNase1 {    
    private String DHS="";    
    private String mutationsIntersectDHS="";
    private String mutationsMappedToDHS="";
    private Mutations mutations;
    private String outputLocation;
    
    public DNase1(String cDHS, String cOutputLocation, Mutations cmut){        
        DHS = cDHS;
        outputLocation = cOutputLocation;
        mutations = cmut;
        mutationsIntersectDHS = cOutputLocation+"MutationsIntersectDHS.bed";
        mutationsMappedToDHS = cOutputLocation+"MutationsMappedToDHS.bed";
        this.overlap();
       
    }
    
//    public static void main(String[] args) {
//        Mutations m = new Mutations("/home/dilmi/DarkMatter/UserInput/test/Mutations-input.bed", "/home/dilmi/DarkMatter/Output/test/");
//        DNase1 d = new DNase1("/home/dilmi/DarkMatter/Cell_type_specific_data/HMEC/DNase1.bed", "/home/dilmi/DarkMatter/Output/test/", m);
//    }
    public void exportToFile(String table, String file) {
        dbConnect.exportTableToFile(table, file , outputLocation);
    }
       
    
    public void overlap(){
        UtilitiesLinux.bedtoolsIntersect(mutations.getMutationsReduced(), DHS, mutationsIntersectDHS, "-c" );
        UtilitiesLinux.bedtoolsIntersect(mutations.getMutationsReduced(), DHS, mutationsMappedToDHS, "-wo" );
        dbConnect.importDataFromFile( "dark_matter_temp.mutations_to_dhs",mutationsIntersectDHS);
        dbConnect.importDataFromFile("dark_matter_temp.mutations_mapped_to_dhs",mutationsMappedToDHS);
    }
    
    
    
}
