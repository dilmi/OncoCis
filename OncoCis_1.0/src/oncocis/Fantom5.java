/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package OncoCis;

/**
 *
 * @author bioinfo
 */
public class Fantom5 {
    private String fantomPromoters = OncoCis.MAIN_FOLDER + "default/Fantom5/TSS_human.bed";
    private String fantomEnhancers = OncoCis.MAIN_FOLDER + "default/Fantom5/permissive_enhancers.bed";    
    private String mutationsIntersectEnhancers="";
    private String mutationsIntersectPromoters="";
    private String mutations;
    private String outputLocation;
    
    public Fantom5(String cOutputLocation, Mutations cmut){
        
        outputLocation = cOutputLocation;
        mutations = cmut.getMutations500bpExtented();
        mutationsIntersectEnhancers = cOutputLocation+"MutationsIntersectFantom5.bed";
        mutationsIntersectPromoters = cOutputLocation+"MutationsIntersectPromoters.bed";
        this.overlap();
             
    } 
    
    public static void main(String[] args) {
        OncoCis.clearDB();
        Preprocessor p=new Preprocessor("/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/Mutations-input.bed", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/");
        Mutations mut=new Mutations("/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/17SamplesOncoCisInput.bed", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/");
        GeneMapping gm = new GeneMapping("/home/bioinfo/Darkmatter-data/DarkMatter/default/regDoms1000kb.bed"
                , mut.getOutputLocation(), mut);
        Fantom5 f=new Fantom5(mut.getOutputLocation(), mut);
    }
     
    
    
    public void overlap(){
        System.out.println("Fantom5");
        UtilitiesLinux.shellCommandExecuter("");
        UtilitiesLinux.bedtoolsIntersect(mutations, fantomPromoters, mutationsIntersectPromoters, "-c" );
        UtilitiesLinux.bedtoolsIntersect(mutationsIntersectPromoters, fantomEnhancers, mutationsIntersectEnhancers, "-c" );
        dbConnect.importDataFromFile( "dark_matter_temp.mutations_to_Fantom",mutationsIntersectEnhancers);
//        dbConnect.execute("UPDATE dark_matter_temp.mutations_to_Fantom\n" +
//            "inner join dark_matter_temp.mutations_to_genes\n" +
//            "ON (mutations_to_Fantom.MutationID=mutations_to_genes.MutationID)\n" +
//            "SET mutations_to_Fantom.Promoter='0'\n" +
//            "where mutations_to_Fantom.Promoter='1' and\n" +
//            " (mutations_to_genes.DistanceToTSS>1000 or \n" +
//            "mutations_to_genes.DistanceToTSS<-10000);");
        dbConnect.execute("UPDATE dark_matter_temp.mutations_to_Fantom\n" +
            "SET mutations_to_Fantom.Enhancer='1'\n" +
            "where mutations_to_Fantom.Enhancer>'1' ");
        dbConnect.execute("UPDATE dark_matter_temp.mutations_to_Fantom\n" +
            "SET mutations_to_Fantom.Promoter='1'\n" +
            "where mutations_to_Fantom.Promoter>'1' ");
        
        //mutations_to_Fantom.Promoter
    }
    
}
