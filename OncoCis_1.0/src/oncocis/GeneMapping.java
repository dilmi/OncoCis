/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package OncoCis;


/**
 *
 * @author Dilmi
 */
public class GeneMapping {
    private String geneMap = "";
    private String mutationsMapd2Genes = "";
    private Mutations mutations;
    
    public GeneMapping(String location, String outputFolder, Mutations m){
        geneMap=location;
        mutationsMapd2Genes = outputFolder+"MutationsMappedToGenes.bed";
        mutations = m;
        this.map2Gene(outputFolder+"temp-genemap.bed");
        this.mapUsingFantom5();
    }
    
    public static void main(String[] args) {
        OncoCis.clearDB();
        Preprocessor p=new Preprocessor("/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/Mutations-input.bed", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/");
        Mutations m=new Mutations("/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/17SamplesOncoCisInput.bed", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/Test/");
        GeneMapping g = new GeneMapping("/home/bioinfo/Darkmatter-data/DarkMatter/default/regDoms1000kb.bed", m.getOutputLocation(), m);
    }
    
    public void map2Gene(String temp){
        UtilitiesLinux.bedtoolsIntersect(mutations.getMutationsReduced(), geneMap, temp , "-wao");
        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==11 { if((($9-$3)>1000000)||(($9-$3)<-1000000)){print $1, $2, $3, $4, \"NONE\", '10000000';}  else {print $1, $2, $3, $4, $8, ($9-$3) ;}}' " + temp + " > "+ mutationsMapd2Genes );
        dbConnect.importDataFromFile("dark_matter_temp.mutations_to_genes", mutationsMapd2Genes);
    
    }
    
    public void mapUsingFantom5(){
        UtilitiesLinux.shellCommandExecuter("bedtools intersect -wao -a "+mutations.getMutationsReduced()+" -b "+OncoCis.MAIN_FOLDER+"default/Fantom5/enhancerToGene.bed"+" > "+mutations.getOutputLocation()+"mutationsMappedToGeneFantom5");
        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF==10 {if($8!=\".\"){print $1,$2,$3,$4,$8,($9-$3);}}' "+mutations.getOutputLocation()+"mutationsMappedToGeneFantom5 > "+mutations.getOutputLocation()+"mutationsMappedToGeneFantom5-reduced");
        dbConnect.importDataFromFile("dark_matter_temp.mutations_to_genes_Fantom5", mutations.getOutputLocation()+"mutationsMappedToGeneFantom5-reduced");
        dbConnect.execute("UPDATE dark_matter_temp.mutations_to_genes\n" +
                        "inner join dark_matter_temp.mutations_to_genes_Fantom5\n" +
                        "ON (mutations_to_genes.MutationID=mutations_to_genes_Fantom5.MutationID)\n" +
                        "SET mutations_to_genes.Gene = mutations_to_genes_Fantom5.Gene,\n" +
                        "mutations_to_genes.DistanceToTSS = mutations_to_genes_Fantom5.DistanceToTSS;");
    }
    
    
}
