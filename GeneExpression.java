/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package OncoCis;

import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author Dilmi
 */
public class GeneExpression {

    private int numberOfCols = 0;
    private String geneExpFile;
    private String sampleToGeneFile;
    private String geneExpMatrix;
    private String outpuFolder;
    private static int numOfMutations;

    public GeneExpression(String cgeneExpFile, String coutputFolder) {
        geneExpFile = cgeneExpFile;
        outpuFolder = coutputFolder;
        numberOfCols =  Preprocessor.getNumOfCols(geneExpFile , outpuFolder);
        sampleToGeneFile = outpuFolder+"sampleToGene.bed";   
        generateSampleToGene();
        filterGeneExpFile();
//        dbConnect.importDataFromFile("geneToSampleToFoldChange", outpuFolder+"t-testforgenes.bed");
 

    }

    public static void main(String[] args) {
        
//        GeneExpression ge = new GeneExpression("/home/bioinfo/Darkmatter-data/DarkMatter/Output/EIP6Y/GeneExpUserInput.bed",
//                "/home/bioinfo/Darkmatter-data/DarkMatter/Output/EIP6Y/");
        dbConnect.connectDriver();
//        dbConnect.importDataFromFile("mutations_to_geneExp", "/home/bioinfo/Darkmatter-data/DarkMatter/Output/EIP6Y/"+"mutationsMappedToGeneExp.bed");
        
        OncoCis.getFinalResults("EIP6Y", 94502, "HMEC",true);
     
    }

    //head -n1 test |  sed 's/\t/\n/g' | wc -l
    public void generateSampleToGene(){        
        dbConnect.exportDataToFile("SELECT distinct Gene,SampleID,'0','0' FROM dark_matter_temp.Final_results WHERE Gene!='NONE' and DHS='1'", "sampleToGene.bed", outpuFolder);
    }
    
    public void filterGeneExpFile() {
        dbConnect.exportDataToFile("SELECT distinct Gene,'0' FROM dark_matter_temp.Final_results WHERE Gene!='NONE' and DHS=1", "geneList.bed", outpuFolder);
        String genes[][] = DataSetReader.readDataSet(outpuFolder + "geneList.bed", 2, "\t");
        String matrix[][] = DataSetReader.readDataSet(geneExpFile, numberOfCols + 1, "\t");
        int count = 0;
        boolean found = false;
        matrix[0][numberOfCols] = "median";
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 0; j < genes.length; j++) {
                if (matrix[i][0].equalsIgnoreCase(genes[j][0])) {
                    
                    genes[j][1] = String.valueOf(Integer.parseInt(genes[j][1]) + 1);
                    if (Integer.parseInt(genes[j][1]) > 1) {
                        matrix[i][numberOfCols] = "0";
                    }
                    else{
                        matrix[i][numberOfCols] = "1";
                        numOfMutations++;
                        count++;
                    }
                    
                    found = true;
                }
            }
            if (!found) {
                matrix[i][numberOfCols] = "0";
            }
            found = false;
        }
//        DataSetWriter.writeToFile2D(outpuFolder+"GeneExpression-edited.bed",matrix);
        System.out.println("FILTERING");
        String geneExp[][] = new String[count+1][numberOfCols+1];
        String geneExpMat[][] = new String [geneExp.length][geneExp[0].length]; 
       
        geneExp[0] = matrix[0];
        geneExpMat[0]= matrix[0];
        int index=1;
        for (int i = 1; i < matrix.length; i++) {
            if(matrix[i][numberOfCols].equals("1")){
                geneExp[index] = matrix[i];
                geneExpMat[index] = matrix[i];
                index++;
            }
        }
        DataSetWriter.writeToFile2D(outpuFolder+"genExp.bed", geneExp);
        
//GENERATE GENE EXP MATRIX
        System.out.println("GENERATE GENE EXP MATRIX");
         String sampleTogeneMatrix[][] = DataSetReader.readDataSet(sampleToGeneFile, 4, "\t");
        
        for (int i = 1; i < geneExpMat.length; i++) {
            for (int j = 1; j < geneExpMat[0].length; j++) {
                geneExpMat[i][j] = "0";
            }
        }
        
        
        for (int i = 1; i < geneExpMat.length; i++) {
            for (int j = 1; j < geneExpMat[0].length; j++) {
                for (int k = 0; k < sampleTogeneMatrix.length; k++) {
                    if (geneExpMat[i][0].equalsIgnoreCase(sampleTogeneMatrix[k][0]) && geneExpMat[0][j].equalsIgnoreCase(sampleTogeneMatrix[k][1])) {
                        geneExpMat[i][j] = "1";
                    }
                }
            }         
        }
        
        DataSetWriter.writeToFile2D(outpuFolder+"genExpMat.bed", geneExpMat);

//GET FOLDCHANGE
        geneExp = DataSetReader.readDataSet(outpuFolder+"genExp.bed",numberOfCols+3,"\t");
        geneExp[0][numberOfCols]="Median/Mean";
        geneExp[0][numberOfCols+1]="std";
        geneExp[0][numberOfCols+2]="numOfSamples";
        
        System.out.println("Calculate foldchange");
        for (int j = 1; j < geneExp.length; j++) {
            ArrayList<Double> expression = new ArrayList<Double>();
            for (int i = 1; i < geneExpMat.length; i++) {
                if (geneExpMat[i][0].equalsIgnoreCase(geneExp[j][0])) {
                    for (int k = 1; k < geneExpMat[1].length; k++) {
                        if (geneExpMat[i][k].equals("0")) {
                            expression.add(Double.parseDouble(geneExp[j][k]));
                        }
                    }
                    Collections.sort(expression);
                    if (!expression.isEmpty()) {
                        if ((expression.size() % 2) == 1) {
                            geneExp[j][numberOfCols] = String.valueOf(expression.get(((expression.size() + 1) / 2) - 1).doubleValue());
                        } else if ((expression.size() % 2) == 0) {
                            geneExp[j][numberOfCols] = String.valueOf((expression.get((expression.size() / 2) - 1).doubleValue() + expression.get(((expression.size() / 2) + 1) - 1).doubleValue()) / 2);
                        }
                    }

                }
            }
        }

        for (int i = 0; i < sampleTogeneMatrix.length; i++) {
            for (int j = 1; j < geneExp.length; j++) {
                if (sampleTogeneMatrix[i][0].equalsIgnoreCase(geneExp[j][0])) {
                    for (int k = 1; k < geneExp[1].length-3; k++) {
                        if (sampleTogeneMatrix[i][1].equalsIgnoreCase(geneExp[0][k])) {
                            if(Double.parseDouble(geneExp[j][numberOfCols])==0){
                                sampleTogeneMatrix[i][2] = String.valueOf((float)(Double.parseDouble(geneExp[j][k]) / 0.01));
                            }
                            else{
                                sampleTogeneMatrix[i][2] = String.valueOf((float)(Double.parseDouble(geneExp[j][k]) / Double.parseDouble(geneExp[j][numberOfCols])));                                
                            }
                            
                        }
                    }
                }
            }
        }
       DataSetWriter.writeToFile2D(outpuFolder+"genesWithFoldChange.bed", sampleTogeneMatrix);
     
//PERFORM T-TEST
       double total = 0;
        double mean = 0;
        double std = 0;
        int n = 0;
        System.out.println("Calculate p-value");
        for (int j = 1; j < geneExp.length; j++) {
            ArrayList<Double> expression = new ArrayList<Double>();
            std = 0;
            mean = 0;
            total = 0;
            n = 0;
            for (int i = 1; i < geneExpMat.length; i++) {
                if (geneExpMat[i][0].equalsIgnoreCase(geneExp[j][0])) {
                    for (int k = 1; k < geneExpMat[1].length-3; k++) {
                        if (geneExpMat[i][k].equals("0")) {
                            expression.add(Double.parseDouble(geneExp[j][k]));
                        }
                    }
                    //Mean
                    for (int k = 0; k < expression.size(); k++) {
                        total += expression.get(k);
                    }
                    mean = total / expression.size();
//                    System.out.println("Mean" + mean);
                    geneExp[j][numberOfCols] = String.valueOf(mean);

                    //Standard Deviation
                    for (int k = 0; k < expression.size(); k++) {
                        std += (expression.get(k) - mean) * (expression.get(k) - mean);
                    }
                    geneExp[j][numberOfCols+1] = String.valueOf(Math.sqrt(std / (expression.size() - 1)));
//                    System.out.println("" + geneExp[j][numberOfCols+1]);
                    //Number of samples
//                    System.out.println("" + (expression.size()));
                    geneExp[j][numberOfCols+2] = String.valueOf(expression.size());
                }
            }
        }
        
        for (int i = 0; i < sampleTogeneMatrix.length; i++) {
            double t = 0;
            int df = 0;
            for (int j = 1; j < geneExp.length; j++) {
                if (sampleTogeneMatrix[i][0].equalsIgnoreCase(geneExp[j][0])) {
                    for (int k = 1; k < geneExp[1].length-3; k++) {
                        if (sampleTogeneMatrix[i][1].equalsIgnoreCase(geneExp[0][k])) {
                            mean = Double.parseDouble(geneExp[j][numberOfCols]);
                            std = Double.parseDouble(geneExp[j][numberOfCols+1]);
                            n = Integer.parseInt(geneExp[j][numberOfCols+2]);
                            t = ((Double.parseDouble(geneExp[j][k]) - mean) / std) * Math.sqrt(n);
                            df = n - 1;
                            sampleTogeneMatrix[i][3] = String.valueOf((float)calcPval(t, df, 2)*numOfMutations);
                        }
                    }
                }
            }
        }
        
        dbConnect.exportDataToFile("SELECT distinct MutationID,Gene,SampleID,'0','0' FROM dark_matter_temp.Final_results WHERE Gene!='NONE' and DHS=1", "mutationIDToSampleToGene.bed", outpuFolder);
        String mutationIDmappedGeneExp[][] = DataSetReader.readDataSet(outpuFolder+"mutationIDToSampleToGene.bed", 5, "\t");
        System.out.println("Mapping to mutation");
        for (int i = 0; i < mutationIDmappedGeneExp.length; i++) {
            for (int j = 0; j < sampleTogeneMatrix.length; j++) {
                if(mutationIDmappedGeneExp[i][1].equals(sampleTogeneMatrix[j][0])&&
                        mutationIDmappedGeneExp[i][2].equals(sampleTogeneMatrix[j][1])){
                    mutationIDmappedGeneExp[i][3]= sampleTogeneMatrix[j][2];
                    mutationIDmappedGeneExp[i][4]= sampleTogeneMatrix[j][3];
                }
            }
        }
    
        DataSetWriter.writeToFile2D(outpuFolder+"mutationsMappedToGeneExp.bed", mutationIDmappedGeneExp);

        dbConnect.importDataFromFile("mutations_to_geneExp", outpuFolder+"mutationsMappedToGeneExp.bed");
 

    }
    //        String cols = "";
//        for (int i = 1; i < numberOfCols+1; i++) {
//            if(i==numberOfCols){
//                cols=cols+"$"+i;
//            }
//            else{
//              cols=cols+"$"+i+",";
//            }  
//        }
//           
//
//        UtilitiesLinux.shellCommandExecuter("awk 'BEGIN { OFS = \"\\t\" } NR>0 && NF=="+(numberOfCols+1)+" { if($"+(numberOfCols+1)+"=1){print "+cols+";}}' "+outpuFolder+"GeneExpression-edited.bed"+" > "+ outpuFolder+"GeneExpression-filtered.bed");

//        geneExpFile = outpuFolder+"GeneExpression-filtered.bed";
    
    public void getGeneExpMatrix() {
        String genes[][] = DataSetReader.readDataSet(sampleToGeneFile, 4, "\t");
        String matrix[][] = DataSetReader.readDataSet(geneExpFile, numberOfCols, "\t");
        
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix[0].length; j++) {
                matrix[i][j] = "0";
            }
        }
        
        
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix[0].length; j++) {
                for (int k = 1; k < genes.length; k++) {
                    if (matrix[i][0].equalsIgnoreCase(genes[k][0]) && matrix[0][j].equalsIgnoreCase(genes[k][1])) {
                        matrix[i][j] = "1";
                    }
                }
            }         
        }

        DataSetWriter.writeToFile2D(outpuFolder +"GeneExpression-Matrix.bed", matrix);
        geneExpMatrix=outpuFolder +"GeneExpression-Matrix.bed";
        
    }

    public void getGeneExp() {
        String matrix[][] = DataSetReader.readDataSet(this.geneExpMatrix, this.numberOfCols + 1, "\t");
        String data[][] = DataSetReader.readDataSet(this.geneExpFile, this.numberOfCols + 1, "\t");
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < data.length; j++) {
                if (matrix[i][2].equalsIgnoreCase(data[j][2])) {
                    matrix[i] = data[j];
                    break;
                }
            }
        }
        DataSetWriter.writeToFile2D("GeneExpression-relevant_genes.bed", matrix);
    }

    public static void getFoldChange() {        
        String geneExpMat[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Non-Dnase1 Regions\\GeneExpression-Matrix.bed", 20, "\t");
        String geneExp[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-QuantileNormalized-duplicates-removed.bed", 21, "\t");
        String sampleTogene[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Non-Dnase1 Regions\\Non-Dnase1-mutations-patient-gene.bed", 4, "\t");
        for (int j = 1; j < geneExp.length; j++) {
            ArrayList<Double> expression = new ArrayList<Double>();
            for (int i = 1; i < geneExpMat.length; i++) {
                if (geneExpMat[i][2].equalsIgnoreCase(geneExp[j][2])) {
                    for (int k = 3; k < geneExpMat[1].length; k++) {
                        if (geneExpMat[i][k].equals("0")) {
                            expression.add(Double.parseDouble(geneExp[j][k]));
                        }
                    }
                    Collections.sort(expression);
                    if (!expression.isEmpty()) {
                        if ((expression.size() % 2) == 1) {
                            geneExp[j][20] = String.valueOf(expression.get(((expression.size() + 1) / 2) - 1).doubleValue());
                        } else if ((expression.size() % 2) == 0) {
                            geneExp[j][20] = String.valueOf((expression.get((expression.size() / 2) - 1).doubleValue() + expression.get(((expression.size() / 2) + 1) - 1).doubleValue()) / 2);
                        }
                    }

                }
            }
        }

        for (int i = 1; i < sampleTogene.length; i++) {
            for (int j = 1; j < geneExp.length - 1; j++) {
                if (sampleTogene[i][2].equalsIgnoreCase(geneExp[j][2])) {
                    for (int k = 3; k < geneExp[1].length; k++) {
                        if (sampleTogene[i][1].equalsIgnoreCase(geneExp[0][k])) {
                            System.out.println("true");
                            sampleTogene[i][3] = String.valueOf(Double.parseDouble(geneExp[j][k]) / Double.parseDouble(geneExp[j][20]));
                        }
                    }
                }
            }
        }
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Breast Cancer Data\\Non-Dnase1 Regions\\genesWithFoldChange.bed", sampleTogene);

    }

    public static void tTest() {
        String geneExpMat[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-matrix-non-empty.bed", 20, "\t");
        String geneExp[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-relevant_genes.bed", 24, "\t");
        String sampleTogene[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\mutations-patient-gene.bed", 4, "\t");
//        String matrix[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\test-matrix.bed", 20, "\t");
//        String geneExpression[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\test-geneExpression.bed", 24, "\t");
//        String genes[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\mutations-patient-gene.bed", 5, "\t");
        double total = 0;
        double mean = 0;
        double std = 0;
        int n = 0;
        for (int j = 1; j < geneExp.length; j++) {
            ArrayList<Double> expression = new ArrayList<Double>();
            std = 0;
            mean = 0;
            total = 0;
            n = 0;
            for (int i = 1; i < geneExpMat.length; i++) {
                if (geneExpMat[i][2].equalsIgnoreCase(geneExp[j][2])) {
                    for (int k = 3; k < geneExpMat[1].length; k++) {
                        if (geneExpMat[i][k].equals("0")) {
                            expression.add(Double.parseDouble(geneExp[j][k]));
                        }
                    }
                    //Mean
                    for (int k = 0; k < expression.size(); k++) {
                        total += expression.get(k);
                    }
                    mean = total / expression.size();
                    System.out.println("Mean" + mean);
                    geneExp[j][20] = String.valueOf(mean);

                    //Standard Deviation
                    for (int k = 0; k < expression.size(); k++) {
                        std += (expression.get(k) - mean) * (expression.get(k) - mean);
                    }
                    geneExp[j][21] = String.valueOf(Math.sqrt(std / (expression.size() - 1)));
                    System.out.println("" + geneExp[j][21]);
                    //Number of samples
                    System.out.println("" + (expression.size()));
                    geneExp[j][22] = String.valueOf(expression.size());
                }
            }
        }
//        genes[0][3]="t";
//        genes[0][4]="df";
//        for (int i = 1; i < genes.length; i++) {
//            for (int j = 1; j < geneExpression.length; j++) {               
//               if(genes[i][2].equalsIgnoreCase(geneExpression[j][2])){
//                   for (int k = 3; k < geneExpression[1].length; k++) {
//                       if(genes[i][1].equalsIgnoreCase(geneExpression[0][k])){
//                           mean=Double.parseDouble(geneExpression[j][20]);
//                           std=Double.parseDouble(geneExpression[j][21]);
//                           n=Integer.parseInt(geneExpression[j][22]);
//                           genes[i][3]=String.valueOf(((Double.parseDouble(geneExpression[j][k])-mean)/std)*Math.sqrt(n));
//                           System.out.println(""+genes[i][3]);
//                           genes[i][4]=String.valueOf(n-1);
//                       }
//                   }                
//                } 
//            }         
//        }          
        sampleTogene[0][3] = "p";
        for (int i = 1; i < sampleTogene.length; i++) {
            double t = 0;
            int df = 0;
            for (int j = 1; j < geneExp.length; j++) {
                if (sampleTogene[i][2].equalsIgnoreCase(geneExp[j][2])) {
                    for (int k = 3; k < geneExp[1].length; k++) {
                        if (sampleTogene[i][1].equalsIgnoreCase(geneExp[0][k])) {
                            mean = Double.parseDouble(geneExp[j][20]);
                            std = Double.parseDouble(geneExp[j][21]);
                            n = Integer.parseInt(geneExp[j][22]);
                            t = ((Double.parseDouble(geneExp[j][k]) - mean) / std) * Math.sqrt(n);
                            df = n - 1;
                            sampleTogene[i][3] = String.valueOf(calcPval(t, df, 2));
                        }
                    }
                }
            }
        }
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\t-testforgenes-sqrt-test.bed", sampleTogene);

    }

    public static double calcPval(double tva, int df, int side) {
        double ttp = TtoP(tva, df);
        if (side == 1) {
            return ttp / 2;
        } else {
            return ttp;
        }
    }

    public static double TtoP(double t, int df) {
        double tsq = t * t;
        double p = 0;
        double abst = Math.abs(t);
        if (df == 1) {
            p = 1 - 2 * Math.atan(abst) / Math.PI;
        } else if (df == 2) {
            p = 1 - abst / Math.sqrt(tsq + 2);
        } else if (df == 3) {
            p = 1 - 2 * (Math.atan(abst / Math.sqrt(3)) + abst * Math.sqrt(3) / (tsq + 3)) / Math.PI;
        } else if (df == 4) {
            p = 1 - abst * (1 + 2 / (tsq + 4)) / Math.sqrt(tsq + 4);
        } else {
            double z = TtoZ(abst, df);
            if (df > 4) {
                p = Norm_p(z);
            } else {
                p = Norm_p(z);
            }
        }
        return p;
    }

    public static double TtoZ(double t, int df) {
        double A9 = df - 0.5;
        double B9 = 48 * A9 * A9;
        double T9 = t * t / df, Z8, P7, B7, z;

        if (T9 >= 0.04) {
            Z8 = A9 * Math.log(1 + T9);
        } else {
            Z8 = A9 * (((1 - T9 * 0.75) * T9 / 3 - 0.5) * T9 + 1) * T9;
        }
        P7 = ((0.4 * Z8 + 3.3) * Z8 + 24) * Z8 + 85.5;
        B7 = 0.8 * Math.pow(Z8, 2) + 100 + B9;
        z = (1 + (-P7 / B7 + Z8 + 3) / B9) * Math.sqrt(Z8);
        return z;
    }

    public static double Norm_p(double z) {
        double absz = Math.abs(z);
        double a1 = 0.0000053830;
        double a2 = 0.0000488906;
        double a3 = 0.0000380036;
        double a4 = 0.0032776263;
        double a5 = 0.0211410061;
        double a6 = 0.0498673470;
        double p = (((((a1 * absz + a2) * absz + a3) * absz + a4) * absz + a5) * absz + a6) * absz + 1;
        p = Math.pow(p, -16);
        return p;
    }

    public static void getBaseline() {
        String genExpression[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-QuantileNormalized-duplicates-removed.bed", 20, "\t");
        String matrix[][] = new String[genExpression.length][genExpression[1].length];
        for (int i = 0; i < matrix.length; i++) {
            if (i == 0) {
                for (int j = 0; j < matrix[1].length; j++) {
                    matrix[i][j] = genExpression[i][j];
                }
            } else {
                matrix[i][0] = genExpression[i][0];
                matrix[i][1] = genExpression[i][1];
                matrix[i][2] = genExpression[i][2];
            }
        }
        double value = 0;
        double median = 0;
        double foldChange = 0;
        int count = 0;

        for (int i = 1; i < genExpression.length; i++) {
            //genExpression[1].length
            for (int j = 3; j < genExpression[1].length; j++) {
                value = Double.parseDouble(genExpression[i][j]);
                ArrayList<Double> expression = new ArrayList<Double>();
                for (int k = 3; k < genExpression[1].length; k++) {
                    if (k != j) {
                        expression.add(Double.parseDouble(genExpression[i][k]));
                    }
                }
                Collections.sort(expression);
                if (expression.size() != 16) {
                    System.out.println("Error");
                }
//                for (int k = 0; k < 16; k++) {
////                    System.out.println(""+expression.get(k));
//                }

                median = (expression.get(7) + expression.get(8)) / 2;
//                System.out.println("meadian" + median);
                foldChange = value / median;
                if (foldChange >= 2 || foldChange <= -2) {
                    count++;
                }
                matrix[i][j] = String.valueOf(foldChange);
            }
        }
        System.out.println("count" + count);
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\Baseline.bed", matrix);

    }

}

//     //GeneExpression-relevant_genes.bed
//        //GeneExpression-matrix-non-empty.bed
////        String originalMatrix[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-Matrix-for-bootstrap.bed", 21, "\t");
////        String geneExpression[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-QuantileNormalized-duplicates-removed.bed", 21, "\t");
//        String originalMatrix[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-matrix-non-empty.bed", 21, "\t");
//        String geneExpression[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-relevant_genes.bed", 21, "\t");
//        String genes[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\mutations-patient-gene.bed", 103, "\t");
//        
//        String bootstrap[][]=new String[1][100];
//        String matrix[][] = new String[originalMatrix.length][originalMatrix[0].length];
//        int count=0;
//        int diff=0;
//        ArrayList <Integer>list = new ArrayList<Integer>();
//        for (int i = 1; i < originalMatrix.length; i++) {
//            list.add(i);
//        }
//        for (int a = 0; a < 1; a++) {
//            System.out.println(""+a);
////            Collections.shuffle(list);            
////            System.out.println("List.size"+list.size());
////            System.out.println("matrix.length"+matrix.length);
//            for (int j = 0; j <list.size() ; j++) {
//                matrix[j+1]=originalMatrix[list.get(j)];
//            }
//            matrix[0]=originalMatrix[0];
//
//            for (int i = 1; i < matrix.length; i++) {
//                for (int j = 1; j < geneExpression.length; j++) {
//                    if(matrix[i][2].equalsIgnoreCase(geneExpression[j][2])){
//                        for (int k = 3; k < matrix[1].length-1; k++) {
//                            if(matrix[i][k].equals("1")){
//                                matrix[i][k]="#N/A";
//                            }
//                            else{
//                                matrix[i][k]=geneExpression[j][k];
//                            }
//                        }
//                        break;
//                    }
//                }
//            }
//            String geneExpression1[][]=matrix;
//            ArrayList<Double> expression= new ArrayList<Double>();        
//            double foldChange=0;        
//            for (int i = 1; i < geneExpression1.length; i++) {
//                expression.clear();
//                for (int j = 3; j < geneExpression1[1].length-1; j++) {
//                    
//                    if(geneExpression1[i][j]!="#N/A"){
//                        try {
//                            expression.add(Double.parseDouble(geneExpression1[i][j]));
//                        } catch (Exception e) {
//                            System.out.println(""+geneExpression1[i][j]);
//                        }                    
//                    }
//                }
//                if(expression.size()>2){
//                    Collections.sort(expression);                
//                    if((expression.size()%2)==1){
//                        geneExpression1[i][20]=String.valueOf(expression.get(((expression.size()+1)/2)-1).doubleValue());
//                    }
//                    else if ((expression.size()%2)==0){                    
//                        geneExpression1[i][20]=String.valueOf((expression.get((expression.size()/2)-1).doubleValue()+expression.get(((expression.size()/2)+1)-1).doubleValue())/2);
//                    }  
//                }
//
//                else if(expression.size()==2){
//                    geneExpression1[i][20]=String.valueOf((expression.get(0).doubleValue()+expression.get(1).doubleValue())/2); 
//                }
//                else{
//                    geneExpression1[i][20]=String.valueOf((expression.get(0).doubleValue()));
//                }
//
//            }
//            //foldchange
//            count=0;
//            diff=0;
//            for (int i = 1; i < matrix.length; i++) {
//                for (int j = 1; j < geneExpression.length; j++) {
//                    if(matrix[i][2].equalsIgnoreCase(geneExpression[j][2])){
//                        for (int k = 3; k < matrix[1].length-1; k++) {
//                            if(matrix[i][k].equals("#N/A")){
////                                matrix[i][k]="#N/A";
////                            }
////                            else{
//                                count++;
////                            if(matrix[i][k].equals("1")){
//                                foldChange=Double.parseDouble(geneExpression[j][k])/Double.parseDouble(geneExpression1[j][20]);
////                                System.out.println(""+foldChange);
////                                matrix[i][k]=geneExpression[j][k];
//                                matrix[i][k]=String.valueOf(foldChange);
//                                if(foldChange>=2||foldChange<=-2){
//                                    diff++;
//                                    System.out.println("true" + diff);
//                                }
//                            }
//                        }
//                        matrix[i][20]=geneExpression1[j][20];                        
//                        break;
//                    }
//                }
//            }
//            System.out.println("Count = " + count);
//            System.out.println("Diff = " + diff);
//            bootstrap[0][a]=String.valueOf(diff);
////            for (int k = 3; k < matrix[1].length-1; k++) {
////            for (int i = 1; i < genes.length ; i++) {
////                if(matrix[0][k].equalsIgnoreCase(genes[i][1])){
////                    for (int j = 1; j < matrix.length; j++) {
////                        if(genes[i][2].equalsIgnoreCase(matrix[j][2])){ 
//////                            System.out.println(""+matrix[j][k]);
//////                            System.out.println(""+matrix[j][20]);
////                                foldChange=Math.log(Double.parseDouble(matrix[j][k])/Double.parseDouble(matrix[j][20]))/Math.log(2.0); 
////                                genes[i][a+3]=String.valueOf(foldChange);
////                                if(foldChange>=2||foldChange<=-2){
////                                    count++;
////                                }
////                                
////                        }                        
////                    }                    
////                }            
////            }   
////        }
////         bootstrap[0][a]=String.valueOf(count);                          
//        }
//        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\bootstrapformutations100.bed", matrix);
//        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\bootstrapforcount100.bed", bootstrap);
//    
//    public static void getMedian(){
//        
////        String geneExpression[][] = DataSetReader.readDataSet("GeneExpression-QuantileNormalized.bed", 21, "\t");
////        String matrix[][] = DataSetReader.readDataSet("GeneExpression-Matrix.bed", 20, "\t");
//        String geneExpression[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-relevant_genes.bed", 21, "\t");
//        String matrix[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-matrix-non-empty.bed", 21, "\t");
//        for (int i = 1; i < matrix.length; i++) {
//            for (int j = 1; j < geneExpression.length; j++) {
//                if(matrix[i][2].equalsIgnoreCase(geneExpression[j][2])){
//                    for (int k = 3; k < matrix[1].length-1; k++) {
//                        if(matrix[i][k].equals("1")){
//                            matrix[i][k]="#N/A";
//                        }
//                        else{
//                            matrix[i][k]=geneExpression[j][k];
//                        }
//                    }
//                    break;
//                }
//            }
//        }
//        
//        geneExpression=matrix;
//       
//        ArrayList<Double> expression= new ArrayList<Double>();        
//        double foldChange=0;        
//        for (int i = 1; i < geneExpression.length; i++) {
//            expression.clear();
//            for (int j = 3; j < geneExpression[1].length-1; j++) {
//                System.out.println("" + geneExpression[i][j]);
//                if(geneExpression[i][j]!="#N/A"){
//                    try {
//                        expression.add(Double.parseDouble(geneExpression[i][j]));
//                    } catch (Exception e) {
//                        System.out.println(""+geneExpression[i][j]);
//                    }                    
//                }
//            }
//            if(expression.size()>2){
//                Collections.sort(expression);                
//                if((expression.size()%2)==1){
//                    geneExpression[i][20]=String.valueOf(expression.get(((expression.size()+1)/2)-1).doubleValue());
//                }
//                else if ((expression.size()%2)==0){                    
//                    geneExpression[i][20]=String.valueOf((expression.get((expression.size()/2)-1).doubleValue()+expression.get(((expression.size()/2)+1)-1).doubleValue())/2);
//                }  
//            }
//            
//            else if(expression.size()==2){
//                geneExpression[i][20]=String.valueOf((expression.get(0).doubleValue()+expression.get(1).doubleValue())/2); 
//            }
//            else{
//                geneExpression[i][20]=String.valueOf((expression.get(0).doubleValue()));
//            }
//                                
//        }
////        
////        for (int i = 1; i < matrix.length; i++) {
////            
////            for (int j = 3; j < matrix[0].length; j++) {
////                if(matrix[i][j].equals("1")){
////                   foldChange=Double.parseDouble(geneExpression[i][j])/Double.parseDouble(geneExpression[i][20]);
////                   matrix[i][j]=String.valueOf(foldChange);
////                }
////                else{
////                   matrix[i][j]="-"; 
////                }
////            }           
////        }
////        
//        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-with-Median-for-t-test2.bed", geneExpression);
////        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Breast Cancer Data\\GeneExpression-with-foldchange", matrix);
//        
//    }
//    
//    public static void getFoldChange(){
//        String genes[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\mutations-patient-gene.bed", 4, "\t");
////        String matrix[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\17BreastCancers-foldChange.bed", 18, "\t");
////        String matrix[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\foldchange2.bed", 18, "\t");
//        String matrix[][] = DataSetReader.readDataSet("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\GeneExpression-with-Median-for-foldchange.bed", 21, "\t");
//        int start=1;
//        int count=0;
//        double sum = 0;
//        String temp = "";
////        for (int i = 0; i < matrix.length; i++) {
////            for (int j = 0; j < matrix[0].length; j++) {
////                if(matrix[i][j].contains(",")&&!matrix[i][j].contains("-")){
////                    count=0;
////                    sum=0;
////                    while(true){
////                        if(matrix[i][j].indexOf(",")!= -1){
////                            sum += Double.parseDouble(matrix[i][j].substring(0, matrix[i][j].indexOf(",")).trim());
////                            temp = matrix[i][j]. substring(matrix[i][j].indexOf(",")+1).trim();
////                            matrix[i][j] = temp;
////                            count++;
////
////                        }
////                        else{
////                            sum+= Double.parseDouble(matrix[i][j]);
////                            count++;
////                            break;
////
////
////                        }
////                    }
////                    matrix[i][j]=String.valueOf(sum/count);
//////                                
////                
////            }
////                else if(matrix[i][j].contains(",")&&matrix[i][j].contains("-")){
////                   matrix[i][j]="-"; 
////                }
////        }
////            
////        }
//       
////                                count=0;
////                                sum=0;
////                                while(true){
////                                    if(matrix[j][a].indexOf(",")!= -1){
////                                        sum += Double.parseDouble(matrix[j][a].substring(0, matrix[j][a].indexOf(",")).trim());
////                                        temp = matrix[j][a]. substring(matrix[j][a].indexOf(",")+1).trim();
////                                        matrix[j][a] = temp;
////                                        count++;
////
////                                    }
////                                    else{
////                                        sum+= Double.parseDouble(matrix[j][a]);
////                                        count++;
////                                        break;
////                                        
////                                        
////                                    }
////                                }
////                                genes[i][3]=String.valueOf(sum/count);
////                                
////                                
////                            }
//        
//        for (int a = 3; a < matrix[1].length; a++) {
//            for (int i = 1; i < genes.length ; i++) {
//                if(matrix[0][a].equalsIgnoreCase(genes[i][1])){
//                    for (int j = 1; j < matrix.length; j++) {
//                        if(genes[i][2].equalsIgnoreCase(matrix[j][2])){                            
//                                genes[i][3]=String.valueOf(Double.parseDouble(matrix[j][a])/Double.parseDouble(matrix[j][20]));                        
//                        }                        
//                    }                    
//                }            
//            }   
//        }
//        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Breast Cancer Data\\Iteration2\\foldChange-for-mutations-new.bed", genes);
//    }

