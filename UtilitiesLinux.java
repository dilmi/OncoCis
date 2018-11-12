/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package OncoCis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

/**
 *
 * @author Dilmi
 */
public class UtilitiesLinux {

    static String outputFolder = OncoCis.MAIN_FOLDER+"temp-data/";

    public static void bedtoolsIntersect(String fileA, String fileB, String outputFile) {
        fromDosCheck(fileA);
        fromDosCheck(fileB);
        Writer output = null;
        Process prc = null;
        String cmd = "intersectBed -a " + fileA + " -b " + fileB + " > " + outputFile;
        try {
            output = new BufferedWriter(new FileWriter(outputFolder + "tempbedtoolsIntersect.sh"));
            output.write(cmd);
            output.close();
            prc = Runtime.getRuntime().exec("chmod 777 " + outputFolder + "tempbedtoolsIntersect.sh");
            prc.waitFor();
            prc.destroy();
            prc = Runtime.getRuntime().exec(outputFolder + "tempbedtoolsIntersect.sh");
            prc.waitFor();
//            Runtime.getRuntime().exec(OncoCis.OUTPUT_FOLDER+"temp.sh");

//            output = new BufferedWriter(new FileWriter(OUTPUTFOLDER+"/temp.sh"));
//            output.write(cmd);
//            output.close();
//            prc = Runtime.getRuntime().exec("chmod 777 " + OUTPUTFOLDER+ "temp.sh");
//            prc.waitFor();
//            Runtime.getRuntime().exec(OUTPUTFOLDER+"/temp.sh");
//            prc.destroy();
//            Runtime.getRuntime().exec(OUTPUTFOLDER+"/temp.sh");
        } catch (IOException ex) {
            System.out.println("" + ex.toString());
        } catch (InterruptedException e) {
            System.out.println("" + e.toString());
        }

//            Runtime.getRuntime().exec("chmod 777 /home/dilmi/data/Breast_Cancer/DarkMatter/temp.sh");
    }

    public static void bedtoolsIntersect(String fileA, String fileB, String outputLocation, String options) {
        fromDosCheck(fileA);
        fromDosCheck(fileB);
        Writer output = null;
        Process prc = null;
        String cmd = "intersectBed " + options + " -a " + fileA + " -b " + fileB + " > " + outputLocation;
        try {
            output = new BufferedWriter(new FileWriter(outputFolder + "tempbedtoolsIntersect1.sh"));
            output.write(cmd);
            output.close();
            prc = Runtime.getRuntime().exec("chmod 777 " + outputFolder + "tempbedtoolsIntersect1.sh");
            prc.waitFor();
            prc.destroy();
            prc = Runtime.getRuntime().exec(outputFolder + "tempbedtoolsIntersect1.sh");
            prc.waitFor();

        } catch (IOException ex) {
            System.out.println("Error" + outputLocation);
        } catch (InterruptedException e) {
            System.out.println("Error" + outputLocation);
        }

    }

    public static void bedtoolsGetFasta(String hg19, String inputFile, String outputLocation) {
        System.out.println("hg19" + hg19);
        System.out.println("Input" + inputFile);
        fromDosCheck(inputFile);
        Writer output = null;
        Process prc = null;
        String cmd = "bedtools getfasta -fi " + hg19 + " -bed " + inputFile + " -fo " + outputLocation;
        try {
            output = new BufferedWriter(new FileWriter(outputFolder + "tempbedtoolsGetFasta.sh"));
            output.write(cmd);
            output.close();
            prc = Runtime.getRuntime().exec("chmod 777 " + outputFolder + "tempbedtoolsGetFasta.sh");
            prc.waitFor();
            prc.destroy();
            prc = Runtime.getRuntime().exec(outputFolder + "tempbedtoolsGetFasta.sh");
            prc.waitFor();
//            output =
//                    new BufferedWriter(new FileWriter("/home/dilmi/data/Breast_Cancer/DarkMatter/temp.sh"));
//            output.write(cmd);
//            output.close();
//            prc = Runtime.getRuntime().exec("chmod 777 /home/dilmi/data/Breast_Cancer/DarkMatter/temp.sh");
//            prc.waitFor();
//
//            Runtime.getRuntime().exec("/home/dilmi/data/Breast_Cancer/DarkMatter/temp.sh");
//            prc.destroy();
//
//            Runtime.getRuntime().exec("/home/dilmi/data/Breast_Cancer/DarkMatter/temp.sh");

        } catch (IOException ex) {
            //Logger.getLogger(PunchGUI.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException e) {
        }

    }

    public static void fimo(String motifFile, String sequenceFile, String outputLocation) {
        fromDosCheck(sequenceFile);
        Writer output = null;
        Process prc = null;
        //fimo  --verbosity 1 --thresh 0.0001 jolma2013.meme temp.fa

        String cmd = "fimo --oc " + outputLocation + " --verbosity 1 --thresh 0.0001 " + motifFile + " " + sequenceFile;
        System.out.println("cmd : " + cmd);
        try {
            output = new BufferedWriter(new FileWriter(outputFolder + "tempfimo.sh"));
            output.write(cmd);
            output.close();
            prc = Runtime.getRuntime().exec("chmod 777 " + outputFolder + "tempfimo.sh");
            prc.waitFor();
            prc.destroy();
            prc = Runtime.getRuntime().exec(outputFolder + "tempfimo.sh");
            prc.waitFor();

        } catch (IOException ex) {
            //Logger.getLogger(PunchGUI.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException e) {
        }

    }

    public static void possum(String sequenceFile, String outputLocation) {
        fromDosCheck(sequenceFile);
        Writer output = null;
        Process prc = null;
        //fimo  --verbosity 1 --thresh 0.0001 jolma2013.meme temp.fa
        /* ./possum-src-JW-edited/possum-src/possum -o Motif-output.bed JASPAR_CORE_2009_EDITED_RED1_4Possum4.dat  NormalSequence.fasta */
        /* grep . Motif-output.bed > Motif-output-edited.bed                                                                            */
        String cmd = "cd /home/bioinfo/Darkmatter-data/DarkMatter/default/possum-src-JW-edited/possum-src/" + "\n" + "./possum -q -p 0 /home/bioinfo/Darkmatter-data/DarkMatter/default/JASPAR-CORE_2014_pfm_vertebrates_4Possum.dat  " + sequenceFile + "> " + outputLocation + "/Motif-output.bed"
                + "\n" + "grep . " + outputLocation + "/Motif-output.bed > " + outputLocation + "/Possum.bed";
        System.out.println("Running possum");
        try {
            output = new BufferedWriter(new FileWriter(outputFolder + "tempPossum.sh"));
            output.write(cmd);
            output.close();
            prc = Runtime.getRuntime().exec("chmod 777 " + outputFolder + "tempPossum.sh");
            prc.waitFor();
            prc.destroy();
            prc = Runtime.getRuntime().exec(outputFolder + "tempPossum.sh");
            prc.waitFor();

        } catch (IOException ex) {
            //Logger.getLogger(PunchGUI.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException e) {
        }

    }

    public static void conservationScore(String inputFile, String outputFile) {
        fromDosCheck(inputFile);
        Writer output = null;
        Process prc = null;
        String cmd = "cd "+OncoCis.MAIN_FOLDER+"default \n" + "while read line ; do  $line; done < " + inputFile + " > " + outputFile;
        try {
            output = new BufferedWriter(new FileWriter(outputFolder + "tempconservationScore.sh"));
            output.write(cmd);
            output.close();
            prc = Runtime.getRuntime().exec("chmod 777 " + outputFolder + "tempconservationScore.sh");
            prc.waitFor();
            prc.destroy();
            prc = Runtime.getRuntime().exec(outputFolder + "tempconservationScore.sh");
            prc.waitFor();

//            
        } catch (IOException ex) {
            //Logger.getLogger(PunchGUI.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException e) {
        }

//            Runtime.getRuntime().exec("chmod 777 /home/dilmi/data/Breast_Cancer/DarkMatter/temp.sh");
    }

    public static void shellCommandExecuter(String cmd) {
//        System.out.println("True");
        Writer output = null;
        Process prc = null;
        //fimo  --verbosity 1 --thresh 0.0001 jolma2013.meme temp.fa
        // 

        try {
//            
            output = new BufferedWriter(new FileWriter(outputFolder + "tempawk.sh"));
            output.write(cmd);
            output.close();
            prc = Runtime.getRuntime().exec("chmod 777 " + outputFolder + "tempawk.sh");
            prc.waitFor();

//            Runtime.getRuntime().exec(OncoCis.OUTPUT_FOLDER+"tempawk.sh");
            prc.destroy();
            prc = Runtime.getRuntime().exec(outputFolder + "tempawk.sh");
            prc.waitFor();

        } catch (IOException ex) {
            System.out.println("Error");
        } catch (InterruptedException e) {
            System.out.println("Error");
        }

    }

    public static void fromDosCheck(String file) {
        Writer output = null;
        Process prc = null;
        String cmd = "fromdos " + file;
        try {
            output = new BufferedWriter(new FileWriter(outputFolder + "temp.sh"));
            output.write(cmd);
            output.close();
            prc = Runtime.getRuntime().exec("chmod 777 " + outputFolder + "temp.sh");
            prc.waitFor();
            prc.destroy();
            prc = Runtime.getRuntime().exec(outputFolder + "temp.sh");
            prc.waitFor();
//            output = new BufferedWriter(new FileWriter(OUTPUTFOLDER+"/temp.sh"));
//            output.write(cmd);
//            output.close();
//            prc = Runtime.getRuntime().exec("chmod 777 " + OUTPUTFOLDER+ "temp.sh");
//            prc.waitFor();
//            Runtime.getRuntime().exec(OUTPUTFOLDER+"/temp.sh");
//            prc.destroy();
//            Runtime.getRuntime().exec(OUTPUTFOLDER+"/temp.sh");

        } catch (IOException ex) {
            //Logger.getLogger(PunchGUI.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException e) {
        }
    }

}
