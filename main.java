import java.io.*;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

public class main {
    public static void main(String[] args) throws IOException {
        String result = "";
        String inputFile = args[0];
        String inputType = args[1];
        if (inputType.equals("1")) {
            File file = new File(inputFile);
            BufferedReader input = new BufferedReader(new FileReader(file));
            String firstLine = input.readLine();
            String[] firstLineArray = firstLine.trim().split("\\s+");
            int k = Integer.parseInt(firstLineArray[0]);
            int t = Integer.parseInt(firstLineArray[1]);
            String[] dnaArray = new String[t];
            for (int i = 0; i < t; i++) {
                dnaArray[i] = input.readLine();
            }
            result = randomizedMotifSearch(k, t, dnaArray);

            String outputFileName = inputFile.replaceAll("[^\\d.]", "");
            outputFileName = "sol_" + outputFileName;
            File outputFile = new File(outputFileName);
            FileWriter writer = new FileWriter(outputFileName);
            writer.write(result);
            writer.close();
        }
        else if(inputType.equals("2")) {
            File file = new File(inputFile);
            BufferedReader input = new BufferedReader(new FileReader(file));
            String firstLine = input.readLine();
            String[] firstLineArray = firstLine.trim().split("\\s+");
            int k = Integer.parseInt(firstLineArray[0]);
            int t = Integer.parseInt(firstLineArray[1]);
            int r = Integer.parseInt(firstLineArray[2]);
            String[] dnaArray = new String[t];
            for (int i = 0; i < t; i++) {
                dnaArray[i] = input.readLine();
            }
            result = gibbsSampling(k, t, r, dnaArray);

            String outputFileName = inputFile.replaceAll("[^\\d.]", "");
            outputFileName = "gibbis_sol_" + outputFileName;
            File outputFile = new File(outputFileName);
            FileWriter writer = new FileWriter(outputFileName);
            writer.write(result);
            writer.close();
        }
        else {
            System.out.println("Incorrect option, program ended");
        }

        String outputFileName = inputFile.replaceAll("[^\\d.]", "");
        outputFileName = "sol_" + outputFileName;
        File outputFile = new File(outputFileName);
        FileWriter writer = new FileWriter(outputFileName);
        writer.write(result);
        writer.close();

    }

    public static String randomizedMotifSearch(int k, int t, String[]dnaArray) {
//        System.out.println(k);
//        System.out.println(t);
//        for (int i = 0; i < dnaArray.length; i++) {
//            System.out.println(dnaArray[i]);
//        }

        //code:
        String finalConsensus = "";
        Random rand = new Random();
        String[] bestMotifs = new String[t];
        float[] bestMotifsValue = new float[t];
        float bestScore = -1;

        for (int i = 0; i < 1500; i++) {
            String[] motifTable = new String[t];
            for (int j = 0; j < t; j++) {
                int rand_start = rand.nextInt(dnaArray[0].length() - k + 1);
                motifTable[j] = dnaArray[j].substring(rand_start, rand_start + k);
            }
            int[] A = new int[k];
            int[] C = new int[k];
            int[] T = new int[k];
            int[] G = new int[k];
            Arrays.fill(A, 1);
            Arrays.fill(C, 1);
            Arrays.fill(T, 1);
            Arrays.fill(G, 1);
            for (int j = 0; j < t; j++) {
                for (int z = 0; z < k; z++) {
                    if (motifTable[j].charAt(z) == 'A') {
                        A[z]++;
                    }
                    else if (motifTable[j].charAt(z) == 'C') {
                        C[z]++;
                    }
                    else if (motifTable[j].charAt(z) == 'T') {
                        T[z]++;
                    }
                    else if (motifTable[j].charAt(z) == 'G') {
                        G[z]++;
                    }
                }
            }
            String[] currentBestMotifs = new String[t];
            float[] currentBestMotifsValue = new float[t];
            Arrays.fill(currentBestMotifsValue, 0);

            for (int j = 0; j < t; j++) {
                for (int z = 0; z < dnaArray[0].length() - k + 1; z++) {
                    float currentValue = motiveGiver(dnaArray[j].substring(z, z + k), A, C, T, G, t);
                    if (currentBestMotifsValue[j] < currentValue) {
                        currentBestMotifsValue[j] = currentValue;
                        currentBestMotifs[j] = dnaArray[j].substring(z, z + k);
                    }
                }
            }
            int currentScore = 0;
            String CurrentConsensus = "";
            for (int z = 0; z < k; z++) {
                int aValue = 0;
                int cValue = 0;
                int tValue = 0;
                int gValue = 0;
                for (int j = 0; j < t; j++) {
                    if (currentBestMotifs[j].charAt(z) == 'A') {
                        aValue++;
                    }
                    else if (currentBestMotifs[j].charAt(z) == 'C') {
                        cValue++;
                    }
                    else if (currentBestMotifs[j].charAt(z) == 'T') {
                        tValue++;
                    }
                    else if (currentBestMotifs[j].charAt(z) == 'G') {
                        gValue++;
                    }
                }
                int maxScore = Math.max(Math.max(tValue, gValue), Math.max(aValue, cValue));
                currentScore += t - maxScore;

                if (maxScore == tValue) {
                    CurrentConsensus += "T";
                }
                else if (maxScore == gValue) {
                    CurrentConsensus += "G";
                }
                else if (maxScore == aValue) {
                    CurrentConsensus += "A";
                }
                else if (maxScore == cValue) {
                    CurrentConsensus += "C";
                }
            }
            if (bestScore == -1 || currentScore < bestScore) {
                bestScore = currentScore;
                bestMotifs = currentBestMotifs;
                finalConsensus = CurrentConsensus;
            }
//            for (int j = 0; j < t; j++) {
//                System.out.println(motifTable[j]);
//            }
        }
        String output = "";
        for (int j = 0; j < t; j++) {
            System.out.println(bestMotifs[j]);
            output += (bestMotifs[j]);
            output += "\n";
        }
        System.out.println("Score: " + bestScore);
        System.out.println("Consensus sequence : " + finalConsensus);
        return output;
    }

    public static float motiveGiver(String text, int[] A, int[] C, int[] T, int[] G, int t) {

//        for (int i = 0; i < text.length(); i++) {
//            System.out.print(A[i]);
//        }
//        System.out.println();
//        for (int i = 0; i < text.length(); i++) {
//            System.out.print(C[i]);
//        }
//        System.out.println();
//        for (int i = 0; i < text.length(); i++) {
//            System.out.print(T[i]);
//        }
//        System.out.println();
//        for (int i = 0; i < text.length(); i++) {
//            System.out.print(G[i]);
//        }
//        System.out.println();
        //System.out.println(text);
        float value = 1;
        for(int i = 0; i < text.length(); i++) {
            if (text.charAt(i) == 'A') {
                value *= ((float)A[i] / (t  * 2));
            }
            else if (text.charAt(i) == 'C') {
                value *= ((float)C[i] / (t  * 2));
            }
            else if (text.charAt(i) == 'T') {
                value *= ((float)G[i] / (t  * 2));
            }
            else if (text.charAt(i) == 'G') {
                value *= ((float)T[i] / (t  * 2));
            }
        }
        return value;
    }

    public static String gibbsSampling(int k, int t, int r, String[]dnaArray) {
        // System.out.println(k);
        // System.out.println(t);
        // System.out.println(r);
//        for (int i = 0; i < dnaArray.length; i++) {
//            System.out.println(dnaArray[i]);
//        }
        //code:
        String finalConsensus = "";

        String[] bestMotifs = new String[t];
        float bestScore = -1;
        for (int i = 0; i < 30; i++) {
            String[] motifTable = new String[t];
            Random rand = new Random();
            for (int j = 0; j < t; j++) {
                int rand_start = rand.nextInt(dnaArray[0].length() - k);
                motifTable[j] = dnaArray[j].substring(rand_start, rand_start + k);
            }

            for (int rCount = 0; rCount < r; rCount++) {
                Random rand2 = new Random();
                int rand_line = rand2.nextInt(t);
                // System.out.println("Randline: " + rand_line);

                int[] A = new int[k];
                int[] C = new int[k];
                int[] T = new int[k];
                int[] G = new int[k];
                Arrays.fill(A, 1);
                Arrays.fill(C, 1);
                Arrays.fill(T, 1);
                Arrays.fill(G, 1);
                for (int j = 0; j < t; j++) {
                    if (j != rand_line) {
                        for (int z = 0; z < k; z++) {
                            if (motifTable[j].charAt(z) == 'A') {
                                A[z]++;
                            } else if (motifTable[j].charAt(z) == 'C') {
                                C[z]++;
                            } else if (motifTable[j].charAt(z) == 'T') {
                                T[z]++;
                            } else if (motifTable[j].charAt(z) == 'G') {
                                G[z]++;
                            }
                        }
                    }
                }

                double[] allValues = new double[dnaArray[0].length() - k];
                String[] allOption = new String[dnaArray[0].length() - k];
                long totalValue = 0;
                for (int q = 0; q < allValues.length; q++) {
                    int currentValue = gibbsMotiveGiver(dnaArray[rand_line].substring(q, q + k), A, C, T, G, t - 1);
                    allValues[q] = currentValue;
                    allOption[q] = dnaArray[rand_line].substring(q, q + k);
                    totalValue += currentValue;
                    //System.out.println(currentValue);
                }

                // System.out.print(totalValue)

                for (int q = 0; q < allValues.length; q++) {
                    allValues[q] = allValues[q] / totalValue;
                    // System.out.println(allValues[q]);
                }

                Random rand3 = new Random();
                double randomPickWithWeight = rand3.nextDouble();
                double randValue = 0.0;
                int counter = 0;
                //System.out.println(randomPickWithWeight);
//                for (int f = 0; f < allValues.length; f++) {
//                    System.out.print(allValues[f] + ", ");
//                }
                while (randValue  < randomPickWithWeight) {
                    randValue += allValues[counter];
                    counter++;
                }
                //System.out.println("counter: " + (counter - 1));
                String newLine = allOption[counter - 1];
                //System.out.println(newLine);


                int oldScore = ScoreGiver(motifTable, k, t);
                String[] motifTableNew = new String[motifTable.length];
                for(int a = 0; a < motifTable.length; a++){
                    motifTableNew[a] = motifTable[a];
                }
                motifTableNew[rand_line] = newLine;
                int newScore = ScoreGiver(motifTableNew, k, t);

                if (newScore < oldScore) {
                    motifTable = motifTableNew;
                }
            }

            int currentScore = 0;
            String CurrentConsensus = "";
            for (int z = 0; z < k; z++) {
                int aValue = 0;
                int cValue = 0;
                int tValue = 0;
                int gValue = 0;
                for (int j = 0; j < t; j++) {
                    if (motifTable[j].charAt(z) == 'A') {
                        aValue++;
                    }
                    else if (motifTable[j].charAt(z) == 'C') {
                        cValue++;
                    }
                    else if (motifTable[j].charAt(z) == 'T') {
                        tValue++;
                    }
                    else if (motifTable[j].charAt(z) == 'G') {
                        gValue++;
                    }
                }
                int maxScore = Math.max(Math.max(tValue, gValue), Math.max(aValue, cValue));
                currentScore += t - maxScore;

                if (maxScore == tValue) {
                    CurrentConsensus += "T";
                }
                else if (maxScore == gValue) {
                    CurrentConsensus += "G";
                }
                else if (maxScore == aValue) {
                    CurrentConsensus += "A";
                }
                else if (maxScore == cValue) {
                    CurrentConsensus += "C";
                }
            }
            //System.out.println(currentScore);
            if (bestScore == -1 || currentScore < bestScore) {
                bestScore = currentScore;
                bestMotifs = motifTable;
                finalConsensus = CurrentConsensus;
            }
        }

        String output = "";
        for (int j = 0; j < t; j++) {
            System.out.println(bestMotifs[j]);
            output += bestMotifs[j];
            output += "\n";
        }
        System.out.println("Score: " + bestScore);
        System.out.println("Consensus sequence : " + finalConsensus);

        return output;
    }

    public static int ScoreGiver(String[] motifTable, int k, int t) {
        int currentScore = 0;
        for (int z = 0; z < k; z++) {
            int aValue = 0;
            int cValue = 0;
            int tValue = 0;
            int gValue = 0;
            for (int j = 0; j < t; j++) {
                if (motifTable[j].charAt(z) == 'A') {
                    aValue++;
                }
                else if (motifTable[j].charAt(z) == 'C') {
                    cValue++;
                }
                else if (motifTable[j].charAt(z) == 'T') {
                    tValue++;
                }
                else if (motifTable[j].charAt(z) == 'G') {
                    gValue++;
                }
            }
            int maxScore = Math.max(Math.max(tValue, gValue), Math.max(aValue, cValue));
            currentScore += t - maxScore;
        }
        return currentScore;
    }
    public static int gibbsMotiveGiver(String text, int[] A, int[] C, int[] T, int[] G, int t) {

//        for (int i = 0; i < text.length(); i++) {
//            System.out.print(A[i]);
//        }
//        System.out.println();
//        for (int i = 0; i < text.length(); i++) {
//            System.out.print(C[i]);
//        }
//        System.out.println();
//        for (int i = 0; i < text.length(); i++) {
//            System.out.print(T[i]);
//        }
//        System.out.println();
//        for (int i = 0; i < text.length(); i++) {
//            System.out.print(G[i]);
//        }
//        System.out.println();
        //System.out.println(text);
        int value = 1;
        for(int i = 0; i < text.length(); i++) {
            if (text.charAt(i) == 'A') {
                value *= A[i];
            }
            else if (text.charAt(i) == 'C') {
                value *= C[i];
            }
            else if (text.charAt(i) == 'T') {
                value *= T[i];
            }
            else if (text.charAt(i) == 'G') {
                value *= G[i];
            }
        }
        return value;
    }
}
