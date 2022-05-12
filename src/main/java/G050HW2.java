import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;


public class G050HW2 {


    private static double r_min = 0.0;
    private static double last_r = 0.0;
    private static int n_iters = 0;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Input reading methods
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    public static Vector strToVector(String str) {
        String[] tokens = str.split(",");
        double[] data = new double[tokens.length];
        for (int i = 0; i < tokens.length; ++i) {
            data[i] = Double.parseDouble(tokens[i]);
        }

        return Vectors.dense(data);
    }

    public static ArrayList<Vector> readVectorsSeq(String filename) throws IOException {
        if (Files.isDirectory(Paths.get(filename))) {
            throw new IllegalArgumentException("readVectorsSeq is meant to read a single file.");
        }

        ArrayList<Vector> result = new ArrayList<>();
        Files.lines(Paths.get(filename))
                .map(str -> strToVector(str))
                .forEach(e -> result.add(e));

        return result;
    }


    /**
     * Implements the weighted variant of kcenterOUT
     * @param P set of points
     * @param W set of weights
     * @param k number of centers
     * @param z number of allowed outliers
     * @param alpha coefficient used in ball generation
     * @return  set of centers
     */
    public static ArrayList<Vector> SeqWeightedOutliers(ArrayList<Vector> P,
                                                        ArrayList<Long> W,
                                                        int k,
                                                        int z,
                                                        double alpha) {
        // compute (min distance between first k + z + 1 points of P)/2

        // kCenterOut

    }

    /**
     * Computes the value of the objective function
     * @param P set of points
     * @param S set of centers
     * @param z outliers
     * @return value of the objective function
     */
    public static double ComputeObjective(ArrayList<Vector> P, ArrayList<Vector> S, int z) {

    }

    public static void main(String[] args) throws IOException {

        if (args.length != 3) {
            throw new IllegalArgumentException("USAGE: file path, number of centers, number of allowed outliers");
        }

        final String path = args[0];
        final int k = Integer.parseInt(args[1]);
        final int z = Integer.parseInt(args[2]);

        ArrayList<Vector> inputPoints = readVectorsSeq(path);
        ArrayList<Long> weights = new ArrayList<>(Collections.nCopies(inputPoints.size(), 0L));

        if (weights.size() != inputPoints.size() || weights.size() ==0)  throw new AssertionError("ERROR: len(W) != len(P)");

        final long startTime = System.currentTimeMillis();
        ArrayList <Vector> solution = SeqWeightedOutliers(inputPoints, weights, k, z, 0);
        final long timeElapsed = System.currentTimeMillis() - startTime;

        double objective = ComputeObjective(inputPoints, solution, z);


        System.out.println("Input size n = " + inputPoints.size());
        System.out.println("Number of centers k = " + k);
        System.out.println("Number of outliers z = " + z);
        System.out.println("Initial guess = " + r_min);
        System.out.println("Final guess = " + last_r);
        System.out.println("Number of guesses = " + n_iters);
        System.out.println("Objective function = " + objective);
        System.out.println("Time of SeqWeightedOutliers = " + timeElapsed);
    }
}