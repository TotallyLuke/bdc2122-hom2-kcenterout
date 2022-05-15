import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;


public class G050HW2 {

    private static double[][] distanceMatrix;

    /**
     * Implements the weighted variant of kcenterOUT
     *
     * @param P     set of points
     * @param W     set of weights
     * @param k     number of centers
     * @param z     number of allowed outliers
     * @param alpha coefficient used in ball generation
     * @return set of centers
     */
    public static ArrayList<Vector> SeqWeightedOutliers(ArrayList<Vector> P,
                                                        ArrayList<Long> W,
                                                        int k,
                                                        int z,
                                                        double alpha) {
        // compute (min distance between first k + z + 1 points of P)/2
        // this is an initial guess for the radius of the ball
        double r_min = Double.MAX_VALUE;
        for (int i = 0; i < k + z + 1; ++i) {
            for (int j = i + 1; j < k + z + 1; ++j) {
                if (distanceMatrix[i][j] < r_min * 2) {
                    r_min = distanceMatrix[i][j] / 2;
                }
            }
        }


        double r = r_min;
        ArrayList<Vector> S = new ArrayList<>(0); // S: set of centers
        int n_iters = 1;
        // kCenterOut
        while (true) {
            List<Integer> Z_indexes = new LinkedList<>(); // Z_indexes: set of indexes of points in Uncovered set
            for (int i = 0; i < P.size(); ++i) { // at the beginning, all points are in Uncovered set
                Z_indexes.add(i);
            }
            S.clear();

            long Wz = 0; // Wz: sum of weights of points in Uncovered set
            for (double w : W) {
                Wz += w;
            }

            while (S.size() < k && Wz > 0) {
                // PREcondition: Z must be greater than 0 (true because Wz > 0)
                //TODO: why here max is not initialized as 0 (as in the pseudocode)? In this way you are also checking twice the first point: here and in the for loop
                //@ftrole newcenter if not initialized then error in line 71, either create an empty Vector or newcenter = P.get(first)
                int first = Z_indexes.get(0);
                Vector newcenter = P.get(first);        // dummy inizialization to avoid error
                double max = getBallWeight(first, r, alpha, Z_indexes, W);   // max ball weight
                int i_max = first;
                for (int i = 1; i < P.size(); ++i) {    // find max ball, assign center and value of it to newcenter and max
                                                        // finds the best center for uncovered points
                    double ball_weight = getBallWeight(i, r,  alpha, Z_indexes, W);
                    if (ball_weight > max) {
                        max = ball_weight;
                        newcenter = P.get(i);
                        i_max = i;
                    }
                }
                S.add(newcenter);

                for (int i = 0; i < Z_indexes.size(); ) { // remove now covered points from Uncovered set, update Wz
                    if (distanceMatrix[i_max][Z_indexes.get(i)] <= (3 + 4 * alpha) * r) {
                        Wz -= W.get(Z_indexes.get(i));
                        Z_indexes.remove(i);
                    } else {
                        ++i;
                    }
                }
            }
            if (Wz <= z) // I've covered all non outliers points
                break;
            // else
            r = 2 * r; // double radius for a new attempt
            ++n_iters;
        }
        System.out.println("Initial guess = " + r_min);
        System.out.println("Final guess = " + r);
        System.out.println("Number of guesses = " + n_iters);
        return S;
    }

    /**
     * Computes the weight of ball centered in c_i of radius (1 + 2*alpha) * r
     *
     * @param i         index of ongoing Point c_i
     * @param r         radius of ball
     * @param alpha     alpha
     * @param Z_indexes indexes in Uncovered set
     * @param W         weights
     * @return ball weight
     */
    private static double getBallWeight(int i, double r, double alpha, List<Integer> Z_indexes, ArrayList<Long> W) {
        double ball_weight = 0.0;
        for (int j : Z_indexes) {
            if (distanceMatrix[i][j] <= (1 + 2*alpha) * r) {
                ball_weight += W.get(j);
            }
        }
        return ball_weight;
    }

    /**
     * Computes the value of the objective function
     *
     * @param P set of points
     * @param S set of centers
     * @param z outliers
     * @return value of the objective function
     */
    public static double ComputeObjective(ArrayList<Vector> P, ArrayList<Vector> S, int z) {
        // no problem if points are overlapping
        ArrayList<Double> distsToCenters = new ArrayList<>(Collections.nCopies(P.size(), Double.MAX_VALUE));

        for(Vector s : S) { // for every center
            final int indS = P.indexOf(s); // compute indexes (in P) of elements in S

            // to keep update distsToCenters with distances to subset of visited centers
            for (int p = 0; p < P.size(); ++p) { // for every point, calculate distance to center
                                                 // since S is a subset of P, we can use precalculated distanceMatrix
                if (distanceMatrix[indS][p] < distsToCenters.get(p)) {
                    distsToCenters.set(p, distanceMatrix[indS][p]);
                }
            }
        }

        distsToCenters.sort(Comparator.naturalOrder());

        if(distsToCenters.get(S.size() - 1) != 0)   // to demonstrate that first s elements are 0
            throw new AssertionError("Every center must be distant 0 from centers");

        return distsToCenters.get(P.size() - 1 - z);
    }

    public static void main(String[] args) throws IOException {

        if (args.length != 3) {
            throw new IllegalArgumentException("USAGE: file path, number of centers, number of allowed outliers");
        }

        final String path = args[0];
        final int k = Integer.parseInt(args[1]);
        final int z = Integer.parseInt(args[2]);

        ArrayList<Vector> inputPoints = readVectorsSeq(path);
        ArrayList<Long> weights = new ArrayList<>(Collections.nCopies(inputPoints.size(), 1L));

        System.out.println("Input size n = " + inputPoints.size());
        System.out.println("Number of centers k = " + k);
        System.out.println("Number of outliers z = " + z);

        final long startTime = System.currentTimeMillis();

        // Distance Matrix calculation
        distanceMatrix = new double[inputPoints.size()][inputPoints.size()];
        for (int i = 0; i < inputPoints.size(); ++i) {
            for (int j = i + 1; j < inputPoints.size(); ++j) {
                distanceMatrix[i][j] = distanceMatrix[j][i] = Math.sqrt(Vectors.sqdist(inputPoints.get(i), inputPoints.get(j)));
            }
        }

        ArrayList<Vector> solution = SeqWeightedOutliers(inputPoints, weights, k, z, 0);
        final long timeElapsed = System.currentTimeMillis() - startTime;

        final double objective = ComputeObjective(inputPoints, solution, z);

        System.out.println("Objective function = " + objective);
        System.out.println("Time of SeqWeightedOutliers = " + timeElapsed);
    }

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
}
