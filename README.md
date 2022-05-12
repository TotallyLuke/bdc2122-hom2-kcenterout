# bdc2122-hom2-kcenterout
Homework 2 for Big Data Computing course

## INTRODUCTION

Homeworks 2 and 3 will focus on the **k-center with z outliers problem**, that is, the robust version of the k-center problem which is useful in the analysis of noisy data (a quite common scenario in big data computing). Given a set P of points and two integers k and z, the problem requires to determine a set S ⊂ P of k centers which **minimize the maximum distance of a point of P-Z from S**, where Z are the z farthest points of P from S. In other words, with respect to the standard k-center problem, in the k-center with z outliers problem, we are allowed to disregard the z farthest points from S in the objective function. Unfortunately, the solution of this problem turns out much harder than the one of the standard k-center problem. The 3-approximation sequential algorithm by Charikar et al. for k-center with z outliers, which we call **kcenterOUT**, is simple to implement but has superlinear complexity (more than quadratic, unless sophisticated data structures are used). The rigorous specification of the problem and the description of kcenterOUT with pseudocode, can be found in this [set of slides](https://esami.elearning.unipd.it/pluginfile.php/319634/mod_page/content/34/PresentationHW2.pdf?time=1651045518761).

The two homeworks will demonstrate that in this case a coreset-based approach can be successfully employed. **In Homework 2 you will implement the 3-approximation sequential algorithm** and will get a **first-hand experience of its inefficiency**. In Homework 3, you will implement a 2-round MapReduce coreset-based algorithm for the problem, where the use of the inefficient 3-approximation is confined to a small coreset computed in parallel through the efficient Farthest-First Traversal.

## REPRESENTATION of POINTS

We will work with points in Euclidean space (real cooordinates) and with the Euclidean L2-distance.

### FOR JAVA USERS

In Spark, points can be represented as instances of the class `org.apache.spark.mllib.linalg.Vector` and can be manipulated through static methods offered by the class `org.apache.spark.mllib.linalg.Vectors`. For example, method `Vectors.dense(x)` transforms an array `x` of `double` into an instance of class `Vector`.  Given two points `x` and `y`, instances of `Vector`, their Euclidean L2-distance can be computed by invoking: `Math.sqrt(Vectors.sqdist(x, y))`. Details on the classes can be found in the Spark Java API.  

#### Warning

Make sure to use the classes from the `org.apache.spark.mllib package`. There are classes with the same name in `org.apache.spark.ml` package which are functionally equivalent, but incompatible with those of the `org.apache.spark.mllib package`.

## TASK for HOMEWORK 2

You must:

1. Develop a method `SeqWeightedOutliers(P,W,k,z,alpha)` which *implements the weighted variant of kcenterOUT* (the 3-approximation algorithm for k-center with z-outliers). The method returns the set of centers `S` computed as specified by the algorithm (as `ArrayList<Vector>`). It is understood that the i-th integer in `W` is the weight of the i-th point in `P`. The method takes as input
    * the set of points `P`,as `ArrayList<Vector>`
    * the set of weights `W`, as `ArrayList<Long>`
    * the number of centers `k`,
    * the number of outliers `z`,
    * and the coefficient `alpha` used by the algorithm.
2. Develop a method `ComputeObjective(P,S,z)` which *computes the value of the objective function* for the set of points `P`, the set of centers `S`, and `z` outliers (the number of centers, which is the size of `S`, is not needed as a parameter). Hint: you may
    1. compute all distances `d(x,S)`, for every `x` in `P`,
    2. sort them,
    3. exclude the `z` largest distances,
    4. and return the largest among the remaining ones.
    Note that in this case we are not using weights!
3. Write a program G050HW2.java, which *receives in input* the following command-line (CLI) arguments:
    * A `path` to a text file containing point set in Euclidean space. Each line of the file contains, separated by commas, the coordinates of a point. Your program should make no assumptions on the number of dimensions!
    * An integer `k` (the number of centers).
    * An integer `z` (the number of allowed outliers).

The program must do the following:

* Read the points in the input file into an `ArrayList<Vector>` called `inputPoints`. To this purpose, you can use the code provided in the file `InputCode.java`, for Java users.
* Create an `ArrayList<Long>` called `weights` of the same cardinality of `inputPoints`, initialized with all 1's. (In this homework we will use unit weights, but in Homework 3 we will need the generality of arbitrary integer weights!).
* Run `SeqWeightedOutliers(inputPoints,weights,k,z,0)` to compute a set of (at most) `k` centers. The output of the method must be saved into an `ArrayList<Vector>` called `solution`.
* Run `ComputeObjective(inputPoints,solution,z)` and save the output in a variable called `objective`.
* Return as output the following quantities: |`P`|, `k`, `z`, the initial guess made by `SeqWeightedOutliers(inputPoints,weights,k,z,0)`, the value objective, and the time (in milliseconds) required by the execution of `SeqWeightedOutliers(inputPoints,weights,k,z,0)`. Use the following [output format](https://esami.elearning.unipd.it/pluginfile.php/319634/mod_page/content/34/OutputFormat.txt)

Test your program using the datasets available [here](https://esami.elearning.unipd.it/mod/page/view.php?id=44241).
**ADD DATASETS and SUGGESTED VALUES OF k and z**

## SUBMISSION INSTRUCTIONS

Each group must submit a single file (G050HW2.java). Only one student per group must submit the files in Moodle Exam using the link provided in the Homework2 section. Make sure that your code is free from compiling/run-time errors and that you use the file/variable names in the homework description, otherwise your score will be penalized.
