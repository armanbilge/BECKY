/*
 * IntegratedMultivariateTraitLikelihood.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.continuous;

import dr.evolution.tree.MultivariateTraitTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.Vector;
import dr.util.Author;
import dr.util.Citation;

import java.util.Arrays;
import java.util.List;

/**
 * A multivariate trait likelihood that analytically integrates out the unobserved trait values at all internal
 * and root nodes
 *
 * @author Marc A. Suchard
 */
public abstract class IntegratedMultivariateTraitLikelihood extends AbstractMultivariateTraitLikelihood {

    public static final double LOG_SQRT_2_PI = 0.5 * Math.log(2 * Math.PI);

    public IntegratedMultivariateTraitLikelihood(String traitName,
                                                 MultivariateTraitTree treeModel,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 CompoundParameter traitParameter,
                                                 List<Integer> missingIndices,
                                                 boolean cacheBranches, boolean scaleByTime, boolean useTreeLength,
                                                 BranchRateModel rateModel, Model samplingDensity,
                                                 boolean reportAsMultivariate,
                                                 boolean reciprocalRates) {

        this(traitName, treeModel, diffusionModel, traitParameter, null, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, samplingDensity, reportAsMultivariate, reciprocalRates);
    }

    public IntegratedMultivariateTraitLikelihood(String traitName,
                                                 MultivariateTraitTree treeModel,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 CompoundParameter traitParameter,
                                                 Parameter deltaParameter,
                                                 List<Integer> missingIndices,
                                                 boolean cacheBranches, boolean scaleByTime, boolean useTreeLength,
                                                 BranchRateModel rateModel, Model samplingDensity,
                                                 boolean reportAsMultivariate,
                                                 boolean reciprocalRates) {
        this(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches,
                scaleByTime, useTreeLength, rateModel, null, samplingDensity, reportAsMultivariate, reciprocalRates);
    }


    private final CacheHelper cacheHelper;

    public IntegratedMultivariateTraitLikelihood(String traitName,
                                                 MultivariateTraitTree treeModel,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 CompoundParameter traitParameter,
                                                 Parameter deltaParameter,
                                                 List<Integer> missingIndices,
                                                 boolean cacheBranches, boolean scaleByTime, boolean useTreeLength,
                                                 BranchRateModel rateModel,
                                                 List<BranchRateModel> driftModels,
                                                 Model samplingDensity,
                                                 boolean reportAsMultivariate,
                                                 boolean reciprocalRates) {

        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, driftModels, samplingDensity, reportAsMultivariate, reciprocalRates);


        // Delegate caches to helper
        meanCache = new double[dim * treeModel.getNodeCount()];
        if (driftModels != null) {
            cacheHelper = new DriftCacheHelper(dim * treeModel.getNodeCount(), cacheBranches); // new DriftCacheHelper ....
        } else {
            cacheHelper = new CacheHelper(dim * treeModel.getNodeCount(), cacheBranches);
        }

        drawnStates = new double[dim * treeModel.getNodeCount()];
        upperPrecisionCache = new double[treeModel.getNodeCount()];
        lowerPrecisionCache = new double[treeModel.getNodeCount()];
        logRemainderDensityCache = new double[treeModel.getNodeCount()];

        if (cacheBranches) {
            storedMeanCache = new double[dim * treeModel.getNodeCount()];
            storedUpperPrecisionCache = new double[treeModel.getNodeCount()];
            storedLowerPrecisionCache = new double[treeModel.getNodeCount()];
            storedLogRemainderDensityCache = new double[treeModel.getNodeCount()];
        }

        // Set up reusable temporary storage
        Ay = new double[dimTrait];
        tmpM = new double[dimTrait][dimTrait];
        tmp2 = new double[dimTrait];

        zeroDimVector = new double[dim];

        missingTraits = new MissingTraits.CompletelyMissing(treeModel, missingIndices, dim);
        setTipDataValuesForAllNodes();

    }

    private void setTipDataValuesForAllNodes() {
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getExternalNode(i);
            setTipDataValuesForNode(node);
        }
        missingTraits.handleMissingTips();
    }

    @SuppressWarnings("unused")
    public double getTotalTreePrecision() {
        getLogLikelihood(); // Do peeling if necessary
        final int rootIndex = treeModel.getRoot().getNumber();
        return lowerPrecisionCache[rootIndex];
    }

    private void setTipDataValuesForNode(NodeRef node) {
        // Set tip data values
        int index = node.getNumber();
        double[] traitValue = traitParameter.getParameter(index).getParameterValues();
        if (traitValue.length < dim) {
            throw new RuntimeException("The trait parameter for the tip with index, " + index + ", is too short");
        }

        cacheHelper.setTipMeans(traitValue, dim, index, node);
//        System.arraycopy(traitValue, 0, meanCache
////                cacheHelper.getMeanCache()
//                , dim * index, dim);
    }

    public double[] getTipDataValues(int index) {
        double[] traitValue = new double[dim];
        System.arraycopy(cacheHelper.getMeanCache(), dim * index, traitValue, 0, dim);
        return traitValue;
    }

    public void setTipDataValuesForNode(int index, double[] traitValue) {
        // Set tip data values
        // cacheHelper.copyToMeanCache(traitValue, dim*index, dim);
        cacheHelper.setTipMeans(traitValue, dim, index);
        //System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
        makeDirty();
    }

    protected String extraInfo() {
        return "\tSample internal node traits: false\n";
    }

    public List<Citation> getCitations() {
        List<Citation> citations = super.getCitations();
        citations.add(
                new Citation(
                        new Author[]{
                                new Author("O", "Pybus"),
                                new Author("P", "Lemey"),
                                new Author("A", "Rambaut"),
                                new Author("MA", "Suchard")
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
        return citations;
    }

    public double getLogDataLikelihood() {
        return getLogLikelihood();
    }

    public abstract boolean getComputeWishartSufficientStatistics();

    public double calculateLogLikelihood() {

        double logLikelihood = 0;
        double[][] traitPrecision = diffusionModel.getPrecisionmatrix();
        double logDetTraitPrecision = Math.log(diffusionModel.getDeterminantPrecisionMatrix());
        double[] conditionalRootMean = tmp2;

        final boolean computeWishartStatistics = getComputeWishartSufficientStatistics();

        if (computeWishartStatistics) {
//            if (wishartStatistics == null) {
            wishartStatistics = new WishartSufficientStatistics(dimTrait);
//            } else {
//                wishartStatistics.clear();
//            }
        }

        // Use dynamic programming to compute conditional likelihoods at each internal node
        postOrderTraverse(treeModel, treeModel.getRoot(), traitPrecision, logDetTraitPrecision, computeWishartStatistics);

        if (DEBUG) {
            System.err.println("mean: " + new Vector(cacheHelper.getMeanCache()));
            System.err.println("correctedMean: " + new Vector(cacheHelper.getCorrectedMeanCache()));
            System.err.println("upre: " + new Vector(upperPrecisionCache));
            System.err.println("lpre: " + new Vector(lowerPrecisionCache));
            System.err.println("cach: " + new Vector(logRemainderDensityCache));
        }

        // Compute the contribution of each datum at the root
        final int rootIndex = treeModel.getRoot().getNumber();

        // Precision scalar of datum conditional on root
        double conditionalRootPrecision = lowerPrecisionCache[rootIndex];

        for (int datum = 0; datum < numData; datum++) {

            double thisLogLikelihood = 0;

            // Get conditional mean of datum conditional on root
            // System.arraycopy(meanCache, rootIndex * dim + datum * dimTrait, conditionalRootMean, 0, dimTrait);
            System.arraycopy(cacheHelper.getMeanCache(), rootIndex * dim + datum * dimTrait, conditionalRootMean, 0, dimTrait);

            if (DEBUG) {
                System.err.println("Datum #" + datum);
                System.err.println("root mean: " + new Vector(conditionalRootMean));
                System.err.println("root prec: " + conditionalRootPrecision);
                System.err.println("diffusion prec: " + new Matrix(traitPrecision));
            }

            // B = root prior precision
            // z = root prior mean
            // A = likelihood precision
            // y = likelihood mean

            // y'Ay
            double yAy = computeWeightedAverageAndSumOfSquares(conditionalRootMean, Ay, traitPrecision, dimTrait,
                    conditionalRootPrecision); // Also fills in Ay

            if (conditionalRootPrecision != 0) {
                thisLogLikelihood += -LOG_SQRT_2_PI * dimTrait
                        + 0.5 * (logDetTraitPrecision + dimTrait * Math.log(conditionalRootPrecision) - yAy);
            }

            if (DEBUG) {
                double[][] T = new double[dimTrait][dimTrait];
                for (int i = 0; i < dimTrait; i++) {
                    for (int j = 0; j < dimTrait; j++) {
                        T[i][j] = traitPrecision[i][j] * conditionalRootPrecision;
                    }
                }
                System.err.println("Conditional root MVN precision = \n" + new Matrix(T));
                System.err.println("Conditional root MVN density = " + MultivariateNormalDistribution.logPdf(
                        conditionalRootMean, new double[dimTrait], T,
                        Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(T)), 1.0));
            }

            if (integrateRoot) {
                // Integrate root trait out against rootPrior
                thisLogLikelihood += integrateLogLikelihoodAtRoot(conditionalRootMean, Ay, tmpM, traitPrecision,
                        conditionalRootPrecision); // Ay is destroyed
            }

            if (DEBUG) {
                System.err.println("yAy = " + yAy);
                System.err.println("logLikelihood (before remainders) = " + thisLogLikelihood +
                        " (should match conditional root MVN density when root not integrated out)");
            }

            logLikelihood += thisLogLikelihood;
        }

        logLikelihood += sumLogRemainders();
        if (DEBUG) {
            System.out.println("logLikelihood is " + logLikelihood);
        }

        if (DEBUG) { // Root trait is univariate!!!
            System.err.println("logLikelihood (final) = " + logLikelihood);
//            checkViaLargeMatrixInversion();
        }

        if (DEBUG_PNAS) {
            checkLogLikelihood(logLikelihood, sumLogRemainders(), conditionalRootMean,
                    conditionalRootPrecision, traitPrecision);
        }

        areStatesRedrawn = false;  // Should redraw internal node states when needed
        return logLikelihood;
    }

    protected void checkLogLikelihood(double loglikelihood, double logRemainders,
                                      double[] conditionalRootMean, double conditionalRootPrecision,
                                      double[][] traitPrecision) {
        // Do nothing; for checking PNAS paper
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == traitParameter) { // A tip value got updated
            if (index > dimTrait * treeModel.getExternalNodeCount()) {
                throw new RuntimeException("Attempting to update an invalid index");
            }
            cacheHelper.setMeanCache(index, traitParameter.getValue(index));
            // meanCache[index] = traitParameter.getValue(index);
            likelihoodKnown = false;
//            if (!cacheBranches) {
//                throw new RuntimeException("Must cache means in IMTL if they are random");
//            }
            // TODO Need better solution.  If tips are random, cacheBranches should be true (to get reset).
            // TODO However, jitter calls setParameterValue() on the tips at initialization
        }
        super.handleVariableChangedEvent(variable, index, type);
    }

    protected static double computeWeightedAverageAndSumOfSquares(double[] y, double[] Ay, double[][] A,
                                                                  int dim, double scale) {
        // returns Ay and yAy
        double yAy = 0;
        for (int i = 0; i < dim; i++) {
            Ay[i] = 0;
            for (int j = 0; j < dim; j++)
                Ay[i] += A[i][j] * y[j] * scale;
            yAy += y[i] * Ay[i];
        }
        return yAy;
    }

    private double sumLogRemainders() {
        double sumLogRemainders = 0;
        for (double r : logRemainderDensityCache)
            sumLogRemainders += r;
        // Could skip leafs
        return sumLogRemainders;
    }

    protected abstract double integrateLogLikelihoodAtRoot(double[] conditionalRootMean,
                                                           double[] marginalRootMean,
                                                           double[][] temporaryStorage,
                                                           double[][] treePrecisionMatrix,
                                                           double conditionalRootPrecision);

    public void makeDirty() {
        super.makeDirty();
        areStatesRedrawn = false;
    }

    void postOrderTraverse(MultivariateTraitTree treeModel, NodeRef node, double[][] precisionMatrix,
                           double logDetPrecisionMatrix, boolean cacheOuterProducts) {

        final int thisNumber = node.getNumber();

        if (treeModel.isExternal(node)) {

            // Fill in precision scalar, traitValues already filled in

            if (missingTraits.isCompletelyMissing(thisNumber)) {
                upperPrecisionCache[thisNumber] = 0;
                lowerPrecisionCache[thisNumber] = 0; // Needed in the pre-order traversal
            } else { // not missing tip trait
                upperPrecisionCache[thisNumber] = 1.0 / getRescaledBranchLengthForPrecision(node);
                lowerPrecisionCache[thisNumber] = Double.POSITIVE_INFINITY;
            }
            return;
        }

        final NodeRef childNode0 = treeModel.getChild(node, 0);
        final NodeRef childNode1 = treeModel.getChild(node, 1);

        postOrderTraverse(treeModel, childNode0, precisionMatrix, logDetPrecisionMatrix, cacheOuterProducts);
        postOrderTraverse(treeModel, childNode1, precisionMatrix, logDetPrecisionMatrix, cacheOuterProducts);

        final int childNumber0 = childNode0.getNumber();
        final int childNumber1 = childNode1.getNumber();
        final int meanOffset0 = dim * childNumber0;
        final int meanOffset1 = dim * childNumber1;
        final int meanThisOffset = dim * thisNumber;

        final double precision0 = upperPrecisionCache[childNumber0];
        final double precision1 = upperPrecisionCache[childNumber1];
        final double totalPrecision = precision0 + precision1;

        lowerPrecisionCache[thisNumber] = totalPrecision;

        // Multiply child0 and child1 densities

        // Delegate this!
        cacheHelper.computeMeanCaches(meanThisOffset, meanOffset0, meanOffset1,
                totalPrecision, precision0, precision1, missingTraits, node, childNode0, childNode1);
//        if (totalPrecision == 0) {
//            System.arraycopy(zeroDimVector, 0, meanCache, meanThisOffset, dim);
//        } else {
//            // Delegate in case either child is partially missing
//            // computeCorrectedWeightedAverage
//            missingTraits.computeWeightedAverage(meanCache,
//                    meanOffset0, precision0,
//                    meanOffset1, precision1,
//                    meanThisOffset, dim);
//        }
        // In this delegation, you can call
        //getShiftForBranchLength(node);

        if (!treeModel.isRoot(node)) {
            // Integrate out trait value at this node
            double thisPrecision = 1.0 / getRescaledBranchLengthForPrecision(node);
            if (Double.isInfinite(thisPrecision)) {
                upperPrecisionCache[thisNumber] = totalPrecision;
            } else {
                upperPrecisionCache[thisNumber] = totalPrecision * thisPrecision / (totalPrecision + thisPrecision);
            }
        }

        // Compute logRemainderDensity

        logRemainderDensityCache[thisNumber] = 0;

        if (precision0 != 0 && precision1 != 0) {

            incrementRemainderDensities(
                    precisionMatrix,
                    logDetPrecisionMatrix, thisNumber, meanThisOffset,
                    meanOffset0,
                    meanOffset1,
                    precision0,
                    precision1,
                    cacheOuterProducts);
        }
    }

    private void incrementRemainderDensities(double[][] precisionMatrix,
                                             double logDetPrecisionMatrix,
                                             int thisIndex,
                                             int thisOffset,
                                             int childOffset0,
                                             int childOffset1,
                                             double precision0,
                                             double precision1,
                                             boolean cacheOuterProducts) {

        final double remainderPrecision = precision0 * precision1 / (precision0 + precision1);

        if (cacheOuterProducts) {
            incrementOuterProducts(thisOffset, childOffset0, childOffset1, precision0, precision1);
        }

        for (int k = 0; k < numData; k++) {

            double childSS0 = 0;
            double childSS1 = 0;
            double crossSS = 0;

            for (int i = 0; i < dimTrait; i++) {

                // In case of no drift, getCorrectedMeanCache() simply returns mean cache
                // final double wChild0i = meanCache[childOffset0 + k * dimTrait + i] * precision0;
                final double wChild0i = cacheHelper.getCorrectedMeanCache()[childOffset0 + k * dimTrait + i] * precision0;
                // final double wChild1i = meanCache[childOffset1 + k * dimTrait + i] * precision1;
                final double wChild1i = cacheHelper.getCorrectedMeanCache()[childOffset1 + k * dimTrait + i] * precision1;

                for (int j = 0; j < dimTrait; j++) {

                    // subtract "correction"
                    //final double child0j = meanCache[childOffset0 + k * dimTrait + j];
                    final double child0j = cacheHelper.getCorrectedMeanCache()[childOffset0 + k * dimTrait + j];
                    // subtract "correction"
                    //final double child1j = meanCache[childOffset1 + k * dimTrait + j];
                    final double child1j = cacheHelper.getCorrectedMeanCache()[childOffset1 + k * dimTrait + j];

                    childSS0 += wChild0i * precisionMatrix[i][j] * child0j;
                    childSS1 += wChild1i * precisionMatrix[i][j] * child1j;

                    // make sure meanCache in following is not "corrected"
                    // crossSS += (wChild0i + wChild1i) * precisionMatrix[i][j] * meanCache[thisOffset + k * dimTrait + j];
                    crossSS += (wChild0i + wChild1i) * precisionMatrix[i][j] * cacheHelper.getMeanCache()[thisOffset + k * dimTrait + j];
                }
            }

            logRemainderDensityCache[thisIndex] +=
                    -dimTrait * LOG_SQRT_2_PI
                            + 0.5 * (dimTrait * Math.log(remainderPrecision) + logDetPrecisionMatrix)
                            - 0.5 * (childSS0 + childSS1 - crossSS);
        }
    }

    private void incrementOuterProducts(int thisOffset,
                                        int childOffset0,
                                        int childOffset1,
                                        double precision0,
                                        double precision1) {

        final double[][] outerProduct = wishartStatistics.getScaleMatrix();

        for (int k = 0; k < numData; k++) {

            for (int i = 0; i < dimTrait; i++) {

                // final double wChild0i = meanCache[childOffset0 + k * dimTrait + i] * precision0;
                // final double wChild1i = meanCache[childOffset1 + k * dimTrait + i] * precision1;
                final double wChild0i = cacheHelper.getCorrectedMeanCache()[childOffset0 + k * dimTrait + i] * precision0;
                final double wChild1i = cacheHelper.getCorrectedMeanCache()[childOffset1 + k * dimTrait + i] * precision1;

                for (int j = 0; j < dimTrait; j++) {

                    //final double child0j = meanCache[childOffset0 + k * dimTrait + j];
                    //final double child1j = meanCache[childOffset1 + k * dimTrait + j];
                    final double child0j = cacheHelper.getCorrectedMeanCache()[childOffset0 + k * dimTrait + j];
                    final double child1j = cacheHelper.getCorrectedMeanCache()[childOffset1 + k * dimTrait + j];

                    outerProduct[i][j] += wChild0i * child0j;
                    outerProduct[i][j] += wChild1i * child1j;

                    //outerProduct[i][j] -= (wChild0i + wChild1i) * meanCache[thisOffset + k * dimTrait + j];
                    outerProduct[i][j] -= (wChild0i + wChild1i) * cacheHelper.getMeanCache()[thisOffset + k * dimTrait + j];
                }
            }
        }
        wishartStatistics.incrementDf(1); // Peeled one node
    }

//    private void computeWeightedMeanCache(int thisOffset,
//                                          int childOffset0,
//                                          int childOffset1,
//                                          double precision0,
//                                          double precision1) {
//
//        final double totalVariance = 1.0 / (precision0 + precision1);
//        for (int i = 0; i < dim; i++) {
//            meanCache[thisOffset + i] = (meanCache[childOffset0 + i] * precision0 +
//                    meanCache[childOffset1 + i] * precision1)
//                    * totalVariance;
//        }
//    }

    protected double[] getRootNodeTrait() {
        return getTraitForNode(treeModel, treeModel.getRoot(), traitName);
    }

    public double[] getTraitForNode(Tree tree, NodeRef node, String traitName) {

//        if (tree != treeModel) {
//            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
//        }

        getLogLikelihood();

        if (!areStatesRedrawn)
            redrawAncestralStates();

        int index = node.getNumber();

        double[] trait = new double[dim];
        System.arraycopy(drawnStates, index * dim, trait, 0, dim);
        return trait;
    }

    public void redrawAncestralStates() {

        double[][] treePrecision = diffusionModel.getPrecisionmatrix();
        double[][] treeVariance = new SymmetricMatrix(treePrecision).inverse().toComponents();

        preOrderTraverseSample(treeModel, treeModel.getRoot(), 0, treePrecision, treeVariance);

        if (DEBUG) {
            System.err.println("all draws = " + new Vector(drawnStates));
        }

        areStatesRedrawn = true;
    }

    public void storeState() {
        super.storeState();

        if (cacheBranches) {
            //     System.arraycopy(meanCache, 0, storedMeanCache, 0, meanCache.length);
            cacheHelper.store();
            System.arraycopy(upperPrecisionCache, 0, storedUpperPrecisionCache, 0, upperPrecisionCache.length);
            System.arraycopy(lowerPrecisionCache, 0, storedLowerPrecisionCache, 0, lowerPrecisionCache.length);
            System.arraycopy(logRemainderDensityCache, 0, storedLogRemainderDensityCache, 0, logRemainderDensityCache.length);
        }
    }

    public void restoreState() {
        super.restoreState();

        if (cacheBranches) {
            double[] tmp;

            cacheHelper.restore();
            //  tmp = storedMeanCache;
            //  storedMeanCache = meanCache;
            //  meanCache = tmp;

            tmp = storedUpperPrecisionCache;
            storedUpperPrecisionCache = upperPrecisionCache;
            upperPrecisionCache = tmp;

            tmp = storedLowerPrecisionCache;
            storedLowerPrecisionCache = lowerPrecisionCache;
            lowerPrecisionCache = tmp;

            tmp = storedLogRemainderDensityCache;
            storedLogRemainderDensityCache = logRemainderDensityCache;
            logRemainderDensityCache = tmp;
        }
    }


    // Computes x^t A y, used many times in these computations

    protected static double computeQuadraticProduct(double[] x, double[][] A, double[] y, int dim) {
        double sum = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                sum += x[i] * A[i][j] * y[j];
            }
        }
        return sum;
    }

    // Computes the weighted average of two vectors, used many times in these computations

    public static void computeWeightedAverage(double[] in0, int offset0, double weight0,
                                              double[] in1, int offset1, double weight1,
                                              double[] out2, int offset2,
                                              int length) {

        final double totalInverseWeight = 1.0 / (weight0 + weight1);
        for (int i = 0; i < length; i++) {
            out2[offset2 + i] = (in0[offset0 + i] * weight0 + in1[offset1 + i] * weight1) * totalInverseWeight;
        }
    }


    protected void computeCorrectedWeightedAverage(int offset0, double weight0, NodeRef childNode0,
                                                   int offset1, double weight1, NodeRef childNode1,
                                                   int offset2,
                                                   int length, NodeRef thisNode) {

        final double totalInverseWeight = 1.0 / (weight0 + weight1);
        // TODO fix
//        final double length0 = getRescaledBranchLength(childNode0);
//        final double length1 = getRescaledBranchLength(childNode1);
//        final double thisLength = getRescaledBranchLength(thisNode);
        double[] shift;
        if (!treeModel.isRoot(thisNode)) {
            shift = getShiftForBranchLength(thisNode);
        } else {
            shift = null;
        }
        double[] shiftChild0 = getShiftForBranchLength(childNode0);
        double[] shiftChild1 = getShiftForBranchLength(childNode1);

        if (treeModel.isExternal(childNode0)) {
            for (int i = 0; i < length; i++) {
                correctedMeanCache[offset0 + i] = meanCache[offset0 + i] - shiftChild0[i];
            }
        }

        if (treeModel.isExternal(childNode1)) {
            for (int i = 0; i < length; i++) {
                correctedMeanCache[offset1 + i] = meanCache[offset1 + i] - shiftChild1[i];
            }
        }

        for (int i = 0; i < length; i++) {

            // meanCache[offset2 + i] = ((meanCache[offset0 + i] -  shiftChild0[i]) * weight0 + (meanCache[offset1 + i] - shiftChild1[i]) * weight1) * totalInverseWeight;
            meanCache[offset2 + i] = (correctedMeanCache[offset0 + i] * weight0 + correctedMeanCache[offset1 + i] * weight1) * totalInverseWeight;
            if (!treeModel.isRoot(thisNode)) {
                correctedMeanCache[offset2 + i] = meanCache[offset2 + i] - shift[i];
            } else {
                correctedMeanCache[offset2 + i] = meanCache[offset2 + i];
            }
        }
    }


    protected abstract double[][] computeMarginalRootMeanAndVariance(double[] conditionalRootMean,
                                                                     double[][] treePrecisionMatrix,
                                                                     double[][] treeVarianceMatrix,
                                                                     double conditionalRootPrecision);


    private void preOrderTraverseSample(MultivariateTraitTree treeModel, NodeRef node, int parentIndex, double[][] treePrecision,
                                        double[][] treeVariance) {

        final int thisIndex = node.getNumber();

        if (treeModel.isRoot(node)) {
            // draw root

            double[] rootMean = new double[dimTrait];
            final int rootIndex = treeModel.getRoot().getNumber();
            double rootPrecision = lowerPrecisionCache[rootIndex];

            for (int datum = 0; datum < numData; datum++) {
                // System.arraycopy(meanCache, thisIndex * dim + datum * dimTrait, rootMean, 0, dimTrait);
                System.arraycopy(cacheHelper.getMeanCache(), thisIndex * dim + datum * dimTrait, rootMean, 0, dimTrait);

                double[][] variance = computeMarginalRootMeanAndVariance(rootMean, treePrecision, treeVariance,
                        rootPrecision);

                double[] draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(rootMean, variance);

                if (DEBUG_PREORDER) {
                    Arrays.fill(draw, 1.0);
                }

                System.arraycopy(draw, 0, drawnStates, rootIndex * dim + datum * dimTrait, dimTrait);

                if (DEBUG) {
                    System.err.println("Root mean: " + new Vector(rootMean));
                    System.err.println("Root var : " + new Matrix(variance));
                    System.err.println("Root draw: " + new Vector(draw));
                }
            }
        } else { // draw conditional on parentState

            if (!missingTraits.isCompletelyMissing(thisIndex)
                    && !missingTraits.isPartiallyMissing(thisIndex)) {

                //System.arraycopy(meanCache, thisIndex * dim, drawnStates, thisIndex * dim, dim);
                System.arraycopy(cacheHelper.getMeanCache(), thisIndex * dim, drawnStates, thisIndex * dim, dim);

            } else {

                if (missingTraits.isPartiallyMissing(thisIndex)) {
                    throw new RuntimeException("Partially missing values are not yet implemented");
                }
                // This code should work for sampling a missing tip trait as well, but needs testing

                // parent trait at drawnStates[parentOffset]
                double precisionToParent = 1.0 / getRescaledBranchLengthForPrecision(node);
                double precisionOfNode = lowerPrecisionCache[thisIndex];
                double totalPrecision = precisionOfNode + precisionToParent;

                double[] mean = Ay; // temporary storage
                double[][] var = tmpM; // temporary storage

                for (int datum = 0; datum < numData; datum++) {

                    int parentOffset = parentIndex * dim + datum * dimTrait;
                    int thisOffset = thisIndex * dim + datum * dimTrait;

                    if (DEBUG) {
                        double[] parentValue = new double[dimTrait];
                        System.arraycopy(drawnStates, parentOffset, parentValue, 0, dimTrait);
                        System.err.println("Parent draw: " + new Vector(parentValue));
                        if (parentValue[0] != drawnStates[parentOffset]) {
                            throw new RuntimeException("Error in setting indices");
                        }
                    }

                    for (int i = 0; i < dimTrait; i++) {
                        mean[i] = (drawnStates[parentOffset + i] * precisionToParent
                                //  + meanCache[thisOffset + i] * precisionOfNode) / totalPrecision;
                                + cacheHelper.getMeanCache()[thisOffset + i] * precisionOfNode) / totalPrecision;
                        for (int j = 0; j < dimTrait; j++) {
                            var[i][j] = treeVariance[i][j] / totalPrecision;
                        }
                    }
                    double[] draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(mean, var);
                    System.arraycopy(draw, 0, drawnStates, thisOffset, dimTrait);

                    if (DEBUG) {
                        System.err.println("Int prec: " + totalPrecision);
                        System.err.println("Int mean: " + new Vector(mean));
                        System.err.println("Int var : " + new Matrix(var));
                        System.err.println("Int draw: " + new Vector(draw));
                        System.err.println("");
                    }
                }
            }
        }

        if (peel() && !treeModel.isExternal(node)) {
            preOrderTraverseSample(treeModel, treeModel.getChild(node, 0), thisIndex, treePrecision, treeVariance);
            preOrderTraverseSample(treeModel, treeModel.getChild(node, 1), thisIndex, treePrecision, treeVariance);
        }
    }

    protected boolean peel() {
        return true;
    }

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LikelihoodColumn(getId())};
    }

    protected boolean areStatesRedrawn = false;

    protected double[] meanCache;
    protected double[] correctedMeanCache;


    class CacheHelper {

        public CacheHelper(int cacheLength, boolean cacheBranches) {
            // modify code later so we can uncomment the following

            meanCache = new double[cacheLength];
            this.cacheBranches = cacheBranches;
            if (cacheBranches) {
                storedMeanCache = new double[cacheLength];
            }

        }

        public double[] getMeanCache() {
            return meanCache;
        }

        public double[] getCorrectedMeanCache() {
            return meanCache;
        }

        public void store() {
            // if (cacheBranches) {
            System.arraycopy(meanCache, 0, storedMeanCache, 0, meanCache.length);
            // }
        }

        public void restore() {
            //  if (cacheBranches) {
            double[] tmp = storedMeanCache;
            storedMeanCache = meanCache;
            meanCache = tmp;
            // }
        }

        public void doFunkyStuff(NodeRef node) {

        }

        protected boolean cacheBranches;
        //  private double[] meanCache;
        //  private double[] storedMeanCache;

        public void computeMeanCaches(int meanThisOffset, int meanOffset0, int meanOffset1,
                                      double precision0, double precision1, MissingTraits missingTraits) {
            //To change body of created methods use File | Settings | File Templates.
        }

        public void computeMeanCaches(int meanThisOffset, int meanOffset0, int meanOffset1,
                                      double totalPrecision, double precision0, double precision1, MissingTraits missingTraits,
                                      NodeRef thisNode, NodeRef node0, NodeRef node1) {
            if (totalPrecision == 0) {
                System.arraycopy(zeroDimVector, 0, meanCache, meanThisOffset, dim);
            } else {
                // Delegate in case either child is partially missing
                // computeCorrectedWeightedAverage
                missingTraits.computeWeightedAverage(meanCache,
                        meanOffset0, precision0,
                        meanOffset1, precision1,
                        meanThisOffset, dim);
            }
        }

        // public void setTipMeans(double[] meanCache, double[] traitValue, int dim, int index, NodeRef node) {
        public void setTipMeans(double[] traitValue, int dim, int index, NodeRef node) {
            System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
        }

        public void setTipMeans(double[] traitValue, int dim, int index) {
            System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
        }

        public void copyToMeanCache(double[] src, int destPos, int length) {
            System.arraycopy(src, 0, meanCache, destPos, length);
        }

        public void setMeanCache(int index, double value) {
            meanCache[index] = value;
        }


    }

    ;


    class DriftCacheHelper extends CacheHelper {

        public DriftCacheHelper(int cacheLength, boolean cacheBranches) {
            super(cacheLength, cacheBranches);
            correctedMeanCache = new double[cacheLength];
        }

        public double[] getCorrectedMeanCache() {
            return correctedMeanCache;
        }

        public void doFunkyStuff(NodeRef node) {
            getShiftForBranchLength(node);
        }

        public void setTipMeans(double[] traitValue, int dim, int index, NodeRef node) {
            System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
            /*
            double[] shift = getShiftForBranchLength(node);
            for (int i = 0; i < dim; i++) {
                correctedMeanCache[dim * index + i] = traitValue[i] - shift[i];
            }
            */
        }


        public void computeMeanCaches(int meanThisOffset, int meanOffset0, int meanOffset1,
                                      double totalPrecision, double precision0, double precision1, MissingTraits missingTraits,
                                      NodeRef thisNode, NodeRef node0, NodeRef node1) {
            if (totalPrecision == 0) {
                System.arraycopy(zeroDimVector, 0, meanCache, meanThisOffset, dim);
            } else {
                // Delegate in case either child is partially missing
                // computeCorrectedWeightedAverage
                computeCorrectedWeightedAverage(
                        meanOffset0, precision0, node0,
                        meanOffset1, precision1, node1,
                        meanThisOffset, dim, thisNode);
            }
        }


        // private double[] correctedMeanCache;
    }


    // protected boolean hasDrift;

    protected double[] upperPrecisionCache;
    protected double[] lowerPrecisionCache;
    private double[] logRemainderDensityCache;

    private double[] storedMeanCache;
    private double[] storedUpperPrecisionCache;
    private double[] storedLowerPrecisionCache;
    private double[] storedLogRemainderDensityCache;

    private double[] drawnStates;

    protected final boolean integrateRoot = true; // Set to false if conditioning on root value (not fully implemented)
    //  protected final boolean integrateRoot = false;
    protected static boolean DEBUG = false;
    protected static boolean DEBUG_PREORDER = false;
    protected static boolean DEBUG_PNAS = false;

    private double[] zeroDimVector;

    protected WishartSufficientStatistics wishartStatistics;

    // Reusable temporary storage
    protected double[] Ay;
    protected double[][] tmpM;
    protected double[] tmp2;

    protected final MissingTraits missingTraits;

}
