/**
 * 
 */
package org.ithinktree.becky.jprimewrappers;

import java.util.Map;

import dr.math.distributions.Distribution;
import se.cbb.jprime.math.Continuous1DPDDependent;
import se.cbb.jprime.math.PRNG;
import se.cbb.jprime.math.RealInterval;
import se.cbb.jprime.mcmc.ChangeInfo;
import se.cbb.jprime.mcmc.Dependent;

/**
 * @author Arman D. Bilge
 *
 */
public class JPrIMEContinuous1DPDDependentWrapperForBEASTDistribution implements
		Continuous1DPDDependent, WrappedBEASTObject<Distribution> {

	private final Distribution distribution;
	
	/**
	 * 
	 */
	public JPrIMEContinuous1DPDDependentWrapperForBEASTDistribution(Distribution d) {
		distribution = d;
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getCDF(double)
	 */
	@Override
	public double getCDF(double d) {
		return distribution.cdf(d);
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getCV()
	 */
	@Override
	public double getCV() {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getDomainInterval()
	 */
	@Override
	public RealInterval getDomainInterval() {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getMean()
	 */
	@Override
	public double getMean() {
		return distribution.mean();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getMedian()
	 */
	@Override
	public double getMedian() {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getMode()
	 */
	@Override
	public double getMode() {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getPDF(double)
	 */
	@Override
	public double getPDF(double d) {
		return distribution.pdf(d);
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getProbability(double, double)
	 */
	@Override
	public double getProbability(double a, double b) {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getQuantile(double)
	 */
	@Override
	public double getQuantile(double d) {
		return distribution.quantile(d);
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getStandardDeviation()
	 */
	@Override
	public double getStandardDeviation() {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#getVariance()
	 */
	@Override
	public double getVariance() {
		return distribution.variance();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#sampleValue(se.cbb.jprime.math.PRNG)
	 */
	@Override
	public double sampleValue(PRNG arg0) {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#setMean(double)
	 */
	@Override
	public void setMean(double arg0) {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#setStandardDeviation(double)
	 */
	@Override
	public void setStandardDeviation(double arg0) {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.Continuous1DPD#setVariance(double)
	 */
	@Override
	public void setVariance(double arg0) {
		throw new UnsupportedOperationException();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.ProbabilityDistribution#getName()
	 */
	@Override
	public String getName() {
		return distribution.toString();
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.ProbabilityDistribution#getNoOfDimensions()
	 */
	@Override
	public int getNoOfDimensions() {
		return 0;
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.math.ProbabilityDistribution#getNoOfParameters()
	 */
	@Override
	public int getNoOfParameters() {
		return 0;
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.mcmc.ProperDependent#cacheAndUpdate(java.util.Map, boolean)
	 */
	@Override
	public void cacheAndUpdate(Map<Dependent, ChangeInfo> arg0, boolean arg1) {
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.mcmc.ProperDependent#clearCache(boolean)
	 */
	@Override
	public void clearCache(boolean arg0) {
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.mcmc.ProperDependent#getParentDependents()
	 */
	@Override
	public Dependent[] getParentDependents() {
		return new Dependent[0];
	}

	/* (non-Javadoc)
	 * @see se.cbb.jprime.mcmc.ProperDependent#restoreCache(boolean)
	 */
	@Override
	public void restoreCache(boolean arg0) {
	}

	/* (non-Javadoc)
	 * @see org.ithinktree.becky.jprimewrappers.WrappedBEASTObject#getWrappedBEASTObject()
	 */
	@Override
	public Distribution getWrappedBEASTObject() {
		return distribution;
	}

}
