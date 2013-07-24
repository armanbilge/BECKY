package org.ithinktree.becky.JPrIMEWrappers;

import se.cbb.jprime.io.SampleDouble;
import se.cbb.jprime.mcmc.DoubleParameter;
import dr.inference.model.Parameter;

public class JPrIMEDoubleParameterWrapperForBEASTParameter extends
		DoubleParameter implements WrappedBEASTObject<Parameter> {

	private final Parameter parameter;
	
	public JPrIMEDoubleParameterWrapperForBEASTParameter(Parameter parameter) {
		super(parameter.getId(), 0);
		this.parameter = parameter;
	}

	@Override
	public Parameter getWrappedBEASTObject() {
		return parameter;
	}

	/**
	 * Caches the current value. May e.g. be used by a <code>Proposer</code>.
	 * @param indices is of no importance.
	 */
	@Override
	public void cache(int[] indices) {
		// Avoid messing with BEAST Internals
	}

	@Override
	public void clearCache() {
		// Avoid messing with BEAST Internals
	}

	@Override
	public void restoreCache() {
		// Avoid messing with BEAST Internals
	}

	@Override
	public String getName() {
		return parameter.getId();
	}


	@Override
	public String getSampleHeader() {
		return getName();
	}

	@Override
	public String getSampleValue(SamplingMode mode) {
		return SampleDouble.toString(getValue());
	}

	@Override
	public Class<?> getSampleType() {
		return SampleDouble.class;
	}

	@Override
	public double getValue(int idx) {
		return getValue();
	}

	@Override
	public void setValue(int idx, double value) {
		// Avoid messing with BEAST Internals
	}
	
	/**
	 * Returns this parameter's current value.
	 * @return the value.
	 */
	public double getValue() {
		return parameter.getParameterValue(0);
	}

	/**
	 * Sets this parameter's current value.
	 * @param value the new value.
	 */
	public void setValue(double value) {
		// Avoid messing with BEAST Internals
	}

	
}
