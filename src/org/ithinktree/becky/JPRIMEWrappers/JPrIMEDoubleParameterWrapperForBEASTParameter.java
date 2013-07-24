package org.ithinktree.becky.JPRIMEWrappers;

import dr.inference.model.Parameter;
import se.cbb.jprime.math.ScaleTransformation;
import se.cbb.jprime.mcmc.DoubleParameter;

public class JPrIMEDoubleParameterWrapperForBEASTParameter extends
		DoubleParameter implements WrappedBEASTObject<Parameter> {

	private final Parameter parameter;
	
	public JPrIMEDoubleParameterWrapperForBEASTParameter(String name,
			double initVal) {
		super(name, initVal);
		// TODO Auto-generated constructor stub
	}

	public JPrIMEDoubleParameterWrapperForBEASTParameter(String name,
			ScaleTransformation scale, double initVal) {
		super(name, scale, initVal);
		// TODO Auto-generated constructor stub
	}

	@Override
	public Parameter getWrappedBEASTObject() {
		return parameter;
	}

}
