package org.ithinktree.becky.JPRIMEWrappers;

import java.util.Arrays;

import se.cbb.jprime.io.SampleDoubleArray;
import se.cbb.jprime.topology.DoubleMap;
import dr.evolution.tree.Tree;

public class JPrIMEDoubleMapWrapperForBEASTTree extends DoubleMap implements
		WrappedBEASTObject<Tree> {

	private final Tree tree;
	
	public JPrIMEDoubleMapWrapperForBEASTTree(Tree tree) {
		super(tree.getId(), 0);
		this.tree = tree;
	}

	@Override
	public Tree getWrappedBEASTObject() {
		return tree;
	}

	public String getName() {
		return tree.getId();
	}

	@Override
	public void setName(String name) {
		// Avoid messing with BEAST internals
	}
	
	@Override
	public Object getAsObject(int x) {
		return new Double(get(x));
	}

	@Override
	public void setAsObject(int x, Object value) {
		set(x, ((Double) value).doubleValue());
	}

	/**
	 * Returns the element of a vertex.
	 * @param x the vertex.
	 * @return the value.
	 */
	public double get(int x) {
		return tree.getBranchLength(tree.getNode(x));
	}
	
	/**
	 * Sets the element of a vertex.
	 * @param x the vertex.
	 * @param val the value.
	 */
	public void set(int x, double val) {
		// Avoid messing with BEAST internals
	}

	@Override
	public int getNoOfSubParameters() {
		return getSize();
	}

	@Override
	public void cache(int[] vertices) {
		// Avoid messing with BEAST internals
	}

	@Override
	public void clearCache() {
		// Avoid messing with BEAST internals
	}

	@Override
	public void restoreCache() {
		// Avoid messing with BEAST internals
	}

	@Override
	public Class<?> getSampleType() {
		return SampleDoubleArray.class;
	}

	@Override
	public String getSampleHeader() {
		return getName();
	}

	@Override
	public String getSampleValue(SamplingMode mode) {
		int nodeCount = tree.getNodeCount();
		double[] values = new double[nodeCount];
		for (int i = 0; i < nodeCount; ++i)
			values[i] = tree.getBranchLength(tree.getNode(i));
		return SampleDoubleArray.toString(values);
	}

	@Override
	public double getValue(int idx) {
		return get(idx);
	}

	@Override
	public void setValue(int idx, double value) {
		set(idx, value);
	}

	@Override
	public int getSize() {
		return tree.getNodeCount();
	}

	@Override
	public String toString() {
		int nodeCount = tree.getNodeCount();
		double[] values = new double[nodeCount];
		for (int i = 0; i < nodeCount; ++i)
			values[i] = tree.getBranchLength(tree.getNode(i));
		return Arrays.toString(values);
	}
	
}
