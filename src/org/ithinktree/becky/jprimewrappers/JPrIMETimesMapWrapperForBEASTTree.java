package org.ithinktree.becky.jprimewrappers;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import se.cbb.jprime.topology.DoubleMap;
import se.cbb.jprime.topology.TimesMap;

public class JPrIMETimesMapWrapperForBEASTTree extends TimesMap implements
		WrappedBEASTObject<Tree> {

	private final Tree tree;
	
	public JPrIMETimesMapWrapperForBEASTTree(Tree tree) {
		super(tree.getId(), new double[0], new double[0]);
		this.tree = tree;
	}

	@Override
	public Tree getWrappedBEASTObject() {
		return tree;
	}

	/**
	 * Returns the absolute time (vertex time) of a vertex.
	 * Identical to get(x).
	 * @param x the vertex.
	 * @return the vertex time.
	 */
	public double getVertexTime(int x) {
		return get(x);
	}
	
	/**
	 * Returns the relative time (arc time) of an arc, indexed by the arc's head x,
	 * i.e. time(p(x))-time(x). If x is the root, the "top time" is returned.
	 * @param x the arc's head.
	 * @return the arc time.
	 */
	public double getArcTime(int x) {
		NodeRef n = tree.getNode(x);
		return tree.isRoot(n) ? tree.getNodeHeight(n) : tree.getBranchLength(n);
	}

	
	/**
	 * Returns the absolute time (vertex time) of a vertex.
	 * Identical to getVertexTime(x).
	 * @param x the vertex.
	 * @return the vertex time.
	 */
	@Override
	public double get(int x) {
		return tree.getNodeHeight(tree.getNode(x));
	}
	
	/**
	 * Unsupported method. Use getVertexTimes() and getArcTimes() instead for low-level manipulation of times.
	 * @param x the vertex.
	 * @param val the new value.
	 */
	@Override
	public void set(int x, double val) {
		throw new UnsupportedOperationException("Cannot set absolute time (vertex time) of a vertex without also " +
				"changing corresponding arc times. Use getVertexTimes() and getArcTimes() instead for low-level " +
				"manipulation of the underlying values.");
	}
	
	/**
	 * Returns the actual vertex times of this map for low-level manipulation.
	 * User must ensure vertex times and arc times are ultrametric and compatible
	 * after changes to its elements. See also sister method <code>getArcTimes()</code>.
	 * @return the internal vertex times.
	 */
	public double[] getVertexTimes() {
		// Avoid messing with BEAST internals (that is what this is for, after all
		throw new UnsupportedOperationException();
	}
	
	/**
	 * Returns the actual arc times of this map for low-level manipulation.
	 * User must ensure vertex times and arc times are ultrametric and compatible
	 * after changes to its elements. See also sister method <code>getVetexTimes()</code>.
	 * @return the internal arc times.
	 */
	public double[] getArcTimes() {
		// Avoid messing with BEAST internals (that is what this is for, after all)
		throw new UnsupportedOperationException();
	}
	
	/**
	 * Returns the total arc time of the entire tree.
	 * @return the total timespan.
	 */
	public double getTotalArcTime() {
		double sum = 0.0;
		for (int i = 0; i < tree.getNodeCount(); ++i) {
			sum += tree.getBranchLength(tree.getNode(i));
		}
		return sum;
	}
	
	@Override
	public void cache(int[] vertices) {
		// Avoiding messing with BEAST internals
	}

	@Override
	public void clearCache() {
		// Avoiding messing with BEAST internals
	}

	@Override
	public void restoreCache() {
		// Avoiding messing with BEAST internals
	}
	
	/**
	 * Returns a map of the arc times for when such a specific need arises.
	 * @return the arc times in a map.
	 */
	public DoubleMap getArcTimesMap() {
		throw new UnsupportedOperationException();
	}

	
}
