/**
 * SimpleGeographicalCophylogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.inference.model.Parameter;

/**
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public class SimpleGeographicalCophylogenyModel extends SimpleCophylogenyModel {

	protected final AbstractMultivariateTraitLikelihood hostGeographicTrait;
	protected final AbstractMultivariateTraitLikelihood symbiontGeographicTrait;
	
	/**
	 * @param duplicationRateParameter
	 * @param hostShiftRateParameter
	 * @param lossRateParameter
	 * @param units
	 */
	public SimpleGeographicalCophylogenyModel(
			Parameter duplicationRateParameter,
			Parameter hostShiftRateParameter, Parameter lossRateParameter,
			Type units,
			AbstractMultivariateTraitLikelihood hostGeographicTrait,
			AbstractMultivariateTraitLikelihood symbiontGeographicTrait) {
		super(duplicationRateParameter, hostShiftRateParameter,
				lossRateParameter, units);
		this.hostGeographicTrait = hostGeographicTrait;
		this.symbiontGeographicTrait = symbiontGeographicTrait;
	}

	public double calculateNodeLogLikelihood(final MutableTree symbiontTree, final NodeRef self,
			final NodeRef child1, final NodeRef child2, final Tree hostTree, final NodeRef selfHost,
			final NodeRef child1Host, final NodeRef child2Host, final BranchRates branchRates) {

//		hostGeographicTrait.
		
		return super.calculateNodeLogLikelihood(symbiontTree, self, child1, child2, hostTree, selfHost, child1Host, child2Host, branchRates);
		
	}
	
}
