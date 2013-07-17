/**
 * SimpleTemporalCophyogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky;

//import dr.evolution.tree.NodeRef;
//import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

/**
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public class SimpleTemporalCophylogenyModel extends SimpleCophylogenyModel {

	public SimpleTemporalCophylogenyModel(Parameter duplicationRateParameter,
			Parameter hostShiftRateParameter, Parameter lossRateParameter,
			Type units) {
		super(duplicationRateParameter, hostShiftRateParameter, lossRateParameter,
				units);
		// TODO Auto-generated constructor stub
	}

//	@Override
//	public double calculateNodeLogLikelihood(NodeRef self, Tree hostTree, NodeRef selfHost,
//			NodeRef child1Host, NodeRef child2Host, double selfBranchTime, double child1BranchTime, double child2BranchTime) {
//		if () {
//			return Double.NEGATIVE_INFINITY;
//		} else {
//			return super.calculateNodeLogLikelihood(hostTree, selfHost, child1Host, child2Host, selfBranchTime, child1BranchTime, child2BranchTime);
//		}
//	}
//	
}
