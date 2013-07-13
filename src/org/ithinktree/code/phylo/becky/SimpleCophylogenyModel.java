/**
 * SimpleCophylogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.code.phylo.becky;

import org.apache.commons.math.util.MathUtils;
import org.ithinktree.code.phylo.becky.CophylogenyModel.Relationship.NodalRelationship;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

/**
 * A simple model for cophylogenetic mappings.
 * 
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public class SimpleCophylogenyModel extends CophylogenyModel {

	private Parameter duplicationRateParameter;
	private Parameter hostShiftRateParameter;
	private Parameter lossRateParameter;
		
	/**
	 * 
	 */
	public SimpleCophylogenyModel(Parameter duplicationRateParameter, Parameter hostShiftRateParameter, Parameter lossRateParameter, Type units) {

		super(SimpleCophylogenyModelParser.SIMPLE_COPHYLOGENY_MODEL, units);
		
		this.duplicationRateParameter = duplicationRateParameter;
		addVariable(duplicationRateParameter);
		duplicationRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
		
		this.hostShiftRateParameter = hostShiftRateParameter;
		addVariable(hostShiftRateParameter);
		hostShiftRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
		
		this.lossRateParameter = lossRateParameter;
		addVariable(lossRateParameter);
		lossRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
		
	}
	
	public double getDuplicationRate() {
		return duplicationRateParameter.getParameterValue(0);
	}
	
	public double getHostShiftRate() {
		return hostShiftRateParameter.getParameterValue(0);
	}
	
	public double getLossRate() {
		return lossRateParameter.getParameterValue(0);
	}
		
	/**
	 * Calculates the probability of a particular cophylogenetic mapping at a
	 * node and its children against this model's current state.
	 * 
	 * @param hostTree host tree
	 * @param selfHost host of the given node
	 * @param child1Host host of its 1st child
	 * @param child2Host host of its 2nd child
	 * @param selfBranchTime branch length/time of the given node
	 * @param child1BranchTime branch length/time of its 1st child
	 * @param child2BranchTime branch length/time of its 2nd child
	 * @return log likelihood of the cophylogenetic mapping at node
	 */
	@Override
	public double calculateNodeLogLikelihood(Tree hostTree,
			NodeRef selfHost, NodeRef child1Host, NodeRef child2Host,
			double selfBranchTime, double child1BranchTime, double child2BranchTime) {
		
		if (selfHost == null) { // Self has no host
			
			// Treat any new colonizations as host-shifts
			double logL = 0.0;
			if (child1Host != null)
				logL += calculateEventLogLikelihood(false, true, 0, child1BranchTime);
			if (child2Host != null)
				logL += calculateEventLogLikelihood(false, true, 0, child2BranchTime);
			return logL;
			
		} else { // Self has a host
			
			// Treat "loss of interest events" also as host shifts
			NodalRelationship child1Relationship = Relationship.determineRelationship(hostTree, selfHost, child1Host);
			NodalRelationship child2Relationship = Relationship.determineRelationship(hostTree, selfHost, child2Host);
			
			if (child1Relationship.relationship == Relationship.DESCENDANT
					&& child2Relationship.relationship == Relationship.DESCENDANT) {

				double logL = 0.0;
				
				NodeRef selfHostChild1 = hostTree.getChild(selfHost, 0);
				if (Relationship.determineRelationship(hostTree, selfHostChild1, child1Host).relationship == Relationship.determineRelationship(hostTree, selfHostChild1, child2Host).relationship) {
					
					// Inherent cospeciation, then a loss and host-shift along one linneage (unknown), and then potentially additional losses along both linneages
					// Because we do not know which linneage experienced the loss and host-shift, we sum over both
					// TODO For now this is an expensive call; maybe it can be fixed?
					logL += Math.log(
							Math.exp(calculateEventLogLikelihood(false, true, 1, child1BranchTime)) +
							Math.exp(calculateEventLogLikelihood(false, true, 1, child2BranchTime))
							);
				}
					
				// Inherent cospeciation, then potential losses along both linneages
				// Note that where there are no losses we have case of simple cospeciation
				logL += calculateEventLogLikelihood(false, false, child1Relationship.generations, child1BranchTime) +
						calculateEventLogLikelihood(false, false, child2Relationship.generations, child2BranchTime);

				return logL;
				
			} else if (child1Relationship.relationship == Relationship.SELF
					&& child2Relationship.relationship == Relationship.SELF) {
				
				// Duplication event
				return calculateEventLogLikelihood(true, false, 0, selfBranchTime);

			} else if (child1Relationship.relationship == Relationship.SELF && (child2Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {

				// Child2 host-shift event
				return calculateEventLogLikelihood(false, true, 0, child2BranchTime);

			} else if (child2Relationship.relationship == Relationship.SELF && (child1Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {
				
				// Child1 host-shift event
				return calculateEventLogLikelihood(false, true, 0, child1BranchTime);
				
			} else if ((child1Relationship.relationship == Relationship.COUSIN || child1Relationship.relationship == Relationship.SISTER)
					&& (child2Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {

				// Inherent cospeciation, then host-shift and loss events along each child lineage
				return calculateEventLogLikelihood(false, true, 1, child1BranchTime)
						+ calculateEventLogLikelihood(false, true, 1, child2BranchTime);

			} else if (child1Relationship.relationship == Relationship.DESCENDANT
					&& (child2Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {

				// Inherent cospeciation, then losses along both lineages and host-shift along child2 lineage
				return calculateEventLogLikelihood(false, true, 1, child2BranchTime)
						+ calculateEventLogLikelihood(false, false, child1Relationship.generations, child1BranchTime);

			} else if ((child1Relationship.relationship == Relationship.COUSIN || child1Relationship.relationship == Relationship.SISTER)
					&& child2Relationship.relationship == Relationship.DESCENDANT) {

				// Inherent cospeciation, then losses along both lineages and host-shift along child1 lineage
				return calculateEventLogLikelihood(false, true, 1, child1BranchTime)
						+ calculateEventLogLikelihood(false, false, child2Relationship.generations, child2BranchTime);

			} else { // Everything else is impossible
				return Double.NEGATIVE_INFINITY;
			}
		}
	}
	
//	/**
//	 * Calculates the log likelihood of a series of cophylogenetic events occurring over a period of time.
//	 * <p/>
//	 * Probability of n events with rate \lambda over t time calculated as follows (poisson process):<p/>
//	 * P(n | \lambda, t) = (\lambda * t) ^ n * e ^ (-\lambda * t) / n! 
//	 * 
//	 * @param nDuplications # of duplication events
//	 * @param nHostShifts # of host shift events
//	 * @param nLosses # of loss events
//	 * @param branchTime branch length/time
//	 * @return log likelihood of 
//	 */
//	private double calculateEventLogLikelihood(int nDuplications, int nHostShifts, int nLosses, double branchTime) {
//		
//		double logL = 0.0;
//		double lambdaxt;
//		
//		// Probably to deal with root?
//		if (branchTime < 0) return logL;
//		
//		if (nDuplications > 0) {
//			lambdaxt = getDuplicationRate() * branchTime;
//			logL += Math.log(Math.pow(lambdaxt, nDuplications));
//			logL -= lambdaxt;
//			logL -= MathUtils.factorialLog(nDuplications);
//		}
//		
//		if (nHostShifts > 0) {
//			lambdaxt = getHostShiftRate() * branchTime;
//			logL += Math.log(Math.pow(lambdaxt, nHostShifts));
//			logL -= lambdaxt;
//			logL -= MathUtils.factorialLog(nHostShifts);
//		}
//		
//		if (nLosses > 0) {
//			lambdaxt = getLossRate() * branchTime;
//			logL += Math.log(Math.pow(lambdaxt, nLosses));
//			logL -= lambdaxt;
//			logL -= MathUtils.factorialLog(nLosses);
//		}
//				
//		return logL;
//		
//	}

	/**
	 * Calculates the log likelihood of a series of cophylogenetic events occurring over a period of time.
	 * <p/>
	 * Probability of n events with rate \lambda over t time calculated as follows:<p/>
	 * P(n | \lambda, t) = (\lambda * t) ^ n * e ^ (-\lambda * t) / n! 
	 * <p/>
	 * This implementation should be more efficient, but cannot handle more than one each of duplication and host-shift events.
	 * 
	 * @param involvesDuplication involves duplication event
	 * @param involvesHostShift involves host-shift event
	 * @param nLosses # of loss events
	 * @param branchTime branch length/time
	 * @return log likelihood of events in time
	 */
	private double calculateEventLogLikelihood(boolean involvesDuplication, boolean involvesHostShift, int nLosses, double branchTime) {
		
		double logL = 0.0;
		double lambdaxt; // Holds lambda * t
		
		if (branchTime < 0) return logL; // To deal with root branch
		
		if (involvesDuplication) {
			lambdaxt = getDuplicationRate() * branchTime;
			logL += Math.log(lambdaxt) - lambdaxt;
		}
		
		if (involvesHostShift) {
			lambdaxt = getHostShiftRate() * branchTime;
			logL += Math.log(lambdaxt) - lambdaxt;
		}
		
		if (nLosses > 0) {
			lambdaxt = getLossRate() * branchTime;
			logL += Math.log(Math.pow(lambdaxt, nLosses));
			logL -= lambdaxt;
			logL -= MathUtils.factorialLog(nLosses);
		}
				
		return logL;
		
	}

	
}
