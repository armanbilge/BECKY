/**
 * SimpleCophylogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.becky;

import java.util.Set;

import org.ithinktree.becky.CophylogenyModel.Utils.NodalRelationship;
import org.ithinktree.becky.CophylogenyModel.Utils.Relationship;
import org.ithinktree.becky.xml.SimpleCophylogenyModelParser;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;

/**
 * A simple model for cophylogenetic mappings.
 * 
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public class SimpleCophylogenyModel extends CophylogenyModel {

	final protected Parameter duplicationRateParameter;
	final protected Parameter hostShiftRateParameter;
	final protected Parameter lossRateParameter;
	private double duplicationRate;
	private double hostShiftRate;
	private double lossRate;
		
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
	
	private void updateVariables() {
		duplicationRate = getDuplicationRate();
		hostShiftRate = getHostShiftRate();
		lossRate = getLossRate();
		overallRate = duplicationRate + hostShiftRate + lossRate;
		dirty = false;
	}
	
	public double getDuplicationRate() {
		return duplicationRateParameter.getParameterValue(0);
	}
	
	protected final double likelihoodDuplicationAtTime(double t) {
		return super.likelihoodEventAtTime(t, duplicationRate);
	}
	
	protected final double likelihoodDuplicationInTime(double t) {
		return likelihoodEventInTime(t, getDuplicationRate());
	}
	
	public double getHostShiftRate() {
		return hostShiftRateParameter.getParameterValue(0);
	}
	
	protected final double likelihoodHostShiftAtTime(double t) {
		return super.likelihoodEventAtTime(t, hostShiftRate);
	}
	
	protected final double likelihoodHostShiftInTime(double t) {
		return likelihoodEventInTime(t, getHostShiftRate());
	}

	public double getLossRate() {
		return lossRateParameter.getParameterValue(0);
	}
	
	protected final double likelihoodLossInTime(double t) {
		return likelihoodEventInTime(t, getLossRate());
	}
	
	protected final double likelihoodHostShiftAndLossInTime(double t) {
		final double twoXoverallRate = 2 * overallRate;
		return (hostShiftRate * (1 - Math.exp(-(twoXoverallRate - lossRate) * t))) / (twoXoverallRate - lossRate) - (hostShiftRate * (1 - Math.exp(twoXoverallRate * t))) / (twoXoverallRate);
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
	public double calculateNodeLogLikelihood(final MutableTree symbiontTree, final NodeRef self,
			final NodeRef child1, final NodeRef child2, final Tree hostTree, final NodeRef selfHost,
			final NodeRef child1Host, final NodeRef child2Host, final BranchRates branchRates) {
		
					if (dirty) updateVariables();
		
					double likelihood = 1.0;
					
					final NodalRelationship child1Relationship = Utils.determineRelationship(hostTree, selfHost, child1Host);
					final NodalRelationship child2Relationship = Utils.determineRelationship(hostTree, selfHost, child2Host);
					final double selfBranchLength = symbiontTree.getBranchLength(self);
					final double selfBranchRate = branchRates.getBranchRate(symbiontTree, self);
					double selfHeight = symbiontTree.getNodeHeight(self);
					final double selfHostHeight = hostTree.getNodeHeight(selfHost);

					// Check if symbiont coexisted temporally with its host
					if ((!hostTree.isRoot(selfHost) && selfHeight >= hostTree.getNodeHeight(hostTree.getParent(selfHost))) ||
							selfHeight < selfHostHeight)
						return Double.NEGATIVE_INFINITY;
										
					if (child1Relationship.relationship == Relationship.DESCENDANT
							&& child2Relationship.relationship == Relationship.DESCENDANT) {
						
						final double child1BranchRate = branchRates.getBranchRate(symbiontTree, child1);
						final double child2BranchRate = branchRates.getBranchRate(symbiontTree, child2);
						
						NodeRef hostChild = hostTree.getChild(selfHost, 0);
						final Utils.NodalRelationship nr1 = Utils.determineRelationship(hostTree, hostChild, child1Host);
						final Utils.NodalRelationship nr2 = Utils.determineRelationship(hostTree, hostChild, child2Host);
						if (nr1.relationship == nr2.relationship || (nr1.relationship == Relationship.SISTER && nr2.relationship == Relationship.COUSIN) || (nr2.relationship == Relationship.SISTER && nr1.relationship == Relationship.COUSIN)) {
							
							if (nr1.relationship != Relationship.DESCENDANT)
								hostChild = hostTree.getChild(selfHost, 1);
							final double hostChildBranchLength = hostTree.getBranchLength(hostChild);
								
							double case1 = 1.0;
							double case2 = 1.0;
							
							// Case 1: duplication, cospeciation, then losses

							// The duplication occurred at the time of their last common ancestor
							case1 *= likelihoodDuplicationAtTime(selfBranchLength * selfBranchRate);
							
							final double potentialLossLength = (selfHeight - selfHostHeight) + hostChildBranchLength;
							case1 *= likelihoodLossInTime(potentialLossLength * child1BranchRate) * likelihoodLossInTime(potentialLossLength * child2BranchRate);
							
							
							// Case 2: cospeciation, then host-shift and loss
							
							case2 *= likelihoodNoEventsInTime(selfBranchLength * selfBranchRate);
							case2 *= likelihoodHostShiftAndLossInTime(hostChildBranchLength * child1BranchRate) +
										likelihoodHostShiftAndLossInTime(hostChildBranchLength * child2BranchRate);
							// UH-OH: Approximating b/c not calculating likelihood losses along host-shift lineages b/c too difficult to integrate over
							// TODO Stop being lazy: can be done by partitioning integral for every point this likelihood changes using start and end times a and b
							// Note that the loss calculation below should hypothetically cover this
							
							likelihood *= case1 + case2;
							
						} else { // Plain old cospeciation
														
							// Check if violates tree validity (parent younger than children)
							if (selfHostHeight < symbiontTree.getNodeHeight(child1) || selfHostHeight < symbiontTree.getNodeHeight(child2))
								return Double.POSITIVE_INFINITY;
							symbiontTree.setNodeHeight(self, selfHostHeight);
							
							likelihood *= likelihoodNoEventsInTime(selfBranchLength * selfBranchRate);
							
						}
							
						// Potential losses along both child lineages
						likelihood *= likelihoodLossesAlongLineages(hostTree, child1Relationship.lostLineages, child1BranchRate);
						likelihood *= likelihoodLossesAlongLineages(hostTree, child2Relationship.lostLineages, child2BranchRate);

						
					} else if (child1Relationship.relationship == Relationship.SELF
							&& child2Relationship.relationship == Relationship.SELF) {
						
						// Duplication event
						likelihood *= likelihoodDuplicationAtTime(selfBranchLength * selfBranchRate);
						
					} else if (child1Relationship.relationship == Relationship.SELF && (child2Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {

						// Child2 host-shift event
						likelihood *= likelihoodHostShiftAtTime(selfBranchLength * selfBranchRate);
						likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child2Host, selfHeight), branchRates.getBranchRate(symbiontTree, child2));
					
					} else if (child2Relationship.relationship == Relationship.SELF && (child1Relationship.relationship == Relationship.COUSIN || child1Relationship.relationship == Relationship.SISTER)) {
						
						// Child1 host-shift event
						likelihood *= likelihoodHostShiftAtTime(selfBranchLength * selfBranchRate);
						likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child1Host, selfHeight), branchRates.getBranchRate(symbiontTree, child1));
						
					} else if ((child1Relationship.relationship == Relationship.COUSIN || child1Relationship.relationship == Relationship.SISTER)
							&& (child2Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {

						// Double host shift event with a loss: no symbionts left on this host lineage
						
						likelihood *= likelihoodHostShiftAtTime(selfBranchLength * selfBranchRate);
						
						final double child1BranchRate = branchRates.getBranchRate(symbiontTree, child1);
						final double child2BranchRate = branchRates.getBranchRate(symbiontTree, child2);
						
						likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child1Host, selfHeight), child1BranchRate) +
								likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child2Host, selfHeight), child2BranchRate);
						
						// Determine the latest time that the host shift could have happened
						final double a = hostTree.getNodeHeight(child2Host);
						final double b = symbiontTree.getNodeHeight(child1);
						final double c = symbiontTree.getNodeHeight(child2);
						final double d = Math.max(a, b);
						double e = Math.max(a, c);
						if (d == e) { // Host speciates first
							likelihood *= 2 * likelihoodHostShiftAndLossInTime((selfHeight - a) * selfBranchRate);
						} else { // Host-shifting lineage speciates first
							boolean isChild1 = Math.max(b, c) == b;
							likelihood *= likelihoodHostShiftAndLossInTime((selfHeight - (isChild1 ? b : c)) * (isChild1 ? child1BranchRate : child2BranchRate));
							e = Math.max(a, (isChild1 ? c : b));
							likelihood *= likelihoodHostShiftAndLossInTime((selfHeight - e) * (e == a ? selfBranchRate : (isChild1 ? child2BranchRate : child1BranchRate)));
						}
						// Note that again I am being lazy about considering losses after the host shift; too difficult to integrate that over time
						
					} else if (child1Relationship.relationship == Relationship.DESCENDANT
							&& (child2Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {

						// Child2 host-shift and child1 losses
						likelihood *= likelihoodHostShiftAtTime(selfBranchLength * selfBranchRate);
						likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child2Host, selfHeight), branchRates.getBranchRate(symbiontTree, child2));
						likelihood *= likelihoodLossesAlongLineages(hostTree, child1Relationship.lostLineages, branchRates.getBranchRate(symbiontTree, child1));

					} else if ((child1Relationship.relationship == Relationship.COUSIN || child1Relationship.relationship == Relationship.SISTER)
							&& child2Relationship.relationship == Relationship.DESCENDANT) {

						// Child1 host-shift and child2 losses 
						likelihood *= likelihoodHostShiftAtTime(selfBranchLength * selfBranchRate);
						likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child1, selfHeight), branchRates.getBranchRate(symbiontTree, child1));
						likelihood *= likelihoodLossesAlongLineages(hostTree, child2Relationship.lostLineages, branchRates.getBranchRate(symbiontTree, child2));

					} else { // Everything else is impossible
						return Double.NEGATIVE_INFINITY;
					}
					
					// Should already be checked above
//					selfHeight = symbiontTree.getNodeHeight(self);
//					if (!symbiontTree.isRoot(self) && 
//							selfHeight > symbiontTree.getNodeHeight(symbiontTree.getParent(self)))
//						return Double.NEGATIVE_INFINITY;
//					for (int i = 0; i < symbiontTree.getChildCount(self); ++i) {
//						if (selfHeight < symbiontTree.getNodeHeight(symbiontTree.getChild(self, i)))
//							return Double.NEGATIVE_INFINITY;
//					}
					
					return Math.log(likelihood);
	}
	
	protected double likelihoodLossesAlongLineages(final Tree tree, final NodeRef[] lineages, double rate) {
		double likelihood = 1.0;
		for (NodeRef n : lineages)
			likelihood *= likelihoodLossInTime(tree.getBranchLength(n) * rate);
		return likelihood;
	}

	@Override
	public double calculateTreeLogLikelihood(Tree arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double calculateTreeLogLikelihood(Tree arg0, Set<Taxon> arg1) {
		// TODO Auto-generated method stub
		return 0;
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

//	/**
//	 * Calculates the log likelihood of a series of cophylogenetic events occurring over a period of time.
//	 * <p/>
//	 * Probability of n events with rate \lambda over t time calculated as follows:<p/>
//	 * P(n | \lambda, t) = (\lambda * t) ^ n * e ^ (-\lambda * t) / n! 
//	 * <p/>
//	 * This implementation should be more efficient, but cannot handle more than one each of duplication and host-shift events.
//	 * 
//	 * @param involvesDuplication involves duplication event
//	 * @param involvesHostShift involves host-shift event
//	 * @param nLosses # of loss events
//	 * @param branchTime branch length/time
//	 * @return log likelihood of events in time
//	 */
//	private double calculateEventLogLikelihood(boolean involvesDuplication, boolean involvesHostShift, int nLosses, double branchTime) {
//		
//		double logL = 0.0;
//		double lambdaxt; // Holds lambda * t
//		
//		if (branchTime < 0) return logL; // To deal with root branch
//		
//		if (involvesDuplication) {
//			lambdaxt = getDuplicationRate() * branchTime;
//			logL += Math.log(lambdaxt) - lambdaxt;
//		}
//		
//		if (involvesHostShift) {
//			lambdaxt = getHostShiftRate() * branchTime;
//			logL += Math.log(lambdaxt) - lambdaxt;
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

	
}
