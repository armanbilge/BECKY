/**
 * SimpleCophylogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.becky;

import java.util.Arrays;
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
	
	protected final double likelihoodHostShiftAtTime(final double t) {
		return likelihoodEventAtTime(t, hostShiftRate);
	}
	
	protected final double likelihoodHostShiftInTime(final double t) {
		return likelihoodEventInTime(t, getHostShiftRate());
	}

	public double getLossRate() {
		return lossRateParameter.getParameterValue(0);
	}
	
	protected final double likelihoodLossInTime(final double t) {
		return likelihoodEventInTime(t, getLossRate());
	}
	
	/**
	 * Evaluates an integral
	 * @param a
	 * @param b
	 * @param t
	 * @param rate
	 * @return
	 */
	protected final double likelihoodHostShiftAndLossInTime(final double a, final double b, final double t, final double rate) {
		final double hostShiftRate = this.hostShiftRate * rate;
		final double lossRate = this.lossRate * rate;
		final double overallRate = this.overallRate * rate;
		return hostShiftRate * lossRate / overallRate * ((Math.exp(-overallRate * a) - Math.exp(-overallRate * b)) / overallRate - Math.exp(-overallRate * t) * (b - a));
	}

	protected final double likelihoodLossesAlongLineages(final Tree tree, final NodeRef[] lineages, double rate) {
		double likelihood = 1.0;
		for (NodeRef n : lineages)
			likelihood *= likelihoodLossInTime(tree.getBranchLength(n) * rate);
		return likelihood;
	}
	
	/**
	 * Discretizes a branch along which a host-shift and loss happened to properly determine the likelihood of these events.
	 * 
	 * @param start the earliest the host-shift could have happened
	 * @param hostShiftStop the latest the host-shift could have happened
	 * @param lossStop the latest the loss could have happened
	 * @param rate the rate on the relevant branch
	 * @param tree the host tree
	 * @param originalLineages lineages lost on original host lineage
	 * @param newHostLineages lineages lost on new host lineage
	 * @return
	 */
	protected final double likelihoodHostShiftAndLossInTime(final double start, final double hostShiftStop, final double lossStop, final double rate, final Tree tree, final NodeRef[] originalLineages, final NodeRef[] newHostLineages) {
		
		final double t = start - lossStop;
		final double stop = start - Math.max(hostShiftStop, lossStop);
		double likelihood = 0;
		double height;
		double subHeight;
		double nextHeight = 0.0;
		double nextSubHeight;
		for (int i = originalLineages.length - 1; nextHeight < stop; --i) {
			height = nextHeight;
			if (i >= 0) nextHeight = Math.min(start - tree.getNodeHeight(tree.getParent(originalLineages[i])), stop);
			else nextHeight = stop;
			nextSubHeight = height;
			for (int j = newHostLineages.length - 1; nextSubHeight < nextHeight; --j) {
				subHeight = nextSubHeight;
				if (j >= 0) nextSubHeight = Math.min(start - tree.getNodeHeight(tree.getParent(newHostLineages[j])), nextHeight);
				else nextSubHeight = nextHeight;
				if (likelihoodHostShiftAndLossInTime(subHeight, nextSubHeight, t, rate) <= 0) throw new RuntimeException(subHeight + " " + nextSubHeight + " " + t + " " + rate);
				likelihood += likelihoodHostShiftAndLossInTime(subHeight, nextSubHeight, t, rate) *
						likelihoodLossesAlongLineages(tree, Arrays.copyOfRange(originalLineages, i+1, originalLineages.length), rate) *
						likelihoodLossesAlongLineages(tree, Arrays.copyOfRange(newHostLineages, 0, j+1), rate);
			}
		}
		
		if (likelihood <= 0 || Double.isNaN(likelihood)) {
			System.out.println(start);
			System.out.println(hostShiftStop);
			System.out.println(lossStop);
			System.out.println(rate);
			System.out.println(tree);
			System.out.println(Arrays.toString(originalLineages));
			System.out.println(Arrays.toString(newHostLineages));
			System.out.println(likelihood);
			System.exit(1);
		}
		
		return likelihood;
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
							
							// Determine along which child lineage the loss(es) happened
							if (nr1.relationship != Relationship.DESCENDANT)
								hostChild = hostTree.getChild(selfHost, 1);
							final double hostChildBranchLength = hostTree.getBranchLength(hostChild);
							final double hostChildHeight = hostTree.getNodeHeight(hostChild);
								
							double case1 = 1.0;
							double case2 = 1.0;
							
							// Case 1: duplication, cospeciation, then losses

							// The duplication occurred at the time of their last common ancestor
							case1 *= likelihoodDuplicationAtTime(selfBranchLength * selfBranchRate);
							
							final double potentialLossLength = (selfHeight - selfHostHeight) + hostChildBranchLength;
							case1 *= likelihoodLossInTime(potentialLossLength * child1BranchRate) * likelihoodLossInTime(potentialLossLength * child2BranchRate);
														
							// Case 2: cospeciation, then host-shift and loss
							
							final double child1Height = symbiontTree.getNodeHeight(child1);
							final double child2Height = symbiontTree.getNodeHeight(child2);
							
							// All possible lineages along which losses may have occurred
							final NodeRef[] child1OriginalHostLineages = Utils.lostLineagesToTime(hostTree, hostChild, selfHeight);
							final NodeRef[] child2OriginalHostLineages = Utils.lostLineagesToTime(hostTree, hostChild, selfHeight);
							final NodeRef[] child1NewHostLineages = child1Relationship.lostLineages; // Utils.lostLineagesToTime(hostTree, child1Host, selfHeight);
							final NodeRef[] child2NewHostLineages = child2Relationship.lostLineages; // Utils.lostLineagesToTime(hostTree, child2Host, selfHeight);
							
							case2 *= likelihoodNoEventsInTime(selfBranchLength * selfBranchRate);
							System.out.println("over here");
							// Sum over two subcases: child1 lineage made host-shift/loss or child2 made host-shift/loss
							case2 *= likelihoodLossesAlongLineages(hostTree, child1NewHostLineages, child1BranchRate) *
										likelihoodHostShiftAndLossInTime(selfHeight, child2Height, hostChildHeight, child2BranchRate, hostTree, child2OriginalHostLineages, child2NewHostLineages) +
										likelihoodLossesAlongLineages(hostTree, child2NewHostLineages, child2BranchRate) +
										likelihoodHostShiftAndLossInTime(selfHeight, child1Height, hostChildHeight, child1BranchRate, hostTree, child1OriginalHostLineages, child1NewHostLineages);

							likelihood *= case1 + case2;
							
						} else { // Plain old cospeciation
														
							// Check if violates tree validity (parent younger than children)
//							if (selfHostHeight < symbiontTree.getNodeHeight(child1) || selfHostHeight < symbiontTree.getNodeHeight(child2))
//								return Double.NEGATIVE_INFINITY;
//							symbiontTree.setNodeHeight(self, selfHostHeight);
							
							likelihood *= likelihoodNoEventsInTime(symbiontTree.getBranchLength(self) * selfBranchRate);
							
							// Potential losses along both child lineages
							likelihood *= likelihoodLossesAlongLineages(hostTree, child1Relationship.lostLineages, child1BranchRate);
							likelihood *= likelihoodLossesAlongLineages(hostTree, child2Relationship.lostLineages, child2BranchRate);							
						}
						
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

						// Double host-shift event with a loss: no symbionts left on this host lineage
						
						likelihood *= likelihoodHostShiftAtTime(selfBranchLength * selfBranchRate);
						
						final double child1Height = symbiontTree.getNodeHeight(child1);
						final double child1BranchRate = branchRates.getBranchRate(symbiontTree, child1);
						final double child2Height = symbiontTree.getNodeHeight(child2);
						final double child2BranchRate = branchRates.getBranchRate(symbiontTree, child2);
						
						final NodeRef[] noLineages = new NodeRef[0];
						final NodeRef[] child1NewHostLineages = Utils.lostLineagesToTime(hostTree, child1Host, selfHeight);
						final NodeRef[] child2NewHostLineages = Utils.lostLineagesToTime(hostTree, child2Host, selfHeight);
						
						// We definitely know the time of the first host-shift
						likelihood *= likelihoodHostShiftAtTime(selfBranchLength * selfBranchRate);
						
						// Case 1: Child2 lineage host-shifted first
						double case1 = likelihoodLossesAlongLineages(hostTree, child1NewHostLineages, child1BranchRate);
						// Case 2: Child1 lineage host-shifted first
						double case2 = likelihoodLossesAlongLineages(hostTree, child2NewHostLineages, child2BranchRate);
						System.out.println("actually over here");
						case1 *= likelihoodHostShiftAndLossInTime(selfHeight, child2Height, selfHostHeight, child2BranchRate, hostTree, noLineages, child2NewHostLineages);
						case2 *= likelihoodHostShiftAndLossInTime(selfHeight, child1Height, selfHostHeight, child1BranchRate, hostTree, noLineages, child1NewHostLineages);
						
//						// Determine the latest time that the second host shift could have happened
//						final double a = hostTree.getNodeHeight(selfHost);
//						final double b = symbiontTree.getNodeHeight(child1);
//						final double c = symbiontTree.getNodeHeight(child2);
//						final double d = Math.max(a, b);
//						double e = Math.max(a, c);
//						if (d == e) { // Host speciates first
//							case1 *= likelihoodHostShiftAndLossInTime(selfHeight, a, a, child1BranchRate, hostTree, noLineages, child1NewHostLineages);
//							case2 *= likelihoodHostShiftAndLossInTime(selfHeight, a, a, child2BranchRate, hostTree, noLineages, child2NewHostLineages);
//						} else { // Host-shifting lineage speciates first
//							boolean isChild1 = Math.max(b, c) == b;
//							if (isChild1)
//								case1 *= likelihoodHostShiftAndLossInTime(selfHeight, a, a, child1BranchRate, hostTree, noLineages, child1NewHostLineages);
//							else
//								case2 *= likelihoodHostShiftAndLossInTime(selfHeight, a, a, child2BranchRate, hostTree, noLineages, child2NewHostLineages);
//							
//							e = Math.max(a, (isChild1 ? c : b));
//							likelihood *= likelihoodHostShiftAndLossInTime((selfHeight - e) * (e == a ? selfBranchRate : (isChild1 ? child2BranchRate : child1BranchRate)));
//						}
//						// Note that again I am being lazy about considering losses after the host shift; too difficult to integrate that over time
						
						likelihood *= case1 + case2;
						
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
						likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child1Host, selfHeight), branchRates.getBranchRate(symbiontTree, child1));
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
					
					if (Math.log(likelihood) == Double.NEGATIVE_INFINITY || Math.log(likelihood) == Double.POSITIVE_INFINITY || Double.isNaN(Math.log(likelihood))) {
						System.out.println(Math.log(likelihood));
						System.out.println(child1Relationship.relationship.toString());
						System.out.println(child2Relationship.relationship.toString());
						System.exit(1);
					}
					
					return Math.log(likelihood);
	}
	
	@Override
	public double calculateTreeLogLikelihood(Tree arg0) {
		throw new UnsupportedOperationException();
	}

	@Override
	public double calculateTreeLogLikelihood(Tree arg0, Set<Taxon> arg1) {
		throw new UnsupportedOperationException();
	}
	
}
