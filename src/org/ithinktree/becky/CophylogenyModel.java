/**
 * CophylogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.becky;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

/**
 * An abstract class to describe a cophylogenetic model.
 * 
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public abstract class CophylogenyModel extends AbstractModel implements Units {

	private Units.Type units;
	protected double overallRate;
	protected boolean dirty = true;
	
	
	/**
	 * 
	 */
	public CophylogenyModel(String name, Units.Type units) {
		super(name);
		setUnits(units);
	}
		
	protected final double likelihoodEventAtTime(double t, double lambda) {
		return lambda * Math.exp(-overallRate * t);
	}
	
	protected final double likelihoodEventInTime(double t, double lambda) {
		return Math.exp((lambda - overallRate) * t) - Math.exp(-overallRate * t);
	}
	
	protected final double likelihoodNoEventsInTime(double t) {
		return Math.exp(overallRate * t);
	}
	
	protected static final double likelihoodNoEventInTime(double t, double lambda) {
		throw new UnsupportedOperationException();
//		return 1.0 - ExponentialDistribution.cdf(t, lambda);
	}
	
	protected static final double likelihoodEventsInTime(double t, double lambda, int k) {
		throw new UnsupportedOperationException();
//		final double tXlambda = t * lambda;
//		return Math.exp(-tXlambda) * Math.pow(tXlambda, k) / MathUtils.factorial(k);
	}
	
	public void setUnits(Units.Type u) {
		units = u;
	}
	
	public Units.Type getUnits() {
		return units;
	}
	
	protected void handleModelChangedEvent(Model model, Object object, int index) {};
	
	@SuppressWarnings("rawtypes")
	protected void handleVariableChangedEvent(Variable variable, int index, ChangeType type) {dirty = true;}
	
	protected void storeState() {};
	protected void restoreState() {};
	protected void acceptState() {};
	
	public static class Utils {
	
		public enum Relationship {
			SELF,
			DESCENDANT,
			ANCESTOR,
			SISTER,
			COUSIN;
		}
		
		public static final class NodalRelationship {
			public final Relationship relationship;
			public final int generations; // The number of generations separating the two members of a relationship
			public final NodeRef[] lostLineages;
			public NodalRelationship(Relationship r, int g) {
				this(r, g, EMPTY_NODE_REF_ARRAY);
			}
			
			public NodalRelationship(final Utils.Relationship r, final int g, final NodeRef[] ll) {
				relationship = r;
				generations = g;
				lostLineages = ll;
			}
			
			/**
			 * Determines the reciprocal of a relationship.
			 * @param r given relationship.
			 * @return the reciprocal relationship
			 */
			public final Utils.NodalRelationship reciprocal() {
				switch (this.relationship) {
				case DESCENDANT: return new NodalRelationship(Utils.Relationship.ANCESTOR, this.generations, this.lostLineages);
				case ANCESTOR: return new NodalRelationship(Utils.Relationship.DESCENDANT, this.generations, this.lostLineages);
				// SELF, SISTER, or COUSIN should return the same
				default: return this;
				}
			}
		}

			/**
			 * Determines the relationship between two given nodes in a given tree.
			 * Does not 
			 * @param tree The tree
			 * @param self 
			 * @param relation The node 
			 * @return
			 */
			public static final NodalRelationship determineRelationship(final Tree tree, final NodeRef self, final NodeRef relation) {
				
				if (self == null || relation == null)
//					return new NodalRelationship(COUSIN, 0);
					throw new IllegalArgumentException();
					
				int selfN = self.getNumber();
				int relationN = relation.getNumber();
				
				if (selfN == relationN)
					return new NodalRelationship(Utils.Relationship.SELF, 0);
							
				if (!tree.isRoot(self) && !tree.isRoot(relation) && tree.getParent(self).getNumber() == tree.getParent(relation).getNumber())
					return new NodalRelationship(Utils.Relationship.SISTER, 0);
				
				List<NodeRef> lostLineages = new ArrayList<NodeRef>();
				
				NodeRef temp = tree.getParent(relation);
				// Arguably g should = 1, but confounds proper counting of loss events
				for (int g = 0; temp != null; ++g, temp = tree.getParent(temp)) {
					lostLineages.addAll(getSisters(tree, temp));
					if (selfN == temp.getNumber()) {
						return new NodalRelationship(Utils.Relationship.DESCENDANT, g, lostLineages.toArray(EMPTY_NODE_REF_ARRAY));
					}
				}
				
				lostLineages.clear();
				temp = tree.getParent(self);
				for (int g = 0; temp != null; ++g, temp = tree.getParent(temp)) {
					lostLineages.addAll(getSisters(tree, temp));
					if (relationN == temp.getNumber()) {
						return new NodalRelationship(Utils.Relationship.ANCESTOR, g, lostLineages.toArray(EMPTY_NODE_REF_ARRAY));
					}
				}
				
				// Otherwise must be a cousin
				return new NodalRelationship(Utils.Relationship.COUSIN, 0);
			}
			
			private static final NodeRef[] EMPTY_NODE_REF_ARRAY = new NodeRef[0];
			private static final List<NodeRef> EMPTY_NODE_LIST = new ArrayList<NodeRef>(0);
			private static final List<NodeRef> getSisters(Tree t, NodeRef n) {
				if (t.isRoot(n)) return EMPTY_NODE_LIST;
				NodeRef p = t.getParent(n);
				int cc = t.getChildCount(p);
				List<NodeRef> sisters = new ArrayList<NodeRef>(cc);
				for (int i = 0; i < cc; ++i) {
					NodeRef c = t.getChild(p, i);
					if (c != n) sisters.add(c);
				}
				return sisters;
			}
			
		public static final NodeRef[] lostLineagesToTime(final Tree t, NodeRef n, final double d) {
			
			if (t.getNodeHeight(n) > d) return EMPTY_NODE_REF_ARRAY;
			List<NodeRef> lostLineages = new ArrayList<NodeRef>();
			while (t.getNodeHeight(n) < d && !t.isRoot(n)) {
				lostLineages.addAll(getSisters(t, n));
				n = t.getParent(n);
			}
//			lostLineages.add(n);
			return lostLineages.toArray(EMPTY_NODE_REF_ARRAY);
			
		}
		
		public static final boolean isTreeValid(Tree t) {
			return isTreeValid(t, null, t.getRoot());
		}
		
		private static final boolean isTreeValid(Tree t, NodeRef p, NodeRef n) {
			
			boolean valid = p != null ? t.getNodeHeight(n) < t.getNodeHeight(p) : true;
			if (!valid) return valid;
			if (!t.isExternal(n)) {
				for (int i = 0; i < t.getChildCount(n); ++i) {
					valid = isTreeValid(t, n, t.getChild(n, i));
					if (!valid) return valid;
				}
			}
			return true;
		}
		
	}

	public abstract double calculateNodeLogLikelihood(final MutableTree symbiontTree, final NodeRef self,
			final NodeRef child1, final NodeRef child2, final Tree hostTree, final NodeRef selfHost,
			final NodeRef child1Host, final NodeRef child2Host, final BranchRates branchRates);

}

