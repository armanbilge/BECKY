/**
 * CophylogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.becky;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.speciation.SpeciationModel;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

/**
 * An abstract class to describe a cophylogenetic model.
 * 
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public abstract class CophylogenyModel extends SpeciationModel {

	protected double overallRate;
	protected boolean dirty = true;

	/**
	 * 
	 */
	public CophylogenyModel(String name, Units.Type units) {
		super(name, units);
	}
	
	public final double getOverallRate() {
	    return overallRate;
	}
	
	protected abstract void updateVariables();
	
	protected final double likelihoodEventAtTime(final double t, final double lambda, final double rate) {
		return lambda * rate * Math.exp(-overallRate * rate * t);
	}
	
	protected final double likelihoodEventInTime(final double t, final double lambda, final double rate) {
	    final double adjustedOverallRate = rate * overallRate;
		return lambda * rate * (1 - Math.exp(-adjustedOverallRate * t)) / adjustedOverallRate;
	}
	
	protected final double likelihoodNoEventsInTime(final double t, final double rate) {
		return Math.exp(-overallRate * rate * t);
	}
	
	@SuppressWarnings("rawtypes")
	protected void handleVariableChangedEvent(Variable variable, int index, ChangeType type) {dirty = true;}
		
	protected void storeState() {dirty = true;}
	protected void restoreState() {dirty = true;}
	
	public static class Utils {
	
		public enum Relationship {
			SELF,
			DESCENDANT,
			ANCESTOR,
			SISTER,
			COUSIN;
		}
		
		public static class NodalRelationship {
			public final Relationship relationship;
			/** The number of generations separating the two members of a relationship */
			public final int generations;
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
					
				final int selfN = self.getNumber();
				final int relationN = relation.getNumber();
				
				if (selfN == relationN)
					return new NodalRelationship(Utils.Relationship.SELF, 0);
							
				if (!tree.isRoot(self) && !tree.isRoot(relation) && tree.getParent(self).getNumber() == tree.getParent(relation).getNumber())
					return new NodalRelationship(Utils.Relationship.SISTER, 0);
				
				List<NodeRef> lostLineages = new ArrayList<NodeRef>();
				
				NodeRef temp = tree.getParent(relation);
				lostLineages.add(relation);
				for (int g = 1; temp != null; ++g, temp = tree.getParent(temp)) {
					if (selfN == temp.getNumber()) {
						return new NodalRelationship(Utils.Relationship.DESCENDANT, g, lostLineages.toArray(EMPTY_NODE_REF_ARRAY));
					}
					lostLineages.addAll(getSisters(tree, temp));
				}
				
				lostLineages.clear();
				
				temp = tree.getParent(self);
				lostLineages.add(self);
				for (int g = 1; temp != null; ++g, temp = tree.getParent(temp)) {
					if (relationN == temp.getNumber()) {
						return new NodalRelationship(Utils.Relationship.ANCESTOR, g, lostLineages.toArray(EMPTY_NODE_REF_ARRAY));
					}
                lostLineages.addAll(getSisters(tree, temp));

				}
				
				// Otherwise must be a cousin
				return new NodalRelationship(Utils.Relationship.COUSIN, 0);
			}
			
			private static final NodeRef[] EMPTY_NODE_REF_ARRAY = new NodeRef[0];
			private static final List<NodeRef> EMPTY_NODE_LIST = new ArrayList<NodeRef>(0);
			public static final List<NodeRef> getSisters(Tree t, NodeRef n) {
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
			
		
			/**
			 * Returns the lost lineages up to d, inclusive.
			 * @param t tree
			 * @param n base of the lineage we are following
			 * @param d inclusive upper limit
			 * @return
			 */
			public static final NodeRef[] lostLineagesToTime(final Tree t, NodeRef n, final double d) {
				
				List<NodeRef> lostLineages = new ArrayList<NodeRef>();
				while (!t.isRoot(n) && t.getNodeHeight(t.getParent(n)) <= d) {
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
		
		public static final int getContemporaneousLineageCount(final Tree tree, final double height, final boolean approachFromPresent) {
			if (approachFromPresent) {
				return getContemporaneousLineageCountApproachingFromPresent(tree, tree.getRoot(), height);
			} else {
				return getContemporaneousLineageCountApproachingFromPast(tree, tree.getRoot(), height);
			}
		}

		public static final int getContemporaneousLineageCount(final Tree tree, final double height) {
			return getContemporaneousLineageCount(tree, height, false);
		}
		
		private static final int getContemporaneousLineageCountApproachingFromPast(final Tree tree, final NodeRef node, final double height) {
			if (height >= tree.getNodeHeight(node) && (tree.isRoot(node) || (tree.getNodeHeight(tree.getParent(node)) > height))) {
				return 1;
			} else {
				int count = 0;
				for (int i = 0; i < tree.getChildCount(node); ++i)
					count += getContemporaneousLineageCountApproachingFromPast(tree, tree.getChild(node, i), height);
				return count;
			}
		}

		private static final int getContemporaneousLineageCountApproachingFromPresent(final Tree tree, final NodeRef node, final double height) {
			if (height > tree.getNodeHeight(node) && (tree.isRoot(node) || (tree.getNodeHeight(tree.getParent(node)) >= height))) {
				return 1;
			} else {
				int count = 0;
				for (int i = 0; i < tree.getChildCount(node); ++i)
					count += getContemporaneousLineageCountApproachingFromPast(tree, tree.getChild(node, i), height);
				return count;
			}
		}
		
//		public static final Set<NodeRef> getContemporaneousLineages(final Tree tree, final double height) {
//			Set<NodeRef> lineages = new HashSet<NodeRef>(tree.getExternalNodeCount());
//			for (int i = 0; i < tree.getNodeCount(); ++i) {
//				NodeRef node = tree.getNode(i);
//				if (isContemporaneous(tree, node, height)) lineages.add(node);
//			}
//			return lineages;
//		}
		
		public static final Set<NodeRef> getContemporaneousLineages(final Tree tree, final double height) {
			Set<NodeRef> lineages = new HashSet<NodeRef>(tree.getExternalNodeCount());
			getContemporaneousLineages(tree, tree.getRoot(), height, lineages);
			return lineages;
		}
	
		private static final void getContemporaneousLineages(final Tree tree, final NodeRef node, final double height, Set<NodeRef> lineages) {
			if (isContemporaneous(tree, node, height)) {
				lineages.add(node);
			} else {
				for (int i = 0; i < tree.getChildCount(node); ++i)
					getContemporaneousLineages(tree, tree.getChild(node, i), height, lineages);
			}
		}

        public static final List<NodeRef> getLineagesInTimeRange(final Tree tree, final double startHeight, final double stopHeight) {
            List<NodeRef> lineages = new ArrayList<NodeRef>(tree.getExternalNodeCount());
            getLineagesInTimeRange(tree, tree.getRoot(), startHeight, stopHeight, lineages);
            return lineages;
        }
    
        private static final void getLineagesInTimeRange(final Tree tree, final NodeRef node, final double startHeight, final double stopHeight, final List<NodeRef> lineages) {
            // requires that startHeight > stopHeight
            if (startHeight >= tree.getNodeHeight(node)) {
                if (tree.isRoot(node) || (tree.getNodeHeight(tree.getParent(node)) > stopHeight)) {
                    lineages.add(node);
                } else {
                    return;
                }
            }
            for (int i = 0; i < tree.getChildCount(node); ++i)
                getLineagesInTimeRange(tree, tree.getChild(node, i), startHeight, stopHeight, lineages);
            
        }

        public static final boolean isContemporaneous(final Tree tree, final NodeRef node, final double height) {
            return height >= tree.getNodeHeight(node) && (tree.isRoot(node) || (tree.getNodeHeight(tree.getParent(node)) > height));
        }
        
	}

	public abstract double calculateNodeLogLikelihood(final MutableTree symbiontTree, final NodeRef self,
			final NodeRef child1, final NodeRef child2, final Tree hostTree, final NodeRef selfHost,
			final NodeRef child1Host, final NodeRef child2Host, final BranchRates branchRates);
	
	public abstract double calculateOriginLogLikelihood(final Tree symbiontTree, final double originHeight, final NodeRef root, final Tree hostTree, final NodeRef originHost, final NodeRef rootHost, final BranchRates branchRates);
	
	public abstract void initialize(final Tree tree);
	
	@Override
	public double calculateTreeLogLikelihood(Tree arg0) {
		throw new UnsupportedOperationException();
	}

	@Override
	public double calculateTreeLogLikelihood(Tree arg0, Set<Taxon> arg1) {
		throw new UnsupportedOperationException();
	}

	
}

