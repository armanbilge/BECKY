/**
 * CophylogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.becky;

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
	
	/**
	 * 
	 */
	public CophylogenyModel(String name, Units.Type units) {
		super(name);
		setUnits(units);
	}
	
	public double calculateNodeLogLikelihood(NodeRef self, Tree hostTree, NodeRef selfHost, NodeRef child1Host, NodeRef child2Host, double selfBranchTime, double child1BranchTime, double child2BranchTime) {
		return this.calculateNodeLogLikelihood(hostTree, selfHost, child1Host, child2Host, selfBranchTime, child1BranchTime, child2BranchTime);
	}
	public abstract double calculateNodeLogLikelihood(Tree hostTree, NodeRef selfHost, NodeRef child1Host, NodeRef child2Host, double selfBranchTime, double child1BranchTime, double child2BranchTime);
	
	public void setUnits(Units.Type u) {
		units = u;
	}
	
	public Units.Type getUnits() {
		return units;
	}
	
	protected void handleModelChangedEvent(Model model, Object object, int index) {};
	
	@SuppressWarnings("rawtypes")
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {}
	
	protected void storeState() {};
	protected void restoreState() {};
	protected void acceptState() {};
	
	public enum Relationship {
		SELF,
		DESCENDANT,
		ANCESTOR,
		SISTER,
		COUSIN;
		
		public static class NodalRelationship {
			public final Relationship relationship;
			public final int generations; // The number of generations separating the two members of a relationship
			public NodalRelationship(Relationship r, int g) {
				relationship = r;
				generations = g;
			}
		}
		
		/**
		 * Determines the reciprocal of a relationship.
		 * @param r given relationship.
		 * @return the reciprocal relationship
		 */
		public static NodalRelationship reciprocal(NodalRelationship r) {
			switch (r.relationship) {
			case DESCENDANT: return new NodalRelationship(ANCESTOR, r.generations);
			case ANCESTOR: return new NodalRelationship(DESCENDANT, r.generations);
			// SELF, SISTER, or COUSIN should return the same
			default: return r;
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
		public static NodalRelationship determineRelationship(Tree tree, NodeRef self, NodeRef relation) {
			
			if (self == null || relation == null)
				return new NodalRelationship(COUSIN, 0);
				
			int selfN = self.getNumber();
			int relationN = relation.getNumber();
			
			if (selfN == relationN)
				return new NodalRelationship(SELF, 0);
						
			if (!tree.isRoot(self) && !tree.isRoot(relation) && tree.getParent(self).getNumber() == tree.getParent(relation).getNumber())
				return new NodalRelationship(SISTER, 0);
			
			NodeRef descendant = tree.getParent(relation);
			// Arguably g should = 1, but confounds proper counting of loss events
			for (int g = 0; descendant != null; g++, descendant = tree.getParent(descendant)) {
				if (selfN == descendant.getNumber()) {
					return new NodalRelationship(DESCENDANT, g);
				}
			}
			
			NodeRef ancestor = tree.getParent(self);
			for (int g = 0; ancestor != null; g++, ancestor = tree.getParent(ancestor)) {
				if (relationN == ancestor.getNumber()) {
					return new NodalRelationship(ANCESTOR, g);
				}
			}
			
			// Otherwise must be a cousin
			// TODO Any potential to collect generational info, e.g. generations to MRCA?
			return new NodalRelationship(COUSIN, 0);
		}
		
	}


}
