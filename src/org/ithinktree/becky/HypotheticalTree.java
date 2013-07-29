/**
 * BECKY - Bayesian Estimation of Coevolution KrYteria
 * 
 * HypotheticalTree.java
 * 
 */
package org.ithinktree.becky;

import dr.evolution.tree.Tree;

/**
 * @author Arman D. Bilge
 *
 */
public interface HypotheticalTree extends Tree {

	public Tree getObservedTree();
	
}
