/**
 * NodeRefProvider.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KryTeria
 * 
 */
package org.ithinktree.becky;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;

/**
 * @author Arman D. Bilge
 *
 */
public class NodeRefProvider implements TreeTraitProvider {
	
	@SuppressWarnings("unused")
	private final Tree tree;
	private final TreeTraitProvider.Helper treeTraits = new Helper();
	
	/**
	 * 
	 */
	public NodeRefProvider(final Tree tree, final String tag) {
		
		this.tree = tree;
		
		treeTraits.addTrait(tag, new TreeTrait.I() {

			@Override
			public String getTraitName() {
				return tag;
			}

			@Override
			public Intent getIntent() {
				return Intent.NODE;
			}

			@Override
			public Integer getTrait(Tree tree, NodeRef node) {
				return node.getNumber();
			}
			
			public String getTraitString(Tree tree, NodeRef node) {
				return Integer.toString(node.getNumber());
			}
			
		});
		
	}

	@SuppressWarnings("rawtypes")
	@Override
	public TreeTrait[] getTreeTraits() {
		return treeTraits.getTreeTraits();
	}

	@SuppressWarnings("rawtypes")
	@Override
	public TreeTrait getTreeTrait(String key) {
		return treeTraits.getTreeTrait(key);
	}

}
