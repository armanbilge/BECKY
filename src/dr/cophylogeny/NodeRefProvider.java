/**
 * NodeRefProvider.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KryTeria
 * 
 */
package dr.cophylogeny;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;

/**
 * @author Arman D. Bilge
 *
 */
public class NodeRefProvider implements TreeTraitProvider {

	public static final String NODE_REF = "nodeRef";
	
	Tree tree;
	TreeTraitProvider.Helper treeTraits = new Helper();
	
	/**
	 * 
	 */
	public NodeRefProvider(Tree tree, final String tag) {
		
		this.tree = tree;
		
		treeTraits.addTrait(NODE_REF, new TreeTrait.I() {

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
