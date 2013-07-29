/**
 * 
 */
package org.ithinktree.becky;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.MutableTreeListener;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;

/**
 * @author armanbilge
 *
 */
@SuppressWarnings("serial")
public class ObservedTreeModel extends TreeModel {

	private final HypotheticalTree hypotheticalTree;
	@SuppressWarnings("unused")
	private Tree observedTree;
	private boolean isDirty = false;
	
	/**
	 * 
	 */
	public ObservedTreeModel(HypotheticalTree hypotheticalTree) {
		super(hypotheticalTree);
		this.hypotheticalTree = hypotheticalTree;
		observedTree = hypotheticalTree.getObservedTree();
		if (hypotheticalTree instanceof MutableTree)
			((MutableTree) hypotheticalTree).addMutableTreeListener(new MutableTreeListener() {
				@Override
				public void treeChanged(Tree tree) {
					isDirty = true;
				}
			});
	}

	@SuppressWarnings("unused")
	private void resolveDirtiness() {
		if (isDirty) {
			observedTree = hypotheticalTree.getObservedTree();
			isDirty = false;
		}
	}

}
