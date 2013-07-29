/**
 * 
 */
package org.ithinktree.becky;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;

/**
 * @author armanbilge
 *
 */
@SuppressWarnings("serial")
public class SimpleHypotheticalTree extends SimpleTree implements
		HypotheticalTree {

	/**
	 * 
	 */
	public SimpleHypotheticalTree() {
		super();
	}

	/**
	 * @param tree
	 */
	public SimpleHypotheticalTree(Tree tree) {
		super(tree);
	}

	/**
	 * @param node
	 */
	public SimpleHypotheticalTree(SimpleNode node) {
		super(node);
	}

	/* (non-Javadoc)
	 * @see org.ithinktree.becky.HypotheticalTree#getObservedTree()
	 */
	@Override
	public Tree getObservedTree() {
		return new SimpleTree(getObservedLineage(getRoot()));
	}
	
	private SimpleNode getObservedLineage(NodeRef nodeRef) {
	
		SimpleNode node = new SimpleNode(this, nodeRef);
		while(node.hasChildren()) node.removeChild(0);
		
		if (isExternal(nodeRef)) {
			
			if (node.getTaxon() == null) return null;
			
		} else {
			
			SimpleNode child1 = getObservedLineage(getChild(nodeRef, 0));
			SimpleNode child2 = getObservedLineage(getChild(nodeRef, 1));
			
			if (child1 != null && child2 != null) {
				// Both lineages observed
				node.addChild(child1);
				node.addChild(child2);
			} else if (child1 == child2) {
				// Both are null, hence this entire lineage was not observed
				return null;
			} else {
				// Just one child lineage is observed
				return child2 == null ? child1 : child2; // Set this node to the continuing child lineage
			}
			
		}
		
		return node;
	}
		
		
}
