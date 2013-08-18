package org.ithinktree.becky.jprimewrappers;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import se.cbb.jprime.topology.RBTree;
import se.cbb.jprime.topology.RootedTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

public class JPrIMERBTreeWrapperForBEASTTree extends RBTree implements WrappedBEASTObject<Tree> {

	private final Tree tree;
	
	public JPrIMERBTreeWrapperForBEASTTree(Tree tree) {
		super(tree.getId(), 1);
		this.tree = tree;
		root = getRoot();
	}

	@Override
	public int getLeftChild(int i) {
		return tree.getChild(tree.getNode(i), 0).getNumber();
	}

	@Override
	public int getRightChild(int i) {
		return tree.getChild(tree.getNode(i), 1).getNumber();
	}

	@Override
	public int getSibling(int i) {
		NodeRef n = tree.getNode(i);
		if (tree.isRoot(n)) return RootedTree.NULL;
		NodeRef p = tree.getParent(n);
		int c1 = tree.getChild(p, 0).getNumber();
		int c2 = tree.getChild(p, 1).getNumber();
		return n.getNumber() == c2 ? c1 : c2;
	}

	@Override
	public List<Integer> getAncestors(int i, boolean b) {
		List<Integer> ancestors = new ArrayList<Integer>();
		if (b) ancestors.add(i);
		for (NodeRef n = tree.getParent(tree.getNode(i)); n != null; n = tree.getParent(n))
			ancestors.add(n.getNumber());
		return ancestors;
	}

	@Override
	public List<Integer> getChildren(int i) {
		NodeRef n = tree.getNode(i);
		int childCount = tree.getChildCount(n);
		List<Integer> children = new ArrayList<Integer>(childCount);
		for (int j = 0; j < childCount; ++j)
			children.add(tree.getChild(n, j).getNumber());
		return children;
	}

	@Override
	public List<Integer> getDescendantLeaves(int i, boolean b) {
		List<Integer> descendants = new ArrayList<Integer>();
		if (!b) descendants.add(i);
		for (NodeRef n : Tree.Utils.getExternalNodes(tree, tree.getNode(i)))
			descendants.add(n.getNumber());
		return descendants;
	}

	@Override
	public List<Integer> getDescendants(int i, boolean b) {
		List<Integer> descendants = new ArrayList<Integer>();
		getDescendants(tree.getNode(i), descendants);
		if (b) descendants.remove(0);
		return descendants;
	}

	private void getDescendants(NodeRef n, List<Integer> descendants) {
		descendants.add(n.getNumber());
		if (tree.isExternal(n)) return;
		for (int i = 0; i < tree.getChildCount(n); ++i) getDescendants(tree.getChild(n, i), descendants);
	}
	
	@Override
	public int getHeight() {
		return getHeight(tree.getRoot().getNumber());
	}

	@Override
	public int getHeight(int i) {
		throw new UnsupportedOperationException();
	}

	@Override
	public int getLCA(int i, int j) {
		return Tree.Utils.getCommonAncestor(tree, tree.getNode(i), tree.getNode(j)).getNumber();
	}

	@Override
	public List<Integer> getLeaves() {
		int externalNodeCount = tree.getExternalNodeCount();
		List<Integer> leaves = new ArrayList<Integer>(externalNodeCount);
		for (int i = 0; i < externalNodeCount; ++i)
			leaves.add(tree.getExternalNode(i).getNumber());
		return leaves;
	}

	@Override
	public int getNoOfAncestors(int i, boolean b) {
		return getAncestors(i, b).size();
	}

	@Override
	public int getNoOfChildren(int i) {
		return tree.getChildCount(tree.getNode(i));
	}

	@Override
	public int getNoOfDescendantLeaves(int i, boolean b) {
		return getDescendantLeaves(i, b).size();
	}

	@Override
	public int getNoOfDescendants(int i, boolean b) {
		return getDescendants(i, b).size();
	}

	@Override
	public int getNoOfLeaves() {
		return tree.getExternalNodeCount();
	}

	@Override
	public int getParent(int i) {
		NodeRef n = tree.getParent(tree.getNode(i));
		return n != null ? n.getNumber() : RootedTree.NULL;
	}

	@Override
	public int getRoot() {
		return tree.getRoot().getNumber();
	}

	@Override
	public boolean isLeaf(int i) {
		return tree.isExternal(tree.getNode(i));
	}

	@Override
	public boolean isRoot(int i) {
		return tree.isRoot(tree.getNode(i));
	}

	@Override
	public List<Integer> getTopologicalOrdering(int i) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean hasArc(int i, int j) {
		return getParent(i) == j;
	}

	@Override
	public boolean hasPath(int i, int j) {
		do {
			j = getParent(j);
			if (i == j) return true;
		} while (j != NULL);
		return false;
	}

	@Override
	public List<Set<Integer>> getComponents() {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getName() {
		return tree.getId();
	}

	@Override
	public int getNoOfVertices() {
		return tree.getNodeCount();
	}

	@Override
	public void setName(String s) {
//		tree.setId(s);
		throw new UnsupportedOperationException(); // Avoid messing with BEAST internals
	}

	@Override
	public String toString() {
		return Tree.Utils.newick(tree);
	}

	@Override
	public Tree getWrappedBEASTObject() {
		return tree;
	}

}
