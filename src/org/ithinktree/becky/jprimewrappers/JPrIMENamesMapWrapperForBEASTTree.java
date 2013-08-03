package org.ithinktree.becky.jprimewrappers;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import se.cbb.jprime.topology.NamesMap;
import se.cbb.jprime.topology.RootedTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

public class JPrIMENamesMapWrapperForBEASTTree extends NamesMap implements
		WrappedBEASTObject<Tree> {

	private final Tree tree;
	
	public JPrIMENamesMapWrapperForBEASTTree(Tree tree) {
		super(tree.getId(), null);
		this.tree = tree;
	}

	@Override
	public Tree getWrappedBEASTObject() {
		return tree;
	}

	/**
	 * Returns the name of this map. Don't confuse this with names
	 * of individual elements (returned by get(x)).
	 * @return the name of the map.
	 */
	@Override
	public String getName() {
		return tree.getId();
	}

	/**
	 * Sets the name of this map. Don't confuse this with names
	 * of individual elements.
	 * @param name the name of the map.
	 */
	@Override
	public void setName(String name) {
		// Avoid messing with BEAST internals
	}
	
	@Override
	public Object getAsObject(int x) {
		return get(x);
	}

	@Override
	public void setAsObject(int x, Object value) {
		set(x, value.toString());
	}

	/**
	 * Returns the name of a vertex/arc.
	 * @param x the vertex/head of arc.
	 * @return the name.
	 */
	public String get(int x) {
		Taxon t = tree.getNodeTaxon(tree.getNode(x));
		return t != null ? t.getId() : "";
	}
	
	/**
	 * Sets the name of a vertex/arc. No check is made for uniqueness.
	 * @param x the vertex/head of arc.
	 * @param val the name.
	 */
	public void set(int x, String val) {
		// Avoid messing with BEAST internals
	}
	
	/**
	 * Swap two vertices numbers given two vertices names.
	 * @param u vertex name.
	 * @param v vertex name.
	 */
	public void swapVertices(String u, String v) {
		// Avoid messing with BEAST internals
	}
	
	/**
	 * Changes the vertex number associated with a given name.
	 * @param val name of the vertex.
	 * @param x vertex number.
	 */
	public void changeVertex(String val, int x) {
		// Avoid messing with BEAST internals
	}
	
	/**
	 * Returns the vertex/arc of a certain name.
	 * @param name the name.
	 * @return the vertex containing the name.
	 */
	public int getVertex(String val) {
		for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
			NodeRef n = tree.getNode(i);
			if (tree.getNodeTaxon(n).getId() == val) return n.getNumber();
		}
		return RootedTree.NULL;
	}
	
	/**
	 * Returns all non-null names of vertices/arcs.
	 * @param excludeBootstrapNames excludes all integer names.
	 * @return the names.
	 */
	public Set<String> getNames(boolean excludeBootstrapNames) {
		Set<String> names = new HashSet<String>(tree.getExternalNodeCount());
		for (Iterator<Taxon> taxa = tree.iterator(); taxa.hasNext(); names.add(taxa.next().getId()));
		return names;
	}
	
	/**
	 * Returns all non-null names of vertices/arcs in a sorted representation.
	 * @param excludeBootstrapNames excludes all integer names.
	 * @return the names.
	 */
	public TreeSet<String> getNamesSorted(boolean excludeBootstrapNames) {
		return new TreeSet<String>(getNames(excludeBootstrapNames));
	}

	
}
