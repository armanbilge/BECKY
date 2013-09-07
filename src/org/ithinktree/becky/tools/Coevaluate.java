/**
 * Coevaluate.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky.tools;

import java.io.FileReader;
import java.io.IOException;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;

/**
 * @author Arman D. Bilge
 *
 */
public class Coevaluate {

	private TaxonList taxa;
	
	/**
	 * @throws ImportException 
	 * @throws IOException 
	 * 
	 */
	public Coevaluate(final String host, final String hostSimulated, final String symbiont, final String symbiontSimulated) throws IOException, ImportException {
		
		FileReader fr = new FileReader(host);
		TreeImporter ti = new NexusImporter(fr);
		Tree tr = ti.importNextTree();
		fr.close();
		taxa = tr;

		final Map<BitSet,Integer> cladeToNodeRef = new HashMap<BitSet,Integer>(tr.getNodeCount());
		treeToBitSet(tr, null, cladeToNodeRef, "nodeRef");
		
		fr = new FileReader(hostSimulated);
		ti = new NexusImporter(fr);
		tr = ti.importNextTree();
		fr.close();
		Map<Integer,BitSet> temp = new HashMap<Integer,BitSet>(tr.getNodeCount());
		treeToBitSet(tr, temp, null, "nodeRef");
		
		// maps Simulated host nodeRefs to host nodeRefs
		final int[] nodeRefToNodeRef = new int[tr.getNodeCount()];
		for (int i = 0; i < nodeRefToNodeRef.length; ++i)
			nodeRefToNodeRef[i] = cladeToNodeRef.get(temp.get(i));
		
		fr = new FileReader(symbiont);
		ti = new NexusImporter(fr);
		final Tree symbiontTree = ti.importNextTree();
		fr.close();
		taxa = symbiontTree;
		
		Map<BitSet,Integer> symbiontMap = new HashMap<BitSet,Integer>(symbiontTree.getNodeCount());
		treeToBitSet(symbiontTree, null, symbiontMap, "host.nodeRef");
		
		fr = new FileReader(symbiontSimulated);
		ti = new NexusImporter(fr);
		final Tree symbiontSimulatedTree = ti.importNextTree();
		fr.close();
		
		Map<BitSet,Integer> symbiontSimulatedMap = new HashMap<BitSet,Integer>(symbiontSimulatedTree.getNodeCount());
		treeToBitSet(symbiontSimulatedTree, null, symbiontSimulatedMap, "host.nodeRef");
		
		int correct = 0;
		for (BitSet bs : symbiontMap.keySet()) {
			if (symbiontMap.get(bs) == nodeRefToNodeRef[symbiontSimulatedMap.get(bs)])
				correct++;
		}
		
		System.out.println("Percent correct: " + correct / (double) tr.getNodeCount() * 100 + "%");
	}

	private void treeToBitSet(final Tree t, final Map<Integer,BitSet> m1, final Map<BitSet,Integer> m2, final String s) {
		treeToBitSet(t, t.getRoot(), m1, m2, s);
	}
	
	private BitSet treeToBitSet(final Tree t, final NodeRef n, final Map<Integer,BitSet> m1, final Map<BitSet,Integer> m2, final String s) {
		final BitSet b = new BitSet(t.getTaxonCount());
		if (t.isExternal(n)) {
			b.set(taxa.getTaxonIndex(t.getNodeTaxon(n).getId()));
		} else {
			for (int i = 0; i < t.getChildCount(n); ++i)
				b.or(treeToBitSet(t, t.getChild(n, i), m1, m2, s));
		}
		Integer i = (Integer) t.getNodeAttribute(n, s);
		if (m1 != null) m1.put(i, b);
		if (m2 != null) m2.put(b, i);
		return b;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			new Coevaluate(args[0], args[1], args[2], args[3]);
		} catch (Exception e) {
			System.out.println("Usage: hostTree hostSimulatedTree symbiontTree symbiontSimulatedTree");
			e.printStackTrace();
			System.exit(1);
		}

	}

}
