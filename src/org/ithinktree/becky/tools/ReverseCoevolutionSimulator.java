package org.ithinktree.becky.tools;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.ithinktree.becky.CophylogenyModel;
import org.ithinktree.becky.SimpleCophylogenyModel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;

public class ReverseCoevolutionSimulator {
	
	public ReverseCoevolutionSimulator() {
		
	}
	
	private static final NodeRef[] EMPTY = new NodeRef[0];
	
	public Tree simulateCoevolution(final Tree hostTree, final Taxa symbiont, final String hostAttributeName, final SimpleCophylogenyModel scm) {
				
		if (scm.getLossRate() > 0.0) {
			throw new RuntimeException("No simulation support for models involving extinction");
		}
		
		final double duplicationRate = scm.getDuplicationRate();
		final double hostSwitchRate = scm.getHostSwitchRate();
		final double origin = hostTree.getNodeHeight(hostTree.getRoot()) + MathUtils.nextExponential(hostTree.getNodeCount() / hostTree.getNodeHeight(hostTree.getRoot()));
		
		@SuppressWarnings("unchecked")
		Set<NodeRef>[] host2symbiont = new Set[hostTree.getNodeCount()];
		
		for (int i = 0; i < host2symbiont.length; ++i) {
			host2symbiont[i] = new HashSet<NodeRef>();
		}
		
		for (int i = 0; i < symbiont.getTaxonCount(); ++i) {
			 NodeRef host = hostTree.getExternalNode(hostTree.getTaxonIndex(((Taxon) symbiont.getTaxonAttribute(i, hostAttributeName)).getId()));
			 SimpleNode node = new SimpleNode();
			 node.setHeight(hostTree.getNodeHeight(host));
			 node.setTaxon(symbiont.getTaxon(i));
			 host2symbiont[host.getNumber()].add(node);
		}
		
		double[] nodeHeights = new double[hostTree.getInternalNodeCount()];
		for (int i = 0; i < nodeHeights.length; ++i) nodeHeights[i] = hostTree.getNodeHeight(hostTree.getInternalNode(i));
		Arrays.sort(nodeHeights);
		
		int n = symbiont.getTaxonCount();
		int i = 0;
//		while (n > 1) {
//			double eventHeight;
//			do {
//				eventHeight = nextYuleSpeciationTime(origin, (n - 1) * (duplicationRate + hostSwitchRate));
//			} while (eventHeight < height);
//			boolean possible = true;
//			for (NodeRef node : CophylogenyModel.Utils.getContemporaneousLineages(hostTree, eventHeight)) {
//				if (host2symbiont[node.getNumber()].isEmpty()) {
//					NodeRef child1 = hostTree.getChild(node, 0);
//					NodeRef child2 = hostTree.getChild(node, 1);
//					if (host2symbiont[child1.getNumber()].size() != host2symbiont[child2.getNumber()].size()) {
//						possible = false;
//						break;
//					}
//				}
//			}
//			if (possible) {
//				for (NodeRef node : CophylogenyModel.Utils.getContemporaneousLineages(hostTree, eventHeight)) {
//					if (host2symbiont[node.getNumber()].isEmpty()) {
//						NodeRef child1 = hostTree.getChild(node, 0);
//						NodeRef child2 = hostTree.getChild(node, 1);
//						NodeRef[] symbionts1 = host2symbiont[child1.getNumber()].toArray(EMPTY);
//						NodeRef[] symbionts2 = host2symbiont[child2.getNumber()].toArray(EMPTY);
//						for (int i = 0; i < symbionts1.length; ++i) {
//							SimpleNode newNode = new SimpleNode();
//							newNode.addChild((SimpleNode) symbionts1[i]);
//							newNode.addChild((SimpleNode) symbionts2[i]);
//							newNode.setHeight(hostTree.getNodeHeight(node));
//							host2symbiont[node.getNumber()].add(newNode);
//						}
//					}
//				}
//				
//			}
//		}
		
		return null;//new SimpleTree();
	}
		
	private final int nextEvent(final double...lambdas) {
		final double lambda = MathUtils.getTotal(lambdas);
		final double[] p = new double[lambdas.length - 1];
		p[0] = lambdas[0] / lambda;
		for (int i = 1; i < p.length; ++i)
			p[i] = lambdas[i] / lambda + p[i - 1];
		final double U = 1 - MathUtils.nextDouble();
		int i;
		for (i = 0; i < p.length && p[i] < U; ++i);
		return i;
	}
	
	private final double nextYuleSpeciationTime(final double origin, final double lambda) {
		return - Math.log(1 - (1-MathUtils.nextDouble()) * (1 - Math.exp(-lambda * origin)) / lambda);
	}
	
}
