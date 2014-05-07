package org.ithinktree.becky.tools;

import java.util.Set;

import org.ithinktree.becky.CophylogenyLikelihood;
import org.ithinktree.becky.SimpleCophylogenyModel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.math.MathUtils;

public class ReverseCoevolutionSimulator {
	
	public ReverseCoevolutionSimulator() {
		
	}
	
	public Tree simulateCoevolution(final Tree hostTree, final Taxa symbiont, final String hostAttributeName) {
		
				
		Set<NodeRef>[] host2symbiont = new Set[hostTree.getNodeCount()];
		
		return null;//new SimpleTree();
	}

	
	public Tree simulateCoevolution(final Tree hostTree, final Taxa symbiont, final String hostAttributeName, final SimpleCophylogenyModel scm, final CophylogenyLikelihood cl) {
		
		// TODO
		
		if (scm.getLossRate() > 0.0) {
			throw new RuntimeException("No simulation support for models involving extinction");
		}
		
		final double duplicationRate = scm.getDuplicationRate();
		final double hostSwitchRate = scm.getHostSwitchRate();
		final double origin = hostTree.getNodeHeight(hostTree.getRoot()) + MathUtils.nextExponential(hostTree.getNodeCount() / hostTree.getNodeHeight(hostTree.getRoot()));
		
		
		
		NodeRef host2symbiont[][] = new NodeRef[hostTree.getNodeCount()][];
		
		return null;//new SimpleTree();
	}
	
	private class EventIndexAndTime {
		public final int index;
		public final double time;
		public EventIndexAndTime(int index, double time) {
			this.index = index;
			this.time = time;
		}
	}
	
	private final EventIndexAndTime nextEvent(final double origin, final double...lambdas) {
		final double lambda = MathUtils.getTotal(lambdas);
		final double[] p = new double[lambdas.length - 1];
		p[0] = lambdas[0] / lambda;
		for (int i = 1; i < p.length; ++i)
			p[i] = lambdas[i] / lambda + p[i - 1];
		final double time = nextYuleSpeciationTime(lambda, origin);
		final double U = 1 - MathUtils.nextDouble();
		int i;
		for (i = 0; i < p.length && p[i] < U; ++i);
		return new EventIndexAndTime(i, time);
	}
	
	private final double nextYuleSpeciationTime(final double lambda, final double origin) {
		
		return - Math.log(1 - (1-MathUtils.nextDouble()) * (1 - Math.exp(-lambda * origin)) / lambda);
		
	}
	
}
