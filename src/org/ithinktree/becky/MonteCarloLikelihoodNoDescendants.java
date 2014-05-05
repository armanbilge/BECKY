package org.ithinktree.becky;

import org.ithinktree.becky.tools.CoevolutionSimulator;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

public class MonteCarloLikelihoodNoDescendants {

	final CoevolutionSimulator simulator = new CoevolutionSimulator();
	final SimpleCophylogenyModel cophylogenyModel; // Eventually should be made generic
	final int iterations;
	
	public MonteCarloLikelihoodNoDescendants(final SimpleCophylogenyModel cophylogenyModel, final int iterations) {
		this.cophylogenyModel = cophylogenyModel;
		this.iterations = iterations;
	}
	
	public double likelihoodNoDescendants(final Tree hostTree, final NodeRef originHost, final double originHeight, final double rate) {
		int noDescendantCount = 0;
		for (int i = 0; i < iterations; ++i) {
			if (simulator.simulateCoevolution(hostTree, originHost, originHeight, rate, this.cophylogenyModel, false, 0.0) == null) ++noDescendantCount;
		}
		return noDescendantCount / iterations;
	}
	
}
