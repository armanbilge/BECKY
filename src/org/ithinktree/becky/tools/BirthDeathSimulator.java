package org.ithinktree.becky.tools;

import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;

public class BirthDeathSimulator {

	private int taxonCount;
	
	public final Tree simulateBirthDeathTree(final double origin, final double birthRate, final double deathRate) {
		taxonCount = 0;
		return new SimpleTree(simulateBirthDeathNode(origin, birthRate, deathRate));
	}
	
	private final SimpleNode simulateBirthDeathNode(final double height, final double birthRate, final double deathRate) {
		
		final SimpleNode node = new SimpleNode();
		
		final SimpleNode child1;
		final SimpleNode child2;
		
		final int nextEvent = nextPoissonEvent(birthRate, deathRate);
		final double eventHeight = height - nextPoissonEventTime(birthRate, deathRate);
		
		if (eventHeight < 0) {
			node.setHeight(0.0);
			node.setTaxon(new Taxon("host" + ++taxonCount));
			return node;
		}
		
		node.setHeight(eventHeight);
		switch (nextEvent) {
		case 0:
			child1 = simulateBirthDeathNode(eventHeight, birthRate, deathRate);
			child2 = simulateBirthDeathNode(eventHeight, birthRate, deathRate);
			break;
		case 1: return null;
		default: throw new RuntimeException("Unknown event: " + nextEvent);
		}
		
		if (child1 != null && child2 != null) {
			// Both lineages survived
			node.addChild(child1);
			node.addChild(child2);
		} else if (child1 == child2) {
			// Both are null, hence this entire lineage was lost
			return null;
		} else {
			// Just one child lineage is null/lost
			return child2 == null ? child1 : child2; // Set this node to the continuing child lineage
		}
		return node;
		
	}
	
	private final double nextPoissonEventTime(final double...lambdas) {
		final double lambda = MathUtils.getTotal(lambdas);
		return MathUtils.nextExponential(lambda);
	}
	
	private final int nextPoissonEvent(final double...lambdas) {
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

	
}
