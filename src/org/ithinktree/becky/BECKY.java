/**
 * BECKY.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 */
package org.ithinktree.becky;

import java.util.HashSet;
import java.util.Set;

import dr.app.plugin.Plugin;
import dr.xml.XMLObjectParser;

/**
 * @author Arman D. Bilge
 *
 */
public class BECKY implements Plugin {

	public Set<XMLObjectParser> getParsers() {
		
		Set<XMLObjectParser> parsers = new HashSet<XMLObjectParser>();
		parsers.add(new SimpleCophylogenyModelParser());
		parsers.add(new CophylogenyLikelihoodParser());
		parsers.add(new HostShiftOperatorParser());
		parsers.add(new CoevolutionSimulatorParser());
		parsers.add(new NodeRefProviderParser());
		return parsers;
	}

}
