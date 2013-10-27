/**
 * BECKY.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 */
package org.ithinktree.becky;

import java.util.HashSet;
import java.util.Set;

import org.ithinktree.becky.jprimewrappers.xml.CophylogenyLikelihoodWrapperForJPrIMEDLTRSModelParser;
import org.ithinktree.becky.xml.CoevolutionSimulatorParser;
import org.ithinktree.becky.xml.CophylogenyLikelihoodParser;
import org.ithinktree.becky.xml.CospeciationOperatorParser;
import org.ithinktree.becky.xml.HostShiftOperatorParser;
import org.ithinktree.becky.xml.NodeRefProviderParser;
import org.ithinktree.becky.xml.PreAnnotatorParser;
import org.ithinktree.becky.xml.SimpleCophylogenyModelParser;

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
		parsers.add(new CospeciationOperatorParser());
		parsers.add(new CoevolutionSimulatorParser());
		parsers.add(new NodeRefProviderParser());
		parsers.add(new CophylogenyLikelihoodWrapperForJPrIMEDLTRSModelParser());
		parsers.add(new PreAnnotatorParser());
		return parsers;
	}

}
