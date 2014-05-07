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
import org.ithinktree.becky.xml.CophylogenySetupParser;
import org.ithinktree.becky.xml.CospeciationOperatorParser;
import org.ithinktree.becky.xml.CospeciationSimulatorParser;
import org.ithinktree.becky.xml.HostSwitchOperatorParser;
import org.ithinktree.becky.xml.HostSwitchingWilsonBaldingParser;
import org.ithinktree.becky.xml.NodeRefProviderParser;
import org.ithinktree.becky.xml.PreAnnotatorParser;
import org.ithinktree.becky.xml.SimpleCophylogenyModelParser;
import org.ithinktree.becky.xml.SimpleStochasticCophylogenyModelParser;
import org.ithinktree.becky.xml.TipHostSwitchOperatorParser;
import org.ithinktree.becky.xml.TugOperatorParser;

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
		parsers.add(new SimpleStochasticCophylogenyModelParser());
		parsers.add(new CophylogenyLikelihoodParser());
		parsers.add(new HostSwitchOperatorParser());
		parsers.add(new TipHostSwitchOperatorParser());
		parsers.add(new CospeciationOperatorParser());
		parsers.add(new CoevolutionSimulatorParser());
		parsers.add(new CospeciationSimulatorParser());
		parsers.add(new CophylogenySetupParser());
		parsers.add(new NodeRefProviderParser());
		parsers.add(new CophylogenyLikelihoodWrapperForJPrIMEDLTRSModelParser());
		parsers.add(new PreAnnotatorParser());
		parsers.add(new TugOperatorParser());
		parsers.add(new HostSwitchingWilsonBaldingParser());
		return parsers;
	}

}
