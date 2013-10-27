/**
 * 
 */
package org.ithinktree.becky;

import dr.inference.model.Parameter;

/**
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public class HostBasedCophylogenyModel extends SimpleCophylogenyModel {

	/**
	 * @param duplicationRateParameter
	 * @param hostSwitchRateParameter
	 * @param lossRateParameter
	 * @param units
	 */
	public HostBasedCophylogenyModel(Parameter duplicationRateParameter,
			Parameter hostSwitchRateParameter, Parameter lossRateParameter,
			Type units) {
		super(duplicationRateParameter, hostSwitchRateParameter,
				lossRateParameter, units);
	}

}
