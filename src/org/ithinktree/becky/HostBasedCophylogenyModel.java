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
	 * @param hostShiftRateParameter
	 * @param lossRateParameter
	 * @param units
	 */
	public HostBasedCophylogenyModel(Parameter duplicationRateParameter,
			Parameter hostShiftRateParameter, Parameter lossRateParameter,
			Type units) {
		super(duplicationRateParameter, hostShiftRateParameter,
				lossRateParameter, units);
		// TODO Auto-generated constructor stub
	}

}
