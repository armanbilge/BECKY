/**
 * GeographicalHostBasedCophylogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky;

import dr.inference.model.Parameter;

/**
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public class GeographicalHostBasedCophylogenyModel extends
		HostBasedCophylogenyModel {

	/**
	 * @param duplicationRateParameter
	 * @param hostShiftRateParameter
	 * @param lossRateParameter
	 * @param units
	 */
	public GeographicalHostBasedCophylogenyModel(
			Parameter duplicationRateParameter,
			Parameter hostShiftRateParameter, Parameter lossRateParameter,
			Type units) {
		super(duplicationRateParameter, hostShiftRateParameter,
				lossRateParameter, units);
		// TODO Auto-generated constructor stub
	}

}
