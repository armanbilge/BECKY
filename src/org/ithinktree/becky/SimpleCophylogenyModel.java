/**
 * SimpleCophylogenyModel.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.becky;

import java.util.Arrays;

import org.ithinktree.becky.CophylogenyModel.Utils.NodalRelationship;
import org.ithinktree.becky.CophylogenyModel.Utils.Relationship;
import org.ithinktree.becky.xml.SimpleCophylogenyModelParser;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

/**
 * A simple model for cophylogenetic mappings.
 * 
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public class SimpleCophylogenyModel extends CophylogenyModel {

    final protected Parameter duplicationRateParameter;
    final protected Parameter hostShiftRateParameter;
    final protected Parameter lossRateParameter;
    private double duplicationRate;
    private double hostShiftRate;
    private double lossRate;
        
    /**
     * 
     */
    public SimpleCophylogenyModel(Parameter duplicationRateParameter, Parameter hostShiftRateParameter, Parameter lossRateParameter, Type units) {

        super(SimpleCophylogenyModelParser.SIMPLE_COPHYLOGENY_MODEL, units);
        
        this.duplicationRateParameter = duplicationRateParameter;
        addVariable(duplicationRateParameter);
        duplicationRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        
        this.hostShiftRateParameter = hostShiftRateParameter;
        addVariable(hostShiftRateParameter);
        hostShiftRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        
        this.lossRateParameter = lossRateParameter;
        addVariable(lossRateParameter);
        lossRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        
    }
    
    protected void updateVariables() {
        duplicationRate = getDuplicationRate();
        hostShiftRate = getHostShiftRate();
        lossRate = getLossRate();
        overallRate = duplicationRate + hostShiftRate + lossRate;
        dirty = false;
    }
    
    protected final double likelihoodEvent(final EventType e, final double t, final double rate) {
        switch (e) {
        case DUPLICATION: return likelihoodDuplicationAtTime(t, rate);
        case HOST_SHIFT: return likelihoodHostShiftAtTime(t, rate);
        case LOSS: return likelihoodLossInTime(t, rate);
        case NO_EVENT: return likelihoodNoEventsInTime(t, rate);
        }
        return 0.0;
    }
    
    protected double getEventRate(EventType e) {
        switch (e) {
        case DUPLICATION: return duplicationRate;
        case HOST_SHIFT: return hostShiftRate;
        case LOSS: return lossRate;
        case NO_EVENT: return 1.0;
        }
        return 0.0;
    }
    
    public double getDuplicationRate() {
        return duplicationRateParameter.getParameterValue(0);
    }
    
    protected final double likelihoodDuplicationAtTime(final double t, final double rate) {
        return likelihoodEventAtTime(t, duplicationRate, rate);
    }
    
    protected final double likelihoodDuplicationInTime(final double t, final double rate) {
        return likelihoodEventInTime(t, duplicationRate, rate);
    }
    
    public double getHostShiftRate() {
        return hostShiftRateParameter.getParameterValue(0);
    }
    
    protected final double likelihoodHostShiftAtTime(final double t, final double rate) {
        return likelihoodEventAtTime(t, hostShiftRate, rate);
    }
    
    protected final double likelihoodHostShiftInTime(final double t, final double rate) {
        return likelihoodEventInTime(t, hostShiftRate, rate);
    }

    public double getLossRate() {
        return lossRateParameter.getParameterValue(0);
    }
    
    protected final double likelihoodLossInTime(final double t, final double rate) {
        return likelihoodEventInTime(t, lossRate, rate);
    }
        
    protected final double likelihoodLineageLoss(final Tree tree, final NodeRef lineage, final double rate, boolean excludeRoot) {
        
        double sum = 0.0;
        final ExtinctionLikelihood[] els = permuteExtinctLineageLikelihoods(tree, lineage, rate);
        final int l = els.length - (excludeRoot ? 1 : 0);
        for (int i = 0; i < l; ++i)
            sum += els[i].getOverallLikelihood();
        return sum;
        
    }
    
    public ExtinctionLikelihood[] permuteExtinctLineageLikelihoods(final Tree tree, final NodeRef lineage, final double rate) {
        
        final double length = tree.getBranchLength(lineage);
        final ExtinctionLikelihood selfEL = new ExtinctionLikelihood(likelihoodLossInTime(length, rate));
        
        if (tree.isExternal(lineage))
            return new ExtinctionLikelihood[]{selfEL};
        
        final ExtinctionLikelihood[] child1ELs = permuteExtinctLineageLikelihoods(tree, tree.getChild(lineage, 0), rate);
        final ExtinctionLikelihood[] child2ELs = permuteExtinctLineageLikelihoods(tree, tree.getChild(lineage, 1), rate);
        
        final ExtinctionLikelihood[] els = new ExtinctionLikelihood[child1ELs.length * child2ELs.length + 1];

        int i = 0;
        final double selfNoEvent = likelihoodNoEventsInTime(length, rate);
        for (ExtinctionLikelihood el1 : child1ELs) {
            for (ExtinctionLikelihood el2 : child2ELs)
                els[i++] = el1.createCombination(el2, selfNoEvent);
        }
        
        els[i] = selfEL;
        
        return els;
    }
    
    protected final double[] emptyDoubleArray = new double[0];
    protected class ExtinctionLikelihood {
        
        private final double[] eventless, extinct;
        
        private ExtinctionLikelihood(double[] eventless, double[] extinct) {
            this.eventless = eventless;
            this.extinct = extinct;
        }
        
        public ExtinctionLikelihood(double d) {
            this(emptyDoubleArray, new double[]{d});
        }
        
        public double getOverallLikelihood() {
            double l = 1.0;
            for (double d : eventless) l *= d;
            for (double d : extinct) l *= d;
            return l;
        }
        
        public ExtinctionLikelihood createCombination(final ExtinctionLikelihood other, final double noEvent) {
            
            final double[] newEventless = new double[eventless.length + other.eventless.length + 1];
            final double[] newExtinct = new double[extinct.length + other.extinct.length];
            System.arraycopy(eventless, 0, newEventless, 0, eventless.length);
            System.arraycopy(other.eventless, 0, newEventless, eventless.length, other.eventless.length);
            newEventless[newEventless.length - 1] = noEvent;
            System.arraycopy(extinct, 0, newExtinct, 0, extinct.length);
            System.arraycopy(other.extinct, 0, newExtinct, extinct.length, other.extinct.length);
            
            return new ExtinctionLikelihood(newEventless, newExtinct);
        }
    }

    
    /**
     * Evaluates an integral
     * @param a
     * @param b
     * @param t
     * @param rate
     * @return
     */
    protected final double likelihoodHostShiftEventAndLossInTime(final double a, final double b, final double e, final double l, final double lambda_event, final double rate, final Tree tree, final NodeRef lostLineage) {
        final double hostShiftRate = this.hostShiftRate * rate;
        final double eventRate = lambda_event * rate;
        final double lossRate = this.lossRate * rate;
        final double overallRate = this.overallRate * rate;
        return hostShiftRate * eventRate * (lossRate / overallRate * Math.exp(-overallRate * e) * (b - a + (Math.exp(-overallRate * (l - a)) - Math.exp(-overallRate * (l - b))) / overallRate) -
                Math.exp(-overallRate * (e + l)) * (Math.exp(overallRate * a) - Math.exp(overallRate * b)) / overallRate * likelihoodLineageLoss(tree, lostLineage, rate, true));
    }

    protected double likelihoodLossesAlongLineages(final Tree tree, final NodeRef[] lineages, double rate) {
        double likelihood = 1.0;
        for (NodeRef n : lineages)
            likelihood *= likelihoodLineageLoss(tree, n, rate, false);
        return likelihood;
    }
    
    /**
     * Discretizes a branch along which a host-shift and loss happened to properly determine the likelihood of these events.
     * 
     * @param start the earliest the host-shift could have happened
     * @param hostShiftStop the latest the host-shift could have happened
     * @param lossStop the latest the loss could have happened
     * @param rate the rate on the relevant branch
     * @param tree the host tree
     * @param originalLineages lineages lost on original host lineage
     * @param newHostLineages lineages lost on new host lineage
     * @return
     */
    protected final double likelihoodHostShiftEventAndLossInTime(final double start, final double hostShiftStop, final double eventStop, final double lossStop, final double eventRate, final double rate, final Tree tree, final NodeRef lostLineage, final NodeRef[] originalLineages, final NodeRef[] newHostLineages) {
        
        final double e = start - eventStop;
        final double l = start - lossStop;
        final double stop = start - Math.max(hostShiftStop, Math.max(eventStop, lossStop));
        double likelihood = 0;
        double height;
        double subHeight;
        double nextHeight = 0.0;
        double nextSubHeight;
        for (int i = originalLineages.length - 1; nextHeight < stop; --i) {
            height = nextHeight;
            if (i >= 0) nextHeight = Math.min(start - tree.getNodeHeight(tree.getParent(originalLineages[i])), stop);
            else nextHeight = stop;
            nextSubHeight = height;
            for (int j = newHostLineages.length - 1; nextSubHeight < nextHeight; --j) {
                subHeight = nextSubHeight;
                if (j >= 0) nextSubHeight = Math.min(start - tree.getNodeHeight(tree.getParent(newHostLineages[j])), nextHeight);
                else nextSubHeight = nextHeight;
                likelihood += likelihoodHostShiftEventAndLossInTime(subHeight, nextSubHeight, e, l, eventRate, rate, tree, lostLineage) *
                        likelihoodLossesAlongLineages(tree, Arrays.copyOfRange(originalLineages, i+1, originalLineages.length), rate) *
                        likelihoodLossesAlongLineages(tree, Arrays.copyOfRange(newHostLineages, 0, j+1), rate) /
                        (Utils.getContemporaneousLineageCount(tree, nextSubHeight) - 1);
            }
        }
                
        return likelihood;
    }
    
    /**
     * Calculates the probability of a particular cophylogenetic mapping at a
     * node and its children against this model's current state.
     * 
     * @param hostTree host tree
     * @param selfHost host of the given node
     * @param child1Host host of its 1st child
     * @param child2Host host of its 2nd child
     * @param selfBranchTime branch length/time of the given node
     * @param child1BranchTime branch length/time of its 1st child
     * @param child2BranchTime branch length/time of its 2nd child
     * @return log likelihood of the cophylogenetic mapping at node
     */
    @Override
    public double calculateNodeLogLikelihood(final MutableTree symbiontTree, final NodeRef self,
            final NodeRef child1, final NodeRef child2, final Tree hostTree, final NodeRef selfHost,
            final NodeRef child1Host, final NodeRef child2Host, final BranchRates branchRates) {
                
                    if (dirty) updateVariables();
        
                    double likelihood = 1.0;
                    double sum;
                    
                    final NodalRelationship child1Relationship = Utils.determineRelationship(hostTree, selfHost, child1Host);
                    final NodalRelationship child2Relationship = Utils.determineRelationship(hostTree, selfHost, child2Host);
                    double selfHeight = symbiontTree.getNodeHeight(self);
                    final double selfHostHeight = hostTree.getNodeHeight(selfHost);
                    final double child1BranchRate = branchRates.getBranchRate(symbiontTree, child1);
                    final double child2BranchRate = branchRates.getBranchRate(symbiontTree, child2);

                    boolean calculatedChild1 = false;
                    boolean calculatedChild2 = false;
                    
                    // Check if symbiont coexisted temporally with its host
                    if ((!hostTree.isRoot(selfHost) && selfHeight >= hostTree.getNodeHeight(hostTree.getParent(selfHost))) ||
                            selfHeight < selfHostHeight)
                        return Double.NEGATIVE_INFINITY;
                                        
                    if (child1Relationship.relationship == Relationship.DESCENDANT
                            && child2Relationship.relationship == Relationship.DESCENDANT) {
                                                
                        NodeRef hostChild = hostTree.getChild(selfHost, 0);
                        final Utils.NodalRelationship nr1 = Utils.determineRelationship(hostTree, hostChild, child1Host);
                        final Utils.NodalRelationship nr2 = Utils.determineRelationship(hostTree, hostChild, child2Host);
                        if (nr1.relationship == nr2.relationship || (nr1.relationship == Relationship.SISTER && nr2.relationship == Relationship.COUSIN) || (nr2.relationship == Relationship.SISTER && nr1.relationship == Relationship.COUSIN)) {
                            
                            // Determine along which child lineage the loss(es) happened
                            if (nr1.relationship != Relationship.DESCENDANT)
                                hostChild = hostTree.getChild(selfHost, 1);
                            final double hostChildBranchLength = hostTree.getBranchLength(hostChild);
                            final double hostChildHeight = hostTree.getNodeHeight(hostChild);
                                
                            
                            final Event[] events = new Event[2];
                            // Case 1: duplication, cospeciation, then losses
                            double case1 = 1.0;

                            // The duplication occurred at the time of their last common ancestor
                            
                            final double potentialLossLength = (selfHeight - selfHostHeight) + hostChildBranchLength;
                            case1 *= likelihoodLossInTime(potentialLossLength, child1BranchRate) * likelihoodLossInTime(potentialLossLength, child2BranchRate)
                                    * likelihoodLineageLoss(hostTree, hostChild, child1BranchRate, true) * likelihoodLineageLoss(hostTree, hostChild, child2BranchRate, true);
                            sum = 0.0;
                            for (Event e : reconstructedEvents[child1.getNumber()])
                                sum += likelihoodEvent(e.event, symbiontTree.getBranchLength(child1), child1BranchRate) * e.partialLikelihood;
                            case1 *= sum;
                            sum = 0.0;
                            for (Event e : reconstructedEvents[child2.getNumber()])
                                sum += likelihoodEvent(e.event, symbiontTree.getBranchLength(child2), child2BranchRate) * e.partialLikelihood;
                            case1 *= sum;
                            
                            events[0] = new Event(EventType.DUPLICATION, case1);
                            
                            // Case 2: cospeciation, then host-shift and loss
                            
                            final double child1Height = symbiontTree.getNodeHeight(child1);
                            final double child2Height = symbiontTree.getNodeHeight(child2);
                            
                            // All possible lineages along which losses may have occurred
                            final NodeRef[] child1OriginalHostLineages = Utils.lostLineagesToTime(hostTree, hostChild, selfHeight);
                            final NodeRef[] child2OriginalHostLineages = Utils.lostLineagesToTime(hostTree, hostChild, selfHeight);
                            final NodeRef[] child1NewHostLineages = child1Relationship.lostLineages; // Utils.lostLineagesToTime(hostTree, child1Host, selfHeight);
                            final NodeRef[] child2NewHostLineages = child2Relationship.lostLineages; // Utils.lostLineagesToTime(hostTree, child2Host, selfHeight);
                            
                            // Sum over two subcases: child1 lineage made host-shift/loss or child2 made host-shift/loss
                            double case2 = 0.0;
                            double sum2;
                            
                            sum = 0.0;
                            for (Event e : reconstructedEvents[child1.getNumber()])
                                sum += likelihoodHostShiftEventAndLossInTime(selfHeight, child1Height, child1Height, hostChildHeight, getEventRate(e.event), child1BranchRate, hostTree, hostChild, child1OriginalHostLineages, child1NewHostLineages) * e.partialLikelihood;
        
                            sum2 = 0.0;
                            for (Event e : reconstructedEvents[child2.getNumber()])
                                sum2 += likelihoodEvent(e.event, symbiontTree.getBranchLength(child2), child2BranchRate) * e.partialLikelihood;

                            case2 += sum * sum2 * likelihoodLossesAlongLineages(hostTree, child1NewHostLineages, child1BranchRate);
                            
                            sum = 0.0;
                            for (Event e : reconstructedEvents[child2.getNumber()])
                                sum += likelihoodHostShiftEventAndLossInTime(selfHeight, child2Height, child2Height, hostChildHeight, getEventRate(e.event), child2BranchRate, hostTree, hostChild, child2OriginalHostLineages, child2NewHostLineages) * e.partialLikelihood;
 
                            sum2 = 0.0;
                            for (Event e : reconstructedEvents[child1.getNumber()])
                                sum2 += likelihoodEvent(e.event, symbiontTree.getBranchLength(child1), child1BranchRate) * e.partialLikelihood;

                            case2 += sum * sum2 * likelihoodLossesAlongLineages(hostTree, child1NewHostLineages, child2BranchRate);
                                                        
                            events[1] = new Event(EventType.NO_EVENT, case2);
                            
                            setReconstructedEvents(self, events);
                            calculatedChild1 = true;
                            calculatedChild2 = true;
                            
                        } else { // Plain old cospeciation

                            // Check if violates tree validity (parent younger than children)
//                          if (selfHostHeight < symbiontTree.getNodeHeight(child1) || selfHostHeight < symbiontTree.getNodeHeight(child2))
//                              return Double.NEGATIVE_INFINITY;
//                          symbiontTree.setNodeHeight(self, selfHostHeight);
                            
                            setReconstructedEvents(self, new Event[]{NO_EVENT});
                            
                            // Potential losses along both child lineages
                            likelihood *= likelihoodLossesAlongLineages(hostTree, child1Relationship.lostLineages, child1BranchRate);
                            likelihood *= likelihoodLossesAlongLineages(hostTree, child2Relationship.lostLineages, child2BranchRate);                           
                        }
                        
                    } else if (child1Relationship.relationship == Relationship.SELF
                            && child2Relationship.relationship == Relationship.SELF) {
                        
                        // Duplication event
                        setReconstructedEvents(self, new Event[]{DUPLICATION});
                        
                    } else if (child1Relationship.relationship == Relationship.SELF && (child2Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {

                        // Child2 host-shift event
                        setReconstructedEvents(self, new Event[]{HOST_SHIFT});
                        likelihood /= (Utils.getContemporaneousLineageCount(hostTree, selfHeight) - 1);
                        likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child2Host, selfHeight), branchRates.getBranchRate(symbiontTree, child2));
                    
                    } else if (child2Relationship.relationship == Relationship.SELF && (child1Relationship.relationship == Relationship.COUSIN || child1Relationship.relationship == Relationship.SISTER)) {
                        
                        // Child1 host-shift event
                        setReconstructedEvents(self, new Event[]{HOST_SHIFT});
                        likelihood /= (Utils.getContemporaneousLineageCount(hostTree, selfHeight) - 1);
                        likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child1Host, selfHeight), branchRates.getBranchRate(symbiontTree, child1));
                        
                    } else if ((child1Relationship.relationship == Relationship.COUSIN || child1Relationship.relationship == Relationship.SISTER)
                            && (child2Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {

                        // Double host-shift event with a loss: no symbionts left on this host lineage
                                                
                        final double child1Height = symbiontTree.getNodeHeight(child1);
                        final double child2Height = symbiontTree.getNodeHeight(child2);
                        
                        final NodeRef[] noLineages = new NodeRef[0];
                        final NodeRef[] child1NewHostLineages = Utils.lostLineagesToTime(hostTree, child1Host, selfHeight);
                        final NodeRef[] child2NewHostLineages = Utils.lostLineagesToTime(hostTree, child2Host, selfHeight);
                        
                        // We definitely know the time of the first host-shift
                        setReconstructedEvents(self, new Event[] {HOST_SHIFT});
                        likelihood /= (Utils.getContemporaneousLineageCount(hostTree, selfHeight) - 1);
                        
                        // Case 1: Child1 lineage host-shifted first
                        double case1 = likelihoodLossesAlongLineages(hostTree, child1NewHostLineages, child1BranchRate);
                        // Case 2: Child2 lineage host-shifted first
                        double case2 = likelihoodLossesAlongLineages(hostTree, child2NewHostLineages, child2BranchRate);

                        sum = 0.0;
                        for (Event e : getReconstructedEvents(child2))
                            sum += likelihoodHostShiftEventAndLossInTime(selfHeight, child2Height, child2Height, selfHostHeight, getEventRate(e.event), child2BranchRate, hostTree, selfHost, noLineages, child2NewHostLineages) *
                                e.partialLikelihood;
                        case1 *= sum;

                        sum = 0.0;
                        for (Event e : getReconstructedEvents(child1))                      
                            sum += likelihoodHostShiftEventAndLossInTime(selfHeight, child1Height, child1Height, selfHostHeight, getEventRate(e.event), child1BranchRate, hostTree, selfHost, noLineages, child1NewHostLineages) *
                                e.partialLikelihood;
                        case2 *= sum;
                        
                        likelihood *= case1 + case2;
                        
                        calculatedChild1 = true;
                        calculatedChild2 = true;
                        
                    } else if (child1Relationship.relationship == Relationship.DESCENDANT
                            && (child2Relationship.relationship == Relationship.COUSIN || child2Relationship.relationship == Relationship.SISTER)) {

                        // Child2 host-shift and child1 losses
                        setReconstructedEvents(self, new Event[]{HOST_SHIFT});
                        likelihood /= (Utils.getContemporaneousLineageCount(hostTree, selfHeight) - 1);
                        likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child2Host, selfHeight), branchRates.getBranchRate(symbiontTree, child2));
                        likelihood *= likelihoodLossesAlongLineages(hostTree, child1Relationship.lostLineages, branchRates.getBranchRate(symbiontTree, child1));

                    } else if ((child1Relationship.relationship == Relationship.COUSIN || child1Relationship.relationship == Relationship.SISTER)
                            && child2Relationship.relationship == Relationship.DESCENDANT) {

                        // Child1 host-shift and child2 losses 
                        setReconstructedEvents(self, new Event[]{HOST_SHIFT});
                        likelihood /= (Utils.getContemporaneousLineageCount(hostTree, selfHeight) - 1);
                        likelihood *= likelihoodLossesAlongLineages(hostTree, Utils.lostLineagesToTime(hostTree, child1Host, selfHeight), branchRates.getBranchRate(symbiontTree, child1));
                        likelihood *= likelihoodLossesAlongLineages(hostTree, child2Relationship.lostLineages, branchRates.getBranchRate(symbiontTree, child2));

                    } else { // Everything else is impossible
                        return Double.NEGATIVE_INFINITY;
                    }
                                                            
                    if (!calculatedChild1) {
                        sum = 0.0;
                        for (Event e : reconstructedEvents[child1.getNumber()])
                            sum += likelihoodEvent(e.event, symbiontTree.getBranchLength(child1), child1BranchRate) * e.partialLikelihood;
                        likelihood *= sum;
                    }
                    
                    if (!calculatedChild2) {
                        sum = 0.0;
                        for (Event e : reconstructedEvents[child2.getNumber()])
                            sum += likelihoodEvent(e.event, symbiontTree.getBranchLength(child2), child2BranchRate) * e.partialLikelihood;
                        likelihood *= sum;
                    }

                    return Math.log(likelihood);
    }

    protected class Event {
        public final EventType event;
        public final double partialLikelihood;
        public Event(final EventType et, final double pl) {
            event = et;
            partialLikelihood = pl;
        }
        public Event(final EventType et) {
            this(et, 1.0);
        }
    }
    
    protected final Event DUPLICATION = new Event(EventType.DUPLICATION);
    protected final Event HOST_SHIFT = new Event(EventType.HOST_SHIFT);
    protected final Event NO_EVENT = new Event(EventType.NO_EVENT);
    
    protected enum EventType {
        DUPLICATION, HOST_SHIFT, LOSS, NO_EVENT;
    }
    
    protected void setReconstructedEvents(NodeRef n, Event[] e) {
        reconstructedEvents[n.getNumber()] = e;
    }
    
    protected Event[] getReconstructedEvents(NodeRef n) {
        return reconstructedEvents[n.getNumber()];
    }
    
    private Event[][] reconstructedEvents;
    public void initialize(final Tree tree) {
        reconstructedEvents = new Event[tree.getNodeCount()][];
        for (int i = 0; i < tree.getExternalNodeCount(); ++i)
            reconstructedEvents[tree.getExternalNode(i).getNumber()] = new Event[]{NO_EVENT};
    }
    
}
