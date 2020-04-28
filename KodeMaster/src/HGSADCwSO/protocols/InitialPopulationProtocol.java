package HGSADCwSO.protocols;

import HGSADCwSO.Individual;

public interface InitialPopulationProtocol {

    public Individual createIndividual();

    public int getNumberOfConstructionHeuristicRestarts();

}
