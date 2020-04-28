package HGSADCwSO.protocols;

import HGSADCwSO.Individual;

import java.util.ArrayList;

public interface ReproductionProtocol {

    public Individual crossover(ArrayList<Individual> parents);

    public int getNumberOfCrossoverRestarts();

}
