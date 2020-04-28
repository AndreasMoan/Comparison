package HGSADCwSO.protocols;

import HGSADCwSO.Genotype;
import HGSADCwSO.Individual;
import HGSADCwSO.Voyage;

import java.util.ArrayList;

public interface FitnessEvaluationProtocol {

    public void evaluate(Individual individual);


    public void updateBiasedFitness(ArrayList<Individual> population);

    public void addDiversityDistance(Individual individual);

    public void removeDiversityDistance(Individual individual);

    public double getHammingDistance(Individual individual1, Individual individual2);

    public void setPenalizedCostIndividual(Individual individual, double durationViolationPenalty, double capacityViolationPenalty);

    public void setPenalizedCostIndividual(Individual individual);

    public double getPenalizedCost(ArrayList<Integer> orderSequence); //har bytta Voyage voyage med ArrayList<Integer> orderSequence

    public double getPenalizedCost(ArrayList<Integer> orderSequence,  double durationViolationPenalty, double capacityViolationPenalty);//har bytta Voyage voyage med ArrayList<Integer> orderSequence

    public double getDurationViolationPenalty();

    public double getCapacityViolationPenalty();

    public double getNumberOfInstallationsPenalty();

    public void setDurationViolationPenalty(double durationViolationPenalty);

    public void setCapacityViolationPenalty(double capacityViolationPenalty);

    public void setNumberOfInstallationsViolationPenalty(double numberOfInstallationsViolationPenalty);

    void setPenalizedCostPopulation(ArrayList<Individual> population);

}
