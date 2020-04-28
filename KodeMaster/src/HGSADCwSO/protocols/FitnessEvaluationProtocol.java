package HGSADCwSO.protocols;

import HGSADCwSO.Individual;

import java.util.ArrayList;

public interface FitnessEvaluationProtocol {

    public void evaluate(Individual individual);

    public void updateBiasedFitness(ArrayList<Individual> population);

    public void addDiversityDistance(Individual individual);

    public void removeDiversityDistance(Individual individual);

    public double getHammingDistance(Individual individual1, Individual individual2);

    public void setPenalizedCostIndividual(Individual individual, double durationViolationPenalty, double capacityViolationPenalty, double deadlineViolationPenalty);

    public void setPenalizedCostIndividual(Individual individual);

    public double getPenalizedCostOfVoyage(ArrayList<Integer> orderSequence);

    public double getPenalizedCostOfVoyage(ArrayList<Integer> orderSequence,  double durationViolationPenalty, double capacityViolationPenalty, double deadlineViolationPenalty);//har bytta Voyage voyage med ArrayList<Integer> orderSequence

    public double getDurationViolationPenalty();

    public double getCapacityViolationPenalty();

    public double getDeadlineViolationPenalty();

    public void setDurationViolationPenalty(double durationViolationPenalty);

    public void setCapacityViolationPenalty(double capacityViolationPenalty);

    public void setDeadlineViolationPenalty(double deadlineViolationPenalty);

    void setPenalizedCostPopulation(ArrayList<Individual> population);




}
