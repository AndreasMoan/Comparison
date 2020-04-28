package HGSADCwSO.implementations.DAG;

import HGSADCwSO.*;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;

import java.util.ArrayList;


public class FitnessEvaluationDAG implements FitnessEvaluationProtocol {

    private ProblemData problemData;

    public FitnessEvaluationDAG(ProblemData problemData){
        this.problemData = problemData;
    }


    public void evaluate(Individual individual) {

        Genotype genotype = individual.getGenotype();

        boolean feasibility = true;
        double fitness = 0;

        int multiplier = (int) problemData.getHeuristicParameterDouble("Number of time periods per hour");

        for (Vessel vessel : problemData.getVessels()){
            DAG graph = new DAG(genotype.getVesselTourChromosome().get(vessel.getNumber()), vessel.getReturnDay()*24*multiplier, this.problemData);
            fitness += graph.getShortestFeasiblePathCost();
            feasibility = graph.getFeasibility() && feasibility;

        }

        System.out.println("Fitness this iteration is equal to:          " + fitness);

        individual.setFitness(fitness);
        individual.setFeasibility(feasibility);
    }


    @Override
    public void updateBiasedFitness(ArrayList<Individual> population) {

    }

    @Override
    public void addDiversityDistance(Individual individual) {

    }

    @Override
    public void removeDiversityDistance(Individual individual) {

    }

    @Override
    public double getHammingDistance(Individual individual1, Individual individual2) {
        return 0;
    }

    @Override
    public void setPenalizedCostIndividual(Individual individual, double durationViolationPenalty, double capacityViolationPenalty) {

    }

    @Override
    public void setPenalizedCostIndividual(Individual individual) {

    }

    @Override
    public double getPenalizedCost(ArrayList<Integer> orderSequence) {
        return 0;
    }

    @Override
    public double getPenalizedCost(ArrayList<Integer> orderSequence, double durationViolationPenalty, double capacityViolationPenalty) {
        return 0;
    }

    @Override
    public double getDurationViolationPenalty() {
        return 0;
    }

    @Override
    public double getCapacityViolationPenalty() {
        return 0;
    }

    @Override
    public double getNumberOfInstallationsPenalty() {
        return 0;
    }

    @Override
    public void setDurationViolationPenalty(double durationViolationPenalty) {

    }

    @Override
    public void setCapacityViolationPenalty(double capacityViolationPenalty) {

    }

    @Override
    public void setNumberOfInstallationsViolationPenalty(double numberOfInstallationsViolationPenalty) {

    }

    @Override
    public void setPenalizedCostPopulation(ArrayList<Individual> population) {

    }
}
