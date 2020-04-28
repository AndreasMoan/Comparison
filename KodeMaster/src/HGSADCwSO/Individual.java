package HGSADCwSO;

import HGSADCwSO.protocols.FitnessEvaluationProtocol;


import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

public class Individual {

    private double penalizedCost, scheduleCost, durationViolation, capacityViolation;
    private Set<Integer> departingVessels;
    private double fitness;

    private boolean feasibility;

    private Genotype genotype;
    private FitnessEvaluationProtocol fitnessEvaluationProtocol;

    public Individual(HashMap<Integer, ArrayList<Integer>>  vesselTourChromosome, FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        this.genotype = new Genotype(vesselTourChromosome);
        this.fitnessEvaluationProtocol = fitnessEvaluationProtocol;

        System.out.println(fitnessEvaluationProtocol);
        this.departingVessels = genotype.getVesselTourChromosome().keySet();
    }

    public Set<Integer> getDepartingVessels() {return departingVessels;}

    public HashMap<Integer, ArrayList<Integer>> getVesselTourChromosome() {
        return genotype.getVesselTourChromosome();
    }

    public void setVesselTourChromosome(HashMap<Integer, ArrayList<Integer>> vesselTour){
        this.genotype.setVesselTourChromosome(vesselTour);
    }

    public Genotype getGenotype() {
        return genotype;
    }

    public double getFitness() {
        return fitness;
    }

    public void setFeasibility(boolean feasibility) {
        this.feasibility = feasibility;
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    public boolean isFeasible(){
        return feasibility;
    }

    public double getPenalizedCost() {
        return penalizedCost;
    }

    public void setPenalizedCost(Double cost) {
        penalizedCost = cost;
    }

    /*
    public void updatePenalizedCostForChromosome(FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        fitnessEvaluationProtocol.setPenalizedCostIndividual(this);
    }
    */

    public double getScheduleCost() {return scheduleCost; }

    //public void setScheduleCost() {this.scheduleCost = fitnessEvaluationProtocol.evaluate(this);}

    public double getDurationViolation() {return durationViolation; }

    public void setDurationViolation(double violation) {this.durationViolation = violation; }

    public double getCapacityViolation() {return capacityViolation; }

    public void setCapacityViolation(double violation) {this.capacityViolation = violation; }
}
