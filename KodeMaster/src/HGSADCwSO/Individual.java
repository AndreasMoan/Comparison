package HGSADCwSO;

import HGSADCwSO.protocols.FitnessEvaluationProtocol;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class Individual {

    private double penalizedCost, scheduleCost, durationViolation, capacityViolation, deadlineViolation;
    private Set<Integer> departingVessels;
    private double fitness;

    private Genotype genotype;
    private FitnessEvaluationProtocol fitnessEvaluationProtocol;

    private boolean feasibility;
    private boolean capacityFeasibility;
    private boolean durationFeasibility;
    private boolean deadlineFeasibility;

    private int costRank;
    private int diversityRank;
    private double biasedFitness;
    private double diversityContribution;

    //private HashMap<Integer, ArrayList<Integer>> vesselTourChromosome;

    public Individual(HashMap<Integer, ArrayList<Integer>> vesselTourChromosome, FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        this.genotype = new Genotype(vesselTourChromosome);
        this.fitnessEvaluationProtocol = fitnessEvaluationProtocol;
        fitnessEvaluationProtocol.evaluate(this);
        setDepartingVessels();
    }

    private void setDepartingVessels(){
        Set<Integer> departingSet = new HashSet<>(); //initializing departingVessels
        HashMap<Integer, ArrayList<Integer>> vesselTourChromosome = new HashMap<>(getVesselTourChromosome());
        for (Integer vessel : vesselTourChromosome.keySet()){
            if (vesselTourChromosome.get(vessel).size()>0){
                departingSet.add(vessel);
            }
        }
        this.departingVessels = new HashSet<>(departingSet);
    }

    public Set<Integer> getDepartingVessels() {
        return departingVessels;
    }

    public HashMap<Integer, ArrayList<Integer>> getVesselTourChromosome() { return genotype.getVesselTourChromosome(); }

    public void setVesselTourChromosome(HashMap<Integer, ArrayList<Integer>> vesselTour){
        this.genotype.setVesselTourChromosome(vesselTour);
        setDepartingVessels();
    }

    public Genotype getGenotype() { return genotype; }

    public double getFitness() {
        return fitness;
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    public void setFeasibility(boolean feasibility) {
        this.feasibility = feasibility;
    }

    public boolean isFeasible() {
        return feasibility;
    }

    public double getPenalizedCost() {
        return penalizedCost;
    }

    public void setPenalizedCost(Double cost) {
        penalizedCost = cost;
    }

    /*public void setDepartingVessels(Set<Integer> departingVessels) {
        this.departingVessels = departingVessels;
    }*/

    public void updatePenalizedCostForChromosome(FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        fitnessEvaluationProtocol.setPenalizedCostIndividual(this);
    }

    public double getScheduleCost() { return scheduleCost; }

    public void setScheduleCost(double cost) { this.scheduleCost = cost; }

    public double getDurationViolation() { return durationViolation; }

    public void setDurationViolation(double violation) { this.durationViolation = violation; }

    public double getCapacityViolation() { return capacityViolation; }

    public void setCapacityViolation(double violation) { this.capacityViolation = violation; }

    public double getDeadlineViolation() { return deadlineViolation; }

    public void setDeadlineViolation(double violation) { this.deadlineViolation = violation; }

    public  boolean getCapacityFeasibility() { return capacityFeasibility; }

    /*
    public boolean isCapacityFeasible() { //TODO - Legg til denne koden i fitnessevaluation
        boolean noCapacityViolation = true;

        for (Integer vessel : vesselTourChromosome.keySet()){
            Vessel vessel1 = problemData.getVesselByNumber().get(vessel);
            int totalLoad = 0;

            for (Integer order : vesselTourChromosome.get(vessel)){
                Order order1 = problemData.getOrdersByNumber().get(order);
                totalLoad += order1.getDemand();
            }
            if (totalLoad > vessel1.getCapacity()){noCapacityViolation = false;}
        }
        this.capacityFeasibility = noCapacityViolation;
        return capacityFeasibility;
    }
     */

    public boolean getDurationFeasibility() { return true;} //TODO return durationFeasibility

    public boolean getDeadlineFeasibility() { return true;}//TODO}return deadlineFeasibility;

    public void setCostRank(int costRank) { this.costRank = costRank; }

    public int getCostRank(){return costRank;}

    public void setBiasedFitness(double biasedFitness) { this.biasedFitness = biasedFitness; }

    public double getBiasedFitness(){return biasedFitness;}

    public void setDiversityContribution(double diversityContribution){ this.diversityContribution = diversityContribution;}

    public double getDiversityContribution() { return diversityContribution; }

    public void setDiversityRank(int diversityRank) { this.diversityRank = diversityRank; }

    public int getDiversityRank() { return diversityRank; }

    /*public double getFactorToIncreaseUpdatedPenalitiesWith(){
        return
    }

    public void setGetFactorToDecreaseUpdatedPenalitiesWith(){
        return
    }*/
}

