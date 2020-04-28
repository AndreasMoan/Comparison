package HGSADCwSO.implementations;

import HGSADCwSO.*;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;
import HGSADCwSO.protocols.SailingLegCalculationsProtocol;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class FitnessEvaluationQuickAndDirty implements FitnessEvaluationProtocol {

    //private Individual individual;
    private ProblemData problemData;
    private SailingLegCalculationsProtocol sailingLegCalculationsProtocol;
    private double value;

    private double nCloseProp;
    protected double nEliteProp;
    private double numberOfOrders;
    private double durationViolationPenalty;
    private double capacityViolationPenalty;
    private double deadlineViolationPenalty;
    private HashMap<Individual, HashMap<Individual, Double>> hammingDistances;

    public FitnessEvaluationQuickAndDirty(ProblemData problemData){
        this.problemData = problemData;
        selectProtocols();
        this.hammingDistances = new HashMap<Individual, HashMap<Individual,Double>>();
        this.nCloseProp = problemData.getHeuristicParameterDouble("Proportion of individuals considered for distance evaluation");
        this.nEliteProp = problemData.getHeuristicParameterDouble("Proportion of elite individuals");
        this.numberOfOrders = problemData.getNumberOfOrders();

        this.capacityViolationPenalty = problemData.getHeuristicParameterDouble("Capacity constraint violation penalty");
        this.durationViolationPenalty = problemData.getHeuristicParameterDouble("Duration constraint violation penalty");
        this.deadlineViolationPenalty = problemData.getHeuristicParameterDouble("Deadline constraint violation penalty");
    }

    public void evaluate(Individual individual) {
        Genotype genotype = individual.getGenotype();
        int nVessels = problemData.getNumberOfVessels();
        double cost = 0;
        for (int i = 0; i < nVessels; i++){
            ArrayList<Integer> route = genotype.getVesselTourChromosome().get(i);
            Vessel vessel = problemData.getVesselByNumber().get(i);
            cost += evaluateRoute(route, vessel);
        }

        individual.setFitness(cost);
        individual.setFeasibility(true);
    }

    public double evaluateRoute(ArrayList<Integer> route, Vessel vessel) {
        double totalConsumption = 0;
        double timeHorizon = 24*vessel.getReturnDay();
        double totalNumberOfHiv = 0;
        double totalDistance = 0;
        Installation departureInstallation = problemData.getInstallationByNumber().get(0);
        Installation destinationInstallation;
        double fuelPrice = problemData.getProblemInstanceParameterDouble("fuel price");

        for (int i : route){
            totalNumberOfHiv += problemData.getOrdersByNumber().get(i).getDemand();
            destinationInstallation = problemData.getOrdersByNumber().get(i).getInstallation();
            totalDistance += problemData.getDistance(departureInstallation, destinationInstallation);
            departureInstallation = destinationInstallation;
        }
        destinationInstallation = problemData.getInstallationByNumber().get(0);
        totalDistance += problemData.getDistance(departureInstallation, destinationInstallation);
        departureInstallation = destinationInstallation;

        double sailingTime = timeHorizon - totalNumberOfHiv/10;
        double sailingSpeed = totalDistance/sailingTime;

        double timeSpent = 0;

        for (int i : route){
            destinationInstallation = problemData.getOrdersByNumber().get(i).getInstallation();
            int weatherState = problemData.getWeatherStateByHour().get((int)timeSpent);
            double distance = problemData.getDistance(departureInstallation, destinationInstallation);

            sailingLegCalculationsProtocol.calculateSailingLeg(distance, sailingSpeed, timeSpent, problemData.getOrdersByNumber().get(i).getDemand());

            timeSpent = sailingLegCalculationsProtocol.getArrivalTime();
            totalConsumption += sailingLegCalculationsProtocol.getFuelConsumption();

            departureInstallation = destinationInstallation;
        }
        return fuelPrice*totalConsumption;
    }

    private void selectProtocols() { selectSailingLegCalculationsProtocol();
    }

    private void selectSailingLegCalculationsProtocol(){
        switch (problemData.getHeuristicParameters().get("Sailing leg calculations protocol")){
            case "quick and dirty":
                sailingLegCalculationsProtocol = new SailingLegCalculationsQuickAndDirty(problemData);
                break;
            default:
                sailingLegCalculationsProtocol = null;
                break;
        }
    }

    public double getValue() {
        return value;
    }


    @Override
    public void updateBiasedFitness(ArrayList<Individual> population) {
        updateDiversityContribution(population);
        updatePenalizedCostRank(population);
        updateDiversityContributionRank(population);
        calculateBiasedFitness(population);
    }

    protected void updateDiversityContribution(ArrayList<Individual> population) {
        int nClose = (int) (population.size() * nCloseProp);
        for (Individual individual : population){
            ArrayList<Individual> nClosest = getnClosest(individual, nClose);
            double distance = 0;
            for (Individual close : nClosest){
                distance += hammingDistances.get(individual).get(close);
            }
            double diversityContribution = distance / nClose;
            individual.setDiversityContribution(diversityContribution);
        }
    }

    protected void updatePenalizedCostRank(ArrayList<Individual> population) {
        Collections.sort(population, Utilities.getPenalizedCostComparator());
        for (int i = 0; i < population.size(); i++){
            Individual individual = population.get(i);
            individual.setCostRank(i+1);
        }
    }

    protected void updateDiversityContributionRank(ArrayList<Individual> population) {
        Collections.sort(population, Utilities.getDiversityContributionComparator());
        for(int i = 0; i < population.size(); i++){
            Individual individual = population.get(i);
            individual.setDiversityRank(i+1);
        }
    }

    protected void calculateBiasedFitness(ArrayList<Individual> population) {
        int nIndividuals = population.size();
        double nElite = nIndividuals * nEliteProp;
        for (Individual individual : population){
            double biasedFitness = individual.getCostRank() + (1 - (nElite/nIndividuals)) * individual.getDiversityRank();
            individual.setBiasedFitness(biasedFitness);
        }
    }

    //The normalized Hamming distance counts the number of orders that are delivered with a different PSVs. Bør man kanskje også telle antall ganger ordre blir levert på ulike dager? Tenker kanskje dette kan droppes i og med at man har deadlines. Da blir jo naturligvis de tidligste fristene besøkt først
    @Override
    public double getHammingDistance(Individual individual1, Individual individual2) {
        return hammingDistances.get(individual1).get(individual2);
    }

    @Override
    public void addDiversityDistance(Individual individual) {
        HashMap<Individual, Double> individualHammingDistances = new HashMap<Individual, Double>();
        for (Individual existingIndividual : hammingDistances.keySet()){
            HashMap<Individual, Double> existingIndividualDistances = hammingDistances.get(existingIndividual);
            double normalizedHammingDistance = getNormalizedHammingDistance(individual, existingIndividual);
            existingIndividualDistances.put(individual, normalizedHammingDistance);
            individualHammingDistances.put(existingIndividual, normalizedHammingDistance);
            hammingDistances.put(existingIndividual, existingIndividualDistances);
        }
        hammingDistances.put(individual, individualHammingDistances);
    }

    @Override
    public void removeDiversityDistance(Individual individual) {
        hammingDistances.remove(individual);
        for (Individual otherIndividual : hammingDistances.keySet()){
            hammingDistances.get(otherIndividual).remove(individual);
        }
    }

    public double getNormalizedHammingDistance(Individual individual1, Individual individual2) { //The normalized Hamming distance counts the number of orders that are delivered with different PSVs.
        HashMap<Integer, ArrayList<Integer>> chromosome1 = individual1.getVesselTourChromosome();
        HashMap<Integer, ArrayList<Integer>> chromosome2 = individual2.getVesselTourChromosome();

        int vesselDifference = 0;
        for (int vessel : chromosome1.keySet()){
            for (int order : chromosome1.get(vessel)){
                if (!chromosome2.get(vessel).contains(order)){
                    vesselDifference++;
                }
            }
        }

        int voyageDifference = 0;
        for (int vessel : chromosome1.keySet())
            for(int index = 0; index < chromosome1.get(vessel).size(); index++){
                if (chromosome1.get(vessel).size() > chromosome2.get(vessel).size()) {
                    if (index < chromosome2.get(vessel).size()-1){
                        if (chromosome2.get(vessel).get(index) != chromosome1.get(vessel).get(index)) {
                            voyageDifference++;
                        }
                    }
                    else if (index >= chromosome2.get(vessel).size()-1) {
                        voyageDifference += chromosome1.get(vessel).size() - chromosome2.get(vessel).size();
                    }
                }
                else {
                    if (chromosome2.get(vessel).get(index) != chromosome1.get(vessel).get(index)) {
                        voyageDifference++;
                    }
                }
            }
        return (vesselDifference + voyageDifference)/ (2*numberOfOrders);
    }

    private ArrayList<Individual> getnClosest(Individual individual, int nClose){
        ArrayList<Individual> nClosest = new ArrayList<Individual>();
        ArrayList<Map.Entry<Individual, Double>> otherIndividuals = new ArrayList<Map.Entry<Individual, Double>>(hammingDistances.get(individual).entrySet());
        Collections.sort(otherIndividuals, Utilities.getMapEntryWithDoubleComparator());
        for (int i = 0; i < nClose; i++){
            nClosest.add(otherIndividuals.get(i).getKey());
        }
        return nClosest;
    }

    @Override
    public void setPenalizedCostIndividual(Individual individual, double durationViolationPenalty, double capacityViolationPenalty, double deadlineViolationPenalty) {

    }

    @Override
    public void setPenalizedCostIndividual(Individual individual) {

    }

    @Override
    public double getPenalizedCostOfVoyage(ArrayList<Integer> orderSequence) {
        return getPenalizedCostOfVoyage(orderSequence, durationViolationPenalty, capacityViolationPenalty, deadlineViolationPenalty);
    }

    @Override
    public double getPenalizedCostOfVoyage(ArrayList<Integer> orderSequence, double durationViolationPenalty, double capacityViolationPenalty, double deadlineViolationPenalty) {
        return 0; //TODO
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
    public double getDeadlineViolationPenalty() {
        return 0;
    }

    @Override
    public void setDurationViolationPenalty(double durationViolationPenalty) {

    }

    @Override
    public void setCapacityViolationPenalty(double capacityViolationPenalty) {

    }

    @Override
    public void setDeadlineViolationPenalty(double deadlineViolationPenalty) {

    }

    @Override
    public void setPenalizedCostPopulation(ArrayList<Individual> population) {

    }

}
