package HGSADCwSO.implementations;

import HGSADCwSO.*;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;
import HGSADCwSO.protocols.SailingLegCalculationsProtocol;

import java.util.ArrayList;

public class FitnessEvaluationQuickAndDirty implements FitnessEvaluationProtocol {

    private ProblemData problemData;
    private SailingLegCalculationsProtocol sailingLegCalculationsProtocol;
    private double value;

    public FitnessEvaluationQuickAndDirty(ProblemData problemData){
        this.problemData = problemData;
        selectProtocols();
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

    private void selectProtocols() {
        selectSailingLegCalculationsProtocol();
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
