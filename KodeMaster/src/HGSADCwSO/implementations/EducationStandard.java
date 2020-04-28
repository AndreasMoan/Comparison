package HGSADCwSO.implementations;

import HGSADCwSO.*;
import HGSADCwSO.protocols.EducationProtocol;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;

import java.lang.reflect.Array;
import java.util.*;

public class EducationStandard implements EducationProtocol {

    protected ProblemData problemData;
    protected FitnessEvaluationProtocol fitnessEvaluationProtocol;
    protected PenaltyAdjustmentProtocol penaltyAdjustmentProtocol;
    //protected GenoToPhenoConverterProtocol genoToPhenoConverter;
    protected boolean isRepair;
    protected int penaltyMultiplier;


    public EducationStandard(ProblemData problemData, FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        this.problemData = problemData;
        this.fitnessEvaluationProtocol = fitnessEvaluationProtocol;
        this.isRepair = false;
        this.penaltyMultiplier = 1;
    }

    @Override
    public void educate(Individual individual) {

    }

    @Override
    public void repairEducate(Individual individual, int penaltyMultiplyer) {

    }
}
/*
    @Override
    public void educate(Individual individual) {
        Random rand = new Random();
        double randomNumber = rand.nextDouble();

        if (randomNumber > problemData.getHeuristicParameterDouble("Education rate")) {
            //Education:
            voyageImprovement(individual);
            patternImprovement(individual);
            voyageImprovement(individual);

            updatePenaltyAdjustmentCounter(individual);
        }
    }


    protected void updatePenaltyAdjustmentCounter(Individual individual) {
        if (!isRepair) {
            penaltyAdjustmentProtocol.countAddedIndividual(individual);
        }
    }


    @Override
    public void repairEducate(Individual individual, int penaltyMultiplier) {
        isRepair = true;
        this.penaltyMultiplier = penaltyMultiplier;

        educate(individual);

        isRepair = false;
        this.penaltyMultiplier = 1;
    }

    protected void voyageImprovement(Individual individual) {
        //HashMap<Vessel, Voyage> giantTour = individual.getGiantTour(); //voyage.installations indikerer hvilke installasjoner som skal besøkes, i rekkefølgen de ligger i
        HashMap<Integer, ArrayList<Integer>> chromosome = individual.getVesselTourChromosome(); //first key is a vessel number, and the value is a set of installations

        for (int vessel : chromosome.keySet()) {
            //Voyage voyage = giantTour.get(vessel)
            ArrayList<Integer> voyage = chromosome.get(vessel);
            if (voyage != null) {
                ArrayList<Integer> improvedVoyage = getImprovedVoyage(voyage);
                chromosome.put(vessel, new ArrayList<Integer>(improvedVoyage));
            }
        }
    }

    private ArrayList<Integer> getImprovedVoyage(ArrayList<Integer> voyage) {
        ArrayList<Integer> orders = new ArrayList<Integer>(voyage);
        ArrayList<Integer> untreatedOrders = new ArrayList<Integer>(voyage);

        while (untreatedOrders.size() > 0) {
            Integer u = Utilities.pickAndRemoveRandomElementFromList(untreatedOrders);
            ArrayList<Integer> neighbours = getNeighbours(u, orders);
            while (neighbours.size() > 0) {
                Integer v = Utilities.pickAndRemoveRandomElementFromList(neighbours);
                orders = doRandomMove(u, v, orders);
            }
        }
        return new ArrayList<Integer>(orders);
    }


    private ArrayList<Integer> doRandomMove(Integer u, Integer v, ArrayList<Integer> orders) {
        ArrayList<Integer> unusedMoves = new ArrayList<Integer>(); //list of unused moves
        Move move = new Move();
        for (int i = 0; i < move.getNumberOfMoves(); i++) {
            unusedMoves.add(i + 1); //moves are 1-indexed, as in Borthen & Loennechen (2016) and Vidal et al (2012)
        }
        while (unusedMoves.size() > 0) { //try moves in random order
            int moveNumber = Utilities.pickAndRemoveRandomElementFromList(unusedMoves);
            ArrayList<Integer> newOrders = move.doMove(u, v, orders, moveNumber);
            double oldVoyagePenalizedCost;
            double newVoyagePenalizedCost;

            if (!isRepair) { //Normal education
                oldVoyagePenalizedCost = fitnessEvaluationProtocol.getPenalizedCost(orders);
                newVoyagePenalizedCost = fitnessEvaluationProtocol.getPenalizedCost(newOrders);
            } else { //repair education
                double durationViolationPenalty = fitnessEvaluationProtocol.getDurationViolationPenalty() * penaltyMultiplier;
                double capacityViolationPenalty = fitnessEvaluationProtocol.getCapacityViolationPenalty() * penaltyMultiplier;

                //her må det legges inn flere dersom vi relaxer noen andre
                oldVoyagePenalizedCost = fitnessEvaluationProtocol.getPenalizedCost(orders, durationViolationPenalty, capacityViolationPenalty);
                newVoyagePenalizedCost = fitnessEvaluationProtocol.getPenalizedCost(newOrders, durationViolationPenalty, capacityViolationPenalty);
            }
            if (newVoyagePenalizedCost < oldVoyagePenalizedCost) {
                return newOrders;
            }
        }
        return orders;
    }

    private ArrayList<Integer> getNeighbours(Integer order, ArrayList<Integer> orders) {
        double granularityThreshold = problemData.getProblemInstanceParameterDouble("Granularity threshold in RI"); //størrelse av nabolaget?
        int numberOfNeighbours = (int) (orders.size() * granularityThreshold);
        ArrayList<Integer> neighbours = new ArrayList<Integer>(orders); //create set of neighbours
        neighbours.remove(order); //remove this order from the set
        //get the distance from all other installations (orders) to this installation (order) as (key,value)-pairs
        ArrayList<Map.Entry<Order, Double>> distances = new ArrayList<Map.Entry<Order, Double>>(problemData.getDistances().get(problemData.getOrdersByNumber().get(order)).entrySet());  //TODO Hvordan fikse distances i problemdata til å bli riktig format?
        //remove (key,value)-pairs that are not neighbours
        ArrayList<Map.Entry<Order, Double>> removeList = new ArrayList<Map.Entry<Order, Double>>();
        for (Map.Entry<Order, Double> distance : distances) { //removing the order considered from the distances-ArrayList
            if (!neighbours.contains(distance.getKey().getNumber())) {
                removeList.add(distance);
            }
        }
        distances.removeAll(removeList);
        //sort the (key,value)-pairs by distance, having the pairs with the highest distance first
        Collections.sort(distances, Collections.reverseOrder(Utilities.getMapEntryWithDoubleComparator()));
        //removes the neighbours with the highest distance until the correct number of neighbours is obtained
        while (neighbours.size() > numberOfNeighbours) {
            neighbours.remove(Integer.valueOf(distances.remove(0).getKey().getNumber()));
        }
        return neighbours;
    }


    private void patternImprovement(Individual individual) {
        /*for all combinations of two voyages departing that day
        Calculate penalized cost of merging the two voyages (inserting each installation at best position in each voyage)
        end-for
        if best merge has reduced cost
        Choose and perform best merge
        end-if

        //TODO - Et individ består av ett kromosom med ordresekvens pr båt?
        Set<Integer> departingVessels = individual.getDepartingVessels();

        if (departingVessels.size()>1){ //if more than one departure next day. Borthen & Loennechen gjør dette for hver departure day
            mergeVoyages(individual);
        }
    }

    private void mergeVoyages(Individual individual) {
        Set<Integer> departingVessels = new HashSet<Integer> (individual.getDepartingVessels());
        if (departingVessels.size() > 1){
            Set<Set<Integer>> allDepartingVesselCombinations = Utilities.cartesianProduct(departingVessels);

            int bestVesselToKeep = -1;
            int bestVesselToRemove = -1;
            ArrayList<Integer> bestNewVoyage = null;
            double bestCostReduction = 0;

            for (Set<Integer> vesselPair : allDepartingVesselCombinations){
                //Randomly select which vessel to remove and which to keep
                Integer vesselNumberToKeep = Utilities.pickAndRemoveRandomElementFromSet(vesselPair);
                Vessel vesselToKeep = problemData.getVesselByNumber(vesselNumberToKeep);
                Integer vesselNumberToRemove = Utilities.pickAndRemoveRandomElementFromSet(vesselPair);
                Vessel vesselToRemove = problemData.getVesselByNumber(vesselNumberToRemove);

                ArrayList<Integer> voyageToMergeInto = individual.getVesselTourChromosome().get(vesselNumberToKeep); //order sequence of the voyage to keep
                ArrayList<Integer> voyageToMove = individual.getVesselTourChromosome().get(vesselNumberToRemove); //order sequence of the voyage to remove

                double currentPenalizedCost = fitnessEvaluationProtocol.getPenalizedCost(voyageToMergeInto) + fitnessEvaluationProtocol.getPenalizedCost(voyageToMove);

                ArrayList<Integer> voyageSeq = new ArrayList<Integer>(voyageToMergeInto);
                //insert each order in voyageToMove into voyageToMergeInto
                for (Integer order : voyageToMergeInto) {
                    VoyageInsertion bestInsertion = getBestInsertionIntoVoyage(order, voyageToMergeInto);
                    int bestPos = bestInsertion.positionInVoyageToInsertInto;
                    voyageSeq.add(bestPos, order);
                }
                ArrayList<Integer> newVoyage = new ArrayList<Integer>(voyageSeq);
                double newPenalizedCost = fitnessEvaluationProtocol.getPenalizedCost(newVoyage);
                double costReduction = currentPenalizedCost-newPenalizedCost;

                if (costReduction > bestCostReduction){
                    bestVesselToKeep = vesselNumberToKeep; //TODO - Hva skjer på denne og den under?
                    bestVesselToRemove = vesselNumberToRemove;
                    bestNewVoyage = newVoyage;
                    bestCostReduction = costReduction;
                }
            }
            //if there exist a merger resulting in reduced cost, change the tour
            if (bestCostReduction > 0){
            //copy current giant tour
                HashMap<Integer, ArrayList<Integer>> vesselTour = Utilities.deepCopyVesselTour(individual.getVesselTourChromosome());

                vesselTour.put(bestVesselToKeep, bestNewVoyage);
                vesselTour.put(bestVesselToRemove, new ArrayList<>());
                individual.setVesselTourChromosome(vesselTour); // erstatter: Genotype newGenotype = new GenotypeHGS(giantTour, problemData.getCustomerInstallations().size(), problemData.getVessels().size());
                individual.updatePenalizedCostForChromosome(vesselTour); // erstatter: individual.setGenotypeAndUpdatePenalizedCost(newGenotype, genoToPhenoConverter, fitnessEvaluationProtocol);
            }
        }
    }

    */