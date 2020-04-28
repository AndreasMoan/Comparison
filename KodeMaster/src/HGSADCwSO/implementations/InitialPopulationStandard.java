package HGSADCwSO.implementations;

import HGSADCwSO.Individual;
import HGSADCwSO.Order;
import HGSADCwSO.ProblemData;
import HGSADCwSO.Vessel;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;
import HGSADCwSO.protocols.InitialPopulationProtocol;

import java.util.*;

public class InitialPopulationStandard implements InitialPopulationProtocol {

    private ProblemData problemData;
    private int numberOfRestarts;
    private FitnessEvaluationProtocol fitnessEvaluationProtocol;

    public InitialPopulationStandard(ProblemData problemData, FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        this.problemData = problemData;
        this.fitnessEvaluationProtocol = fitnessEvaluationProtocol;
        System.out.println("1: " + fitnessEvaluationProtocol);
        numberOfRestarts = 0;
    }

    public Individual createIndividual(){
        HashMap<Integer, Set<Integer>> vesselOrderChromsome = createVesselOrderChromosome();
        HashMap<Integer, ArrayList<Integer>> vesselTourChromosome = createVesselTourChromosome(vesselOrderChromsome);
        Individual kid = new Individual(vesselTourChromosome, fitnessEvaluationProtocol);
        return kid;
    }

    public int getNumberOfConstructionHeuristicRestarts() {
        return 0;
    }


    private HashMap<Integer, Set<Integer>> createVesselOrderChromosome() {

        HashMap<Integer, Set<Integer>> vesselOrderChromosome = new HashMap<Integer, Set<Integer>>();
        HashMap<Integer, Order> orders = problemData.getOrdersByNumber();

        for (int vessel = 0; vessel < problemData.getNumberOfVessels(); vessel++){
            vesselOrderChromosome.put(vessel, new HashSet<>());
        }

        int n = 1;
        while (n < problemData.getNumberOfOrders() && orders.get(n).getDay() == 0) {
            int randomVessel = new Random().nextInt(problemData.getNumberOfVessels());
            vesselOrderChromosome.get(randomVessel).add(n);
            n++;
        }
        return vesselOrderChromosome;
    }

    private HashMap<Integer, ArrayList<Integer>> createVesselTourChromosome(HashMap<Integer,Set<Integer>> vesselOrderChromosome) {

        HashMap<Integer, ArrayList<Integer>> vesselTourChromosome = new HashMap<Integer, ArrayList<Integer>>();

        for (int i = 0; i < problemData.getNumberOfVessels(); i++) {

            vesselTourChromosome.put(i, new ArrayList<Integer>());
            for (int j : vesselOrderChromosome.get(i)) {
                vesselTourChromosome.get(i).add(j);
            }
            if (vesselOrderChromosome.get(i).size() > 0) {
                Collections.shuffle(vesselTourChromosome.get(i));
            }
        }
        System.out.println(vesselTourChromosome);
        return vesselTourChromosome;
    }
}
