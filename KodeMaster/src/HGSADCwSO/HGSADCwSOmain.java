package HGSADCwSO;


import javax.sound.midi.Soundbank;
import java.sql.SQLOutput;
import java.util.ArrayList;

public class HGSADCwSOmain {

    private ArrayList<Individual> feasiblePopulation, infeasiblePopulation;
    private Individual bestFeasibleIndividual;
    private long startTime, stopTime;
    private ProblemData problemData;
    private Process process;

    private double bestCost = Double.POSITIVE_INFINITY;

    private String[] args;

    private int iteration;

    private double repairChance;
    private double educationChance;

    public static int COST_EDUCATIONS = 0;
    public static int INFEASIBLE_EDUCATIONS = 0;

    public static void main(String[] args) {
        HGSADCwSOmain main = new HGSADCwSOmain();
        main.initialize(args);
        main.fullEvolutionaryRun();
    }


    private void initialize(String[] changeParameters) {
        startTime = System.nanoTime();
        //io = new IO(inputFileName); //TODO
        //this.args = changeParameter;
        problemData = HackInitProblemData.hack();
    }

    private void fullEvolutionaryRun(){

        process = new Process(problemData);
        feasiblePopulation = new ArrayList<Individual>();
        infeasiblePopulation = new ArrayList<Individual>();
        bestFeasibleIndividual = null;

        iteration = 1;
        problemData.printProblemData(); //TODO
        System.out.println("Creating initial population...");
        createInitialPopulation();

        process.updateIterationsSinceImprovementCounter(true);

        runEvolutionaryLoop();
        terminate();
    }

    private void createInitialPopulation() {
        int initialPopulationSize = 100; //TODO
        for (int i = 0; i < initialPopulationSize; i++){
            Individual kid = process.createIndividual();
            process.educate(kid);
            if (! kid.isFeasible()) {
                process.repair(kid, repairChance);
            }

            addToSubpopulation(kid);
        }
    }

    private void runEvolutionaryLoop() {
        process.recordRunStatistics(0, feasiblePopulation, infeasiblePopulation, bestFeasibleIndividual);
        while (!stoppingCriterion()) {
            System.out.println("Iteration " + iteration + "                                                          Best cost thus far: " + bestCost);
            evolve();
            process.recordRunStatistics(iteration, feasiblePopulation, infeasiblePopulation, bestFeasibleIndividual);
            iteration++;
        }
    }

    private boolean stoppingCriterion() {
        //TODO
        return false;
    }

    private void evolve() {
        ArrayList<Individual> parents = process.selectParents(feasiblePopulation, infeasiblePopulation);
        Individual kid = process.mate(parents);
        process.educate(kid);
        process.repair(kid);
        boolean isImprovingSolution = addToSubpopulation(kid);
        if (kid.getFitness() < bestCost) {
            bestCost = kid.getFitness();
        }
        process.updateIterationsSinceImprovementCounter(isImprovingSolution);
        process.adjustPenaltyParameters(feasiblePopulation, infeasiblePopulation);
        /*
        if (process.isDiversifyIteration()) {
            diversify(feasiblePopulation, infeasiblePopulation);
        }
         */
    }

    private void diversify(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) {

        System.out.println("Diversifying...");

        genocide(feasiblePopulation, 2.0/3.0);
        genocide(infeasiblePopulation, 2.0/3.0);

        System.out.println("Breeding new population...");
    }

    private void genocide(ArrayList<Individual> population, double proportionToKill){

        System.out.println(population.size());

        population.sort(Utilities.getFitnessComparator());
        int numberOfIndividualsToKill = (int) Math.round(population.size()*proportionToKill);
        ArrayList<Individual> individualsToBeKilled = new ArrayList<>(population.subList(0,numberOfIndividualsToKill));
        for (Individual individual : individualsToBeKilled) {
            removeFromSubpopulation(individual, population);
        }
    }

    private void removeFromSubpopulation(Individual individual, ArrayList<Individual> population) {
        population.remove(individual);
    }

    public boolean addToSubpopulation(Individual kid) {
        boolean isImprovingSolution = false;
        process.evaluate(kid);
        System.out.println(kid.getGenotype().getVesselTourChromosome());

        if (kid.isFeasible()) {
            if ((bestFeasibleIndividual == null) || (kid.getFitness() < bestFeasibleIndividual.getFitness())) {
                bestFeasibleIndividual = kid;
                isImprovingSolution = true;
            }
            else {
                isImprovingSolution = false;
            }
            feasiblePopulation.add(kid);
        }
        else {
            isImprovingSolution = false;
            infeasiblePopulation.add(kid);
        }
        process.addDiversityDistance(kid);
        process.updateBiasedFitness(feasiblePopulation, infeasiblePopulation);

        if (kid.isFeasible()) {
            checkSubpopulationSize(feasiblePopulation, infeasiblePopulation);
        }
        else {
            checkSubpopulationSize(infeasiblePopulation, feasiblePopulation);
        }
        return isImprovingSolution;
    }

    private void checkSubpopulationSize(ArrayList<Individual> subpopulation, ArrayList<Individual> otherSubpopulation) {
        int maxPopulationSize = problemData.getHeuristicParameterInt("Population size")
                + problemData.getHeuristicParameterInt("Number of offspring in a generation");

        if (subpopulation.size() + otherSubpopulation.size() >= maxPopulationSize) {
            genocide(subpopulation, 3.0/4.0);
            genocide(otherSubpopulation, 3.0/4.0);
        }
    }

    private void terminate() {
        //TODO
    }


}
