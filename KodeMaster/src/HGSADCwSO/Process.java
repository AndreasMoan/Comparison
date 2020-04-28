package HGSADCwSO;

import HGSADCwSO.implementations.*;
import HGSADCwSO.implementations.DAG.FitnessEvaluationDAG;
import HGSADCwSO.protocols.EducationProtocol;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;
import HGSADCwSO.protocols.InitialPopulationProtocol;
import HGSADCwSO.protocols.ParentSelectionProtocol;
import HGSADCwSO.protocols.ReproductionProtocol;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

public class Process {

    private ProblemData problemData;

    private InitialPopulationProtocol initialPopulationProtocol;
    private ParentSelectionProtocol parentSelectionProtocol;
    private EducationProtocol educationProtocol;
    private FitnessEvaluationProtocol fitnessEvaluationProtocol;
    private ReproductionProtocol reproductionProtocol;
    private int processIteration;
    private HashMap<Integer, ArrayList<Individual>> feasibleSubPopulationByIteration = new HashMap<Integer, ArrayList<Individual>>();
    private HashMap<Integer, ArrayList<Individual>> infeasibleSubPopulationByIteration = new HashMap<Integer, ArrayList<Individual>>();
    private HashMap<Integer, Individual> bestFeasibleIndividualByIteration = new HashMap<Integer, Individual>();

    public Process(ProblemData problemData) {
        this.problemData = problemData;
        selectProtocols();
    }

    public Individual createIndividual() {
        return initialPopulationProtocol.createIndividual();
    }

    public ArrayList<Individual> selectParents(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) {
        ArrayList<Individual> entirePopulation = Utilities.getAllElements(feasiblePopulation, infeasiblePopulation);
        return parentSelectionProtocol.selectParents(entirePopulation);
    }

    public Individual mate(ArrayList<Individual> parents) {
        Individual kid = reproductionProtocol.crossover(parents);
        return kid;
    }

    public void repair(Individual individual, double probability) {
        if (!individual.isFeasible()) {
            double randomDouble = new Random().nextDouble();
            if (randomDouble < probability) {
                int penaltyMultiplier = 10;
                educationProtocol.repairEducate(individual, penaltyMultiplier);
                if (!individual.isFeasible()) {
                    penaltyMultiplier = 100;
                    educationProtocol.repairEducate(individual, penaltyMultiplier);
                }
            }
        }
    }

    public void evaluate(Individual individual) {
        fitnessEvaluationProtocol.evaluate(individual);
    }

    public void repair(Individual individual) {
        double probability = problemData.getHeuristicParameterDouble("Repair rate");
        repair(individual, probability);
    }

    public void educate(Individual individual) {
        educationProtocol.educate(individual);
    }

    private void selectProtocols() {
        selectFitnessEvaluationProtocol();
        selectInitialPopulationProtocol();
        selectParentSelectionProtocol();
        selectEducationProtocol();
        selectReproductionProtocol();
    }

    private void selectInitialPopulationProtocol() {
        switch (problemData.getHeuristicParameters().get("Initial population protocol")) {
            case "standard":
                initialPopulationProtocol = new InitialPopulationStandard(problemData, fitnessEvaluationProtocol);
                break;
            default:
                initialPopulationProtocol = null;
                break;
        }
    }

    private void selectParentSelectionProtocol(){
        switch (problemData.getHeuristicParameters().get("Parents selection protocol")) {
            case "binary tournament":
                parentSelectionProtocol = new ParentsSelectionBinaryTournament();
                break;
            default:
                parentSelectionProtocol = null;
                break;
        }
    }

    private void selectEducationProtocol(){
        switch (problemData.getHeuristicParameters().get("Education protocol")) {
            case "cost":
                educationProtocol = new EducationStandard(problemData, fitnessEvaluationProtocol);
                break;
            default:
                educationProtocol = null;
                break;
        }
    }

    private void selectFitnessEvaluationProtocol() {
        switch (problemData.getHeuristicParameters().get("Fitness evaluation protocol")) {
            case "standard":
                fitnessEvaluationProtocol = new FitnessEvaluationQuickAndDirty(problemData);
                break;
            case "dag":
                fitnessEvaluationProtocol = new FitnessEvaluationDAG(problemData);
                break;
            default:
                fitnessEvaluationProtocol = null;
                break;
        }
    }

    private void selectReproductionProtocol() {
        switch (problemData.getHeuristicParameters().get("Reproduction protocol")) {
            case "standard":
                reproductionProtocol = new ReproductionStandard(problemData, fitnessEvaluationProtocol);
                break;
            default:
                reproductionProtocol = null;
                break;
        }
    }

    public void updateIterationsSinceImprovementCounter(boolean isImprovingSolution) {
    }

    public void adjustPenaltyParameters(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) {
    }

    public boolean isDiversifyIteration() {
        return false;
    }

    public void recordRunStatistics(int iteration, ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation, Individual bestFeasibleIndividual) {
        processIteration = iteration;
        feasibleSubPopulationByIteration.put(iteration, feasiblePopulation);
        infeasibleSubPopulationByIteration.put(iteration, infeasiblePopulation);
        bestFeasibleIndividualByIteration.put(iteration, bestFeasibleIndividual);

        // System.out.println("Iteration: " + iteration + " Feasible pop size: " + feasiblePopulation.size() + " Infeasible pop size: " + infeasiblePopulation.size());
    }

    public void addDiversityDistance(Individual kid) {
    }

    public void updateBiasedFitness(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) {
    }


}
