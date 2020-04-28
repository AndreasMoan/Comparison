package HGSADCwSO.implementations;

import HGSADCwSO.Individual;
import HGSADCwSO.ProblemData;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;

import java.util.ArrayList;

public class PenaltyAdjustmentProtocol {

    private double solutionsSincePenaltyAdjustment;

    private double capacityFeasibleSolutions, durationFeasibleSolutions, deadlineFeasibleSolutions;

    private double targetFeasibleProportion;

    private double factorToIncreasePenaltiesWith, factorToDecreasePenaltiesWith;


    public PenaltyAdjustmentProtocol(double targetFeasibleProportion, ProblemData problemData) {
        this.targetFeasibleProportion = targetFeasibleProportion;
        this.solutionsSincePenaltyAdjustment = 0;
        this.capacityFeasibleSolutions = 0;
        this.durationFeasibleSolutions = 0;
        this.deadlineFeasibleSolutions = 0;

        this.factorToIncreasePenaltiesWith = problemData.getHeuristicParameterDouble("Factor to increase penalties with");
        this.factorToDecreasePenaltiesWith = problemData.getHeuristicParameterDouble("Factor to decrease penalties with");
    }

    public void countAddedIndividual(Individual individual) {
        solutionsSincePenaltyAdjustment++;

        /*if (individual.isCapacityFeasible()){ //TODO flytt til stedet vi vil ha den. ligger foreløpig kommentert ut i Individualklassen
            capacityFeasibleSolutions++;
        }
        if (individual.isDurationFeasible()){ //TODO fiks funksjonen i individual-klassen
            durationFeasibleSolutions++;
        }
        if (individual.isDeadlineFeasible()){ //TODO fiks funksjonen i individual-klassen
            deadlineFeasibleSolutions++;
        }
            //System.out.println("Number of solutions since penalty adjustment: " + solutionsSincePenaltyAdjustment);*/
    }


    public void adjustPenalties(ArrayList<Individual> entirePopulation, FitnessEvaluationProtocol fitnessProtocol){
        if (solutionsSincePenaltyAdjustment >= 100){
            System.out.println("Adjusting penalty parameters...");

            adjustCapacityPenalty(fitnessProtocol);
            adjustDurationPenalty(fitnessProtocol);
            adjustDeadlinePenalty(fitnessProtocol);

            System.out.println("New duration penalty: " + fitnessProtocol.getDurationViolationPenalty());
            System.out.println("New capacity penalty: " + fitnessProtocol.getCapacityViolationPenalty());
            System.out.println("New deadline penalty: " + fitnessProtocol.getDeadlineViolationPenalty());

            //Resetting counting variables
            solutionsSincePenaltyAdjustment = 0;
            capacityFeasibleSolutions = 0;
            durationFeasibleSolutions = 0;
            deadlineFeasibleSolutions = 0;

            fitnessProtocol.setPenalizedCostPopulation(entirePopulation);
        }
    }


    private void adjustCapacityPenalty(FitnessEvaluationProtocol fitnessProtocol){
        double capacityFeasibleProportion = capacityFeasibleSolutions/solutionsSincePenaltyAdjustment;
        if (capacityFeasibleProportion <= targetFeasibleProportion - 0.05){
            double currentPenalty = fitnessProtocol.getCapacityViolationPenalty();
            fitnessProtocol.setCapacityViolationPenalty(currentPenalty * factorToIncreasePenaltiesWith);
        }
        else if (capacityFeasibleProportion >= targetFeasibleProportion +0.05){
            double currentPenalty = fitnessProtocol.getCapacityViolationPenalty();
            fitnessProtocol.setCapacityViolationPenalty(currentPenalty * factorToDecreasePenaltiesWith);
        }
    }

    private void adjustDurationPenalty(FitnessEvaluationProtocol fitnessProtocol){
        double durationFeasibleProportion = durationFeasibleSolutions/solutionsSincePenaltyAdjustment;
        if (durationFeasibleProportion <= targetFeasibleProportion - 0.05){
            double currentPenalty = fitnessProtocol.getDurationViolationPenalty();
            fitnessProtocol.setDurationViolationPenalty(currentPenalty * factorToIncreasePenaltiesWith);
        }
        else if (durationFeasibleProportion >= targetFeasibleProportion +0.05){
            double currentPenalty = fitnessProtocol.getDurationViolationPenalty();
            fitnessProtocol.setDurationViolationPenalty(currentPenalty * factorToDecreasePenaltiesWith);
        }
    }

    private void adjustDeadlinePenalty(FitnessEvaluationProtocol fitnessProtocol){
        double deadlineFeasibleProportion = deadlineFeasibleSolutions/solutionsSincePenaltyAdjustment;
        if (deadlineFeasibleProportion <= targetFeasibleProportion - 0.05){
            double currentPenalty = fitnessProtocol.getDeadlineViolationPenalty();
            fitnessProtocol.setDeadlineViolationPenalty(currentPenalty * factorToIncreasePenaltiesWith);
        }
        else if (deadlineFeasibleProportion >= targetFeasibleProportion +0.05){
            double currentPenalty = fitnessProtocol.getDeadlineViolationPenalty();
            fitnessProtocol.setDurationViolationPenalty(currentPenalty * factorToDecreasePenaltiesWith);
        }
    }


}
