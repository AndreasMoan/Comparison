package HGSADCwSO.implementations;

import HGSADCwSO.Individual;
import HGSADCwSO.Utilities;
import HGSADCwSO.protocols.ParentSelectionProtocol;

import java.util.ArrayList;

public class ParentsSelectionBinaryTournament implements ParentSelectionProtocol {

    @Override
    public ArrayList<Individual> selectParents(ArrayList<Individual> population) {
        ArrayList<Individual> parents = new ArrayList<Individual>();
        Individual father = pickParent(population);
        Individual mother = pickParent(population);
        parents.add(father);
        parents.add(mother);
        return parents;
    }

    private Individual pickParent(ArrayList<Individual> population) {
        Individual parentCandidate1 = Utilities.pickRandomElementFromList(population);
        Individual parentCandidate2 = Utilities.pickRandomElementFromList(population);
        if (parentCandidate1.getFitness() <= parentCandidate2.getFitness()) {
            return parentCandidate1;
        }
        else {
            return parentCandidate2;
        }
    }

}
