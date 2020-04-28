package HGSADCwSO;

import java.util.ArrayList;
import java.util.HashMap;

public class Genotype {


    private HashMap<Integer, ArrayList<Integer>> vesselTourChromosome;

    public Genotype(HashMap<Integer, ArrayList<Integer>> vesselTourChromosome) {

        this.vesselTourChromosome = vesselTourChromosome;

    }

    public void setVesselTourChromosome(HashMap<Integer, ArrayList<Integer>> vesselTourChromosome) {
        this.vesselTourChromosome = vesselTourChromosome;
    }

    public HashMap<Integer, ArrayList<Integer>> getVesselTourChromosome() {
        return vesselTourChromosome;
    }
}
