package HGSADCwSO.implementations;

import java.util.ArrayList;


public class VoyageInsertion {

    public int vessel, positionInVoyageToInsertInto, orderNumber;
    public double insertionCost;

    public VoyageInsertion(int vessel, int orderNumber, int positionInVoyageToInsertInto, double insertionCost){

        this.orderNumber = orderNumber;
        this.positionInVoyageToInsertInto = positionInVoyageToInsertInto;
        this.insertionCost = insertionCost;
        this.vessel = vessel;
    }

    public int getVessel() { return vessel; }

    public int getOrderNumber() { return orderNumber; }

    public int getPositionInVoyageToInsertInto() { return positionInVoyageToInsertInto; }

    public double getInsertionCost() { return insertionCost; }

    public static double getTotalInsertionCost(ArrayList<VoyageInsertion> insertions){
        double sum = 0;
        for (VoyageInsertion insertion : insertions){
            sum += insertion.insertionCost;
        }
        return sum;
    }
}
