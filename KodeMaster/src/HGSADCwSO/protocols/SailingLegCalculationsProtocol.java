package HGSADCwSO.protocols;

public interface SailingLegCalculationsProtocol {

    public void calculateSailingLeg(Double distance, Double speed, Double startTime, Integer demand);

    public double getArrivalTime();

    public double getFuelConsumption();

}
