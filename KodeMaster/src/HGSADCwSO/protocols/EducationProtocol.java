package HGSADCwSO.protocols;

import HGSADCwSO.Individual;

public interface EducationProtocol {

    public void educate(Individual individual);

    public void repairEducate(Individual individual, int penaltyMultiplyer);
}
