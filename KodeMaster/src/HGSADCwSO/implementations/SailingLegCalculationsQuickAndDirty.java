package HGSADCwSO.implementations;

import HGSADCwSO.ProblemData;
import HGSADCwSO.Utilities;
import HGSADCwSO.protocols.SailingLegCalculationsProtocol;

public class SailingLegCalculationsQuickAndDirty implements SailingLegCalculationsProtocol {

    private double arrivalTime;
    private double fuelConsumption;
    private ProblemData problemData;


    public SailingLegCalculationsQuickAndDirty(ProblemData problemData) {
        this.problemData = problemData;
    }

    public void calculateSailingLeg(Double distance, Double speed, Double startTime, Integer demand) {

        //Sailing:

        double time = startTime;
        double distanceLeft = distance;
        double consumption = 0.0;
        while (distanceLeft - speed > 0) {
            if (time - Math.floor(time) > 0) {
                consumption += getTimeConsumption(speed, time, time - Math.floor(time));
                distanceLeft -= (1 - (time - Math.floor(time))) * speed;
                time = Math.ceil(time);
            }
            else {
                consumption += getTimeConsumption(speed, time);
                distanceLeft -= speed;
                time++;
            }
        }
        if (distanceLeft > 0) {
            double timeLeft = distance / speed;
            consumption += getTimeConsumption(speed, time, timeLeft);
            time += timeLeft;
        }

        //Servicing:

        double consumptionFactor = problemData.getProblemInstanceParameterDouble("Servicing consumption");
        double servicingTimeLeft = (double)demand/10;

        while (servicingTimeLeft -1 > 0) {
            if (time - Math.floor(time) > 0) {
                consumption += consumptionFactor * (time - Math.floor(time)) * problemData.getWeatherImpactByHour((int) Math.floor(time));
                servicingTimeLeft -= (1 - (time - Math.floor(time)));
                time = Math.ceil(time);
            } else {
                consumption += consumptionFactor * problemData.getWeatherImpactByHour((int) Math.floor(time));
                servicingTimeLeft--;
                time++;
            }
        }
        if (servicingTimeLeft > 0){
            consumption += consumptionFactor*servicingTimeLeft*problemData.getWeatherImpactByHour((int)Math.floor(time));
            time += servicingTimeLeft;
        }


        this.fuelConsumption = consumption;
        this.arrivalTime = time;
    }

    public double getTimeConsumption(Double speed, Double timePeriod, Double duration) {
        double consumptionPerHour = (0.8125*speed*speed-13.00*speed+72.75);
        double weatherImpact = problemData.getWeatherImpactByHour((int)Math.floor(timePeriod));
        return consumptionPerHour*duration*weatherImpact;
    }

    public double getTimeConsumption(Double speed, Double timePeriod) {
        return getTimeConsumption(speed,timePeriod,1.0);
    }

    public double getArrivalTime(){
        return arrivalTime;
    }

    public double getFuelConsumption(){
        return fuelConsumption;
    }
}
