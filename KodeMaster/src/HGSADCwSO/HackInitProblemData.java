package HGSADCwSO;

import com.sun.org.apache.xpath.internal.operations.Or;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class HackInitProblemData {

    public static ProblemData hack(){

        HashMap<String, String> problemInstanceParameters = new HashMap<String, String>();

        problemInstanceParameters.put("fuel price","10,0");
        problemInstanceParameters.put("Servicing consumption","0,4");
        problemInstanceParameters.put("Idling consumption","0,25");
        problemInstanceParameters.put("Max speed","14,0");
        problemInstanceParameters.put("Min speed","7,0");
        problemInstanceParameters.put("Impact on sailing from weather state 2","2,0");
        problemInstanceParameters.put("Impact on sailing from weather state 3","3,0");
        problemInstanceParameters.put("Time per hiv","0,1");
        problemInstanceParameters.put("","");
        problemInstanceParameters.put("","");

        HashMap<String, String> heuristicInstanceParameters = new HashMap<String, String>();

        heuristicInstanceParameters.put("Initial population protocol","standard");
        heuristicInstanceParameters.put("Parents selection protocol","binary tournament");
        heuristicInstanceParameters.put("Education protocol","cost");
        heuristicInstanceParameters.put("Reproduction protocol","standard");
        heuristicInstanceParameters.put("Fitness evaluation protocol","dag");
        heuristicInstanceParameters.put("Population size","25");
        heuristicInstanceParameters.put("Number of offspring in a generation","75");
        heuristicInstanceParameters.put("Education rate", "0.75");
        heuristicInstanceParameters.put("Repair rate", "0.75");
        heuristicInstanceParameters.put("Number of time periods per hour", "1");
        heuristicInstanceParameters.put("Sailing leg calculations protocol","quick and dirty");
        heuristicInstanceParameters.put("Cost penalty per excessive time period","10");
        heuristicInstanceParameters.put("","");
        heuristicInstanceParameters.put("","");
        heuristicInstanceParameters.put("","");
        heuristicInstanceParameters.put("","");



        ArrayList<Vessel> vessels = new ArrayList<Vessel>();

        vessels.add(new Vessel("PSV1", 100, 7, 14, 0, 0, 3));
        vessels.add(new Vessel("PSV2", 100, 7, 14, 0, 1, 2));

        ArrayList<Installation> installations = new ArrayList<Installation>();

        installations.add(new Installation("DEP", 0,23, 0));
        installations.add(new Installation("TRO", 7,19, 1));
        installations.add(new Installation("TRB", 7,19, 2));
        installations.add(new Installation("TRC", 7,19, 3));
        installations.add(new Installation("CPR", 0,23, 4));
        installations.add(new Installation("SEN", 0,23, 5));
        installations.add(new Installation("SDO", 0,23, 6));
        installations.add(new Installation("SEQ", 0,23, 7));
        installations.add(new Installation("OSE", 0,23, 8));
        installations.add(new Installation("OSB", 0,23, 9));
        installations.add(new Installation("OSC", 0,23, 10));
        installations.add(new Installation("OSO", 0,23, 12));
        installations.add(new Installation("SSC", 0,23, 13));
        installations.add(new Installation("OSS", 0,23, 14));
        installations.add(new Installation("DSD", 0,23, 15));
        installations.add(new Installation("KVB", 0,23, 16));
        installations.add(new Installation("VMO", 0,23, 17));
        installations.add(new Installation("WEL", 0,23, 18));
        installations.add(new Installation("VFB", 0,23, 19));
        installations.add(new Installation("WEP", 0,23, 20));
        installations.add(new Installation("HUL", 0,23, 21));
        installations.add(new Installation("STA", 7,19, 22));
        installations.add(new Installation("STB", 0,23, 23));
        installations.add(new Installation("STC", 0,23, 24));
        installations.add(new Installation("GFA", 0,23, 25));
        installations.add(new Installation("GFB", 0,23, 26));
        installations.add(new Installation("GFC", 0,23, 27));
        installations.add(new Installation("SOD", 0,23, 28));

        ArrayList<Order> orders = new ArrayList<Order>();

        orders.add(new Order(30,0, installations.get(4),1));
        orders.add(new Order(30,0, installations.get(16),2));
        orders.add(new Order(30,0, installations.get(11),3));
        orders.add(new Order(30,0, installations.get(17),4));
        orders.add(new Order(30,0, installations.get(20),5));
        orders.add(new Order(30,0, installations.get(21),6));
        orders.add(new Order(30,0, installations.get(26),7));
        orders.add(new Order(40,0, installations.get(9), 8));
        orders.add(new Order(30,0, installations.get(1), 9));

        HashMap<Installation, HashMap<Installation, Double>> distances = new HashMap<Installation, HashMap<Installation, Double>>();

        ArrayList<ArrayList<Double>> distanceMatrix = new ArrayList<ArrayList<Double>>();
        distanceMatrix.add(0, new ArrayList<Double>(Arrays.asList(0.00, 43.47, 47.25, 43.76, 44.47, 44.65, 43.21, 41.73, 71.65, 71.65, 70.54, 64.49, 64.49, 75.08, 89.05, 76.84, 81.18, 81.18, 64.88, 71.58, 71.58, 97.52, 97.69, 96.75, 87.12, 86.97, 85.01, 37.67)));
        distanceMatrix.add(1, new ArrayList<Double>(Arrays.asList(43.47, 0.00, 10.14, 14.85, 6.83, 19.08, 12.96, 15.09, 28.27, 28.27, 28.11, 23.53, 23.53, 31.64, 46.71, 44.07, 47.01, 47.01, 25.81, 33.86, 33.86, 65.74, 64.90, 65.90, 55.13, 55.70, 54.09, 15.85)));
        distanceMatrix.add(2, new ArrayList<Double>(Arrays.asList(47.25, 10.14, 0.00, 7.23, 3.70, 11.07, 5.96, 8.76, 26.56, 26.56, 23.78, 17.26, 17.26, 31.44, 48.84, 34.32, 37.58, 37.58, 17.91, 25.37, 25.37, 56.04, 55.36, 56.08, 45.38, 45.86, 44.19, 11.98)));
        distanceMatrix.add(3, new ArrayList<Double>(Arrays.asList(43.76, 14.85, 7.23, 0.00, 8.42, 4.25, 1.90, 2.13, 33.24, 33.24, 29.62, 22.43, 22.43, 38.37, 56.00, 34.04, 37.99, 37.99, 21.65, 27.85, 27.85, 55.52, 55.24, 55.20, 44.86, 45.04, 43.21, 6.26)));
        distanceMatrix.add(4, new ArrayList<Double>(Arrays.asList(44.47, 6.83, 3.70, 8.42, 0.00, 12.65, 6.62, 9.18, 28.05, 28.05, 26.12, 20.13, 20.13, 32.47, 49.18, 38.02, 41.27, 41.27, 21.27, 28.91, 28.91, 59.74, 59.06, 59.76, 49.07, 49.54, 47.86, 11.26)));
        distanceMatrix.add(5, new ArrayList<Double>(Arrays.asList(44.65, 19.08, 11.07, 4.25, 12.65, 0.00, 6.12, 4.46, 36.01, 36.01, 31.73, 24.25, 24.25, 41.40, 59.35, 32.27, 36.53, 36.53, 22.62, 27.84, 27.84, 53.38, 53.31, 52.86, 42.80, 42.83, 40.93, 7.36)));
        distanceMatrix.add(6, new ArrayList<Double>(Arrays.asList(43.21, 12.96, 5.96, 1.90, 6.62, 6.12, 0.00, 2.81, 32.39, 32.39, 29.15, 22.16, 22.16, 37.37, 54.78, 35.24, 39.06, 39.06, 21.80, 28.39, 28.39, 56.82, 56.46, 56.56, 46.14, 46.38, 44.57, 6.32)));
        distanceMatrix.add(7, new ArrayList<Double>(Arrays.asList(41.73, 15.09, 8.76, 2.13, 9.18, 4.46, 2.81, 0.00, 35.11, 35.11, 31.66, 24.52, 24.52, 40.15, 57.59, 35.78, 39.82, 39.82, 23.78, 29.93, 29.93, 57.16, 56.95, 56.76, 46.52, 46.65, 44.79, 4.14)));
        distanceMatrix.add(8, new ArrayList<Double>(Arrays.asList(71.65, 28.27, 26.56, 33.24, 28.05, 36.01, 32.39, 35.11, 0.00, 0.00, 7.36, 13.61, 13.61, 6.07, 24.69, 36.68, 36.48, 36.48, 18.15, 22.79, 22.79, 54.28, 52.24, 55.61, 45.50, 46.92, 46.17, 38.54)));
        distanceMatrix.add(9, new ArrayList<Double>(Arrays.asList(71.65, 28.27, 26.56, 33.24, 28.05, 36.01, 32.39, 35.11, 0.00, 0.00, 7.36, 13.61, 13.61, 6.07, 24.69, 36.68, 36.48, 36.48, 18.15, 22.79, 22.79, 54.28, 52.24, 55.61, 45.50, 46.92, 46.17, 38.54)));
        distanceMatrix.add(10, new ArrayList<Double>(Arrays.asList(70.54, 28.11, 23.78, 29.62, 26.12, 31.73, 29.15, 31.66, 7.36, 7.36, 0.00, 7.64, 7.64, 13.23, 31.52, 29.33, 29.29, 29.29, 11.38, 15.43, 15.43, 47.41, 45.51, 48.62, 38.35, 39.72, 38.90, 35.46)));
        distanceMatrix.add(11, new ArrayList<Double>(Arrays.asList(64.49, 23.53, 17.26, 22.43, 20.13, 24.25, 22.16, 24.52, 13.61, 13.61, 7.64, 0.00, 0.00, 19.67, 38.31, 25.54, 26.75, 26.75, 4.95, 12.20, 12.20, 45.65, 44.18, 46.45, 35.74, 36.83, 35.69, 28.46)));
        distanceMatrix.add(12, new ArrayList<Double>(Arrays.asList(64.49, 23.53, 17.26, 22.43, 20.13, 24.25, 22.16, 24.52, 13.61, 13.61, 7.64, 0.00, 0.00, 19.67, 38.31, 25.54, 26.75, 26.75, 4.95, 12.20, 12.20, 45.65, 44.18, 46.45, 35.74, 36.83, 35.69, 28.46)));
        distanceMatrix.add(13, new ArrayList<Double>(Arrays.asList(75.08, 31.64, 31.44, 38.37, 32.47, 41.40, 37.37, 40.15, 6.07, 6.07, 13.23, 19.67, 19.67, 0.00, 18.65, 42.33, 41.81, 41.81, 24.22, 28.54, 28.54, 59.09, 56.90, 60.58, 50.73, 52.22, 51.58, 43.35)));
        distanceMatrix.add(14, new ArrayList<Double>(Arrays.asList(89.05, 46.71, 48.84, 56.00, 49.18, 59.35, 54.78, 57.59, 24.69, 24.69, 31.52, 38.31, 38.31, 18.65, 0.00, 59.62, 58.32, 58.32, 42.77, 46.27, 46.27, 73.96, 71.41, 75.82, 66.83, 68.48, 68.13, 60.40)));
        distanceMatrix.add(15, new ArrayList<Double>(Arrays.asList(76.84, 44.07, 34.32, 34.04, 38.02, 32.27, 35.24, 35.78, 36.68, 36.68, 29.33, 25.54, 25.54, 42.33, 59.62, 0.00, 4.99, 4.99, 20.82, 13.93, 13.93, 21.73, 21.22, 21.84, 11.06, 11.69, 10.27, 39.52)));
        distanceMatrix.add(16, new ArrayList<Double>(Arrays.asList(81.18, 47.01, 37.58, 37.99, 41.27, 36.53, 39.06, 39.82, 36.48, 36.48, 29.29, 26.75, 26.75, 41.81, 58.32, 4.99, 0.00, 0.00, 22.41, 14.57, 14.57, 19.00, 17.89, 19.71, 9.09, 10.44, 9.83, 43.68)));
        distanceMatrix.add(17, new ArrayList<Double>(Arrays.asList(81.18, 47.01, 37.58, 37.99, 41.27, 36.53, 39.06, 39.82, 36.48, 36.48, 29.29, 26.75, 26.75, 41.81, 58.32, 4.99, 0.00, 0.00, 22.41, 14.57, 14.57, 19.00, 17.89, 19.71, 9.09, 10.44, 9.83, 43.68)));
        distanceMatrix.add(18, new ArrayList<Double>(Arrays.asList(64.88, 25.81, 17.91, 21.65, 21.27, 22.62, 21.80, 23.78, 18.15, 18.15, 11.38, 4.95, 4.95, 24.22, 42.77, 20.82, 22.41, 22.41, 0.00, 8.19, 8.19, 41.41, 40.12, 42.06, 31.27, 32.27, 31.04, 27.89)));
        distanceMatrix.add(19, new ArrayList<Double>(Arrays.asList(71.58, 33.86, 25.37, 27.85, 28.91, 27.84, 28.39, 29.93, 22.79, 22.79, 15.43, 12.20, 12.20, 28.54, 46.27, 13.93, 14.57, 14.57, 8.19, 0.00, 0.00, 33.45, 32.03, 34.27, 23.60, 24.78, 23.76, 34.06)));
        distanceMatrix.add(20, new ArrayList<Double>(Arrays.asList(71.58, 33.86, 25.37, 27.85, 28.91, 27.84, 28.39, 29.93, 22.79, 22.79, 15.43, 12.20, 12.20, 28.54, 46.27, 13.93, 14.57, 14.57, 8.19, 0.00, 0.00, 33.45, 32.03, 34.27, 23.60, 24.78, 23.76, 34.06)));
        distanceMatrix.add(21, new ArrayList<Double>(Arrays.asList(97.52, 65.74, 56.04, 55.52, 59.74, 53.38, 56.82, 57.16, 54.28, 54.28, 47.41, 45.65, 45.65, 59.09, 73.96, 21.73, 19.00, 19.00, 41.41, 33.45, 33.45, 0.00, 3.13, 2.80, 10.69, 10.56, 12.52, 60.72)));
        distanceMatrix.add(22, new ArrayList<Double>(Arrays.asList(97.69, 64.90, 55.36, 55.24, 59.06, 53.31, 56.46, 56.95, 52.24, 52.24, 45.51, 44.18, 44.18, 56.90, 71.41, 21.22, 17.89, 17.89, 40.12, 32.03, 32.03, 3.13, 0.00, 5.88, 10.58, 11.00, 13.03, 60.62)));
        distanceMatrix.add(23, new ArrayList<Double>(Arrays.asList(96.75, 65.90, 56.08, 55.20, 59.76, 52.86, 56.56, 56.76, 55.61, 55.61, 48.62, 46.45, 46.45, 60.58, 75.82, 21.84, 19.71, 19.71, 42.06, 34.27, 34.27, 2.80, 5.88, 0.00, 10.84, 10.22, 11.99, 60.22)));
        distanceMatrix.add(24, new ArrayList<Double>(Arrays.asList(87.12, 55.13, 45.38, 44.86, 49.07, 42.80, 46.14, 46.52, 45.50, 45.50, 38.35, 35.74, 35.74, 50.73, 66.83, 11.06, 9.09, 9.09, 31.27, 23.60, 23.60, 10.69, 10.58, 10.84, 0.00, 1.89, 3.17, 50.13)));
        distanceMatrix.add(25, new ArrayList<Double>(Arrays.asList(86.97, 55.70, 45.86, 45.04, 49.54, 42.83, 46.38, 46.65, 46.92, 46.92, 39.72, 36.83, 36.83, 52.22, 68.48, 11.69, 10.44, 10.44, 32.27, 24.78, 24.78, 10.56, 11.00, 10.22, 1.89, 0.00, 2.03, 50.18)));
        distanceMatrix.add(26, new ArrayList<Double>(Arrays.asList(85.01, 54.09, 44.19, 43.21, 47.86, 40.93, 44.57, 44.79, 46.17, 46.17, 38.90, 35.69, 35.69, 51.58, 68.13, 10.27, 9.83, 9.83, 31.04, 23.76, 23.76, 12.52, 13.03, 11.99, 3.17, 2.03, 0.00, 48.28)));
        distanceMatrix.add(27, new ArrayList<Double>(Arrays.asList(37.67, 15.85, 11.98, 6.26, 11.26, 7.36, 6.32, 4.14, 38.54, 38.54, 35.46, 28.46, 28.46, 43.35, 60.40, 39.52, 43.68, 43.68, 27.89, 34.06, 34.06, 60.72, 60.62, 60.22, 50.13, 50.18, 48.28, 0.00)));

        for (int i=0; i<28; i++) {
            distances.put(installations.get(i), new HashMap<Installation, Double>());
            for (int j=0; j<28; j++) {
                distances.get(installations.get(i)).put(installations.get(j), distanceMatrix.get(i).get(j));
            }
        }

        System.out.println(distances.get(installations.get(0)).get(installations.get(1)));

        ArrayList<Integer> weatherStatesByHour = new ArrayList<Integer>(Arrays.asList(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2));

        HashMap<Integer, Double> weatherImpactByState = new HashMap<Integer, Double>();
        weatherImpactByState.put(0, 1.0);
        weatherImpactByState.put(1, 1.2);
        weatherImpactByState.put(2, 1.3);
        weatherImpactByState.put(3, 2.0);

        return new ProblemData(problemInstanceParameters, heuristicInstanceParameters, installations, vessels,distances, orders, weatherStatesByHour, weatherImpactByState);


    }
}
