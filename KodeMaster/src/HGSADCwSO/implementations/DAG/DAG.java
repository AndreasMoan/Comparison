package HGSADCwSO.implementations.DAG;

import HGSADCwSO.Order;
import HGSADCwSO.ProblemData;

import java.util.*;

import static java.lang.Math.*;

public class DAG {

    private double shortestPathCost;
    private double shortestFeasiblePathCost;
    private boolean feasibility;

    private ArrayList<Integer> voyage;
    private ArrayList<Integer> voyageWithDepot;
    private ProblemData problemData;

    private HashMap<Integer, HashMap<Integer, Node>> graph;

    private double timePeriodLength;
    private double maxSpeed;
    private double maxSpeedWS2;
    private double maxSpeedWS3;
    private double minSpeed;
    private int vesselReturnTime;
    private double multiplier;
    private double timePerHiv;
    private double speedImpactWS2;
    private double speedImpactWS3;
    private double idlingConsumption;
    private double servicingConsumption;

    private HashMap<Integer, Order> orderByNumber;


    public DAG(ArrayList<Integer> voyage, int vesselReturnTime, ProblemData problemData) {
        this.voyage = voyage;
        // System.out.println(voyage);
        this.problemData = problemData;
        this.graph = new HashMap<Integer, HashMap<Integer, Node>>();

        createVoyageWithDepot();

        this.maxSpeed = this.problemData.getProblemInstanceParameterDouble("Max speed");
        this.minSpeed = this.problemData.getProblemInstanceParameterDouble("Min speed");
        this.speedImpactWS2 = problemData.getProblemInstanceParameterDouble("Impact on sailing from weather state 2");
        this.speedImpactWS3 = problemData.getProblemInstanceParameterDouble("Impact on sailing from weather state 3");
        this.idlingConsumption = problemData.getProblemInstanceParameterDouble("Idling consumption");
        this.servicingConsumption = problemData.getProblemInstanceParameterDouble("Servicing consumption");
        this.maxSpeedWS2 = maxSpeed - speedImpactWS2;
        this.maxSpeedWS3 = maxSpeed - speedImpactWS3;
        this.timePerHiv = problemData.getProblemInstanceParameterDouble("Time per hiv");
        this.vesselReturnTime = vesselReturnTime;
        this.multiplier = this.problemData.getHeuristicParameterDouble("Number of time periods per hour");
        this.timePeriodLength = 1/multiplier;
        

        this.orderByNumber = problemData.getOrdersByNumber();

        if (voyage.size() == 0) {
            shortestPathCost = 0;
            shortestFeasiblePathCost = 0;
        }
        else {
            buildGraph();
            doDijkstra();
            setCost();
            setFeasibility();
        }
    }

    private void createVoyageWithDepot(){
        voyageWithDepot = new ArrayList<Integer>(voyage);
        voyageWithDepot.add(0,0);
        voyageWithDepot.add(0);
    }


    //--------------------------------------- BUILD GRAPH ----------------------------------------



    private void buildGraph(){

        for (int i = 0; i < voyageWithDepot.size(); i++) {
            graph.put(i, new HashMap<Integer, Node>());
        }

        graph.get(0).put(0, new Node(0,0));
        graph.get(0).get(0).setBestCost(0);

        for (int j = 0; j < voyageWithDepot.size() - 1; j++ ) {

            // System.out.println("------------------- 1 --------------------");

            for (int time : graph.get(j).keySet()){

                // System.out.println("-------------------------------- 2 -------------------------------");
                buildEdgesFromNode(graph.get(j).get(time), voyageWithDepot.get(j + 1), j);
            }
        }
    }

    private void buildEdgesFromNode(Node node, int destinationOrderNumber, int legNumber) {


        double distance = problemData.getDistanceByIndex(problemData.getInstallationNumberByOrderNumber(node.getOrderNumber()),problemData.getInstallationNumberByOrderNumber(destinationOrderNumber));

        int nodeStartTime = node.getTime();
        double realStartTime = convertNodeTimeToRealTime(nodeStartTime);

        // System.out.println("From node Iteration  -  order number: " +  node.getOrderNumber());


        int earliestTheoreticalEndTime = nodeStartTime + (int) ceil((distance/maxSpeed + problemData.getDemandByOrderNumber(destinationOrderNumber)*timePerHiv)*multiplier);
        int latestTheoreticalEndTime = nodeStartTime + (int) floor((distance/minSpeed + problemData.getDemandByOrderNumber(destinationOrderNumber)*timePerHiv)*multiplier);

        // System.out.println("ETET : " + earliestTheoreticalEndTime);
        // System.out.println("LTET : " + latestTheoreticalEndTime);

        int finServicingTime = earliestTheoreticalEndTime;

        boolean continue_ = true;
        while (/*continue_ && */ finServicingTime <= latestTheoreticalEndTime ) {

            // System.out.println("Edge creation iteration, ");
            
            finServicingTime = getEarliestFeasibleSercivingFinishingTime(finServicingTime , destinationOrderNumber);

            // System.out.println("Fin servicing time: " + finServicingTime);

            double[] serviceInfo = servicingCalculations(finServicingTime, destinationOrderNumber);

            double servicingCost = serviceInfo[0];
            double finIdlingTime = serviceInfo[1];

            if (!isArrivalPossible(realStartTime, distance, finIdlingTime)) {
                finServicingTime++;
                continue;
            }

            double[] idlingNeededInfo = isIdlingNeeded(realStartTime, distance, finIdlingTime);

            continue_ = (!(idlingNeededInfo[0] >= 0));

            double[] idlingInfo = idlingCalculations(idlingNeededInfo[1], finIdlingTime);

            double idlingCost = idlingInfo[0];
            double finSailingTime = idlingInfo[1];

            double[] timeInAllWeatherStates = getTimeInAllWS(realStartTime, convertNodeTimeToRealTime(finServicingTime));

            double adjustedAverageSpeed = calculateAdjustedAverageSpeed(timeInAllWeatherStates, distance);

            double sailingCost = sailingCalculations(timeInAllWeatherStates, adjustedAverageSpeed);

            Node childNode;

            if (isNodeInGraph(legNumber + 1, finServicingTime)) {
                childNode = graph.get(legNumber + 1).get(finServicingTime);
            }
            else {
                childNode = new Node(finServicingTime, destinationOrderNumber);
                graph.get(legNumber + 1).put(finServicingTime, childNode);
            }

            double edgeCost = sailingCost + idlingCost + servicingCost;

            // System.out.println("Fin servicing time: " + finServicingTime);

            // System.out.println("EDGE COST IS EQUAL TO: " + edgeCost);

            Edge currEdge = new Edge(node, childNode, edgeCost);

            childNode.addParentEdge(currEdge);
            node.addChildEdge(currEdge);

            finServicingTime++;
        }
    }
    
    private double convertNodeTimeToRealTime(int nodeTime) {
        return nodeTime/multiplier;
    }
    
    private int convertRealTimeToNodeTimeFloor(double time) {
        return (int) floor(time/multiplier);
    }

    private boolean isServicePossible(int finServicingNodeTime, int destinationOrder) {

        // System.out.println("Service possible check");

        double realTime = convertNodeTimeToRealTime(finServicingNodeTime);

        double servicingTimeLeft = problemData.getDemandByOrderNumber(destinationOrder) * timePerHiv;

        while (servicingTimeLeft - 1 / problemData.getWeatherImpactByHour((int) realTime) > 0) {
            if (problemData.getWeatherStateByHour().get((int) realTime) == 3) {
                return false;
            }
            if (realTime % 1 > 0) {
                realTime = floor(realTime);
                servicingTimeLeft -= realTime % 1 / problemData.getWeatherImpactByHour((int) floor(realTime));
            } else {
                realTime--;
                servicingTimeLeft -= 1 / problemData.getWeatherImpactByHour((int) floor(realTime));
            }
        }

        return true;
    }

    private int getEarliestFeasibleSercivingFinishingTime(int finServicingTime, int destinationOrderNumber) {
        int earliestFeasibleSercivingFinishingTime = finServicingTime;
        if (!isServicePossible(finServicingTime, destinationOrderNumber)) {
            earliestFeasibleSercivingFinishingTime = getEarliestFeasibleSercivingFinishingTime(finServicingTime +1, destinationOrderNumber);
        }
        return earliestFeasibleSercivingFinishingTime;
    } //TODO fix bug: How do we deal with long periods of bad weather? How long waiting do we accept?

    private double[] servicingCalculations(int time, int order) {
        
        double realTime = convertNodeTimeToRealTime(time);
        double servicingTimeLeft = problemData.getDemandByOrderNumber(order)*timePerHiv;
        
        double consumption = 0;

        while (servicingTimeLeft -1/problemData.getWeatherImpactByHour((int)floor(realTime)) > 0) {
            if (realTime%1 > 0) {
                consumption += servicingConsumption * (realTime%1 * problemData.getWeatherImpactByHour((int) floor(realTime)));
                servicingTimeLeft -= (1 - (realTime%1) / problemData.getWeatherImpactByHour((int) floor(realTime)));
                realTime = floor(realTime);
            } else {
                consumption += servicingConsumption * problemData.getWeatherImpactByHour((int) realTime - 1);
                servicingTimeLeft -= 1 / problemData.getWeatherImpactByHour((int) realTime - 1);
                realTime--;
            }
        }
        if (servicingTimeLeft > 0){
            consumption += servicingConsumption*servicingTimeLeft*problemData.getWeatherImpactByHour((int) realTime - 1);
            realTime -= servicingTimeLeft;
        }

        return new double[] {consumption, realTime};
    }

    private boolean isArrivalPossible(double startTime, double distance, double finIdlingRealTime) {

        double[] tiws = getTimeInAllWS(startTime,finIdlingRealTime);

        double maxDistance = tiws[0]*maxSpeed + tiws[1]*maxSpeedWS2 + tiws[2]*maxSpeedWS3;

        if (maxDistance >= distance) {
            return true;
        }

        return false;
    }


    private double[] isIdlingNeeded(double startTime, double distance, double finIdlingRealTime) {

        double longestSailingTime = distance/minSpeed;

        if (longestSailingTime >= finIdlingRealTime - startTime) {
            return new double[] {-1.0 , 0};
        }
        return new double[] {1, finIdlingRealTime - longestSailingTime - startTime};
    }

    private double[] idlingCalculations(double idlingTime, double finIdlingTime) {
        return new double [] {idlingTime*idlingConsumption , finIdlingTime - idlingTime};
    }

    private double[] getTimeInAllWS(double t1, double t2) {
        double tiws3 = getTimeInWS(t1, t2, 3);
        double tiws2 = getTimeInWS(t1, t2, 2);
        double tiw01 = t2 - t1 - tiws3 - tiws2;
        return new double[] {tiw01, tiws2, tiws3};
    }

    private double getTimeInWS(double t1, double t2, int weatherState){
        double time = 0;
        time += isWeatherState(weatherState, (int) t1) ? 1 - t1%1 : 0;
        time += isWeatherState(weatherState, (int) t2) ? t2%1 : 0;
        for (int i = (int) ceil(t1); i < (int) t2; i++ ){
            time += isWeatherState(weatherState, i) ? 1 : 0;
        }
        return time;
    }

    public boolean isWeatherState(int weatherState, int hour) {
        return problemData.getWeatherStateByHour().get(hour) == weatherState;
    }

    private double calculateAdjustedAverageSpeed(double[] timeInAllWeatherStates, double distance) {
        double durationWS3 = timeInAllWeatherStates[2];
        double durationWS2 = timeInAllWeatherStates[1];
        double durationWS01 = timeInAllWeatherStates[0];
        double duration = durationWS01 + durationWS2 + durationWS3;

        double speed = distance/(duration);

        if (durationWS01 == 0 && speed > (durationWS2*maxSpeedWS2 + durationWS3*maxSpeedWS3)/(durationWS2 + durationWS3)) {
            speed = maxSpeed+1;
        }

        else {
            if (speed > maxSpeedWS3) {
                speed += (durationWS3 * (speed - maxSpeedWS3)) / (durationWS01 + durationWS2);
            }

            if (speed > maxSpeedWS3) {
                speed += (durationWS2 * (speed - maxSpeedWS2)) / (durationWS01);
            }
        }

        return speed;
    }

    private double sailingCalculations(double[] timeInAllWeatherStates, double adjustedAverageSpeed) {
        double cost = 0;
        cost += timeInAllWeatherStates[0]*consumptionFunction(adjustedAverageSpeed);
        cost += (adjustedAverageSpeed < maxSpeedWS2) ? timeInAllWeatherStates[1]*consumptionFunction(adjustedAverageSpeed + speedImpactWS2) : timeInAllWeatherStates[1]*consumptionFunction(maxSpeed);
        cost += (adjustedAverageSpeed < maxSpeedWS3) ? timeInAllWeatherStates[2]*consumptionFunction(adjustedAverageSpeed + speedImpactWS3) : timeInAllWeatherStates[2]*consumptionFunction(maxSpeed);
        return cost;
    }

    private double consumptionFunction(double speed) {
        return (0.8125*speed*speed-13.00*speed+72.75);
    }

    private boolean isNodeInGraph(int columnNumber, int finServicingTime) {
        return graph.get(columnNumber).containsKey(finServicingTime);
    }

    //--------------------------------------- DIJKSTRA ----------------------------------------

    private void doDijkstra(){

        // System.out.println("================================== DIJKSTRA ======================================");

        for (int i = 0; i < graph.size(); i++) {
            // System.out.println("1 ------");
            for (Map.Entry<Integer, Node> pair : graph.get(i).entrySet()){
                // System.out.println("2");
                // System.out.println(pair.getValue().getChildEdges().size());

                expand(pair.getValue());
            }
        }
        for (Map.Entry<Integer, Node> entry : graph.get(graph.size()-1).entrySet()){
            double nodeTime = entry.getValue().getTime();
            double nodeCost = entry.getValue().getBestCost();
            if (nodeTime > vesselReturnTime) {
                entry.getValue().setBestCost(nodeCost + (nodeTime - vesselReturnTime)*problemData.getHeuristicParameterDouble("Cost penalty per excessive time period") );
                entry.getValue().setFeasibility(false);
                }
            else {
                entry.getValue().setFeasibility(true);
            }
        }
    }

    private void expand(Node node) {
        for (Edge childEdge : node.getChildEdges()){
            // System.out.println("3 !!!!");
            if (childEdge.getChildNode().getBestCost() > node.getBestCost() + childEdge.getCost()) {
                childEdge.getChildNode().setBestCost(node.getBestCost() + childEdge.getCost());
                childEdge.getChildNode().setBestParentEdge(childEdge);
                // System.out.println("Parent node cost: " + node.getBestCost() + ", Edge cost: " + childEdge.getCost() + ", Child node cost: " + childEdge.getChildNode().getBestCost());
            }
        }
    }


    //--------------------------------------- SET COST ----------------------------------------



    private void setCost() {

        for (int n = 0; n < graph.size(); n++) {
            // System.out.println("Size of column number " + n + " is : " + graph.get(n).size());
            for ( int k : graph.get(n).keySet()) {
                // System.out.println("Node time: " + graph.get(n).get(k).getTime() + " | Node order number; " + graph.get(n).get(k).getOrderNumber() + " | Node cost: " + graph.get(n).get(k).getBestCost());
            }
        }

        double leastFeasibleCost = Double.POSITIVE_INFINITY;
        double leastInfeasibleCost = Double.POSITIVE_INFINITY;

        if (graph.get(graph.size()-1).containsKey(vesselReturnTime)) {
            leastFeasibleCost = graph.get(graph.size()-1).get(vesselReturnTime).getBestCost();
        }

        for ( Map.Entry<Integer, Node> entry : graph.get(graph.size()-1).entrySet()){

            boolean nodeFeasibility = entry.getValue().getFeasibility();
            double nodeCost = entry.getValue().getBestCost();

            if (nodeFeasibility) {
                if (nodeCost < leastFeasibleCost) {
                    leastFeasibleCost = nodeCost;
                }
            }
            else {
                if (nodeCost < leastInfeasibleCost) {
                    leastInfeasibleCost = nodeCost;
                }
            }
        }

        // System.out.println("LFC : " + leastFeasibleCost);
        // System.out.println("LIC : " + leastInfeasibleCost);

        shortestPathCost = max(leastFeasibleCost, leastInfeasibleCost);
        shortestFeasiblePathCost = leastFeasibleCost;
    }

    //--------------------------------------- SET FEASIBILITY ----------------------------------------

    public void setFeasibility() {
        this.feasibility = shortestFeasiblePathCost < Double.POSITIVE_INFINITY;
    }


    //--------------------------------------- GET  ----------------------------------------

    public double getShortestPathCost(){
        return shortestPathCost;
    }

    public double getShortestFeasiblePathCost() {
        return shortestFeasiblePathCost;
    }

    public boolean getFeasibility() {
        return feasibility;
    }
}
