package HGSADCwSO.implementations.DAG;

import java.util.ArrayList;

public class Edge {

    private Node parentNode;
    private Node childNode;

    private double cost;


    /*
    private double startSailingTime;
    private double startIdlingTime;
    private double startServicingTime;
    private double endTime;

    private double speed;

     */

    public Edge(Node parentNode, Node childNode, double cost){
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.cost = cost;
    }

    public double getCost() {
        return cost;
    }

    public Node getChildNode() {
        return childNode;
    }

    public Node getParentNode() {
        return parentNode;
    }
}
