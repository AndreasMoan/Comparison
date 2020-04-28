package HGSADCwSO.implementations;

import java.util.ArrayList;
import java.util.Collections;

public class Move {

    private int numberOfMoves = 7;

    public ArrayList<Integer> doMove(Integer u, Integer v, ArrayList<Integer> orders, int moveNumber){
        ArrayList<Integer> newOrders;
        Integer x = getSuccessor(u, orders);
        Integer y = getSuccessor(v, orders);

        switch (moveNumber){
            case 1: newOrders = move1(u, v, orders);
                break;
            case 2: newOrders = move2(u,v, x, orders);
                break;
            case 3: newOrders = move3(u, v, x, orders);
                break;
            case 4: newOrders = move4(u,v, orders);
                break;
            case 5: newOrders = move5(u, v, x, orders);
                break;
            case 6: newOrders = move6(u,v, x, y, orders);
                break;
            case 7: newOrders = move7(v, x, orders);
                break;
            default: newOrders = orders;
            System.out.println("Illegal move number");
                break;
        }
        return newOrders;
    }

    public int getNumberOfMoves() {return numberOfMoves;}

    //returns the successor of a in orders, null if a is the last element
    private Integer getSuccessor(Integer a, ArrayList<Integer> orders) {
        int indexA = orders.indexOf(a);
        if(indexA==orders.size()-1){ //at end of list. Hence, no successor
            return null;
        }
        else {
            return orders.get(indexA+1);
        }
    }


    //x is the successor of u, y is the successor of v

    //remove u and place it after v
    private ArrayList<Integer> move1 (Integer u, Integer v, ArrayList<Integer> orders){
        ArrayList<Integer> newOrders = new ArrayList<>(orders);
        newOrders.remove(u);
        int indexV = newOrders.indexOf(v);
        newOrders.add(indexV+1, u);
        return newOrders;
    }

    //remove u and x and place u and x after v
    private ArrayList<Integer> move2(Integer u, Integer v, Integer x, ArrayList<Integer> orders){
        ArrayList<Integer> newOrders = new ArrayList<> (orders);
        if (x!=null){ //if x is null, then u is the last element of the list of orders and no changes are performed
            newOrders.remove(u);
            newOrders.remove(x);
            int indexV = newOrders.indexOf(v);
            newOrders.add(indexV+1, u);
            newOrders.add(indexV+2, x);
        }
        return newOrders;
    }

    //remove u and x and place x and u after v (opposite sequence of move 2)
    private ArrayList<Integer> move3(Integer u, Integer v, Integer x, ArrayList<Integer> orders){
        ArrayList<Integer> newOrders = new ArrayList<Integer> (orders);
        if (x!=null){ //if x is null, then u is the last element of the list of orders and no changes are performed
            newOrders.remove(u);
            newOrders.remove(x);
            int indexV = newOrders.indexOf(v);
            newOrders.add(indexV+1, x);
            newOrders.add(indexV+2, u);
        }
        return newOrders;
    }

    //swap the position of u and v
    private ArrayList<Integer> move4(Integer u, Integer v, ArrayList<Integer> orders){
        ArrayList<Integer> newOrders = new ArrayList<Integer>(orders);
        int indexU = newOrders.indexOf(u);
        int indexV = newOrders.indexOf(v);
        Collections.swap(newOrders, indexU, indexV);
        return newOrders;
    }

    //swap the position of u and x with v
    private ArrayList<Integer> move5 (Integer u, Integer v, Integer x, ArrayList<Integer> orders){
        ArrayList<Integer> newOrders = new ArrayList<Integer>(orders);
        if (x!=null){ //if x is null, then u is the last element of the list of orders and no changes are performed
            //swap place of u and v
            int indexU = newOrders.indexOf(u);
            int indexV = newOrders.indexOf(v);
            Collections.swap(newOrders, indexU, indexV);
            //remove x and insert it after u (previous index of v)
            newOrders.remove(x);
            indexU = newOrders.indexOf(u);
            newOrders.add(indexU+1, x);
        }
        return newOrders;
    }

    //swap u and x with v and y
    private ArrayList<Integer> move6 (Integer u, Integer v, Integer x, Integer y, ArrayList<Integer> orders){
        ArrayList<Integer> newOrders = new ArrayList<Integer>(orders);
        if (x!= null && y!= null) {
            int indexU = newOrders.indexOf(u);
            int indexV = newOrders.indexOf(v);
            int indexX = newOrders.indexOf(x);
            int indexY = newOrders.indexOf(y);
            Collections.swap(newOrders, indexU, indexV);
            Collections.swap(newOrders, indexX, indexY);
        }
        return newOrders;
    }

    //swap the position of x and v
    private ArrayList<Integer> move7(Integer v, Integer x, ArrayList<Integer> orders){
        ArrayList<Integer> newOrders = new ArrayList<Integer>(orders);
        if(x!=null){
            int indexX = newOrders.indexOf(x);
            int indexV = newOrders.indexOf(v);
            Collections.swap(newOrders, indexX, indexV);
        }
        return newOrders;
    }

}
