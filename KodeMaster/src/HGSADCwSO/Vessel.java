package HGSADCwSO;

public class Vessel {

    private String name;
    private int capacity, minSpeed, maxSpeed, timeCharterCost, number, returnDay;

    public Vessel(String name, int capacity, int minSpeed, int maxSpeed, int timeCharterCost, int number, int returnDay) {
        this.name = name;
        this.capacity = capacity;
        this.minSpeed = minSpeed;
        this.maxSpeed = maxSpeed;
        this.timeCharterCost = timeCharterCost;
        this.returnDay = returnDay;
        this.number = number;
    }

    public String getName() {
        return name;
    }

    public int getCapacity() {
        return capacity;
    }

    public int getMinSpeed() {
        return minSpeed;
    }

    public int getMaxSpeed() {
        return maxSpeed;
    }

    public int getTimeCharterCost() {
        return timeCharterCost;
    }

    public int getReturnDay() {
        return returnDay;
    }

    public int getNumber() {
        return number;
    }

}
