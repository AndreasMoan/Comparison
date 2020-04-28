package HGSADCwSO;

public class Order {

    private int demand;
    private int day;
    private int number;
    private Installation installation;
    private int deadline;

    public Order(int demand, int day, Installation installation, int number, int deadline) {
        this.demand = demand;
        this.day = day;
        this.installation = installation;
        this.number = number;
        this.deadline = deadline;
    }

    public Installation getInstallation() {
        return installation;
    }

    public int getDay() {
        return day;
    }

    public int getDemand() {
        return demand;
    }

    public int getNumber() {
        return number;
    }

    public int getDeadline() { return deadline; }
}
