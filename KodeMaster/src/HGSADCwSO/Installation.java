package HGSADCwSO;

public class Installation {

    private String name;
    private int number;
    private double demand, openingHour, closingHour;

    public Installation(String name, double openingHour, double closingHour, int number) {
        this.name = name;
        this.openingHour = openingHour;
        this.closingHour = closingHour;
        this.number = number;
    }

    public String getName() {
        return name;
    }

    public double getOpeningHour() {
        return openingHour;
    }

    public double getClosingHour() {
        return closingHour;
    }

    public int getNumber() {
        return number;
    }
}
