package HGSADCwSO;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.util.*;

public class Utilities {

    public static <T> T pickRandomElementFromList(List<T> list) {
        int randomIndex = new Random().nextInt(list.size());
        return list.get(randomIndex);
    }

    public static <T> T pickAndRemoveRandomElementFromList(ArrayList<T> list) {
        T element = pickRandomElementFromList(list);
        list.remove(element);
        return element;
    }

    public static <T> T pickRandomElementFromSet(Set<T> set) {
        int randomIndex = new Random().nextInt(set.size());
        int counter = 0;
        for (T t : set) {
            if (counter == randomIndex){
                return t;
            }
            counter++;
        }
        return null;
    }

    public static <T> T pickAndRemoveRandomElementFromSet(Set<T> set) {
        int randomIndex = new Random().nextInt(set.size());
        int counter = 0;
        for (T t : set) {
            if (counter == randomIndex){
                set.remove(t);
                return t;
            }
            counter++;
        }
        return null;
    }

    public static double parseDouble(String commaSeparatedDouble) {
        DecimalFormat df = new DecimalFormat();
        DecimalFormatSymbols symbols = new DecimalFormatSymbols();
        symbols.setDecimalSeparator(',');
        df.setDecimalFormatSymbols(symbols);
        try {
            return df.parse(commaSeparatedDouble).doubleValue();
        } catch (ParseException e) {
            e.printStackTrace();
            System.out.println("Something went wrong when parsing doubles");
            return -1.0;
        }
    }

    public static <T> ArrayList<T> getAllElements(ArrayList<T> list1, ArrayList<T> list2) {
        ArrayList<T> allElements = new ArrayList<T>();
        allElements.addAll(list1);
        allElements.addAll(list2);
        return allElements;
    }

    public static <T> Set<Set<T>> cartesianProduct(Set<T> set) { //returns cartesian product of itself, i.e. set x set, excluding those pairs where both elements in a pair are equal
        HashSet<Set<T>> cartProduct = new HashSet<>();
        for (T element : set) {
            for (T element2 : set){
                if (element != element2){
                    HashSet<T> pair = new HashSet<T>();
                    pair.add(element);
                    pair.add(element2);
                    cartProduct.add(pair);
                }
            }
        }
        return cartProduct;
    }

    public static HashMap<Integer, ArrayList<Integer>> deepCopyVesselTour(HashMap<Integer, ArrayList<Integer>> vesselTour) {

        HashMap<Integer, ArrayList<Integer>> vesselTourCopy = new HashMap<>();

        for (Integer vessel : vesselTour.keySet()){
            ArrayList<Integer> visitSequence= new ArrayList<Integer>(vesselTour.get(vessel));
            vesselTourCopy.put(vessel, visitSequence);
        }
        return vesselTourCopy;
    }

    public static Comparator<Individual> getFitnessComparator() {
        return new Comparator<Individual>() {
            @Override
            public int compare(Individual individual1, Individual individual2) {
                if (individual1.getFitness() < individual2.getFitness()){
                    return  1;
                }
                else if (individual1.getFitness() > individual2.getFitness()) {
                    return  -1;
                }
                else {
                    return 0;
                }
            }
        };
    }

    public static Comparator<Order> getDeadlineComparator() { //testet! Den gir tidligste deadline først
        return new Comparator<Order>() {
            @Override
            public int compare(Order order1, Order order2) {
                if (order1.getDeadline() > order2.getDeadline()){
                    return  1;
                }
                else if (order1.getDeadline() < order2.getDeadline()) {
                    return  -1;
                }
                else {
                    return 0;
                }
            }
        };
    }


    public static <K> Comparator <Map.Entry<K, Double>> getMapEntryWithDoubleComparator() {
        return new Comparator<Map.Entry<K, Double>>() {
            public int compare(Map.Entry<K, Double> o1, Map.Entry<K, Double> o2) {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        };
    }

    //sorts so the elements with the lowest biased fitness are first in the list //TODO test
    public static Comparator<Individual> getBiasedFitnessComparator() {
        return new Comparator<Individual>() {
            public int compare(Individual ind1, Individual ind2) {
                if (ind1.getBiasedFitness() < ind2.getBiasedFitness()) {
                    return -1;
                }
                else if (ind1.getBiasedFitness() > ind2.getBiasedFitness()) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        };
    }

    //sorts so the elements with the penalized costs are first in the list //TODO test
    public static Comparator<Individual> getPenalizedCostComparator() {
        return new Comparator<Individual>() {
            public int compare(Individual ind1, Individual ind2) {
                if (ind1.getPenalizedCost() < ind2.getPenalizedCost()) {
                    return -1;
                }
                else if (ind1.getPenalizedCost() > ind2.getPenalizedCost()) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        };
    }

    //sorts so the elements with the diversity contribution are first in the list //TODO test
    public static Comparator<Individual> getDiversityContributionComparator() {
        return new Comparator<Individual>() {
            public int compare(Individual ind1, Individual ind2) {
                if (ind1.getDiversityContribution() < ind2.getDiversityContribution()) {
                    return 1;
                }
                else if (ind1.getDiversityContribution() > ind2.getDiversityContribution()) {
                    return -1;
                }
                return 0;
            }
        };
    }

}


