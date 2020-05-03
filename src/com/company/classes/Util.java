package com.company.classes;

public class Util {

    public static double rms(double[] d) {
        double sum = 0;
        for (int i=0; i < d.length; i++)
            sum += d[i] * d[i];
        return Math.sqrt(sum / d.length);
    }

    public static double min(double[] d) {
        double min = d[0];
        for (int i=1; i < d.length; i++)
            if (d[i] < min)
                min = d[i];
        return min;
    }

    public static double max(double[] d) {
        double max = d[0];
        for (int i=1; i < d.length; i++)
            if (d[i] > max)
                max = d[i];
        return max;
    }

    public static double threshold(double value, double min, double max) {
        if (value < min)
            return min;
        if (value > max)
            return max;
        return value;
    }

} // class Util