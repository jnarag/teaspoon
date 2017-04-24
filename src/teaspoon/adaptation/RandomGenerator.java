package teaspoon.adaptation;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import cern.jet.random.Poisson;
import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

/**
 * Created with IntelliJ IDEA.
 * User: jayna
 * Date: 15/08/2013
 * Time: 11:40
 * To change this template use File | Settings | File Templates.
 */
public class RandomGenerator {

    public static java.util.Random random;
    public static RandomEngine engine;
    public static Uniform uniform;


    public RandomGenerator() {

        random = new Random();
        engine = new MersenneTwister((int)System.currentTimeMillis( ) % Integer.MAX_VALUE);
        uniform = new Uniform(engine);

    }

    public int nextInt() {

        return engine.nextInt();
    }

    public int nextInt(int upper) {

        return uniform.nextIntFromTo(0, upper);
    }

    public double nextDouble() {

        return uniform.nextDouble();
    }

    public int nextInt(int lower, int upper) {

        return uniform.nextIntFromTo(lower, upper);
    }

    public long nextLong() {

        return engine.nextLong();
    }

    public int nextPoissonInt(double expectedValue) {

        Poisson poisson = new Poisson(expectedValue, engine);
        return poisson.nextInt();

    }

    public void shuffle(Object[] a) {

        List<Object> l = Arrays.asList(a);
        Collections.shuffle(l);
        a = l.toArray();

    }


}
