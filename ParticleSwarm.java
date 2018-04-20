/**
 * Created by hyma on 2/4/18.
 */
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;

public class ParticleSwarm {
    private double[] fitnessValues;
    private double[] solutionVector;
    private double[] bestSolutionVector;
    private double fitness;
    private double bestFitness;
    private double[] velocity;

    public int[] neighbours;

    public ParticleSwarm(int numOfGames, int numOfFeatures) {
        fitnessValues = new double[numOfGames];
        solutionVector = new double[numOfFeatures];
        bestSolutionVector = new double[numOfFeatures];
        bestFitness = Double.NEGATIVE_INFINITY;
        velocity = new double[numOfFeatures];
    }

    public void setSolutionVector(double[] sv) {

        solutionVector = sv.clone();
    }

    public double[] getSolutionVector() {

        return solutionVector;
    }

    public double[] getBestSolutionVector() {

        return bestSolutionVector;
    }

    public double[] getVelocity() {

        return velocity;
    }

    public void setVelocity(double[] vel) {

        velocity = vel.clone();
    }

    public void addFitnessValue(double fitness, int index) {

        fitnessValues[index] = fitness;
    }

    public void calculateAverageFitness() {
        Arrays.sort(fitnessValues);
        if (fitnessValues.length % 2 == 0) {
            fitness = (fitnessValues[fitnessValues.length / 2 - 1] +
                    fitnessValues[fitnessValues.length / 2]) / 2;
        } else {
            fitness = fitnessValues[fitnessValues.length / 2 - 1];
        }

        if (fitness > bestFitness) {
            bestFitness = fitness;
            bestSolutionVector = solutionVector.clone();
        }
    }

    public double getFitness() {

        return fitness;
    }

    public double getBestFitness() {

        return bestFitness;
    }

    private static final int NUMOFPARTICLES = 25;
    private static final int TOPOLOGYWIDTH = (int)Math.sqrt(NUMOFPARTICLES);
    private static final int NUMOFGAMES = 7;
    private static final int NUMOFITER = 30;

    private static final int NUMOFFEATURES = Heuristics.NUM_FEATURES;
    public static final double INERTIA = 0.72;
    public static final double ACCELERATION1 = 1.42;
    public static final double ACCELERATION2 = 1.42;
    public static final double VMAX = 0.5;

    public ParticleSwarm[] particles;

    public void Optimization() {
        final ExecutorService pool;
        particles = new ParticleSwarm[NUMOFPARTICLES];
            Random rand = new Random();
            for (int i = 0; i < NUMOFPARTICLES; i++) {
                particles[i] = new ParticleSwarm(NUMOFGAMES, NUMOFFEATURES);
                particles[i].neighbours = findNeighbours(i).clone();
                double[] sv = new double[NUMOFFEATURES];
                for (int j = 0; j < NUMOFFEATURES; j++) {
                    sv[j] = rand.nextDouble()*2 - 1;
                }
                particles[i].setSolutionVector(sv);
            }
    }

    public void runPSO() {
        System.out.println("Starting");
        long start = System.currentTimeMillis();
        for (int i = 0; i < NUMOFITER; i++) {

            try {
                for (int particle = 0; particle < NUMOFPARTICLES; particle++) {
                    for (int game = 0; game < NUMOFGAMES; game++) {
                        int index = particle*NUMOFGAMES + game;
                    }
                    particles[particle].calculateAverageFitness();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
            updateParticles();

        }

        double bestParticleFitnessScore = particles[0].getBestFitness();
        for (int i = 1; i < particles.length; i++) {
            if (particles[i].getBestFitness() > bestParticleFitnessScore) {
                bestParticleFitnessScore = particles[i].getBestFitness();
            }
        }
    }

    class Handler implements Callable<Double> {
        private ParticleSwarm particleSwarm;
        public int[] utilityValues;

        /* Handler(ParticleSwarm particleSwarm) {
            this.particleSwarm = particleSwarm;
            this.utilityValues = new int[9];
        } */

        public Double call() {
            Heuristics heuristics = Heuristics.getInstance();
            State s = new State();
            PlayerSkeleton player = new PlayerSkeleton();
            heuristics.setWeights(particleSwarm.getSolutionVector());
            // Reset all the values for the
            while(!s.hasLost()) {
                s.makeMove(player.pickMove(s,s.legalMoves()));
            }
            return particleSwarm.calculateFitnessForGame(utilityValues, s.getRowsCleared());
        }
    }

    public static double calculateFitnessForGame(int[] uv, int rowsCleared) {
        double fitness = rowsCleared;
        if (uv[1] != 0) {
            fitness += ((double) uv[1] - ((double)uv[0] / (double) uv[8])) / uv[1] * 500;

        }
        if (uv[3] != 0) {
            fitness += ((double) uv[3] - ((double)uv[2] / (double)uv[8])) / uv[3] * 500;

        }
        if (uv[5] != 0) {
            fitness += ((double) uv[5] - ((double)uv[4] / (double) uv[8])) / uv[5] * 500;

        }
        if (uv[7] != 0) {
            fitness += ((double) uv[7] - ((double)uv[6] / (double)uv[8])) / uv[7] * 500;

        }
        return rowsCleared;
    }

    public void updateParticles() {
        Random rand = new Random();
        double r1 = rand.nextDouble();
        double r2 = rand.nextDouble();

        for (int i = 0; i < NUMOFPARTICLES; i++) {
            double[] bestSolution = particles[i].getBestSolutionVector();
            double[] currentSolution = particles[i].getSolutionVector();

            double bestNeighbouringFitness = particles[particles[i].neighbours[0]].getBestFitness();
            int bestNeighbourIndex = particles[i].neighbours[0];
            for (int j = 1; j < particles[i].neighbours.length; j++) {
                if (particles[particles[i].neighbours[j]].getBestFitness() > bestNeighbouringFitness) {
                    bestNeighbouringFitness = particles[particles[i].neighbours[j]].getBestFitness();
                    bestNeighbourIndex = particles[i].neighbours[j];
                }
            }
            double[] bestNeighbouringSolutionVector = particles[bestNeighbourIndex].getBestSolutionVector();
            double[] newVelocity = particles[i].getVelocity().clone();
            for (int j = 0; j < NUMOFFEATURES; j++) {
                newVelocity[j] = newVelocity[j]*INERTIA +
                        ACCELERATION1*r1*(bestSolution[j] - currentSolution[j]) +
                        ACCELERATION2*r2*(bestNeighbouringSolutionVector[j] - currentSolution[j]);
            }

            double[] newSolution = particles[i].getSolutionVector().clone();
            for (int j = 0; j < NUMOFFEATURES; j++) {
                if (newVelocity[j] > VMAX) {
                    newSolution[j] += VMAX;
                } else if (newVelocity[j] < -VMAX) {
                    newSolution[j] -= VMAX;
                } else {
                    newSolution[j] += newVelocity[j];
                }
            }

            particles[i].setSolutionVector(newSolution);
            particles[i].setVelocity(newVelocity);
        }
    }

    public int[] findNeighbours(int index) {
        int[] neighbours;
        if (NUMOFPARTICLES == 1) {
            neighbours = new int[0];
            return neighbours;
        }

        if (index == 0) {
            neighbours = new int[2];
            neighbours[0] = 1;
            neighbours[1] = TOPOLOGYWIDTH;
        } else if (index == TOPOLOGYWIDTH-1) {
            neighbours = new int[2];
            neighbours[0] = TOPOLOGYWIDTH - 2;
            neighbours[1] = index + TOPOLOGYWIDTH;
        } else if (index == NUMOFPARTICLES - TOPOLOGYWIDTH) {
            neighbours = new int[2];
            neighbours[0] = NUMOFPARTICLES - TOPOLOGYWIDTH + 1;
            neighbours[1] = index - TOPOLOGYWIDTH;
        } else if (index == NUMOFPARTICLES - 1) {
            neighbours = new int[2];
            neighbours[0] = NUMOFPARTICLES - 1;
            neighbours[1] = index - TOPOLOGYWIDTH;
        } else if (index < TOPOLOGYWIDTH) {
            neighbours = new int[3];
            neighbours[0] = index - 1;
            neighbours[1] = index + 1;
            neighbours[2] = index + TOPOLOGYWIDTH;
        } else if (index > NUMOFPARTICLES - TOPOLOGYWIDTH) {
            neighbours = new int[3];
            neighbours[0] = index - 1;
            neighbours[1] = index + 1;
            neighbours[2] = index - TOPOLOGYWIDTH;
        } else if (index % TOPOLOGYWIDTH == 0) {
            neighbours = new int[3];
            neighbours[0] = index - TOPOLOGYWIDTH;
            neighbours[1] = index + TOPOLOGYWIDTH;
            neighbours[2] = index + 1;
        } else if ((index + 1) % TOPOLOGYWIDTH == 0) {
            neighbours = new int[3];
            neighbours[0] = index - TOPOLOGYWIDTH;
            neighbours[1] = index + TOPOLOGYWIDTH;
            neighbours[2] = index - 1;
        } else {
            neighbours = new int[4];
            neighbours[0] = index - TOPOLOGYWIDTH;
            neighbours[1] = index + TOPOLOGYWIDTH;
            neighbours[2] = index - 1;
            neighbours[3] = index + 1;
        }
        return neighbours;
    }
}