/**
 * Composite Jacket Thermal Resistance - Java Solution
 * Time Complexity: O(Q * log(MAX_D0)) for Type 4 queries
 * Space Complexity: O(N) where N = 6 layers
 */

import java.io.*;
import java.util.*;

class Layer {
    double d;      // thickness (meters)
    double k;      // base conductivity (W/m·K)
    double mu;     // moisture coefficient
    double c;      // compression coefficient
    double W;      // accumulated moisture
    double C;      // accumulated compression
    
    public Layer(double d, double k, double mu, double c) {
        this.d = d;
        this.k = k;
        this.mu = mu;
        this.c = c;
        this.W = 0.0;
        this.C = 0.0;
    }
    
    // Calculate effective conductivity
    public double kEff(double beta) {
        return k * (1.0 + mu * W) * Math.exp(beta * C);
    }
    
    // Calculate thermal resistance of this layer
    public double thermalResistance(double beta) {
        if (d < 1e-9) return 0.0;
        return d / kEff(beta);
    }
}

class ThermalJacketSystem {
    private static final double EPS = 1e-9;
    private static final double T_S = 37.0;      // Body temperature (°C)
    private static final double T_EXT = -13.0;   // Ambient temperature (°C)
    private static final double Q_MAX = 20.0;    // Maximum heat flux (W/m²)
    private static final double MOISTURE_RATE = 0.01;  // Moisture per hour
    
    private Layer[] layers;
    private double beta;
    
    public ThermalJacketSystem() {
        layers = new Layer[7];
    }
    
    // Calculate total thermal resistance excluding a specific layer
    private double calculateRthWithoutLayer(int excludeLayer) {
        double Rth = 0.0;
        for (int i = 0; i < 7; i++) {
            if (i != excludeLayer) {
                Rth += layers[i].thermalResistance(beta);
            }
        }
        return Rth;
    }
    
    // Calculate total thermal resistance
    private double calculateRth() {
        double Rth = 0.0;
        for (int i = 0; i < 7; i++) {
            Rth += layers[i].thermalResistance(beta);
        }
        return Rth;
    }
    
    // Calculate heat flux for a given d0 (layer 0 thickness)
    private double heatFlux(double d0) {
        // Calculate R_th with layer 0 having thickness d0
        double Rth = calculateRthWithoutLayer(1);  // layer index 1 is layer 0
        
        // Add layer 0's contribution
        if (d0 > EPS) {
            double k0Eff = layers[1].kEff(beta);
            Rth += d0 / k0Eff;
        }
        
        // q = ΔT / R_th
        double deltaT = T_S - T_EXT;  // 50°C
        if (Rth < 1e-12) return 1e18;
        return deltaT / Rth;
    }
    
    // Read initial layer properties and beta
    public void readInput(Scanner sc) {
        // Read 7 layers (index -1 to 5)
        for (int i = 0; i < 7; i++) {
            double d = sc.nextDouble();
            double k = sc.nextDouble();
            double mu = sc.nextDouble();
            double c = sc.nextDouble();
            layers[i] = new Layer(d, k, mu, c);
        }
        
        // Read beta
        beta = sc.nextDouble();
    }
    
    // Type 1: Simulate t hours of snowfall
    private void queryType1EnvironmentalExposure(double t) {
        // Add moisture to layer 5 (index 6)
        layers[6].W += t * MOISTURE_RATE;
        
        // Propagate moisture inward (layer 5 to layer -1)
        for (int i = 6; i > 0; i--) {
            Layer layerI = layers[i];
            Layer layerPrev = layers[i - 1];
            // W_(i-1) += W_i × μ_i
            layerPrev.W += layerI.W * layerI.mu;
        }
    }
    
    // Type 2: Apply compression stress x to layer i
    private void queryType2MechanicalStress(int i, double x) {
        // Convert layer index (-1 to 5) to array index (0 to 6)
        int layerIdx = i + 1;
        
        // Add compression to layer i
        layers[layerIdx].C += x;
        
        // Propagate compression outward (layer i to layer 5)
        for (int idx = layerIdx; idx < 6; idx++) {
            Layer layerCurr = layers[idx];
            Layer layerNext = layers[idx + 1];
            // C_(i+1) += C_i × c_i
            layerNext.C += layerCurr.C * layerCurr.c;
        }
    }
    
    // Type 3: Replace all properties of layer i
    private void queryType3ReplaceLayer(int i, double d, double k, double mu, double c) {
        // Convert layer index (-1 to 5) to array index (0 to 6)
        int layerIdx = i + 1;
        
        // Replace properties
        layers[layerIdx].d = d;
        layers[layerIdx].k = k;
        layers[layerIdx].mu = mu;
        layers[layerIdx].c = c;
        
        // Reset accumulated values
        layers[layerIdx].W = 0.0;
        layers[layerIdx].C = 0.0;
    }
    
    // Type 4: Calculate minimum d0 such that q ≤ 20
    private double queryType4MinimumFoamThickness() {
        // Binary search on d0
        double left = 0.0;
        double right = 1e9;
        
        // Check if even d0 = 0 satisfies the constraint
        if (heatFlux(0.0) <= Q_MAX + EPS) {
            return 0.0;
        }
        
        // Check if even large d0 cannot satisfy
        if (heatFlux(right) > Q_MAX + EPS) {
            return 1e18;
        }
        
        // Binary search
        int iterations = 0;
        final int MAX_ITERATIONS = 100;
        
        while (right - left > 1e-8 && iterations < MAX_ITERATIONS) {
            double mid = (left + right) / 2.0;
            double q = heatFlux(mid);
            
            if (q <= Q_MAX + 1e-10) {
                right = mid;
            } else {
                left = mid;
            }
            
            iterations++;
        }
        
        return right;
    }
    
    // Type 5: Check if any finite d0 can satisfy the constraint
    private String queryType5FeasibilityCheck() {
        // Check with d0 = 0
        double qMin = heatFlux(0.0);
        if (qMin <= Q_MAX + EPS) {
            return "POSSIBLE";
        }
        
        // Check with very large d0
        double qWithLargeD0 = heatFlux(1e9);
        if (qWithLargeD0 <= Q_MAX + EPS) {
            return "POSSIBLE";
        }
        
        // Check the trend
        double q1 = heatFlux(1e6);
        double q2 = heatFlux(1e7);
        
        if (q2 < q1 && q2 <= Q_MAX + EPS) {
            return "POSSIBLE";
        }
        
        // If heat flux is still too high with very large d0
        if (qWithLargeD0 > Q_MAX - EPS) {
            return "IMPOSSIBLE";
        }
        
        return "POSSIBLE";
    }
    
    // Process all queries
    public void processQueries(Scanner sc, PrintWriter out) {
        int Q = sc.nextInt();
        
        for (int q = 0; q < Q; q++) {
            int queryType = sc.nextInt();
            
            if (queryType == 1) {
                // Environmental exposure
                double t = sc.nextDouble();
                queryType1EnvironmentalExposure(t);
            }
            else if (queryType == 2) {
                // Mechanical stress
                int i = sc.nextInt();
                double x = sc.nextDouble();
                queryType2MechanicalStress(i, x);
            }
            else if (queryType == 3) {
                // Replace layer
                int i = sc.nextInt();
                double d = sc.nextDouble();
                double k = sc.nextDouble();
                double mu = sc.nextDouble();
                double c = sc.nextDouble();
                queryType3ReplaceLayer(i, d, k, mu, c);
            }
            else if (queryType == 4) {
                // Minimum foam thickness
                double d0 = queryType4MinimumFoamThickness();
                out.printf("%.10f\n", d0);
            }
            else if (queryType == 5) {
                // Feasibility check
                String result = queryType5FeasibilityCheck();
                out.println(result);
            }
        }
    }
}

public class Solution {
    public static void main(String[] args) throws IOException {
        // Fast I/O
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        Scanner sc = new Scanner(br);
        PrintWriter out = new PrintWriter(new BufferedOutputStream(System.out));
        
        ThermalJacketSystem system = new ThermalJacketSystem();
        system.readInput(sc);
        system.processQueries(sc, out);
        
        out.flush();
        out.close();
        sc.close();
    }
}

/*
COMPILATION:
javac Solution.java

RUN:
java Solution < input.txt

SAMPLE INPUT:
0.02 0.04 0.1 0.2
0.05 0.03 0.15 0.25
0.0 0.035 0.12 0.22
0.0 0.038 0.11 0.21
0.015 0.045 0.13 0.23
0.018 0.042 0.14 0.24
0.001 0.25 0.05 0.1
0.5
5
4
1 10
4
2 3 5.0
4

KEY FEATURES:
1. Fast I/O with BufferedReader and PrintWriter
2. Efficient memory management
3. Clean object-oriented design
4. Binary search with optimal precision
5. Proper floating-point handling

COMPLEXITY:
- Time: O(Q * log(MAX_D0) * N)
- Space: O(N) where N = 7

OPTIMIZATIONS:
1. BufferedOutputStream for faster output
2. Direct mathematical calculations
3. Minimal object creation
4. Efficient loop structures
5. Early termination in binary search
*/
