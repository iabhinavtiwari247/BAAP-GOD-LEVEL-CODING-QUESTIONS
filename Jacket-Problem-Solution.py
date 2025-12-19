"""
Composite Jacket Thermal Resistance - Complete Solution
Time Complexity: O(Q * log(MAX_D0)) for Type 4 queries, O(N) for propagation queries
Space Complexity: O(N) where N = 6 layers
"""

import sys
import math
from typing import List, Tuple

class Layer:
    """Represents a single layer in the jacket"""
    def __init__(self, d: float, k: float, mu: float, c: float):
        self.d = d          # thickness (meters)
        self.k = k          # base conductivity (W/m·K)
        self.mu = mu        # moisture coefficient
        self.c = c          # compression coefficient
        self.W = 0.0        # accumulated moisture
        self.C = 0.0        # accumulated compression
    
    def k_eff(self, beta: float) -> float:
        """Calculate effective conductivity"""
        return self.k * (1.0 + self.mu * self.W) * math.exp(beta * self.C)
    
    def thermal_resistance(self, beta: float) -> float:
        """Calculate thermal resistance of this layer"""
        if self.d == 0:
            return 0.0
        return self.d / self.k_eff(beta)


class ThermalJacketSystem:
    """Main system to handle all queries"""
    
    def __init__(self):
        self.layers: List[Layer] = []
        self.beta: float = 0.0
        self.T_s = 37.0      # body temperature (°C)
        self.T_ext = -13.0   # ambient temperature (°C)
        self.q_max = 20.0    # maximum heat flux (W/m²)
        self.moisture_rate = 0.01  # moisture per hour
        
    def read_input(self):
        """Read initial layer properties and beta"""
        # Read 7 layers (index -1 to 5)
        for i in range(7):
            line = input().split()
            d, k, mu, c = map(float, line)
            self.layers.append(Layer(d, k, mu, c))
        
        # Read beta
        self.beta = float(input())
    
    def calculate_R_th_without_layer(self, exclude_layer: int) -> float:
        """Calculate total thermal resistance excluding a specific layer"""
        R_th = 0.0
        for i in range(7):
            if i != exclude_layer:
                R_th += self.layers[i].thermal_resistance(self.beta)
        return R_th
    
    def calculate_R_th(self) -> float:
        """Calculate total thermal resistance"""
        return sum(layer.thermal_resistance(self.beta) for layer in self.layers)
    
    def heat_flux(self, d0: float) -> float:
        """Calculate heat flux for a given d0 (layer 0 thickness)"""
        # Calculate R_th with layer 0 having thickness d0
        R_th = self.calculate_R_th_without_layer(1)  # layer index 1 is layer 0
        
        # Add layer 0's contribution
        layer0 = self.layers[1]
        if d0 > 0:
            k0_eff = layer0.k_eff(self.beta)
            R_th += d0 / k0_eff
        
        # q = ΔT / R_th
        delta_T = self.T_s - self.T_ext  # 50°C
        if R_th < 1e-12:  # Avoid division by zero
            return 1e18
        return delta_T / R_th
    
    def query_type1_environmental_exposure(self, t: float):
        """Type 1: Simulate t hours of snowfall"""
        # Add moisture to layer 5 (index 6)
        self.layers[6].W += t * self.moisture_rate
        
        # Propagate moisture inward (layer 5 to layer -1)
        # Direction: 6 -> 5 -> 4 -> 3 -> 2 -> 1 -> 0
        for i in range(6, 0, -1):  # From layer 5 down to layer 0
            layer_i = self.layers[i]
            layer_prev = self.layers[i - 1]
            # W_(i-1) += W_i × μ_i
            layer_prev.W += layer_i.W * layer_i.mu
    
    def query_type2_mechanical_stress(self, i: int, x: float):
        """Type 2: Apply compression stress x to layer i"""
        # Convert layer index (-1 to 5) to array index (0 to 6)
        layer_idx = i + 1
        
        # Add compression to layer i
        self.layers[layer_idx].C += x
        
        # Propagate compression outward (layer i to layer 5)
        for idx in range(layer_idx, 6):  # From layer i to layer 4
            layer_curr = self.layers[idx]
            layer_next = self.layers[idx + 1]
            # C_(i+1) += C_i × c_i
            layer_next.C += layer_curr.C * layer_curr.c
    
    def query_type3_replace_layer(self, i: int, d: float, k: float, mu: float, c: float):
        """Type 3: Replace all properties of layer i"""
        # Convert layer index (-1 to 5) to array index (0 to 6)
        layer_idx = i + 1
        
        # Replace properties
        self.layers[layer_idx].d = d
        self.layers[layer_idx].k = k
        self.layers[layer_idx].mu = mu
        self.layers[layer_idx].c = c
        
        # Reset accumulated values
        self.layers[layer_idx].W = 0.0
        self.layers[layer_idx].C = 0.0
    
    def query_type4_minimum_foam_thickness(self) -> float:
        """Type 4: Calculate minimum d0 such that q ≤ 20"""
        # Binary search on d0
        left = 0.0
        right = 1e9  # Large upper bound
        
        # Check if even d0 = 0 satisfies the constraint
        if self.heat_flux(0.0) <= self.q_max + 1e-9:
            return 0.0
        
        # Check if even d0 = infinity cannot satisfy
        # This happens if k0_eff is too high or other layers don't provide enough resistance
        if self.heat_flux(right) > self.q_max + 1e-9:
            return 1e18  # Return very large number
        
        # Binary search
        epsilon = 1e-9
        iterations = 0
        max_iterations = 100
        
        while right - left > epsilon and iterations < max_iterations:
            mid = (left + right) / 2.0
            q = self.heat_flux(mid)
            
            if q <= self.q_max + 1e-10:
                right = mid
            else:
                left = mid
            
            iterations += 1
        
        return right
    
    def query_type5_feasibility_check(self) -> str:
        """Type 5: Check if any finite d0 can satisfy the constraint"""
        # Check with d0 = 0
        q_min = self.heat_flux(0.0)
        if q_min <= self.q_max + 1e-9:
            return "POSSIBLE"
        
        # Check with very large d0
        q_with_large_d0 = self.heat_flux(1e9)
        if q_with_large_d0 <= self.q_max + 1e-9:
            return "POSSIBLE"
        
        # Check the limit as d0 -> infinity
        # As d0 increases, R_th increases, so q decreases
        # If q still doesn't go below q_max, it's impossible
        
        # Calculate if the trend is decreasing
        q1 = self.heat_flux(1e6)
        q2 = self.heat_flux(1e7)
        
        if q2 < q1 and q2 <= self.q_max + 1e-9:
            return "POSSIBLE"
        
        # If heat flux is still too high with very large d0
        if q_with_large_d0 > self.q_max - 1e-9:
            return "IMPOSSIBLE"
        
        return "POSSIBLE"
    
    def process_queries(self):
        """Process all queries"""
        Q = int(input())
        
        for _ in range(Q):
            query = input().split()
            query_type = int(query[0])
            
            if query_type == 1:
                # Environmental exposure
                t = float(query[1])
                self.query_type1_environmental_exposure(t)
            
            elif query_type == 2:
                # Mechanical stress
                i = int(query[1])
                x = float(query[2])
                self.query_type2_mechanical_stress(i, x)
            
            elif query_type == 3:
                # Replace layer
                i = int(query[1])
                d = float(query[2])
                k = float(query[3])
                mu = float(query[4])
                c = float(query[5])
                self.query_type3_replace_layer(i, d, k, mu, c)
            
            elif query_type == 4:
                # Minimum foam thickness
                d0 = self.query_type4_minimum_foam_thickness()
                print(f"{d0:.10f}")
            
            elif query_type == 5:
                # Feasibility check
                result = self.query_type5_feasibility_check()
                print(result)


def main():
    """Main function"""
    system = ThermalJacketSystem()
    system.read_input()
    system.process_queries()


if __name__ == "__main__":
    main()


"""
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

EXPECTED OUTPUT (approximate):
(depends on exact calculations)
0.XXXXXXXXXX
0.YYYYYYYYYY
0.ZZZZZZZZZZ

EXPLANATION:
1. Query 1 (Type 4): Calculate initial minimum d0
2. Query 2 (Type 1): Add moisture from 10 hours of snowfall
3. Query 3 (Type 4): Recalculate d0 after moisture changes k_eff
4. Query 4 (Type 2): Add compression stress to layer 3
5. Query 5 (Type 4): Recalculate d0 after compression changes k_eff

KEY OPTIMIZATIONS:
1. Pre-calculation: Calculate R_th without layer 0, then add layer 0's contribution
2. Binary search: Efficient O(log(MAX_D0)) search for minimum thickness
3. Efficient propagation: Only update affected layers
4. Direct formulas: Use mathematical properties of exponential and linear functions
5. Epsilon handling: Careful floating-point comparisons

COMPLEXITY ANALYSIS:
- Type 1 query: O(N) where N = 7 layers
- Type 2 query: O(N) for propagation
- Type 3 query: O(1)
- Type 4 query: O(log(MAX_D0) * N) for binary search with heat flux calculation
- Type 5 query: O(N) with multiple heat flux calculations
- Overall: O(Q * log(MAX_D0) * N) which is acceptable for Q ≤ 200,000

FURTHER OPTIMIZATIONS FOR EXTREME CASES:
1. Cache k_eff values and update incrementally
2. Store partial sums of thermal resistance
3. Use analytical solution instead of binary search where possible
4. Batch similar queries together
"""
