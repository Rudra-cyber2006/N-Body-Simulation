#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;
const double G = 6.67430e-11; // Gravitational constant

struct Body {
    std::string name;
    double mass;
    double x, y;      // Position
    double vx, vy;    // Velocity
    double fx, fy;    // Accumulated force
};

void computeForces(std::vector<Body>& bodies) {
    int n = bodies.size();
    
    // Reset forces to zero
    for (int i = 0; i < n; ++i) {
        bodies[i].fx = bodies[i].fy = 0.0;
    }

    // O(N^2) Brute Force Calculation
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;

            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double distSq = dx*dx + dy*dy;
            double dist = std::sqrt(distSq);

            if (dist > 0) {
                // Newton's Law of Universal Gravitation
                double force = (G * bodies[i].mass * bodies[j].mass) / distSq;
                
                bodies[i].fx += force * (dx / dist);
                bodies[i].fy += force * (dy / dist);
            }
        }
    }
}

void updateBodies(std::vector<Body>& bodies, double dt) {
    // Forward Euler Integration
    for (auto& body : bodies) {
        

        body.x += body.vx * dt + (body.fx / body.mass) * dt * dt/2;
        body.y += body.vy * dt + (body.fy / body.mass) * dt * dt/2;

        
        body.vx += (body.fx / body.mass) * dt/2;
        body.vy += (body.fy / body.mass) * dt/2;
        
        

    }
}

vector<Body> generateGalaxy(int numStars) {
    std::vector<Body> bodies;
    
    // 1. The Supermassive Black Hole at the center
    bodies.push_back({"BlackHole", 1.0e34, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

    // Setup C++ random number generators
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Gaussian distribution for positions (dense center, sparse edges)
    // 1.0e11 is roughly the distance from Earth to Sun, so this makes a massive galaxy
    std::normal_distribution<double> distPos(0.0, 5.0e11); 
    
    // Uniform distribution for small star masses
    std::uniform_real_distribution<double> distMass(1.0e22, 1.0e26);

    std::uniform_real_distribution<double> distDivisor(0.5, sqrt(sqrt(2)));


    for (int i = 1; i < numStars; ++i) {
        double x = distPos(gen);
        double y = distPos(gen);
        double mass = distMass(gen);

        // Calculate distance from center
        double r = std::sqrt(x*x + y*y);

        // Calculate orbital velocity to keep it in a stable spin
        // v = sqrt(G * Mass_Center / r)
        double v = 0.0;
        if (r > 0) {
            v = std::sqrt(G * bodies[0].mass / r);
        }

        // To make it spin, the velocity vector must be perpendicular to the position vector
        // If position is (x,y), a perpendicular vector is (-y, x)

        double randDivisor = distDivisor(gen);
        double randDivisor2 = distDivisor(gen);
        
        double vx = -v * (y / r) * randDivisor;
        double vy = v * (x / r) * randDivisor2;

        bodies.push_back({"Star_" + std::to_string(i), mass, x, y, vx, vy, 0.0, 0.0});
    }

    return bodies;
}



int main() {
    vector<Body> bodies = generateGalaxy(100);

    double dt = 3600.0; 
    

    int steps = 100 * 365 * 24; 

    std::ofstream outFile("orbits.csv");

    for(int i = 0; i < bodies.size(); i++) {
        outFile << "X" << i << ",Y" << i;
        if (i != bodies.size() - 1) outFile << ",";
    }
    outFile << "\n";

    computeForces(bodies);

    for (int i = 0; i < steps; i++) {
        updateBodies(bodies, dt);
        computeForces(bodies);                      

        for (auto& body : bodies) {
            body.vx += (body.fx / body.mass) * dt/2;
            body.vy += (body.fy / body.mass) * dt/2;
        }

        // 3. The I/O Shield: Only save the data once every 3 hours (3 steps)
        if (i % 3 == 0) {
            for(int j = 0; j < bodies.size(); j++) {
                outFile << (bodies[j].x - bodies[0].x) << "," << (bodies[j].y - bodies[0].y);
                if (j != bodies.size() - 1) outFile << ",";
            }
            outFile << "\n";
            
        }
        if(i%10000==0){
            cout << "Saving Coordinates...." << " " << i/10000 << "\n";
        }

    }

    outFile.close();
    std::cout << "Simulation complete. Data written to orbits.csv\n";
    return 0;
}
