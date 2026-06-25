#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

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
        // Update velocity: v = v + (F/m) * dt
        body.vx += (body.fx / body.mass) * dt;
        body.vy += (body.fy / body.mass) * dt;
        
        // Update position: x = x + v * dt
        body.x += body.vx * dt;
        body.y += body.vy * dt;

    }
}

int main() {
    // Initializing the Sun, Earth, and Moon with real-world approximations
    vector<Body> bodies = {
        {"Sun", 1.989e30, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        // Earth is ~149.6 million km from Sun, moving at ~29.78 km/s
        {"Earth", 5.972e24, 1.496e11, 0.0, 0.0, 29780.0, 0.0, 0.0},
        // Moon is ~384,400 km from Earth, moving at Earth's vy + 1.022 km/s
        {"Moon", 7.342e22, 1.496e11 + 3.844e8, 0.0, 0.0, 29780.0 + 1022.0, 0.0, 0.0}
    };

    double dt = 86400.0; // Time step: 24 hours per tick
    int steps = 1000 * 24 * 365; // total time is dt*steps

    std::ofstream outFile("orbits.csv");
    outFile << "SunX,SunY,EarthX,EarthY,MoonX,MoonY\n";

    for (int i = 0; i < steps; i++) {
        computeForces(bodies);
        updateBodies(bodies, dt);
        for(int i=1;i<bodies.size();i++){
            bodies[i].x -= bodies[0].x;
            bodies[i].y -= bodies[0].y;
        }
        bodies[0].x = 0.0;
        bodies[0].y = 0.0;                      // moving relative to the first body

        // Write coordinates to CSV
        outFile << bodies[0].x << "," << bodies[0].y << ","
                << bodies[1].x << "," << bodies[1].y << ","
                << bodies[2].x << "," << bodies[2].y << "\n";
    }

    outFile.close();
    std::cout << "Simulation complete. Data written to orbits.csv\n";
    
    return 0;
}