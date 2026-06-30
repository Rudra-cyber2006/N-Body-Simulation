#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;

// The universal gravitational constant. 
const double G = 6.67430e-11; 

// The Barnes-Hut accuracy dial. 
// 0.0 means perfect accuracy but slow (O(N^2)). 0.5 is the standard sweet spot.
const double THETA = 0.5;

// The "anti-explosion" factor. Smooths out gravity if stars get dangerously close.
const double SOFTENING = 1.0e9; 

struct Body {
    std::string name;
    double mass;
    double x, y;      
    double vx, vy;    
    double fx, fy;    
};

struct TreeNode {
    double cx, cy, width, height;
    double x, y, mass;
    Body* body;
    TreeNode *tl, *tr, *bl, *br;

    TreeNode() = default;

    TreeNode(double center_x, double center_y, double w, double h) {
        cx = center_x; cy = center_y; width = w; height = h;
        x = 0.0; y = 0.0; mass = 0.0; body = nullptr;
        tl = tr = bl = br = nullptr;
    }

    void insert(Body* b, int depth = 0);
};

// --- THE BUMP ALLOCATOR ---
// Pre-allocates memory to bypass the incredibly slow 'new'/'delete' cycle.
std::vector<TreeNode> nodePool(500000); 
int nodeCount = 0;

TreeNode* allocateNode(double cx, double cy, double w, double h) {
    if (nodeCount >= nodePool.size()) {
        nodePool.resize(nodePool.size() * 2); // Failsafe growth
    }
    TreeNode* node = &nodePool[nodeCount++];
    *node = TreeNode(cx, cy, w, h);
    return node;
}
// ------------------------------

void TreeNode::insert(Body* newBody, int depth) {
    // The Infinite Zoom Shield: Stops dividing if stars are perfectly overlapping.
    if (depth > 50) return; 

    if (this->body == nullptr && this->tl == nullptr) {
        this->body = newBody;
        this->x = newBody->x;
        this->y = newBody->y;
        this->mass = newBody->mass;
        return;
    }

    if (this->tl == nullptr) {
        double halfW = width / 2.0;
        double halfH = height / 2.0;
        double qW = width / 4.0;
        double qH = height / 4.0;

        tl = allocateNode(cx - qW, cy + qH, halfW, halfH);
        tr = allocateNode(cx + qW, cy + qH, halfW, halfH);
        bl = allocateNode(cx - qW, cy - qH, halfW, halfH);
        br = allocateNode(cx + qW, cy - qH, halfW, halfH);

        Body* existingBody = this->body;
        this->body = nullptr;

        // Route existing body to sub-quadrant
        if (existingBody->x < cx) {
            if (existingBody->y < cy) bl->insert(existingBody, depth + 1);
            else tl->insert(existingBody, depth + 1);
        } else {
            if (existingBody->y < cy) br->insert(existingBody, depth + 1);
            else tr->insert(existingBody, depth + 1);
        }
    }

    // Update internal node's Center of Mass
    double newTotalMass = this->mass + newBody->mass;
    this->x = ((this->x * this->mass) + (newBody->x * newBody->mass)) / newTotalMass;
    this->y = ((this->y * this->mass) + (newBody->y * newBody->mass)) / newTotalMass;
    this->mass = newTotalMass;

    // Route new body to sub-quadrant
    if (newBody->x < cx) {
        if (newBody->y < cy) bl->insert(newBody, depth + 1);
        else tl->insert(newBody, depth + 1);
    } else {
        if (newBody->y < cy) br->insert(newBody, depth + 1);
        else tr->insert(newBody, depth + 1);
    }
}

void calculateTreeForce(Body& body, TreeNode* node) {
    if (node == nullptr || node->mass == 0.0) return;
    if (node->body == &body) return;

    double dx = node->x - body.x;
    double dy = node->y - body.y;
    double distSq = dx*dx + dy*dy + (SOFTENING * SOFTENING);
    double s = std::max(node->width, node->height);

    // Accept this node as a single approximation if EITHER:
    //  (a) it's a leaf containing exactly one other body — nothing further to recurse into, OR
    //  (b) it's an internal node that's "far enough" per the theta criterion
    bool isLeaf = (node->tl == nullptr);

    if (isLeaf || (s * s) < (THETA * THETA * distSq)) {
        double dist = std::sqrt(distSq);
        double force = (G * body.mass * node->mass) / distSq;
        body.fx += force * (dx / dist);
        body.fy += force * (dy / dist);
    } else {
        calculateTreeForce(body, node->tl);
        calculateTreeForce(body, node->tr);
        calculateTreeForce(body, node->bl);
        calculateTreeForce(body, node->br);
    }
}

void computeForcesBarnesHut(std::vector<Body>& bodies) {
    if (bodies.empty()) return;

    double minX = bodies[0].x, maxX = bodies[0].x;
    double minY = bodies[0].y, maxY = bodies[0].y;

    for (auto& b : bodies) {
        b.fx = 0.0;
        b.fy = 0.0;
        if (b.x < minX) minX = b.x;
        if (b.x > maxX) maxX = b.x;
        if (b.y < minY) minY = b.y;
        if (b.y > maxY) maxY = b.y;
    }

    double width = maxX - minX;
    double height = maxY - minY;
    double padding = std::max(width, height) * 0.02 + 1.0; 
    
    double centerX = (minX + maxX) / 2.0;
    double centerY = (minY + maxY) / 2.0;

    // Instant O(1) memory recycle for the new frame
    nodeCount = 0; 
    TreeNode* root = allocateNode(centerX, centerY, width + padding, height + padding);
    
    for (auto& b : bodies) {
        root->insert(&b); 
    }

    for (auto& b : bodies) {
        calculateTreeForce(b, root);
    }
}

void updateBodies(std::vector<Body>& bodies, double dt) {
    // Velocity Verlet - Step 1
    for (auto& body : bodies) {
        body.x += body.vx * dt + (body.fx / body.mass) * dt * dt/2;
        body.y += body.vy * dt + (body.fy / body.mass) * dt * dt/2;
        
        body.vx += (body.fx / body.mass) * dt/2;
        body.vy += (body.fy / body.mass) * dt/2;
    }
}

vector<Body> generateGalaxy(int numStars) {
    std::vector<Body> bodies;
    
    // Supermassive Black Hole
    bodies.push_back({"BlackHole", 1.0e34, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distPos(0.0, 5.0e11); 
    std::uniform_real_distribution<double> distMass(1.0e22, 1.0e26);
    std::uniform_real_distribution<double> distDivisor(0.5, sqrt(sqrt(2)));

    for (int i = 1; i < numStars; ++i) {
        double x, y, r;
        
        // The Danger Zone Fix: Do not spawn stars too close to the Black Hole!
        // This prevents infinite gravity slingshots that blow up the bounding box.
        // If the bounding box blows up then it will take more computation time in comparision to the trivial O(N^2)
        do {
            x = distPos(gen);
            y = distPos(gen);
            r = std::sqrt(x*x + y*y);
        } while (r < 1.5e11); 

        double mass = distMass(gen);
        double v = 0.0;
        if (r > 0) v = std::sqrt(G * bodies[0].mass / r);

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

    // Time Compression Fix: 4 hours per step (safe because Velocity Verlet is highly stable)
    double dt = 14400.0; 
    
    // Simulate 100 years. 6 steps = 1 simulated day.
    int steps = 100 * 365 * 6; 

    std::ofstream outFile("orbits.csv");

    for(int i = 0; i < bodies.size(); i++) {
        outFile << "X" << i << ",Y" << i;
        if (i != bodies.size() - 1) outFile << ",";
    }
    outFile << "\n";

    computeForcesBarnesHut(bodies);

    for (int i = 0; i < steps; i++) {
        
        updateBodies(bodies, dt);
        computeForcesBarnesHut(bodies);                      

        // Velocity Verlet - Step 2
        for (auto& body : bodies) {
            body.vx += (body.fx / body.mass) * dt/2;
            body.vy += (body.fy / body.mass) * dt/2;
        }

        // Hard Drive Saver: Log coordinates to CSV exactly once per simulated day
        if (i % 6 == 0) {
            for(int j = 0; j < bodies.size(); j++) {
                outFile << (bodies[j].x - bodies[0].x) << "," << (bodies[j].y - bodies[0].y);
                if (j != bodies.size() - 1) outFile << ",";
            }
            outFile << "\n";
        }
        
        // Console output: 2190 steps (365 * 6) is exactly one year
        if(i % 2190 == 0){
            cout << "Saving Coordinates.... Year: " << i / 2190 << "\n";
        }
    }

    outFile.close();
    std::cout << "Simulation complete. Data written to orbits.csv\n";
    return 0;
}