#include <vector>
#include<cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>



// Define the Particle structure
struct Particle {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
};

// Declare global variables (avoid multiple definitions)
extern std::vector<Particle> particles;


void initialize_particles(std::vector<Particle> &particles, int num_particles) {
    particles.resize(num_particles);
    if (num_particles == 2) {
        // Sun-Earth system
        particles[0] = {1.989e30, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        particles[1] = {5.972e24, 1.496e11, 0, 0, 0, 29780, 0, 0, 0, 0};
    } else if (num_particles == 3) {
        // Sun-Earth-Moon system
        particles[0] = {1.989e30, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        particles[1] = {5.972e24, 1.496e11, 0, 0, 0, 29780, 0, 0, 0, 0};
        particles[2] = {7.347e22, 1.496e11 + 3.844e8, 0, 0, 0, 29780 + 1022, 0, 0, 0, 0};
    } else {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(1.0, 2.0); // for testing purpose, small random scale
        
        for (auto &p : particles) {
            p.mass = 1.0 + dist(gen); 
            p.x = dist(gen);
            p.y = dist(gen);
            p.z = dist(gen);
            p.vx = p.vy = p.vz = 0.0;
            p.fx = p.fy = p.fz = 0.0;
        }
    
    }
}



const double G = 6.67430e-11;
const double SOFTENING = 1e-9;

void computeForces(std::vector<Particle> &particles) {
    // Reset forces to zero
    for (auto &p : particles) {
        p.fx = 0;
        p.fy = 0;
        p.fz = 0;
    }

    // Calculate forces between each pair of particles
    for (size_t i = 0; i < particles.size(); i++) {
        for (size_t j = 0; j < particles.size(); j++) {
            if(i == j) continue;

            Particle &p1 = particles[i];
            Particle &p2 = particles[j];

            // Compute distance components
            double dx = p2.x - p1.x;
            double dy = p2.y - p1.y;
            double dz = p2.z - p1.z;
            double r2 = dx * dx + dy * dy + dz * dz + SOFTENING;
            double r = std::sqrt(r2);

            
            double forceMagnitude = (G * p1.mass * p2.mass) / r2;
            double fx = forceMagnitude * (dx / r);
            double fy = forceMagnitude * (dy / r);
            double fz = forceMagnitude * (dz / r);

            // Apply Newton's Third Law (equal and opposite forces)
            p1.fx += fx;
            p1.fy += fy;
            p1.fz += fz;

            p2.fx -= fx;
            p2.fy -= fy;
            p2.fz -= fz;
            
        }
    }
}

void integrate(std::vector<Particle> &particles, double dt) {
    for (auto &p : particles) {
        double ax = p.fx / p.mass;
        double ay = p.fy / p.mass;
        double az = p.fz / p.mass;
        
        p.vx += ax * dt;
        p.vy += ay * dt;
        p.vz += az * dt;
        
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
    }
}

void log_state(const std::vector<Particle> &particles, std::ofstream &outfile) {
    outfile << particles.size();
    for (const auto &p : particles) {
        outfile << "\t" << p.mass << "\t" << p.x << "\t" << p.y << "\t" << p.z
                << "\t" << p.vx << "\t" << p.vy << "\t" << p.vz
                << "\t" << p.fx << "\t" << p.fy << "\t" << p.fz;
    }
    outfile << "\n";
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <num_particles> <time_step> <iterations> <log_interval>" << std::endl;
        return 1;
    }
    
    int num_particles = std::stoi(argv[1]);
    double dt = std::stod(argv[2]);
    int iterations = std::stoi(argv[3]);
    int log_interval = std::stoi(argv[4]);
    
    // start the timer

    auto start = std::chrono::high_resolution_clock::now();
 

    std::vector<Particle> particles(num_particles);
    initialize_particles(particles, num_particles);
    
    std::ofstream outfile("solar.tsv");
    
    for (int step = 0; step < iterations; step++) {
        computeForces(particles);
        integrate(particles, dt);

        // Dump state at every `log_interval` steps
        if (step % log_interval == 0) {
            log_state(particles, outfile);
            // std::cout << "Logged iteration: " << step << std::endl;  // Optional console update
        }
    }
    
    outfile.close();


    // End the timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Output the duration
    std::cout << "Task took " << duration.count() << " milliseconds.\n";

    

    return 0;
}





