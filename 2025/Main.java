import javax.swing.*;
import javax.swing.Timer;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

// NOTE: Removed delta times and tied all gravity directly to the Java.util timer - 3/7/2025

// For Pi Day 2025 (3/14/2025): The variable collisions tracks the collisions between the spheres and the right wall!
// To calculate Pi with this, we can use the formula pi = collisions / sqrt(m2 / m1) (usually take the floor of this number)
// See 3b1b's video for more details: https://www.youtube.com/watch?v=6dTyOl1fmDo
// Note: Collisions will break after a certain point, don't make the masses too high, I haven't coded this for super high precision
public class Main {
    public static double collisions = 0;
    public static void main(String[] args) {
        JFrame f = new JFrame();
        
        OrbitalSystem os = new OrbitalSystem(0.01, 1000); // Speed = how many times to update per second, up to 1000
        // os.addPlanet(100, 100, 10, 10);
        // os.addPlanet(200, 200, 10, 10);
        // os.addPlanet(200, 400, 10, 20);
        // os.planets.get(0).velocity[0] += 0.025;
        // os.planets.get(0).velocity[1] += 0.025;
        // os.planets.get(0).velocity[0] += 8;
        // os.planets.get(1).velocity[0] -= 8;
        


        // Testing solution to planets phasing through each other
        os.addPlanet(200, 200, 1, 10);
        os.addPlanet(150, 200, 64, 10);
        // os.planets.get(1).velocity[1] += -0.5;
        os.planets.get(1).velocity[0] += 0.05;



        // Testing RigidLines
        // os.lines.add(new RigidLine(100, 100, 150, 150));
        // os.lines.add(new RigidLine(200, 100, 250, 150));
        // os.addPlanet(125, 50, 20, 10);
        // os.planets.get(0).velocity[1] += 0.08;


        // Barriers around the screen
        os.lines.add(new RigidLine(0, 0, 500, 0));
        // os.lines.add(new RigidLine(0, 0, 0, 500));
        os.lines.add(new RigidLine(500, 0, 500, 500));
        os.lines.add(new RigidLine(0, 500, 500, 500));


        // os.lines.add(new RigidLine(100, 100, 200, 120));
        // os.lines.add(new RigidLine(80, 200, 100, 100));
        // os.lines.add(new RigidLine(200, 120, 180, 220));
        // os.lines.add(new RigidLine(100, 100, 200, 120));

        // os.addLine(175, 225, 225, 175);
        // os.lines.add(new RigidLine(275, 225, 350, 175));
        

        // for (int i = 0; i < 10; i++) {
        //     double massAndRadius = Math.random() * 25;
        //     os.addPlanet(Math.random() * 500, Math.random() * 500, massAndRadius, massAndRadius);
        //     os.planets.get(i).velocity[0] = 0.05 * (Math.random() - 0.5);
        //     os.planets.get(i).velocity[1] = 0.05 * (Math.random() - 0.5);
        // }

        // os.addPlanet(200, 200, 50, 10);
        // os.addPlanet(200, 100, 1, 3);
        // os.planets.get(1).velocity[0] = 0.0707106781187; // Velocity required for circular orbit calculated using v = sqrt(GM/r) where m = mass of the body being orbited (50 here)
        // os.planets.get(0).fixed = true;



        // os.addPlanet(200, 200, 10, 25);
        // os.addPlanet(200, 203, 3, 3);



        
        JPanel mainPanel = new JPanel() {
            @Override
            public void paintComponent(Graphics g) {
                Graphics2D g2d = (Graphics2D) g;
                // os.updateSystem();

                ArrayList<Planet> planets = new ArrayList<Planet>();
                for (int i = 0; i < os.planets.size(); i++) {
                    planets.add(os.planets.get(i).clone());
                }

                ArrayList<RigidLine> lines = new ArrayList<RigidLine>();
                for (int i = 0; i < os.lines.size(); i++) {
                    lines.add(os.lines.get(i)); // not cloning for now because I'm lazy
                }
                // planets.get(0).position[0] = 200;
                // planets.get(0).position[1] = 200;

                for (Planet p : planets) {
                    int radius = (int) (p.radius + 0.5);
                    int x = (int) (p.position[0] - radius + 0.5);
                    int y = (int) (p.position[1] - radius + 0.5);
                    g2d.drawOval(x, y, radius * 2, radius * 2);
                }

                for (RigidLine l : lines) {
                    int x1 = (int) (l.position[0][0] + 0.5);
                    int y1 = (int) (l.position[0][1] + 0.5);
                    int x2 = (int) (l.position[1][0] + 0.5);
                    int y2 = (int) (l.position[1][1] + 0.5);
                    g2d.drawLine(x1, y1, x2, y2);
                }

                double sumP = 0;
                double sumKE = 0;
                for (int i = 0; i < planets.size(); i++) {
                    Planet p = planets.get(i);
                    double[] momentum = {p.velocity[0] * p.mass, p.velocity[1] * p.mass};
                    double[] ke = {Math.pow(p.velocity[0], 2) * p.mass * (0.5), Math.pow(p.velocity[1], 2) * p.mass * (0.5)};
                    sumP += momentum[0];
                    sumP += momentum[1];
                    sumKE += ke[0];
                    sumKE += ke[1];
                }

                // System.out.println(sumKE);
                // System.out.println(sumP);
                // System.out.println(planets.get(0).dist(planets.get(1)));
                // System.out.println(os.planets.get(0).velocity[1]);
                System.out.println(collisions);
            }
        };
        
        mainPanel.setBounds(new Rectangle(500, 500));
        
        f.add(mainPanel);
        f.pack();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setSize(500, 500);
        f.setVisible(true);
        
        Timer updater = new Timer(1, new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                f.repaint();
            }
        });
        updater.start();
        
    }
}

class OrbitalSystem {
    public ArrayList<Planet> planets = new ArrayList<Planet>();
    public ArrayList<RigidLine> lines = new ArrayList<RigidLine>();
    public double gravitationalConstant = 0.01;

    // Delta timers are currently deprecated, may add back later
    // public long startTime = System.currentTimeMillis(); // May update to nanoTime for more precision
    // public long lastUpdateTime = startTime;
    // public double deltaTime = startTime;
    
    public double speed = 1;

    // Variables for preventing clipping (one object inside another) by pushing objects inside each other out of each other
    public boolean preventClipping = true;

    // Fraction of sum of radii of two objects,
    // when the distance between two objects is less than radiiSum * clippingPreventionThreshold, clipping prevention will be applied
    // and the objects will be pushed away from each other
    // Recommended to keep below 1, but conservation of momentum will be conserved regardless of the value
    public double clippingPreventionThreshold = 0.1;

    // Strength of clippingPrevention -- Force applied to each object when clippingPreventionThreshold is met
    public double clippingPreventionStrength = 0.1;
    
    public OrbitalSystem(double gravitationalConstant) {
        this.gravitationalConstant = gravitationalConstant;
        this.startTimer();
    }

    public OrbitalSystem(double gravitationalConstant, double speed) {
        this.gravitationalConstant = gravitationalConstant;
        this.speed = speed;
        this.startTimer();
    }

    public void startTimer() {
        new java.util.Timer().schedule(new TimerTask() {
            @Override
            public void run() {
                updateSystem();
            }
        }, 0, (int) Math.max(1000 / speed, 1));
    }
    
    public void addPlanet(double x, double y, double mass, double radius) {
        Planet newPlanet = new Planet(x, y, mass, radius, this);
        planets.add(newPlanet);
    }
    
    public void addLine(double x1, double y1, double x2, double y2) {
        lines.add(new RigidLine(x1, y1, x2, y2));
    }
    
    public void addPlanet(Planet newPlanet) {
        planets.add(newPlanet);
    }
    
    public void updateSystem() {
        for (int i = 0; i < planets.size(); i++) {
            updateAccelerationOfPlanet(i);
        }

        for (int i = 0; i < lines.size(); i++) {
            updateLine(i);
        }
        
        for (int i = 0; i < planets.size(); i++) {
            updatePlanet(i);
        }
    }
    
    // Updates the position of the given planet
    public void updatePlanet(int index) {
        Planet p = planets.get(index);
        p.updatePosition();
    }

    // Updates the position of the given planet
    public void updateAccelerationOfPlanet(int index) {
        Planet p = planets.get(index);
        for (int i = 0; i < planets.size(); i++) {
            if (i != index) {
                Planet other = planets.get(i);
                p.updateAcceleration(other);
            }
        }
    }

    public void updateLine(int index) {
        RigidLine l = lines.get(index);
        for (int i = 0; i < planets.size(); i++) {
            RigidBodySphere other = planets.get(i).rigidBody;
            l.updateAcceleration(other);
        }
    }
     
    // Updates the delta time, deprecated
    // public void updateDeltaTime() {
    //     long currentTime = System.currentTimeMillis();
    //     deltaTime = (double) (currentTime - lastUpdateTime) / 1000;
    //     lastUpdateTime = System.currentTimeMillis();
    // }
}

class Planet {
    public double mass = 1; // default 1
    public double radius = 10; // default 10
    public double[] position = {0, 0};
    public double[] velocity = {0, 0};
    public double[] acceleration = {0, 0};
    public OrbitalSystem parent = null;
    public RigidBodySphere rigidBody = null;
    public boolean fixed = false; // May remove later or keep in, just keeps the position and velocity of the planet fixed no matter its acceleration
    
    // Precondition: x, y, and mass must not be null
    public Planet(double x, double y, double mass, double radius, OrbitalSystem parent) {
        this.position[0] = x;
        this.position[1] = y;
        this.mass = mass;
        this.radius = radius;
        this.parent = parent;
        this.rigidBody = new RigidBodySphere(position, velocity, acceleration, mass, radius, parent);
        this.resetVelocity();
        this.resetAcceleration();
    }

    // Returns the total acceleration to be applied on this object as a sum of all forces affecting it (including collisions)
    public void updateAcceleration(Planet other) {
        if (this.rigidBody.isColliding(other.rigidBody)) {
            this.rigidBody.collision(other.rigidBody);
            Main.collisions += 0.5;
        } else {
            // this.rigidBody.gravity(other.rigidBody);
        }
    }
    
    // Updates the position of the planet and resets its acceleration
    public void updatePosition() {
        if (fixed) {
            this.resetVelocity();
            this.resetAcceleration();
        }

        this.velocity[0] += this.acceleration[0];
        this.velocity[1] += this.acceleration[1];
        
        this.position[0] += this.velocity[0]/*  * deltaTime */;
        this.position[1] += this.velocity[1]/*  * deltaTime */;
        
        this.resetAcceleration();
    }
    
    // Returns distance between two planets
    public double dist(Planet other) {
        double x1 = this.position[0];
        double y1 = this.position[1];
        double x2 = other.position[0];
        double y2 = other.position[1];
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
    }
    
    // Sets the values for velocity to 0
    public void resetVelocity() {
        this.velocity[0] = 0;
        this.velocity[1] = 0;
    }
    
    // Sets the values for acceleration to 0
    public void resetAcceleration() {
        this.acceleration[0] = 0;
        this.acceleration[1] = 0;
    }

    @Override
    public Planet clone() {
        Planet newPlanet = new Planet(this.position[0], this.position[1], this.mass, this.radius, this.parent);
        double[] newVelocity = {this.velocity[0], this.velocity[1]};
        double[] newAcceleration = {this.acceleration[0], this.acceleration[1]};
        newPlanet.velocity = newVelocity;
        newPlanet.acceleration = newAcceleration;
        return newPlanet;
    }
}

class RigidBodySphere {
    public double mass = 0;
    public double radius = 0;
    public double[] position = {0, 0};
    public double[] velocity = {0, 0};
    public double[] acceleration = {0, 0};
    public OrbitalSystem parent = null;

    public RigidBodySphere(double x, double y, double mass, double radius, OrbitalSystem parent) {
        this.position[0] = x;
        this.position[1] = y;
        this.mass = mass;
        this.radius = radius;
        this.parent = parent;
    }

    public RigidBodySphere(double[] position, double mass, double radius, OrbitalSystem parent) {
        this.position = position;
        this.mass = mass;
        this.radius = radius;
        this.parent = parent;
    }

    public RigidBodySphere(double[] position, double[] velocity, double mass, double radius, OrbitalSystem parent) {
        this.position = position;
        this.velocity = velocity;
        this.mass = mass;
        this.radius = radius;
        this.parent = parent;
    }
    
    public RigidBodySphere(double[] position, double[] velocity, double[] acceleration, double mass, double radius, OrbitalSystem parent) {
        this.position = position;
        this.velocity = velocity;
        this.acceleration = acceleration;
        this.mass = mass;
        this.radius = radius;
        this.parent = parent;
    }

    // Applies the acceleration due to collision applied on this RigidBody from another RigidBody
    public void collision(RigidBodySphere other) {
        double m1 = this.mass;
        double m2 = other.mass;

        double[] dv = this.collisionDV(other);

        this.acceleration[0] += dv[0];
        this.acceleration[1] += dv[1];
    }

    public boolean isColliding(RigidBodySphere other) {
        double dist = this.dist(other);
        double radiiSum = this.radius + other.radius;

        // Checking if the spheres are touching/inside one another (ayo)
        if (dist > radiiSum) {
            return false;
        }
        
        double x1 = this.position[0];
        double y1 = this.position[1];
        double x2 = other.position[0];
        double y2 = other.position[1];

        double thisVelocityX = this.velocity[0];
        double thisVelocityY = this.velocity[1];
        double otherVelocityX = other.velocity[0];
        double otherVelocityY = other.velocity[1];

        double[] closestPoint = this.findClosestPoint(other);

        
        double[] distanceVector = {(x2 - x1), (y2 - y1)}; // Not sure if I need to normalize this vector but I don't think it matters in this case
        // double[] distanceVector = {(x2 - closestPoint[0]), (y2 - closestPoint[1])}; // Not currently working
        
        // NOTE: Not sure if I should be accounting for velocity next frame here, like in the isColliding() method of RigidLine, or not
        // For now, I am and I'll see if there are any problems   EDIT: I think there are problems, removing for now
        // double[] relativeVelocity = {(otherVelocityX + other.acceleration[0]) - (thisVelocityX + this.acceleration[0]), (otherVelocityY + other.acceleration[1]) - (thisVelocityY + this.acceleration[1])};
        double[] relativeVelocity = {otherVelocityX - thisVelocityX, otherVelocityY - thisVelocityY};

        double dot = this.dot(distanceVector, relativeVelocity);
        // Not sure if this is needed/useful, it just checks if the distance between the two objects is decreasing (if other is moving into this)
        if (dot > 0) {
            return false;
        }
        
        // This is just a check to see if other is stuck in this and applying a force if it is, to separate them eventually
        if (parent.preventClipping && dist <= radiiSum - (radiiSum * parent.clippingPreventionThreshold)) {
            double[] normal = other.normal(this);
            normal[0] *= parent.clippingPreventionStrength;
            normal[1] *= parent.clippingPreventionStrength;
            this.acceleration[0] -= normal[0] / this.mass;
            this.acceleration[1] -= normal[1] / this.mass;

            other.acceleration[0] += normal[0] / other.mass;
            other.acceleration[1] += normal[1] / other.mass;
        }

        return true;
    }

    // Finds the closest point on the line perpendicular to the velocity of other
    public double[] findClosestPoint(RigidBodySphere other) {
        double x1 = other.position[0];
        double y1 = other.position[1];
        double x2 = x1 + other.velocity[0];
        double y2 = y1 + other.velocity[1];

        double pm = (y2 - y1) / (x2 - x1); // Perpendicular slope = negative reciprocal of slope
        double m = -1 / pm; // Slope of this line

        // Accounting for infinite or zero slope
        if (m == 0) {
            double[] result = {x1, this.position[1]};
            return result;
        } else if (!Double.isFinite(m)) {
            double[] result = {this.position[0], y1};
            return result;
        }

        double b = -(m * x1) + y1; // yint of this line
        // double pb = -(pm * x1); // yint of perpendicular line

        double h = x1;
        double k = y1;

        // Equation of perpendicular line passing through point p is y = pm(x - h) + k
        // pm(x - h) + k = mx + b, solve to find intersection x-value
        double leftX = pm;
        double leftConstants = pm * -h + k;
        double rightX = m;
        double rightConstants = b;

        double totalX = rightX - leftX;
        double totalConstants = leftConstants - rightConstants;

        double x = totalConstants / totalX;
        double y = m * x + b;
    
        double[] result = {x, y};
        return result;
    }

    public double[] collisionDV(RigidBodySphere other) {
        double[] normal = this.normal(other);
        double m1 = this.mass;
        double m2 = other.mass;

        double thisVelocityX = this.velocity[0];
        double thisVelocityY = this.velocity[1];
        double otherVelocityX = other.velocity[0];
        double otherVelocityY = other.velocity[1];
        
        double[] distanceVector = other.normal(this); // I think I have to use other.normal(this) because the y space is flipped for Java graphics space -- not sure
        double[] relativeMomentum = {(otherVelocityX * m2) - (thisVelocityX * m1), (otherVelocityY * m2) - (thisVelocityY * m1)};
        // We need to find the dv that keeps total kinetic energy and total momentum equivalent
        // See this link for more details: https://physics.stackexchange.com/questions/81959/perfect-elastic-collision-and-velocity-transfer
        double[] v1 = {(this.velocity[0] * (this.mass - other.mass) + (2 * other.mass * other.velocity[0])) / (other.mass + this.mass), (this.velocity[1] * (this.mass - other.mass) + (2 * other.mass * other.velocity[1])) / (other.mass + this.mass)};
        // double[] v2 = {other.velocity[0] * (other.mass - this.mass) + (2 * this.mass * this.velocity[0]), other.velocity[1] * (other.mass - this.mass) + (2 * this.mass * this.velocity[1])};

        // double dot = this.dot(distanceVector, relativeMomentum);
        
        // normal[0] *= -dot;
        // normal[1] *= -dot;
        // normal[0] /= m1;
        // normal[1] /= m1;
        // return normal;

        double[] dv = {v1[0] - this.velocity[0], v1[1] - this.velocity[1]};
        return dv;
    }

    public double[] normal(RigidBodySphere other) {
        double dist = this.dist(other);
        double[] normal = {(this.position[0] - other.position[0]) / dist, (this.position[1] - other.position[1]) / dist};
        // normal[0] *= this.radius;
        // normal[1] *= this.radius;
        return normal;
    }

    // Returns acceleration due to gravity exerted on this planet from another planet p as a vector
    public void gravity(RigidBodySphere other) {
        double dist = this.dist(other);
        double m1 = this.mass;
        double m2 = other.mass;
        double g = this.parent.gravitationalConstant;
        double speed = parent.speed;
        
        double forceMagnitude = g * ((m1 * m2) / (Math.pow(dist, 2)));
        double[] normal = this.normal(other);
        double[] force = {normal[0] * forceMagnitude, normal[1] * forceMagnitude};
        double[] acceleration = {force[0] / m1, force[1] / m1};
        this.acceleration[0] -= acceleration[0];
        this.acceleration[1] -= acceleration[1];
    }

    // Returns distance between two rigid bodies
    public double dist(RigidBodySphere other) {
        double x1 = this.position[0];
        double y1 = this.position[1];
        double x2 = other.position[0];
        double y2 = other.position[1];
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
    }

    // Returns the dot product of two 2D vectors
    public double dot(double[] v1, double[] v2) {
        return (v1[0] * v2[0]) + (v1[1] * v2[1]);
    }

    public double[] rotate(double[] v, double radians) {
        double sin = Math.sin(radians);
        double cos = Math.cos(radians);
        double x = v[0];
        double y = v[1];
        double[] result = {(x * cos) - (y * sin), (x * sin) + (y * cos)};
        return result;
    }
}

// A line with collision, stays at a fixed point and does not have mass but can have velocity and acceleration
// Right now, I'm checking collision in this class, NOT the RigidBodySphere class
class RigidLine {
    public double[][] position = {{0, 0}, {0, 0}}; // First row is the first point, second row is the second point
    
    // All of the following are not in use for now
    public double[][] velocity = {{0, 0}, {0, 0}}; // First row is the first point, second v is the second point
    public double[][] acceleration = {{0, 0}, {0, 0}}; // First row is the first point, second row is the second point -- Not reset every frame (for now)
    public double[] angularVelocity = {0, 0}; // First row is the first point, second row is the second point
    public double[] angularAcceleration = {0, 0}; // First row is the first point, second row is the second point -- Not reset every frame (for now)
    public RigidLine(double x1, double y1, double x2, double y2) {
        double[][] temp = {{x1, y1}, {x2, y2}};
        this.position = temp;
    }

    public void updateAcceleration(RigidBodySphere other) {
        if (this.isColliding(other)) {
            this.collision(other);
            Main.collisions += 1;
        }
    }

    public void collision(RigidBodySphere other) {
        double[] normal = this.normal(other);
        // double[] closestPoint = this.findClosestPoint(other.position);

        // double x1 = closestPoint[0];
        // double y1 = closestPoint[1];
        // double x2 = other.position[0];
        // double y2 = other.position[1];

        double otherVelocityX = other.velocity[0];
        double otherVelocityY = other.velocity[1];
        
        double[] relativeVelocity = {otherVelocityX, otherVelocityY};

        double dot = other.dot(normal, relativeVelocity);

        // Multiplying by -dot here because y-values are flipped in the Java coordinate space (I think)
        normal[0] *= -dot;
        normal[1] *= -dot;

        double[] acceleration = {normal[0], normal[1]};
        
        // Multiplying by two for bouncy collisions
        other.acceleration[0] += acceleration[0] * 2;
        other.acceleration[1] += acceleration[1] * 2;
    }

    public double[] normal(RigidBodySphere other) {
        double[] closestPoint = this.findClosestPoint(other.position);
        double dist = this.dist(closestPoint, other.position);
        double[] normal = {(other.position[0] - closestPoint[0]) / dist, (other.position[1] - closestPoint[1]) / dist};
        return normal;
    }


    public boolean isColliding(RigidBodySphere other) {
        double[] closestPoint = this.findClosestPoint(other.position);
        double dist = this.dist(closestPoint, other.position);
        double radiiSum = other.radius;

        // Checking if the spheres are touching/inside one another (ayo)
        if (dist > radiiSum) {
            // System.out.println(closestPoint[0] + ", " + closestPoint[1]);
            return false;
        }

        
        double x1 = closestPoint[0];
        double y1 = closestPoint[1];
        double x2 = other.position[0];
        double y2 = other.position[1];

        double otherVelocityX = other.velocity[0];
        double otherVelocityY = other.velocity[1];
        
        double[] distanceVector = {x2 - x1, y2 - y1};

        // Accounting for velocity next frame by adding acceleration to prevent double collisions on corners because double collisions will greatly
        // increase the velocity of an object, more than intended
        double[] relativeVelocity = {otherVelocityX + other.acceleration[0], otherVelocityY + other.acceleration[1]};

        double dot = other.dot(distanceVector, relativeVelocity);
        // Not sure if this is needed/useful, it just checks if the distance between the two objects is decreasing (if other is moving into this)
        if (dot > 0) {
            return false;
        }

        
        return true;
    }


    public double dist(double[] p1, double[] p2) {
        return Math.sqrt(Math.pow(p2[0] - p1[0], 2) + Math.pow(p2[1] - p1[1], 2));
    }


    // Finds the closest point on this line to the given point p
    public double[] findClosestPoint(double[] p) {
        double x1 = position[0][0];
        double y1 = position[0][1];
        double x2 = position[1][0];
        double y2 = position[1][1];

        double m = (y2 - y1) / (x2 - x1); // Slope of this line
        double pm = -1 / m; // Perpendicular slope = negative reciprocal of slope

        // Accounting for infinite or zero slope
        if (m == 0) {
            double[] result = {p[0], y1};
            return result;
        } else if (!Double.isFinite(m)) {
            double[] result = {x1, p[1]};
            return result;
        }

        double b = -(m * x1) + y1; // yint of this line
        // double pb = -(pm * x1); // yint of perpendicular line

        double h = p[0];
        double k = p[1];

        // Equation of perpendicular line passing through point p is y = pm(x - h) + k
        // pm(x - h) + k = mx + b, solve to find intersection x-value
        double leftX = pm;
        double leftConstants = pm * -h + k;
        double rightX = m;
        double rightConstants = b;

        double totalX = rightX - leftX;
        double totalConstants = leftConstants - rightConstants;

        double x = totalConstants / totalX;
        double y = m * x + b;
    
        
        
        
        // Checking if the sphere is hitting where the line segment is, not where the infinite line is
        double[] leftright = {0, 0};
        double[] topbottom = {0, 0};
        
        if (this.position[0][0] <= this.position[1][0]) {
            leftright[0] = this.position[0][0];
            leftright[1] = this.position[1][0];
        } else {
            leftright[0] = this.position[1][0];
            leftright[1] = this.position[0][0];
        }
        
        if (this.position[0][1] <= this.position[1][1]) {
            topbottom[0] = this.position[0][1];
            topbottom[1] = this.position[1][1];
        } else {
            topbottom[0] = this.position[1][1];
            topbottom[1] = this.position[0][1];
        }
        
        double left = leftright[0];
        double right = leftright[1];
        double top = topbottom[0];
        double bottom = topbottom[1];
        
        // Checking if the closest point found is on this line segment, and if it isn't it sets the closest point to the closest endpoint instead
        if (x < left) x = left;
        if (x > right) x = right;
        if (y < top) y = top;
        if (y > bottom) y = bottom;


        double[] result = {x, y};
        return result;
    }
}
