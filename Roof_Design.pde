import java.util.Arrays;
import java.util.Comparator;
import java.util.Collections;

Graph base = new Graph();
float r = 60;
boolean center = true, PickedFocus=false, showBalls = false;
boolean showWall = true, showRoof = true, showFence = true;

void setup() {
  myFace = loadImage("./data/pic.jpg");
  RodFace = loadImage("./data/rod_pic.jpg");
  textureMode(NORMAL);  
  size(600, 600, P3D);
  noSmooth();
  base.declare();
  base.loadPts2D("./data/pts2D4");
}

float dz = -500; // distance to camera. Manipulated with wheel or when 
float rx = -0.06 * TWO_PI, ry = -0.04 * TWO_PI;    // view angles manipulated when space pressed but not mouse
float w = 600; // half-size of the cube containing the scene displayed

void draw() {
  background(255);
  displayHeader();
  hint(ENABLE_DEPTH_TEST);
  pushMatrix();
    camera();
    translate(width / 2, height / 2, dz); // puts origin of model at screen center and moves forward/away by dz
    lights();
    rotateX(rx); rotateY(ry); // rotates the model around the new origin (center of screen)
    rotateX(PI / 2); // rotates frame around X to make X and Y basis vectors parallel to the floor
    noStroke();
    //showFrame(50);
    fill(grey); pushMatrix(); translate(0, 0, - w / 2 - 1.5); box(w, w, 1); popMatrix(); // draws floor as thin plate
    noFill(); strokeWeight(1); stroke(black); showBlock(w, w, w, 0, 0, 0 ,0);
    noStroke(); 
    pushMatrix();
      translate(-width / 2, -height / 2, - w / 2 - 0.5); //the points is obtained in twoD with no translation, so we translate them back
      scale(600.0 / 800.0, 600.0 / 800.0, 600.0 / 800.0);
      doPick();
      base.createAllDarts();
      base.createDartLoopTree();
      //base.showBoundedArea();
      base.drawPoints();
      base.drawAllDarts(); //<>//
      if(showWall) base.drawWall();
      if(showFence) base.drawFence();
      if(showRoof) base.drawRoof();
      //base.testGetFacePolygon();
    popMatrix();
  popMatrix();
}