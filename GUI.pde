void keyPressed(){
  switch(key){
    case '-': r -= 1; break;
    case '+': r += 1; break;
    case 's': base.swingDart(); break;
    case 'u': base.unswingDart(); break;
    case 'o': base.oppositeDart(); break;
    case 'n': base.nextDart(); break;
    case 'p': base.preDart(); break;
    case 'r': showRoof = !showRoof; break;
    case 'w': showWall = !showWall; break;
    case 'f': showFence = !showFence; break;
    case 'E': base.empty();break;
    case 'L': base.loadPts("./data/pts3D");break;
    case 'S': base.savePts("./data/pts3D");break;
  }
}

void mousePressed(){
  //pt pickedPoint = pick( mouseX, mouseY ); 
  //base.pickClosest(pickedPoint);
  if(!keyPressed){
    if(mouseButton == LEFT){
      base.setPicked();
      base.setPickedDart(Of);
    }
  }
  else{
    switch(key){
      case 'a': base.addVert(Of); break;
       
    }
  }
}

void mouseMoved(){
  if(keyPressed && key == ' ') {rx -= PI * (mouseY - pmouseY) / height; ry += PI * (mouseX - pmouseX) / width;};
  if(keyPressed && key == 'z') dz += (float)(mouseY - pmouseY) * 5; // approach view (same as wheel)
}

void mouseDragged(){
  if(!keyPressed){
    if(mouseButton == LEFT) base.setPickedTo(Of);
    
  }
}

// function 
void displayR(){
  noFill();
  stroke(red);
  ellipse(0, 0, r, r);
  fill(red);
  text("r = " + r, r / 2 , r / 2);
}