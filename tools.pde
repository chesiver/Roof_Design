// ************************************ IMAGES & VIDEO 
int pictureCounter=0, frameCounter=0;
Boolean filming=false, change=false;
PImage myFace, RodFace; // picture of author's face, should be: data/pic.jpg in sketch folder
void snapPicture() {saveFrame("PICTURES/P"+nf(pictureCounter++,3)+".jpg"); }

// ******************************************COLORS 
color black=#000000, white=#FFFFFF, // set more colors using Menu >  Tools > Color Selector
   red=#FF0000, green=#00FF01, blue=#0300FF, yellow=#FEFF00, cyan=#00FDFF, magenta=#FF00FB,
   grey=#818181, orange=#FFA600, brown=#B46005, metal=#B5CCDE, dgreen=#157901;
void pen(color c, float w) {stroke(c); strokeWeight(w);}

// ******************************** TEXT , TITLE, and USER's GUIDE
void scribeHeader(String S, int i) { text(S,10,20+i*20); noFill();} // writes black at line i
void scribeHeaderRight(String S) {fill(0); text(S,width - 100, 20); noFill();} // writes black on screen top, right-aligned
void scribeHeaderRightDown(String S) {fill(0); text(S,width - 250, 20); noFill();} // 

void displayHeader(){
  scribeHeaderRight("Yidong Liu"); 
  scribeHeaderRightDown("Rodrigo Borela Valente");
  image(myFace, width-100, 25, myFace.width / 4, myFace.height / 4); 
  image(RodFace, width-200, 25, RodFace.width / 2, RodFace.height / 2); 
}