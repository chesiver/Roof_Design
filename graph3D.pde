class Graph{
  int nv = 0, nvMax = 1000;
  int pv = 0, pp = 0, pd = 0;
  pt[] G = new pt[nvMax];
  int[] D = new int[nvMax];
  int nd = 0, ndMax = nvMax*2; 
  int[] V = new int[ndMax];
  int[] S = new int[ndMax];
  int[] C = new int[ndMax];
  final float o = 5;
  
  boolean isColoring = false;
  int nf = 0, nfMax = 100;
  int[] F = new int[nfMax];
  int[] Fdepth = new int[nfMax];
  IntList[] NestedFaces = new IntList[nfMax];
  
  Graph(){}
  void declare() { for(int i=0; i<nvMax; i++) G[i]=P(); for(int i = 0; i < nfMax; ++i) NestedFaces[i] = new IntList();}            // creates all points, MUST BE DONE AT INITALIZATION
  void empty() { nv = 0; nd = 0; pv = 0;}        // empties the graph
  
  int addVert(pt P) { G[nv].setTo(P); pv=nv; nv++; return nv-1; }    // adds a vertex at position P
  int addVert(float x, float y) {return addVert(P(x,y));}
  
  void savePts(String fn) 
  {
    String [] inppts = new String [nv+1];
    int s=0;
    inppts[s++]=str(nv);
    for (int i=0; i<nv; i++) {inppts[s++]=str(G[i].x)+","+str(G[i].y)+","+str(G[i].z);}
    saveStrings(fn,inppts);
  };
  
  void loadPts(String fn) 
  {
    println("loading: "+fn); 
    String [] ss = loadStrings(fn);
    int s=0;
    nv = int(ss[s++]); print("nv="+nv);
    for(int k=0; k<nv; k++) {int i=k+s; float [] xy = float(split(ss[i],",")); G[k].setTo(xy[0],xy[1],xy[2]);}
    pv=0;
  }; 
  void loadPts2D(String fn) 
  {
    println("loading: "+fn); 
    String [] ss = loadStrings(fn);
    int s=0;   int comma;   float x, y;
    nv = int(ss[s++]); print("nv="+nv);
    for(int k=0; k<nv; k++) {
      int i=k+s; 
      comma=ss[i].indexOf(',');   
      x=float(ss[i].substring(0, comma));
      y=float(ss[i].substring(comma+1, ss[i].length()));
      G[k].setTo(x,y, 0);
      };
    pv=0;
  };   
  //UI
  int setIdOfVertexWithClosestScreenProjectionTo(pt M){ // sets pp to the index of the vertex that projects closest to the mouse 
    pp=0; 
    for(int i = 1; i < nv; ++i) if(d(M, ToScreen(G[i])) <= d(M,ToScreen(G[pp]))) pp = i; 
    return pp;
  }
  void pickClosest(pt M) { // for picking a vertex with the mouse
    pv=0; 
    for (int i = 1; i < nv; ++i) if(d(M, G[i]) <= d(M, G[pv])) pv=i; 
  }
  void setPicked() {pv = pp;}
  void setPickedTo(pt Q) { G[pv].setTo(Q);}      // moves selected point (index p) by amount mouse moved recently
  void setPickedDart(pt Q){
    pd = 0;
    for(int i = 0; i < nd; ++i){
      if (D[V[pd]] < 0 || d(Q, G[V[i]]) < d(Q, G[V[pd]])) pd = i;
    }
  }
  
  //display
  void drawBalls(float r) {for(int v = 0; v < nv; ++v) show(G[v],r);}
  void showPicked() {show(G[pv], 13);}
  void showPicked(float r) {show(G[pv], r);}
  void drawPoints(){
    pen(blue, 8);
    beginShape(POINTS); 
      for(int v=0; v<nv; v++) v(G[v]); 
    endShape(); 
  }
  
//************************************************************************
//**** Function For Alpha-Shape and Dart Structure
//************************************************************************
  
  //create dart
  int insertSwingRing(int d, int v){
    //V[d] = v;
    if(D[v] < 0){
      D[v] = d;
      S[d] = - v - 1;
    }
    else if(S[D[v]] < 0){
      int p = D[v];
      S[p] = d; S[d] = - v - 1;
    }
    else{
      int p = D[v];
      int dz = d % 2 == 0 ? d + 1: d - 1;
      while(S[p] >= 0){
        int pz = p % 2 == 0 ? p + 1 : p - 1;
        int qz = S[p] % 2 == 0 ? S[p] + 1 : S[p] - 1;
        if(isClockWise(V(G[v], G[V[pz]]), V(G[v], G[V[dz]]), V(G[v], G[V[qz]]))){
          S[d] = S[p]; S[p] = d;
          break;
        }
        p = S[p];
      }
      if(S[p] < 0){
        S[d] = D[v]; D[v] = d;
      }
    }
    return d;
  }
  void createAllDarts(){
    for(int i = 0; i < nv; ++i){
      D[i] = - i - 1;
    }
    nd = 0;
    for(int i = 0; i < nv; ++i){
        for(int j = i + 1; j < nv; ++j){
          boolean flag1 = true, flag2 = true;
          for(int k = 0; k < nv; ++k){
            if(k == i || k == j) continue;
            if(isInDiskOnRight(G[k], G[i], G[j], r)) flag1 = false;
            if(isInDiskOnRight(G[k], G[j], G[i], r)) flag2 = false;
          }
          if(flag1 && flag2){
            V[nd] = i; V[nd + 1] = j;
            insertSwingRing(nd, i);
            insertSwingRing(nd + 1, j);
            C[nd] = 0; C[nd + 1] = 0;
            nd += 2;
          }
          else if(flag1){
            V[nd] = i; V[nd + 1] = j;
            insertSwingRing(nd, i);
            insertSwingRing(nd + 1, j);
            C[nd] = 0; C[nd + 1] = 1;
            nd += 2;
          }
          else if(flag2){
            V[nd] = i; V[nd + 1] = j;
            insertSwingRing(nd, i);
            insertSwingRing(nd + 1, j);
            C[nd] = 1; C[nd + 1] = 0;
            nd += 2;
          }
          else{
            
          }
        }
      }
  }
  
  //dart semantic
  int s(int d){
    if(S[d] < 0) return D[- S[d] - 1];
    else return S[d];
  }
  void swingDart(){
    pd = s(pd);
  }
  int u(int d){
    int p = d;
    while(s(p) != d) p = s(p);
    return p;
  }
  void unswingDart(){
    pd = u(pd);
  }
  int o(int d){
    if(d % 2 == 0) return d + 1;
    else return d - 1;
  }
  void oppositeDart(){
    pd = o(pd);
    //println("pd: " + pd);
  }
  int n(int d){
    return u(o(d));
  }
  void nextDart(){
    pd = n(pd);
  }
  int p(int d){
    return o(s(d));
  }
  void preDart(){
    pd = p(pd);
  }
  //green loop
  boolean isRayIntersectDart(int d, int d2){
    int v1 = V[d], v2 = V[o(d)];
    int v3 = V[d2], v4 = V[o(d2)];
    if(det2(V(G[v1], G[v2]), V(G[v1], G[v3])) * det2(V(G[v1], G[v2]), V(G[v1], G[v4])) < 0){
      vec normal = U(R(V(G[v3], G[v4])));
      if(dot(V(G[v1], G[v3]), normal) > 0) normal = M(normal);
      if(dot(V(G[v1], G[v2]), normal) < 0) return true;
    }
    return false;
  }
  int nRayIntersectOfDartLoop(int d, int d2){
    //println("d " + d + " d2 " + d2);
    int sum = 0;
    int p = d2;
    do{
      if(isRayIntersectDart(d, p)) ++sum;
      p = n(p);
    }while(p != d2);
    //println("d: " + d + " " + "d2 " + d2 + " " + "sum: " + sum);
    return sum;
  }
  boolean isDartLoopInLDartloop(int d, int d2){
    return nRayIntersectOfDartLoop(d, d2) % 2 == 1;
  }
  void getLoopDart(){
    nf = 0;
    boolean[] mark = new boolean[nd];
    for(int i = 0; i < nd; ++i) mark[i] = false;
    while(true){
      boolean flag = false;
      for(int i = 0; i < nd; ++i){
        if(C[i] == 0) continue;
        if(!mark[i]){
          mark[i] = true;
          flag = true;
          F[nf++] = i;
          int p = n(i);
          while(p != i){
            mark[p] = true;
            p = n(p);
          }
          break;
        }
      }
      if(!flag) break;
    }
  }
  void dfSearch(int[][] Map, boolean[] isVisited,int i, IntList L){
    if(!isVisited[i]){
      isVisited[i] = true;
      for(int j = 0; j < Map.length; ++j){
        if(Map[i][j] == 1) dfSearch(Map, isVisited, j, L);
      }
      L.append(i);
    }
  }
  void createDartLoopTree(){
    getLoopDart();
    int[][] Map = new int[nf][nf];
    int[] Out = new int[nf];
    for(int i = 0; i < nf ; ++i){
      for(int j = i + 1; j < nf; ++j){
        if(isDartLoopInLDartloop(F[j], F[i])){
          Map[j][i] = 1;
          ++Out[i];
        }
        else if(isDartLoopInLDartloop(F[i], F[j])){
          Map[i][j] = 1;
          ++Out[j];
        }
      }
    }
    boolean[] isVisited = new boolean[nf];
    IntList L = new IntList();
    for(int i = 0; i < nf; ++i){
      if(Out[i] == 0){
        dfSearch(Map, isVisited, i, L);
      }
    }
    for(int i = 0; i < nf; ++i){
      NestedFaces[i].clear();
    }
    for(int i = 0; i < L.size(); ++i){
      Fdepth[i] = 0;
      for(int j = i - 1; j >= 0; --j){
        if(Map[i][j] == 1){
          NestedFaces[j].append(i);
          Fdepth[i] = Fdepth[j] + 1;
          break;
        }
      }
      
    }
    for(int i = 0; i < nf; ++i){
      if(NestedFaces[i].size() != 0){
        //println("NestedFaces " + i + ": " + join(nf(NestedFaces[i].array(), 0), ", "));
      }
    }
  }
  
  //green loop test
  void showBoundedArea(){
    //boolean[] mark = new boolean[nf];
    for(int i = 0; i < nv; ++i){
      G[i].z = 0;
    }
    for(int i = 0; i < nf; ++i){
      if(Fdepth[i] % 2 == 0) fill(yellow);
      else fill(white);
      int d = F[i], p = F[i];
      beginShape();
        do{
          vertex(G[V[p]]);
          p = n(p);
        }while(p != d);
      endShape(CLOSE);
    }
  }
  
  //display dart
  void drawArrow(pt P, pt Q, color c){
    stroke(c);
    pen(c, 1); fill(c);
    vec PQuo = R(U(P, Q));
    pt A = P(P, o, PQuo);
    pt B = P(Q, o, PQuo);
    pt S = L(A, 0.05, B);
    pt E = L(A, 0.95, B);
    arrow(S, E);
  }
  void drawAllDarts(){
    for(int i = 0; i < nd; i += 2){
      if(C[i] == 0) drawArrow(G[V[i]], G[V[i + 1]], red);
      else if(C[i] == 1) drawArrow(G[V[i]], G[V[i + 1]],green);
      if(C[i + 1] == 0) drawArrow(G[V[i + 1]], G[V[i]], red);
      else if(C[i + 1] == 1) drawArrow(G[V[i + 1]], G[V[i]],green);
    }
    strokeWeight(20);
    //println("pd: " + pd + " v: " + V[pd] + " D[v] " + D[V[pd]]);
    drawArrow(G[V[pd]], G[V[o(pd)]], black);
  }

//************************************************************************
//**** Function For Fences and Walls
//************************************************************************

  pt Centroid(){
    pt C = P(); // will collect sum of points before division
    for (int i = 0; i < nv; ++i) C.add(G[i]); 
    return P(1.0 / nv, C); // returns divided sum
  }
  
  float wallHeight = 100;
  vec wallVector = V(0, 0, wallHeight);
  void drawWall(){
    noStroke(); fill(cyan);
    for(int i = 0; i < nf; ++i){
      //if(Fdepth[i] % 2 == 0){
        int sd = F[i], p = sd;
        do{
          beginShape();
            v(G[V[p]]);
            v(G[V[o(p)]]);
            v(P(G[V[o(p)]], wallVector));
            v(P(G[V[p]], wallVector));
          endShape(CLOSE);
          p = n(p);
        }while(p != sd);
      //}
    }
  }
  float fenceHeight = 50;
  vec fenceVector = V(0, 0, fenceHeight);
  void drawFence(){
    noStroke(); fill(magenta);
    for(int i = 0; i < nd; ++i){
      if(C[i] == 0 && C[o(i)] == 0){
        beginShape();
            v(G[V[i]]);
            v(G[V[o(i)]]);
            v(P(G[V[o(i)]], fenceVector));
            v(P(G[V[i]], fenceVector));
        endShape(CLOSE);
      }
    }
  }
  
//**************************************************************************
// MEDIAL AXIS TRANSFORM
//**************************************************************************
  pt[] getFacePolyloop(int iFace){
    ArrayList<pt> polyloop = new ArrayList<pt>();
    int sd = F[iFace], p = sd;
    do{
      polyloop.add(G[V[p]]);
      p = n(p);
    }while(p != sd);
    return polyloop.toArray(new pt[0]);
  }
  void testGetFacePolygon(){
    pt[] polyloop = getFacePolyloop(0);
    for(int i = 0; i < polyloop.length / 2; ++i){
      show(polyloop[i], V(0, 0, 100));
    }
  }
  
  vec bisector(pt A, pt B, pt C, pt D) {
    vec bisector = V(0, 0, 0); // initialize bisector
    vec uAB = U(A,B);
    vec uCD = U(C,D);
    // if vectors are oriented in opposite directions
    bisector = V((uCD.x-uAB.x)*0.5,(uCD.y-uAB.y)*0.5, 0);
    return bisector;
  }
  // intersection between two lines
  pt intersection(pt A, pt B, pt C, pt D) {
    float xa=A.x, ya=A.y, xb=B.x, yb = B.y;
    float xc=C.x, yc=C.y, xd=D.x, yd = D.y;
    float m1=(yb-ya)/(xb-xa);
    float m2=(yd-yc)/(xd-xc);
    float xi=((yc-ya)+(m1*xa - m2*xc))/(m1-m2);
    float yi = (m1*(xi-xa))+ya;
    return P(xi, yi, 0);
  }
  // circle tangent to two edges from a point P tangent to edge AB
  pt [] tangentArc(pt T1, pt A, pt B, pt C, pt D) {
    pt[] points= new pt[3];
    vec normalAB = R(U(A,B));
    vec normalCD = R(U(C,D));
    pt I = intersection(A, B, C, D);
    vec BISECT = bisector(A, B, C, D);
    pt O = intersection(I, P(I.x+BISECT.x,I.y+BISECT.y), T1, P(T1.x+normalAB.x,T1.y+normalAB.y));
    pt T2 = intersection(O,P(O.x-normalCD.x,O.y-normalCD.y),C,D);
    points[0] = T1;
    points[1] = O;
    points[2] = T2;
    return points;
  }
  // circle at concave vertex and tangent to another edge
  pt [] concaveArc(pt CONCAVE, pt A, pt B, pt C, pt D, pt E, pt F, float n) {
    pt [] points= new pt [3];
    vec normalEF = R(U(E,F));
    vec BISECT = bisector(A, B, C, D);
    pt T2 = L(E, n, F);
    pt O = intersection(CONCAVE, P(CONCAVE.x-BISECT.x,CONCAVE.y-BISECT.y), T2, P(T2.x+normalEF.x,T2.y+normalEF.y));
    points[0] = CONCAVE;
    points[1] = O;
    points[2] = T2;
    return points;
  }
  
  int next(pt[] polyloop, int i){ return i == polyloop.length - 1 ? 0 : i + 1; }
  int pre(pt[] polyloop, int i){ return i == 0 ? polyloop.length - 1 : i - 1; }
  boolean isOnBisectorOnIndex(pt[] polyloop, int i, pt S){
    final float epsilon=0.005;
    vec V1 = U(polyloop[i], polyloop[next(polyloop, i)]), V2 = U(polyloop[i], polyloop[pre(polyloop, i)]);
    vec V3 = U(polyloop[i], S);
    if(abs(det2(V1, V3) - det2(V3, V2)) < epsilon && dot(V3, V(0.5, V1, 0.5, V2)) > 0) return true;
    else return false;
  }
  
  pt[] medialAxisTransform(pt[] polyloop) {
    pt[] roofPoints = new pt[polyloop.length];
    for(int i = 0; i < polyloop.length; ++i) {
      float farthest_dist = MIN_FLOAT; pt f = P();
      for(int j = 0; j < polyloop.length; ++j) {
          if(j == i) continue;
          // concave vertices
          if (j == next(polyloop, i) && det2(V(polyloop[i], polyloop[next(polyloop, i)]),V(polyloop[j], polyloop[next(polyloop, j)]))<0) {
            for (int k = 0; k < polyloop.length; ++k) {
              if(k == i || k == j) continue;
              for (int n = 1; n < 50; ++n){
                pt [] trialArc = concaveArc(polyloop[next(polyloop, i)], polyloop[i], polyloop[next(polyloop, i)], 
                  polyloop[j], polyloop[next(polyloop, j)], polyloop[k], polyloop[next(polyloop, k)], 0.02 * n);
                float radius = d(trialArc[0], trialArc[1]);
                // check whether point is inside the domain
                if (det2(V(polyloop[k], polyloop[next(polyloop, k)]), V(L(polyloop[k], 0.02 * n, polyloop[next(polyloop, k)]),trialArc[1])) < 0) continue;
                if (det2(V(polyloop[i], polyloop[next(polyloop, i)]), V(L(polyloop[i], 0.02 * n, polyloop[next(polyloop, i)]), trialArc[1])) < 0) continue;
                // check if circle does not intercept the polyloop
                boolean doesNotintercept=true;
                for (int p = 0; p < polyloop.length; ++p) {
                  if(p == i || p == j) continue;
                  // check if any vertex is in the interior of the circle
                  if (d(polyloop[p],trialArc[1]) < radius) doesNotintercept=false;
                  // check if circle does not intercept any edges
                  if (projectsBetween(trialArc[1], polyloop[p], polyloop[next(polyloop, p)]) 
                    && disToLine(trialArc[1], polyloop[p], polyloop[next(polyloop, p)]) < radius) doesNotintercept=false;
                }
                if (doesNotintercept){
                  if(farthest_dist < d(trialArc[1], polyloop[i])){
                    farthest_dist = d(trialArc[1], polyloop[i]);
                    f = trialArc[1];
                  }
                }
              }
            }
          }
          else{
            for (int n = 1; n < 50; ++n){
              // obtain tangent arc
              pt [] trialArc = tangentArc(L(polyloop[i], 0.02 * n, polyloop[next(polyloop, i)]), polyloop[i], polyloop[next(polyloop, i)], polyloop[j], polyloop[next(polyloop, j)]);
              float radius = d(trialArc[0], trialArc[1]);
              // check whether point is inside the domain ************TEMPORARY***********
              if(det2(V(polyloop[i], polyloop[next(polyloop, i)]), V(L(polyloop[i], 0.05 * n, polyloop[next(polyloop, i)]), trialArc[1])) < 0) continue;
              // check whether tangent point is on the edge bounds
              boolean notInBounds = true;
              if(min(polyloop[j].x, polyloop[next(polyloop, j)].x) <= trialArc[2].x && trialArc[2].x <= max(polyloop[j].x, polyloop[next(polyloop, j)].x))  notInBounds=false;
              if(notInBounds)  continue;
              // check if circle does not intercept the polyloop
              boolean doesNotintercept = true;
              for(int k = 0; k < polyloop.length; ++k){
                if(k == i || k == j) continue;
                  // check if any vertex is in the interior of the circle
                if(d(polyloop[k], trialArc[1]) < radius) doesNotintercept=false;
                // check if circle does not intercept any edges
                if(projectsBetween(trialArc[1], polyloop[k], polyloop[next(polyloop, k)]) && 
                  disToLine(trialArc[1], polyloop[k], polyloop[next(polyloop, k)]) < radius) doesNotintercept = false;
              }
              if (doesNotintercept) {
                if(isOnBisectorOnIndex(polyloop, i, trialArc[1])){
                  if(farthest_dist < d(trialArc[1], polyloop[i])){
                    farthest_dist = d(trialArc[1], polyloop[i]);
                    f = trialArc[1];
                  }
                }
              }
            }
          }
      }
      roofPoints[i] = f;
    }
    return roofPoints;
  }
  
  pt[] getMergeRoofPoints(pt[] polyloop, pt[] roofPoints){
    ArrayList<pt> mergedRoofPoints = new ArrayList<pt>();
    int[] mergeMap = new int[polyloop.length];
    float eps = 20.0;
    for(int i = 0; i < polyloop.length; ++i) mergeMap[i] = i;
    for(int i = 0; i < polyloop.length - 1; ++i){
      if(mergeMap[i] != i) continue;
      if(roofPoints[i].x == 0 && roofPoints[i].y == 0) continue;
      for(int j = i + 1; j < polyloop.length; ++j){
        if(d(roofPoints[i], roofPoints[j]) < eps){
          mergeMap[j] = i;
        }
      }
    }
    for(int i = 0; i < polyloop.length; ++i){
      if(roofPoints[i].x == 0 && roofPoints[i].y == 0) continue;
      if(mergeMap[i] == i){
        mergedRoofPoints.add(roofPoints[i]);
      }
    }
    return mergedRoofPoints.toArray(new pt[0]);
  }
  
  class mergedRoofPoints_comparator implements Comparator<pt>{
    pt centroid = new pt();
    public mergedRoofPoints_comparator(pt _centroid){
      centroid.x = _centroid.x;
      centroid.y = _centroid.y;
    }
    @Override
    public int compare(pt a, pt b){
      float ax = a.x - this.centroid.x, ay = a.y - this.centroid.y;
      float bx = b.x - this.centroid.x, by = b.y - this.centroid.y;
      float sa = atan2(ay, ax), sb = atan2(by, bx);
      if(sa < sb) return -1;
      else if(sa == sb) return 0;
      else return 1;
    }
  }
  void reOrderMergedRoofPoints(pt[] mergedRoofPoints){
    pt centroid = P();
    for(int i = 0; i < mergedRoofPoints.length; ++i){
      centroid.x += mergedRoofPoints[i].x;
      centroid.y += mergedRoofPoints[i].y;
    }
    centroid.x /= mergedRoofPoints.length; centroid.y /= mergedRoofPoints.length;
    Arrays.sort(mergedRoofPoints, new mergedRoofPoints_comparator(centroid));
  }
  
  int[] findClosestRoofPoints(pt[] polyloop, pt[] mergedRoofPoints){
    int[] closestRoofPointsMap = new int[polyloop.length];
    for(int i = 0; i < polyloop.length; ++i){
      closestRoofPointsMap[i] = 0;
      float closesDist = d(polyloop[i], mergedRoofPoints[0]);
      for(int j = 1; j < mergedRoofPoints.length; ++j){
        if(d(polyloop[i], mergedRoofPoints[j]) < closesDist){
          closestRoofPointsMap[i] = j;
          closesDist = d(mergedRoofPoints[j], polyloop[i]);
        }
      }
    }
    for(int i = 0; i < polyloop.length; ++i){
      //println("closestRoofPoints " + i + " " + closestRoofPointsMap[i]);
    }
    return closestRoofPointsMap;
  }
  
  void drawRoof(int iFace){
    pt[] polyloop = getFacePolyloop(iFace);
    //if(area(polyloop) < 0) Collections.reverse(Arrays.asList(polyloop));
    pt[] roofPoints = medialAxisTransform(polyloop);
    pt[] mergedRoofPoints = getMergeRoofPoints(polyloop, roofPoints); //<>//
    pen(blue, 5); noFill();
    for(int i = 0; i < mergedRoofPoints.length; ++i){
      show(mergedRoofPoints[i], 10);
    }
    reOrderMergedRoofPoints(mergedRoofPoints);
    int[] closestRoofPointsMap = findClosestRoofPoints(polyloop, mergedRoofPoints);

    noStroke(); fill(metal);
    vec roofUpVector = A(V(0, 0, 100), wallVector);
    for(int i = 0; i < polyloop.length; ++i){
      int p = closestRoofPointsMap[i];
      int q = closestRoofPointsMap[next(polyloop, i)];
      beginShape();
        v(P(polyloop[i], wallVector));
        v(P(mergedRoofPoints[p], roofUpVector));
        v(P(mergedRoofPoints[q], roofUpVector));
        v(P(polyloop[next(polyloop, i)], wallVector));
      endShape();
    }
  }
  void drawRingRoof(int iFace){
    ArrayList<pt> polyloopList = new ArrayList<pt>();
    pt[] outerPolyloop = getFacePolyloop(iFace);
    polyloopList.addAll(Arrays.asList(outerPolyloop));
    for(int i = 0; i < NestedFaces[iFace].size(); ++i){
      pt[] innerPolyloop = getFacePolyloop(NestedFaces[iFace].get(i));
      if(area(innerPolyloop) > 0) Collections.reverse(Arrays.asList(innerPolyloop));
      polyloopList.addAll(Arrays.asList(innerPolyloop));
    }
    pt[] polyloop = polyloopList.toArray(new pt[0]);
    //if(area(polyloop) < 0) Collections.reverse(Arrays.asList(polyloop));
    pt[] roofPoints = medialAxisTransform(polyloop);
    pt[] mergedRoofPoints = getMergeRoofPoints(polyloop, roofPoints);
    pen(blue, 5); noFill();
    for(int i = 0; i < mergedRoofPoints.length; ++i){
      show(mergedRoofPoints[i], 10);
    }
    reOrderMergedRoofPoints(mergedRoofPoints);
    int[] closestRoofPointsMap = findClosestRoofPoints(polyloop, mergedRoofPoints);

    noStroke(); fill(metal);
    vec roofUpVector = A(V(0, 0, 100), wallVector);
    for(int i = 0; i < polyloop.length; ++i){
      int p = closestRoofPointsMap[i];
      int q = closestRoofPointsMap[next(polyloop, i)];
      beginShape();
        v(P(polyloop[i], wallVector));
        v(P(mergedRoofPoints[p], roofUpVector));
        v(P(mergedRoofPoints[q], roofUpVector));
        v(P(polyloop[next(polyloop, i)], wallVector));
      endShape();
    }
  }
  void drawRoof(){
    for(int i = 0; i < nf; ++i){
      //println("i: " + i + " " + "Fdepth: " + Fdepth[i] + " " + "NestedFaces[i].size() " + NestedFaces[i].size());
      if(Fdepth[i] % 2 == 0){
        if(NestedFaces[i].size() == 0) drawRoof(i);
        else drawRingRoof(i);
      }
    }
  }
}