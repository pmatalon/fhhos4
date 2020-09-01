N = 32;
Width = 24;
Height = 16;
h = Width/N;

Point(1) = {0,          0, 0, h}; // bottom-left corner
Point(2) = {Width,      0, 0, h}; // bottom-right corner
Point(3) = {Width, Height, 0, h}; // top-right corner
Point(4) = {0,     Height, 0, h}; // top-left corner

Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

Line Loop(1) = {1, 2, 3, 4};

s = 2; // size of the hole
holeId = 2;

xShift = 1;
yShift = 1;

nHolesPerRow = 6;
nHolesPerCol = 4;

For j In {1:nHolesPerCol}
	For i In {1:nHolesPerRow}
		// square hole
		holeId = 10*((j-1)*nHolesPerRow + i);
		
		x = xShift + (i-1)*2*s;
		y = yShift + (j-1)*2*s;
		
		bottomLeftCorner  = holeId;
		bottomRightCorner = holeId+1;
		topRightCorner    = holeId+2;
		topLeftCorner     = holeId+3;
		
		Point(bottomLeftCorner)  = {x,     y,     0, h}; // bottom-left corner
		Point(bottomRightCorner) = {x + s, y,     0, h}; // bottom-right corner
		Point(topRightCorner)    = {x + s, y + s, 0, h}; // top-right corner
		Point(topLeftCorner)     = {x    , y + s, 0, h}; // top-left corner
		
		bottom = holeId+4;
		right  = holeId+5;
		top    = holeId+6;
		left   = holeId+7;
		
		Line(bottom) = {bottomLeftCorner, bottomRightCorner}; // bottom
		Line(right)  = {bottomRightCorner, topRightCorner};   // right
		Line(top)    = {topRightCorner, topLeftCorner};       // top
		Line(left)   = {topLeftCorner, bottomLeftCorner};     // left
		
		Line Loop(holeId) = {bottom, right, top, left};
	EndFor
EndFor

Plane Surface(100) = {1, 10,20,30,40,50,60, 70,80,90,100,110,120, 130,140,150,160,170,180, 190,200,210,220,230,240};

Physical Line("externalBoundary") = {1};
Physical Line("holes")            = {10,20,30,40,50,60, 70,80,90,100,110,120, 130,140,150,160,170,180, 190,200,210,220,230,240};

Physical Surface("domain") = {100};
