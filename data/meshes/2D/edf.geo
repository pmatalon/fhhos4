N = 32;
domainWidth = 32;
domainHeight = 15;
h = domainWidth/N;

//---------------//
// Bottom stripe //
//---------------//
stripeHeight = 5;
Point(1) = {          0,            0, 0, h}; // bottom-left corner
Point(2) = {domainWidth,            0, 0, h}; // bottom-right corner
Point(3) = {domainWidth, stripeHeight, 0, h}; // top-right corner
Point(4) = {          0, stripeHeight, 0, h}; // top-left corner

Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

Line Loop(100) = {1, 2, 3, 4};
Plane Surface(200) = {100}; // Loop Ids
Physical Surface("bottomStripe") = {200};

//---------------//
// Top rectangle //
//---------------//
rectangleOriginX = 4;
rectangleOriginY = 10;
rectangleWidth = 24;
Point(5) = {rectangleOriginX,                rectangleOriginY, 0, h}; // bottom-left corner
Point(6) = {rectangleOriginX+rectangleWidth, rectangleOriginY, 0, h}; // bottom-right corner
Point(7) = {rectangleOriginX+rectangleWidth, domainHeight,     0, h}; // top-right corner
Point(8) = {rectangleOriginX,                domainHeight,     0, h}; // top-left corner

Line(5) = {5, 6}; // bottom
Line(6) = {6, 7}; // right
Line(7) = {7, 8}; // top
Line(8) = {8, 5}; // left

Line Loop(101) = {5, 6, 7, 8};
Plane Surface(201) = {101}; // Loop Ids
Physical Surface("topRectangle") = {201};

//-------------------//
// Left little piece //
//-------------------//
pieceHeight = 1;
pieceWidth = 3;
pieceOriginX = 12;
pieceOriginY = rectangleOriginY-pieceHeight;
Point(9)  = {pieceOriginX,            pieceOriginY,             0, h}; // bottom-left corner
Point(10) = {pieceOriginX+pieceWidth, pieceOriginY,             0, h}; // bottom-right corner
Point(11) = {pieceOriginX+pieceWidth, pieceOriginY+pieceHeight, 0, h}; // top-right corner
Point(12) = {pieceOriginX,            pieceOriginY+pieceHeight, 0, h}; // top-left corner

Line(9)  = { 9, 10}; // bottom
Line(10) = {10, 11}; // right
Line(11) = {11, 12}; // top
Line(12) = {12,  9}; // left

Line Loop(102) = {9, 10, 11, 12};
Plane Surface(202) = {102}; // Loop Ids
Physical Surface("leftLittlePiece") = {202};

//--------------------//
// Right little piece //
//--------------------//
pieceOriginX += 5;
Point(13) = {pieceOriginX,            pieceOriginY,             0, h}; // bottom-left corner
Point(14) = {pieceOriginX+pieceWidth, pieceOriginY,             0, h}; // bottom-right corner
Point(15) = {pieceOriginX+pieceWidth, pieceOriginY+pieceHeight, 0, h}; // top-right corner
Point(16) = {pieceOriginX,            pieceOriginY+pieceHeight, 0, h}; // top-left corner

Line(13) = {13, 14}; // bottom
Line(14) = {14, 15}; // right
Line(15) = {15, 16}; // top
Line(16) = {16, 13}; // left

Line Loop(103) = {13, 14, 15, 16};
Plane Surface(203) = {103}; // Loop Ids
Physical Surface("rightLittlePiece") = {203};

//---------------------------//
// Weird shape in the middle //
//---------------------------//
Point(17) = {          0, domainHeight, 0, h}; // top-left corner of the domain
Point(18) = {domainWidth, domainHeight, 0, h}; // top-right corner of the domain

Line(17) = {17,  4}; // left
Line(19) = { 3, 18}; // right     {bottom pt, top pt}
Line(20) = { 8, 17}; // top-left  {right pt, left pt}
Line(21) = {18,  7}; // top-right {right pt, left pt}
Line(22) = { 6, 15}; // bottom-right of the top rectangle {right pt, left pt}
Line(23) = {16, 11}; // between the little pieces {right pt, left pt}
Line(24) = {12,  5}; // bottom-left of the top rectangle {right pt, left pt}

Line Loop(104) = {-3, 19, 21, -6, 22, -14, -13, -16, 23, -10, -9, -12, 24, -8, 20, 17};
Plane Surface(204) = {104}; // Loop Ids
Physical Surface("weirdShapeInTheMiddle") = {204};



Physical Curve("externalBoundary") = {1, 2, 19, 21, 7, 20, 17, 4}; // Line Ids
