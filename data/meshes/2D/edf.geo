N = 32;
domainWidth = 32;
domainHeight = 15;
h = domainWidth/N;

//---------------//
// Bottom stripe //
//---------------//
stripeHeight = 5;

stripeBL = 1;
stripeBR = 2;
stripeTR = 3;
stripeTL = 4;
Point(stripeBL) = {          0,            0, 0, h}; // bottom-left corner
Point(stripeBR) = {domainWidth,            0, 0, h}; // bottom-right corner
Point(stripeTR) = {domainWidth, stripeHeight, 0, h}; // top-right corner
Point(stripeTL) = {          0, stripeHeight, 0, h}; // top-left corner

stripeB = 1;
stripeR = 2;
stripeT = 3;
stripeL = 4;
Line(stripeB) = {stripeBL, stripeBR}; // bottom
Line(stripeR) = {stripeBR, stripeTR}; // right
Line(stripeT) = {stripeTR, stripeTL}; // top
Line(stripeL) = {stripeTL, stripeBL}; // left

Line Loop(100) = {stripeB, stripeR, stripeT, stripeL}; // Line Ids
Plane Surface(200) = {100}; // Loop Ids
Physical Surface("bottomStripe") = {200};

//---------------//
// Top rectangle //
//---------------//
rectangleOriginX = 4;
rectangleOriginY = 10;
rectangleWidth = 24;

rectangleBL = 5;
rectangleBR = 6;
rectangleTR = 7;
rectangleTL = 8;
Point(rectangleBL) = {rectangleOriginX,                rectangleOriginY, 0, h}; // bottom-left corner
Point(rectangleBR) = {rectangleOriginX+rectangleWidth, rectangleOriginY, 0, h}; // bottom-right corner
Point(rectangleTR) = {rectangleOriginX+rectangleWidth, domainHeight,     0, h}; // top-right corner
Point(rectangleTL) = {rectangleOriginX,                domainHeight,     0, h}; // top-left corner

rectangleB = 5;
rectangleR = 6;
rectangleT = 7;
rectangleL = 8;
Line(rectangleB) = {rectangleBL, rectangleBR}; // bottom
Line(rectangleR) = {rectangleBR, rectangleTR}; // right
Line(rectangleT) = {rectangleTR, rectangleTL}; // top
Line(rectangleL) = {rectangleTL, rectangleBL}; // left

Line Loop(101) = {rectangleB, rectangleR, rectangleT, rectangleL}; // Line Ids
Plane Surface(201) = {101}; // Loop Ids
Physical Surface("topRectangle") = {201};

//-------------------//
// Left little piece //
//-------------------//
pieceHeight = 1;
pieceWidth = 3;
pieceOriginX = 12;
pieceOriginY = rectangleOriginY-pieceHeight;

lpieceBL =  9;
lpieceBR = 10;
lpieceTR = 11;
lpieceTL = 12;
Point(lpieceBL) = {pieceOriginX,            pieceOriginY,             0, h}; // bottom-left corner
Point(lpieceBR) = {pieceOriginX+pieceWidth, pieceOriginY,             0, h}; // bottom-right corner
Point(lpieceTR) = {pieceOriginX+pieceWidth, pieceOriginY+pieceHeight, 0, h}; // top-right corner
Point(lpieceTL) = {pieceOriginX,            pieceOriginY+pieceHeight, 0, h}; // top-left corner

lpieceB =  9;
lpieceR = 10;
lpieceT = 11;
lpieceL = 12;
Line(lpieceB) = {lpieceBL, lpieceBR}; // bottom
Line(lpieceR) = {lpieceBR, lpieceTR}; // right
Line(lpieceT) = {lpieceTR, lpieceTL}; // top
Line(lpieceL) = {lpieceTL, lpieceBL}; // left

Line Loop(102) = {lpieceB, lpieceR, lpieceT, lpieceL}; // Line Ids
Plane Surface(202) = {102}; // Loop Ids
Physical Surface("leftLittlePiece") = {202};

//--------------------//
// Right little piece //
//--------------------//
pieceOriginX += 5;

rpieceBL = 13;
rpieceBR = 14;
rpieceTR = 15;
rpieceTL = 16;
Point(rpieceBL) = {pieceOriginX,            pieceOriginY,             0, h}; // bottom-left corner
Point(rpieceBR) = {pieceOriginX+pieceWidth, pieceOriginY,             0, h}; // bottom-right corner
Point(rpieceTR) = {pieceOriginX+pieceWidth, pieceOriginY+pieceHeight, 0, h}; // top-right corner
Point(rpieceTL) = {pieceOriginX,            pieceOriginY+pieceHeight, 0, h}; // top-left corner

rpieceB = 13;
rpieceR = 14;
rpieceT = 15;
rpieceL = 16;
Line(rpieceB) = {rpieceBL, rpieceBR}; // bottom
Line(rpieceR) = {rpieceBR, rpieceTR}; // right
Line(rpieceT) = {rpieceTR, rpieceTL}; // top
Line(rpieceL) = {rpieceTL, rpieceBL}; // left

Line Loop(103) = {rpieceB, rpieceR, rpieceT, rpieceL}; // Line Ids
Plane Surface(203) = {103}; // Loop Ids
Physical Surface("rightLittlePiece") = {203};

//---------------------------//
// Weird shape in the middle //
//---------------------------//
domainTL = 17;
domainTR = 18;
Point(domainTL) = {          0, domainHeight, 0, h}; // top-left corner of the domain
Point(domainTR) = {domainWidth, domainHeight, 0, h}; // top-right corner of the domain

verticalTL = 17;
verticalTR = 19;
horizontalTL = 20;
horizontalTR = 21;
rectangleHorizBR = 22;
edgeBetweenPieces = 23;
rectangleHorizBL = 24;
Line(verticalTL)        = {domainTL,    stripeTL};    // left      {top pt, bottom pt}
Line(verticalTR)        = {stripeTR,    domainTR};    // right     {bottom pt, top pt}
Line(horizontalTL)      = {rectangleTL, domainTL};    // top-left  {right pt, left pt}
Line(horizontalTR)      = {domainTR,    rectangleTR}; // top-right {right pt, left pt}
Line(rectangleHorizBR)  = {rectangleBR, rpieceTR};    // bottom-right of the top rectangle {right pt, left pt}
Line(edgeBetweenPieces) = {rpieceTL,    lpieceTR};    // between the little pieces {right pt, left pt}
Line(rectangleHorizBL)  = {lpieceTL,    rectangleBL}; // bottom-left of the top rectangle {right pt, left pt}

Line Loop(104) = {-stripeT, verticalTR, horizontalTR, -rectangleR, rectangleHorizBR, -rpieceR,-rpieceB,-rpieceL, edgeBetweenPieces, -lpieceR,-lpieceB,-lpieceL, rectangleHorizBL, -rectangleL, horizontalTL, verticalTL}; // Line Ids
Plane Surface(204) = {104}; // Loop Ids
Physical Surface("weirdShapeInTheMiddle") = {204};


Physical Curve("externalBoundary") = {stripeB, stripeR, verticalTR, horizontalTR, rectangleT, horizontalTL, verticalTL, stripeL}; // Line Ids
