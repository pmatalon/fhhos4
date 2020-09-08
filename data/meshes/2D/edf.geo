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

//-------------------//
// Left little piece //
//-------------------//
pieceHeight = 1;
pieceWidth = 3;
pieceOriginX = 12;
pieceOriginY = 9;

lpieceBL = 5;
lpieceBR = 6;
lpieceTR = 7;
lpieceTL = 8;
Point(lpieceBL) = {pieceOriginX,            pieceOriginY,             0, h}; // bottom-left corner
Point(lpieceBR) = {pieceOriginX+pieceWidth, pieceOriginY,             0, h}; // bottom-right corner
Point(lpieceTR) = {pieceOriginX+pieceWidth, pieceOriginY+pieceHeight, 0, h}; // top-right corner
Point(lpieceTL) = {pieceOriginX,            pieceOriginY+pieceHeight, 0, h}; // top-left corner

lpieceB = 5;
lpieceR = 6;
lpieceT = 7;
lpieceL = 8;
Line(lpieceB) = {lpieceBL, lpieceBR}; // bottom
Line(lpieceR) = {lpieceBR, lpieceTR}; // right
Line(lpieceT) = {lpieceTR, lpieceTL}; // top
Line(lpieceL) = {lpieceTL, lpieceBL}; // left

Line Loop(101) = {lpieceB, lpieceR, lpieceT, lpieceL}; // Line Ids
Plane Surface(201) = {101}; // Loop Ids
Physical Surface("leftLittlePiece") = {201};

//--------------------//
// Right little piece //
//--------------------//
pieceOriginX += 5;

rpieceBL =  9;
rpieceBR = 10;
rpieceTR = 11;
rpieceTL = 12;
Point(rpieceBL) = {pieceOriginX,            pieceOriginY,             0, h}; // bottom-left corner
Point(rpieceBR) = {pieceOriginX+pieceWidth, pieceOriginY,             0, h}; // bottom-right corner
Point(rpieceTR) = {pieceOriginX+pieceWidth, pieceOriginY+pieceHeight, 0, h}; // top-right corner
Point(rpieceTL) = {pieceOriginX,            pieceOriginY+pieceHeight, 0, h}; // top-left corner

rpieceB =  9;
rpieceR = 10;
rpieceT = 11;
rpieceL = 12;
Line(rpieceB) = {rpieceBL, rpieceBR}; // bottom
Line(rpieceR) = {rpieceBR, rpieceTR}; // right
Line(rpieceT) = {rpieceTR, rpieceTL}; // top
Line(rpieceL) = {rpieceTL, rpieceBL}; // left

Line Loop(102) = {rpieceB, rpieceR, rpieceT, rpieceL}; // Line Ids
Plane Surface(202) = {102}; // Loop Ids
Physical Surface("rightLittlePiece") = {202};

//---------------//
// Top rectangle //
//---------------//
rectangleOriginX = 4;
rectangleOriginY = pieceOriginY+pieceHeight;
rectangleWidth = 24;

rectangleBL = 13;
rectangleBR = 14;
rectangleTR = 15;
rectangleTL = 16;
Point(rectangleBL) = {rectangleOriginX,                rectangleOriginY, 0, h}; // bottom-left corner
Point(rectangleBR) = {rectangleOriginX+rectangleWidth, rectangleOriginY, 0, h}; // bottom-right corner
Point(rectangleTR) = {rectangleOriginX+rectangleWidth, domainHeight,     0, h}; // top-right corner
Point(rectangleTL) = {rectangleOriginX,                domainHeight,     0, h}; // top-left corner

rectangleHorizBL = 13;
edgeBetweenPieces = 14;
rectangleHorizBR = 15;
rectangleR = 16;
rectangleT = 17;
rectangleL = 18;
Line(rectangleHorizBL)  = {rectangleBL, lpieceTL}; // bottom-left of the top rectangle
Line(edgeBetweenPieces) = {lpieceTR, rpieceTL};    // between the little pieces
Line(rectangleHorizBR)  = {rpieceTR, rectangleBR}; // bottom-right of the top rectangle
Line(rectangleR) = {rectangleBR, rectangleTR}; // right
Line(rectangleT) = {rectangleTR, rectangleTL}; // top
Line(rectangleL) = {rectangleTL, rectangleBL}; // left

Line Loop(103) = {rectangleHorizBL, -lpieceT, edgeBetweenPieces, -rpieceT, rectangleHorizBR, rectangleR, rectangleT, rectangleL}; // Line Ids
Plane Surface(203) = {103}; // Loop Ids
Physical Surface("topRectangle") = {203};


//---------------------------//
// Weird shape in the middle //
//---------------------------//
domainTL = 17;
domainTR = 18;
Point(domainTL) = {          0, domainHeight, 0, h}; // top-left corner of the domain
Point(domainTR) = {domainWidth, domainHeight, 0, h}; // top-right corner of the domain

verticalTL = 19;
verticalTR = 20;
horizontalTL = 21;
horizontalTR = 22;
Line(verticalTL)        = {domainTL,    stripeTL};    // left      {top pt, bottom pt}
Line(verticalTR)        = {stripeTR,    domainTR};    // right     {bottom pt, top pt}
Line(horizontalTL)      = {rectangleTL, domainTL};    // top-left  {right pt, left pt}
Line(horizontalTR)      = {domainTR,    rectangleTR}; // top-right {right pt, left pt}

Line Loop(104) = {-stripeT, verticalTR, horizontalTR, -rectangleR, -rectangleHorizBR, -rpieceR,-rpieceB,-rpieceL, -edgeBetweenPieces, -lpieceR,-lpieceB,-lpieceL, -rectangleHorizBL, -rectangleL, horizontalTL, verticalTL}; // Line Ids
Plane Surface(204) = {104}; // Loop Ids
Physical Surface("weirdShapeInTheMiddle") = {204};


Physical Curve("externalBoundary") = {stripeB, stripeR, verticalTR, horizontalTR, rectangleT, horizontalTL, verticalTL, stripeL}; // Line Ids
