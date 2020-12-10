#pragma once
#include "ParallelLoop.h"
#include "../Mesh/Element.h"
using namespace std;

template<int Dim>
using ElementParallelLoop = ParallelLoop<Element<Dim>*, CoeffsChunk>;
template<int Dim>
using FaceParallelLoop = ParallelLoop<Face<Dim>*, CoeffsChunk>;