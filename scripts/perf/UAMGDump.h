#pragma once
// Dump of the inputs of the setup of U-AMG, to run the setup alone in uamg_harness.cpp (seconds instead of a 2-minute
// rebuild of Program.cpp plus the assembly). Not part of the build: include it in UncondensedAMG.h while measuring,
// then revert it:
//
//     #include "../../../../scripts/perf/UAMGDump.h"   // in UncondensedAMG.h, after the other includes
//     ...
//     void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
//     {
//         UAMG_DUMP_INPUTS();   // first line: dumps and exits if FHHOS4_UAMG_DUMP is set
//
// then, from build/: FHHOS4_UAMG_DUMP=<existing dir> ./bin/fhhos4 <arguments of the run>
// The directory receives S, A_T_T, A_T_F, A_F_F (binary, see WriteSparse), cellInterpOfOne, faceInterpOfOne (see
// WriteVector) and params.txt (the parameters of the multigrid). The program exits after the dump.
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include "Utils/Types.h"

inline void WriteSparse(const SparseMatrix& M0, const std::string& path)
{
	SparseMatrix M = M0;
	M.makeCompressed();
	std::ofstream f(path, std::ios::binary);
	int64_t h[3] = { (int64_t)M.rows(), (int64_t)M.cols(), (int64_t)M.nonZeros() };
	f.write((const char*)h, sizeof(h));
	f.write((const char*)M.outerIndexPtr(), (M.outerSize() + 1) * sizeof(SparseMatrixIndex));
	f.write((const char*)M.innerIndexPtr(), M.nonZeros() * sizeof(SparseMatrixIndex));
	f.write((const char*)M.valuePtr(), M.nonZeros() * sizeof(double));
	if (!f)
	{
		std::cerr << "cannot write " << path << std::endl;
		std::exit(1);
	}
}

inline SparseMatrix ReadSparse(const std::string& path)
{
	std::ifstream f(path, std::ios::binary);
	if (!f)
	{
		std::cerr << "cannot read " << path << std::endl;
		std::exit(1);
	}
	int64_t h[3];
	f.read((char*)h, sizeof(h));
	SparseMatrix M(h[0], h[1]);
	M.resizeNonZeros(h[2]);
	f.read((char*)M.outerIndexPtr(), (M.outerSize() + 1) * sizeof(SparseMatrixIndex));
	f.read((char*)M.innerIndexPtr(), h[2] * sizeof(SparseMatrixIndex));
	f.read((char*)M.valuePtr(), h[2] * sizeof(double));
	return M;
}

inline void WriteVector(const Vector& v, const std::string& path)
{
	std::ofstream f(path, std::ios::binary);
	int64_t n = v.rows();
	f.write((const char*)&n, sizeof(n));
	f.write((const char*)v.data(), n * sizeof(double));
	if (!f)
	{
		std::cerr << "cannot write " << path << std::endl;
		std::exit(1);
	}
}

inline Vector ReadVector(const std::string& path)
{
	std::ifstream f(path, std::ios::binary);
	if (!f)
	{
		std::cerr << "cannot read " << path << std::endl;
		std::exit(1);
	}
	int64_t n;
	f.read((char*)&n, sizeof(n));
	Vector v(n);
	f.read((char*)v.data(), n * sizeof(double));
	return v;
}

// "key value" lines
inline std::map<std::string, std::string> ReadParams(const std::string& path)
{
	std::map<std::string, std::string> params;
	std::ifstream f(path);
	std::string key, value;
	while (f >> key >> value)
		params[key] = value;
	return params;
}

// In UncondensedAMG::Setup(A, A_T_T, A_T_F, A_F_F, cellInterpOfOne, faceInterpOfOne)
#define UAMG_DUMP_INPUTS() \
	if (const char* dumpDir = getenv("FHHOS4_UAMG_DUMP")) \
	{ \
		std::string d = dumpDir; \
		WriteSparse(A, d + "/S.bin"); \
		WriteSparse(A_T_T, d + "/A_T_T.bin"); \
		WriteSparse(A_T_F, d + "/A_T_F.bin"); \
		WriteSparse(A_F_F, d + "/A_F_F.bin"); \
		WriteVector(cellInterpOfOne, d + "/cellInterpOfOne.bin"); \
		WriteVector(faceInterpOfOne, d + "/faceInterpOfOne.bin"); \
		std::ofstream p(d + "/params.txt"); \
		p << "dim " << _dim << "\ndegree " << _degree << "\ncellBS " << _cellBlockSize << "\nfaceBS " << _faceBlockSize \
		  << "\nstrong " << _strongCouplingThreshold << "\nfaceProlong " << (unsigned)_faceProlong \
		  << "\ncoarseningProlong " << (unsigned)_coarseningProlong << "\nmgProlong " << (unsigned)_multigridProlong \
		  << "\nnLevels " << this->NumberOfLevels() << "\nmaxSizeCoarsest " << this->MatrixMaxSizeForCoarsestLevel \
		  << "\ncycle " << this->Cycle << "\nwLoops " << this->WLoops << "\ngalerkin " << this->UseGalerkinOperator \
		  << "\npreSmoother " << this->PreSmootherCode << "\npostSmoother " << this->PostSmootherCode \
		  << "\npreIt " << this->PreSmoothingIterations << "\npostIt " << this->PostSmoothingIterations \
		  << "\nomega " << this->RelaxationParameter << "\nblockSize " << this->BlockSizeForBlockSmoothers \
		  << "\nclcCoeff " << this->CoarseLevelChangeSmoothingCoeff << "\nclcOp " << this->CoarseLevelChangeSmoothingOperator \
		  << "\nHP_CS " << (unsigned)this->HP_CS << "\nH_CS " << (unsigned)this->H_CS << "\nP_CS " << (unsigned)this->P_CS \
		  << "\nfaceCoarsening " << (unsigned)this->FaceCoarseningStgy << "\nbdryFaceCollapsing " << (unsigned)this->BdryFaceCollapsing \
		  << "\ncoarseningFactor " << this->CoarseningFactor << "\ncoarsePolyDegree " << this->CoarsePolyDegree \
		  << "\nnMeshes " << this->NumberOfMeshes << "\nmanageAniso " << this->ManageAnisotropy << std::endl; \
		std::cout << "U-AMG inputs dumped to " << d << std::endl; \
		std::exit(0); \
	}
