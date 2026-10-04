// Runs the setup of U-AMG alone, on the inputs dumped by UAMG_DUMP_INPUTS() (see UAMGDump.h): iterating on the setup
// then takes a 30-second compilation and a few seconds per run, instead of a rebuild of Program.cpp and the assembly
// (178 s sequential at k=2 n=16). The setup code is that of src/ (header-only), so the instrumentation of Prof.h works.
// Compile it with build_uamg_harness.sh.
//
// Usage: uamg_harness <dump dir> <threads> [<save dir>]   setup time; with a save dir, writes the operator and the
//                                                         prolongation of each level (levelN_A.bin, levelN_P.bin)
//        uamg_harness compare <save dir A> <save dir B>   compares the levels saved by two runs (structure, differences)
#include <chrono>
#include <cmath>
#include <set>
#include "Solver/Multigrid/Multigrid.h"
#include "Solver/Multigrid/UncondensedAMG/UncondensedAMG.h"
#include "UAMGDump.h"

class HarnessAMG : public UncondensedAMG
{
public:
	using UncondensedAMG::UncondensedAMG;
	Level* Fine() { return this->_fineLevel; }
};

class LevelAccess : public Level
{
public:
	static const SparseMatrix& GetP(Level* l) { return static_cast<LevelAccess*>(l)->P; }
};

// Same structure? Then max |x - y| (relative to sqrt(|y_ii y_jj|) if scaleByDiagonal) and max relative difference
void Compare(const std::string& what, const SparseMatrix& X, const SparseMatrix& Y, bool scaleByDiagonal)
{
	bool sameStructure = X.rows() == Y.rows() && X.cols() == Y.cols() && X.nonZeros() == Y.nonZeros()
		&& std::equal(X.outerIndexPtr(), X.outerIndexPtr() + X.rows() + 1, Y.outerIndexPtr())
		&& std::equal(X.innerIndexPtr(), X.innerIndexPtr() + X.nonZeros(), Y.innerIndexPtr());
	cout << what << ": " << X.rows() << " x " << X.cols() << ", nnz " << X.nonZeros() << " vs " << Y.nonZeros() << (sameStructure ? ", same structure" : ", DIFFERENT STRUCTURE");
	if (!sameStructure)
	{
		cout << endl;
		return;
	}
	double maxDiff = 0, maxScaled = 0, maxRel = 0;
	BigNumber identical = 0;
	for (BigNumber i = 0; i < X.rows(); ++i)
	{
		for (SparseMatrixIndex p = X.outerIndexPtr()[i]; p < X.outerIndexPtr()[i + 1]; ++p)
		{
			double x = X.valuePtr()[p], y = Y.valuePtr()[p];
			BigNumber j = X.innerIndexPtr()[p];
			double d = abs(x - y);
			identical += d == 0;
			maxDiff = max(maxDiff, d);
			if (scaleByDiagonal)
				maxScaled = max(maxScaled, d / sqrt(abs(Y.coeff(i, i) * Y.coeff(j, j))));
			if (y != 0)
				maxRel = max(maxRel, d / abs(y));
		}
	}
	cout << ", identical " << identical << ", max |diff| " << maxDiff;
	if (scaleByDiagonal)
		cout << ", max |diff|/sqrt(|d_i d_j|) " << maxScaled;
	cout << ", max rel " << maxRel << endl;
}

int main(int argc, char** argv)
{
	if (argc < 3 || (string(argv[1]) == "compare" && argc < 4))
	{
		cerr << "Usage: uamg_harness <dump dir> <threads> [<save dir>]" << endl << "       uamg_harness compare <save dir A> <save dir B>" << endl;
		return 1;
	}

	if (string(argv[1]) == "compare")
	{
		string a = argv[2], b = argv[3];
		for (int n = 0; ; n++)
		{
			string f = "/level" + to_string(n) + "_A.bin";
			if (!ifstream(a + f) || !ifstream(b + f))
				break;
			Compare("level " + to_string(n) + " A", ReadSparse(a + f), ReadSparse(b + f), true);
			string g = "/level" + to_string(n) + "_P.bin";
			if (ifstream(a + g) && ifstream(b + g))
				Compare("level " + to_string(n) + " P", ReadSparse(a + g), ReadSparse(b + g), false);
		}
		return 0;
	}

	string dir = argv[1];
	Parallelism::SetNThreads(atoi(argv[2]));
	auto p = ReadParams(dir + "/params.txt");
	auto I = [&](const string& k) { return stoi(p.at(k)); };
	auto D = [&](const string& k) { return stod(p.at(k)); };

	Utils::ProgramArgs.Solver.MG.ManageAnisotropy = I("manageAniso");

	SparseMatrix S = ReadSparse(dir + "/S.bin");
	SparseMatrix A_T_T = ReadSparse(dir + "/A_T_T.bin");
	SparseMatrix A_T_F = ReadSparse(dir + "/A_T_F.bin");
	SparseMatrix A_F_F = ReadSparse(dir + "/A_F_F.bin");

	HarnessAMG mg(I("dim"), I("degree"), I("cellBS"), I("faceBS"), D("strong"), (UAMGFaceProlongation)I("faceProlong"),
		(UAMGProlongation)I("coarseningProlong"), (UAMGProlongation)I("mgProlong"), I("nLevels"));
	mg.MatrixMaxSizeForCoarsestLevel = I("maxSizeCoarsest");
	mg.Cycle = p.at("cycle")[0];
	mg.WLoops = I("wLoops");
	mg.UseGalerkinOperator = I("galerkin");
	mg.PreSmootherCode = p.at("preSmoother");
	mg.PostSmootherCode = p.at("postSmoother");
	mg.PreSmoothingIterations = I("preIt");
	mg.PostSmoothingIterations = I("postIt");
	mg.RelaxationParameter = D("omega");
	mg.BlockSizeForBlockSmoothers = I("blockSize");
	mg.CoarseLevelChangeSmoothingCoeff = I("clcCoeff");
	mg.CoarseLevelChangeSmoothingOperator = p.at("clcOp")[0];
	mg.HP_CS = (HP_CoarsStgy)I("HP_CS");
	mg.H_CS = (H_CoarsStgy)I("H_CS");
	mg.P_CS = (P_CoarsStgy)I("P_CS");
	mg.FaceCoarseningStgy = (FaceCoarseningStrategy)I("faceCoarsening");
	mg.BdryFaceCollapsing = (FaceCollapsing)I("bdryFaceCollapsing");
	mg.CoarseningFactor = D("coarseningFactor");
	mg.CoarsePolyDegree = I("coarsePolyDegree");
	mg.NumberOfMeshes = I("nMeshes");

	auto start = chrono::steady_clock::now();
	mg.Setup(S, A_T_T, A_T_F, A_F_F);
	cout << "HARNESS setup " << chrono::duration<double>(chrono::steady_clock::now() - start).count() << " s, " << mg.NumberOfLevels() << " levels" << endl;

	if (argc > 3)
	{
		string out = argv[3];
		int n = 0;
		for (Level* l = mg.Fine(); l; l = l->CoarserLevel, n++)
		{
			WriteSparse(*l->OperatorMatrix, out + "/level" + to_string(n) + "_A.bin");
			if (l->CoarserLevel)
				WriteSparse(LevelAccess::GetP(l), out + "/level" + to_string(n) + "_P.bin");
		}
		cout << "Saved " << n << " levels to " << out << endl;
	}
	return 0;
}
