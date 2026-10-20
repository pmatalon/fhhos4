#pragma once
#include <memory>
#include "StaticCondensation.h"
#include "../Utils/DistributedMatrix.h"
using namespace std;

// The static condensation of a symmetric hybrid system distributed by rows over the processes (see StaticCondensation):
// each process owns cells and faces; A_TT and A_TF are the rows of its cells, the columns of A_TF numbered in Faces (its
// faces, then the other faces of its cells: the ghosts); A_FF the rows of its faces, in full (both triangles), with the
// same columns. Each process condenses its cells (StaticCondensation); the contributions of a cell to the faces of other
// processes (their rows of A, their right-hand side) are added to their owners'.
class DistributedStaticCondensation
{
private:
	shared_ptr<const IndexMap> _faces;
	int _faceBlockSize;
	StaticCondensation _local; // the cells of the process: its A_FT = A_TF^T has the rows of all its local faces

public:
	DistributedStaticCondensation(const SparseMatrix& A_TT, const SparseMatrix& A_TF, shared_ptr<const IndexMap> faces, int cellBlockSize, int faceBlockSize)
		: _faces(faces), _faceBlockSize(faceBlockSize), _local(A_TT, A_TF, nullptr, cellBlockSize)
	{
		assert(A_TF.cols() == faces->NLocal() * faceBlockSize);
	}

	DistributedStaticCondensation(const DistributedStaticCondensation&) = delete;
	DistributedStaticCondensation& operator=(const DistributedStaticCondensation&) = delete;

	// A = A_FF - A_FT A_TT^{-1} A_TF, distributed by rows (the rows of the faces of each process), from the rows of A_FF of
	// the faces of the process, in full. Collective.
	DistributedMatrix CondensedMatrix(const SparseMatrix& A_FF) const
	{
		const IndexMap& faces = *_faces;
		BigNumber nLocal = faces.NLocal() * _faceBlockSize;
		assert(A_FF.rows() == faces.NOwned * _faceBlockSize && A_FF.cols() == nLocal);
		// A_FF - A_FT W on all the local faces (A_FF: 0 on the ghost faces), with the contributions of the cells of the
		// process only
		// The condensation of the cells of the process in the pattern of A_FF (StaticCondensation::CondensedMatrix()), whose
		// rows of the ghosts are empty: only these are completed
		SparseMatrix localRows = _local.CondensedMatrix(A_FF, false);
		return DistributedMatrix::AssembleFromLocalRows(_faces, localRows, _faceBlockSize);
	}

	// A_TT^{-1} b_T (the cells of the process)
	Vector SolveCells(const Vector& b_T) const
	{
		return _local.SolveCells(b_T);
	}

	// b = b_F - A_FT A_TT^{-1} b_T (the faces of the process), from invA_TT_b_T = A_TT^{-1} b_T. Collective.
	Vector CondensedRHS(const Vector& b_F, const Vector& invA_TT_b_T) const
	{
		const IndexMap& faces = *_faces;
		BigNumber nOwned = faces.NOwned * _faceBlockSize;
		Vector extendedB_F = Vector::Zero(faces.NLocal() * _faceBlockSize);
		extendedB_F.head(nOwned) = b_F;
		Vector b = _local.CondensedRHS(extendedB_F, invA_TT_b_T);
		Vector owned = b.head(nOwned);
		faces.AddGhostsToOwners(b.data() + nOwned, owned.data(), _faceBlockSize);
		return owned;
	}

	// x_T = A_TT^{-1} b_T - A_TT^{-1} A_TF x_F (the cells of the process), from invA_TT_b_T = A_TT^{-1} b_T and x_F on the
	// faces of the process. Collective.
	Vector CellUnknowns(const Vector& invA_TT_b_T, const Vector& x_F) const
	{
		return _local.CellUnknowns(invA_TT_b_T, _faces->WithGhosts(x_F, _faceBlockSize));
	}
};
