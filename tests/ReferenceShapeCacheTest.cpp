// Unit tests of the matrices that the reference shapes compute once and store, for the whole process, for each basis
// and each diffusion tensor (src/Geometry/ReferenceShape.h, ReferenceCartesianShape.h). They were keyed by the address
// of the basis or of the tensor: a basis or a tensor allocated at the address of a deleted one (from a previous run in
// the same process) got its matrices. In the test suite, the face basis of a k=1 run got the 1x1 mass matrix of the
// face basis of a k=0 run: wrong discretization, U-AMG stalled at a convergence rate of 0.92. Here the objects are
// built at the same address on purpose (placement new), instead of relying on the allocator to reuse it.
#include <gtest/gtest.h>
#include <new>
#include "Utils/Utils.h"
#include "Utils/MatlabScript.h"
#include "FunctionalBasis/Legendre/LegendreBasis.h"
#include "Geometry/2D/Segment.h"
#include "TestCases/Diffusion/Tensor.h"

// The quadrature rules are allocated for each run of a program (Program_*::Execute())
class ReferenceShapeCache : public ::testing::Test
{
protected:
	void SetUp() override { GaussLegendre::Init(); }
	void TearDown() override { GaussLegendre::Free(); }
};

TEST_F(ReferenceShapeCache, BasisAtTheAddressOfADeletedOne)
{
	ReferenceCartesianShape<1>* refSegment = Segment::InitReferenceShape();
	alignas(LegendreBasis1D) unsigned char memory[sizeof(LegendreBasis1D)];

	LegendreBasis1D* p0 = new (memory) LegendreBasis1D(0);
	refSegment->ComputeAndStoreMassMatrix(p0);
	refSegment->ComputeAndStoreStiffnessMatrices(p0);
	refSegment->ComputeAndStoreIntegralVector(p0);
	p0->~LegendreBasis1D();

	LegendreBasis1D* p1 = new (memory) LegendreBasis1D(1);
	ASSERT_EQ((void*)p1, (void*)p0);
	refSegment->ComputeAndStoreMassMatrix(p1);
	refSegment->ComputeAndStoreStiffnessMatrices(p1);
	refSegment->ComputeAndStoreIntegralVector(p1);
	EXPECT_EQ(refSegment->StoredMassMatrix(p1).rows(), 2);
	EXPECT_EQ(refSegment->StoredStiffnessMatrices(p1).tt.rows(), 2);
	EXPECT_EQ(refSegment->StoredIntegralVector(p1).rows(), 2);
	p1->~LegendreBasis1D();
}

TEST_F(ReferenceShapeCache, TensorAtTheAddressOfADeletedOne)
{
	ReferenceCartesianShape<1>* refSegment = Segment::InitReferenceShape();
	LegendreBasis1D basis(1);
	BasisFunction<1>* phi1 = basis.LocalFunctions()[1]; // degree 1: non-zero gradient
	alignas(Tensor<1>) unsigned char memory[sizeof(Tensor<1>)];

	Tensor<1>* K1 = new (memory) Tensor<1>(1.0);
	refSegment->ComputeAndStoreReconstructStiffnessMatrix(*K1, &basis);
	double withK1 = refSegment->ReconstructStiffnessTerm(*K1, phi1, phi1);
	K1->~Tensor<1>();

	Tensor<1>* K2 = new (memory) Tensor<1>(2.0);
	ASSERT_EQ((void*)K2, (void*)K1);
	refSegment->ComputeAndStoreReconstructStiffnessMatrix(*K2, &basis);
	EXPECT_DOUBLE_EQ(refSegment->ReconstructStiffnessTerm(*K2, phi1, phi1), 2 * withK1);
	K2->~Tensor<1>();
}
