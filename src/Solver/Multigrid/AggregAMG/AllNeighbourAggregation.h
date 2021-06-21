#pragma once
#include "../../../Utils/Utils.h"
using namespace std;

template <class E, class A>
class AllNeighbourAggregation
{
public:
	vector<A> Perform(vector<E>& elements, bool& coarsestPossibleMeshReached)
	{
		vector<A> aggregates;
		aggregates.reserve(elements.size()); // we need to be sure that the vector won't be resized
		
		BigNumber remainingElements = elements.size();
		vector<E*> orderedList = OrderElementsByPriority(elements);

		// Element aggregation
		for (E* elem : orderedList)
		{
			if (elem->CoarseElement)
				continue;

			// New aggregate
			aggregates.emplace_back(aggregates.size());
			A* aggregate = &aggregates.back();

			// Add elem
			AddToAggregate(*elem, *aggregate);
			remainingElements--;

			for (int i = 0; i < elem->Neighbours.size(); i++)
			{
				const pair<E*, double>& p = elem->Neighbours[i];
				E* n = p.first;
				if (n->CoarseElement)
					continue;

				AddToAggregate(*n, *aggregate);
				remainingElements--;
			}

			if (remainingElements == 0)
				break;
		}

		double coarseningRatio = (double)elements.size() / (double)aggregates.size();
		if (coarseningRatio <= 1.01)
			coarsestPossibleMeshReached = true;

		//cout << "Singletons: " << nSingletons << ", pairs: " << nPairs << endl;
		return aggregates;
	}
private:
	vector<E*> OrderElementsByPriority(vector<E>& elements)
	{
		vector<E*> aggregationOrdering(elements.size());
		BigNumber orderNumber = 0;

		vector<bool> isNumbered(elements.size());
		for (BigNumber i = 0; i < elements.size(); ++i)
			isNumbered[i] = false;

		// Find one element that has the minimum value of NElementsIAmStrongNeighbourOf
		E* startingElement = &elements[0];
		int min = elements[0].NElementsIAmStrongNeighbourOf;
		for (E& e : elements)
		{
			if (e.NElementsIAmStrongNeighbourOf < min)
			{
				startingElement = &e;
				min = e.NElementsIAmStrongNeighbourOf;
			}
		}

		// Number this element
		aggregationOrdering[orderNumber++] = startingElement;
		isNumbered[startingElement->Number] = true;

		// Number the neighbours
		for (BigNumber i = 0; i < elements.size(); i++)
			NumberNeighbours(aggregationOrdering[i], aggregationOrdering, isNumbered, orderNumber);
		assert(orderNumber == elements.size());
		
		return aggregationOrdering;
	}

	void NumberNeighbours(E* e, vector<E*>& aggregationOrdering, vector<bool>& isNumbered, BigNumber& orderNumber)
	{
		// Get unnumbered neighbours
		vector<E*> unnumberedNeighbours;
		for (pair<E*, double>& p : e->Neighbours)
		{
			if (!isNumbered[p.first->Number])
				unnumberedNeighbours.push_back(p.first);
		}
		if (unnumberedNeighbours.empty())
			return;

		// Sort neighbours by ascending NElementsIAmStrongNeighbourOf
		sort(unnumberedNeighbours.begin(), unnumberedNeighbours.end(), [](E* n1, E* n2) { return n1->NElementsIAmStrongNeighbourOf < n2->NElementsIAmStrongNeighbourOf; });
		
		// Number neighbours
		for (E* n : unnumberedNeighbours)
		{
			aggregationOrdering[orderNumber++] = n;
			assert(!isNumbered[n->Number]);
			isNumbered[n->Number] = true;
		}
	}

	void AddToAggregate(E& e, A& aggregate)
	{
		aggregate.FineElements.push_back(&e);
		e.CoarseElement = &aggregate;

		for (E* n : e.StrongNeighbours)
		{
			if (!n->CoarseElement)
				n->NElementsIAmStrongNeighbourOf--;
		}
	}
};