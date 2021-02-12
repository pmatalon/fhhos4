#pragma once
#include "../Utils/Utils.h"
using namespace std;

template <class E, class A>
class PairwiseAggregation
{
public:
	vector<A> Perform(vector<E>& elements, bool& coarsestPossibleMeshReached)
	{
		vector<A> aggregates;
		aggregates.reserve(elements.size()); // we need to be sure that the vector won't be resized

		int nSingletons = 0;
		int nPairs = 0;

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

			E* neighbour = StrongestAvailableNeighbour(*elem);
			if (neighbour)
			{
				AddToAggregate(*neighbour, *aggregate);
				remainingElements--;
				nPairs++;
			}
			else
				nSingletons++;

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
		BigNumber currentElement = 0;
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
			isNumbered[n->Number] = true;
		}
	}

	E* StrongestAvailableNeighbour(const E& e)
	{
		// !!!!!! Be very careful in changing this !!!!!!
		// A bad choice of neighbour may result, after several aggregation pass, at elements that have only one neighbour,
		// which leads to a very slow coarsening.
		double min = 100000000;
		E* chosenOne = nullptr;
		for (E* n : e.StrongNeighbours)
		{
			if (!n->CoarseElement && n->NElementsIAmStrongNeighbourOf < min)
			{
				chosenOne = n;
				min = n->NElementsIAmStrongNeighbourOf;
			}
		}
		return chosenOne;

		/*auto it = find_if(e.StrongNeighbours.begin(), e.StrongNeighbours.end(), [](E* n) { return !n->CoarseElement; });
		if (it != e.StrongNeighbours.end())
			return *it;
		return nullptr;*/

		/*E* firstAvailableStrongNeighbour = nullptr;
		auto it = find_if(e.StrongNeighbours.begin(), e.StrongNeighbours.end(), [](E* n) { return !n->CoarseElement; });
		if (it != e.StrongNeighbours.end())
			firstAvailableStrongNeighbour = *it;
		else
			return nullptr;

		vector<E*> candidates;
		candidates.push_back(firstAvailableStrongNeighbour);
		double bestAvailableCoupling = 1;
		for (int i = 0; i < e.Neighbours.size(); i++)
		{
			const pair<E*, double>& p = e.Neighbours[i];
			E* n = p.first;
			if (n == firstAvailableStrongNeighbour)
			{
				bestAvailableCoupling = p.second;
				continue;
			}
			if (bestAvailableCoupling != 1 && !n->CoarseElement)
			{
				if (abs(p.second - bestAvailableCoupling) < Utils::Eps * bestAvailableCoupling)
					candidates.push_back(n);
			}
		}

		double min = firstAvailableStrongNeighbour->NElementsIAmStrongNeighbourOf;
		E* chosenOne = firstAvailableStrongNeighbour;
		for (E* n : candidates)
		{
			if (n->NElementsIAmStrongNeighbourOf < min)
			{
				chosenOne = n;
				min = n->NElementsIAmStrongNeighbourOf;
			}
		}
		return chosenOne;*/

		/*double bestCoupling = e.Neighbours[0].second;

		for (int i = 0; i< e.Neighbours.size(); i++)
		{
			const pair<E*, double>& p = e.Neighbours[i];
			E* n = p.first;
			if (n->CoarseElement)
				continue;

			if (i == 0)
				return n;

			double coupling = p.second;
			if (abs(coupling - bestCoupling) < Utils::Eps * bestCoupling)
				return n;
			else
				return nullptr;
		}
		return nullptr;*/
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