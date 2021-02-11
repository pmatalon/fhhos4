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
		E* nextElement = &elements[0];
		int min = elements[0].NElementsIAmStrongNeighbourOf;
		for (E& e : elements)
		{
			if (e.NElementsIAmStrongNeighbourOf < min)
			{
				nextElement = &e;
				min = e.NElementsIAmStrongNeighbourOf;
			}
		}

		// Element aggregation
		while (nextElement)
		{
			E* elem = nextElement;
			nextElement = nullptr;
			int currentMin = elem->NElementsIAmStrongNeighbourOf;

			assert(!elem->CoarseElement);

			// New aggregate
			aggregates.emplace_back(aggregates.size());
			A* aggregate = &aggregates.back();

			// Add elem
			AddToAggregate(*elem, *aggregate);
			remainingElements--;

			E* neighbour = StrongestAvailableNeighbour(*elem);
			if (neighbour)
			{
				// Add neighbour?
				//bool strongConnectionIsReciprocal = find(neighbour->StrongNeighbours.begin(), neighbour->StrongNeighbours.end(), elem) != neighbour->StrongNeighbours.end();
				//if (strongConnectionIsReciprocal)
				//{
				AddToAggregate(*neighbour, *aggregate);
				remainingElements--;
				nPairs++;
				//}
				//else
					//nSingletons++;
			}
			else
				nSingletons++;

			nextElement = NextElementInTheNeighbourhood(*aggregate);

			if (!nextElement && remainingElements > 0)
			{
				// Browse the previous aggregates in the reverse order to find the next element
				BigNumber i = aggregate->Number - 1;
				while (!nextElement && i >= 0)
					//for (BigNumber i = aggregate->Number - 1; i >= 0; i--)
				{
					const A& previousAgg = aggregates[i];
					nextElement = NextElementInTheNeighbourhood(previousAgg);
					i--;
				}

				if (!nextElement)
				{
					auto it = find_if(elements.begin(), elements.end(), [](const E& e) { return !e.CoarseElement; });
					if (it != elements.end())
						nextElement = &*it;
				}
				/*min = 100000;
				for (E& e : elements)
				{
					if (!e.CoarseElement && e.NElementsIAmStrongNeighbourOf < min)
					{
						nextElement = &e;
						min = e.NElementsIAmStrongNeighbourOf;
					}
				}*/
			}
		}

		//cout << "Singletons: " << nSingletons << ", pairs: " << nPairs << endl;
		return aggregates;
	}
private:
	E* StrongestAvailableNeighbour(const E& e)
	{
		auto it = find_if(e.StrongNeighbours.begin(), e.StrongNeighbours.end(), [](E* n) { return !n->CoarseElement; });
		if (it != e.StrongNeighbours.end())
			return *it;
		return nullptr;
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

	E* NextElementInTheNeighbourhood(const A& aggregate)
	{
		E* nextElement = nullptr;
		int min = 1000000000;
		for (E* e : aggregate.FineElements)
		{
			for (pair<E*, double>& p : e->Neighbours)
			{
				E* n = p.first;
				if (!n->CoarseElement && n->NElementsIAmStrongNeighbourOf < min)
				{
					nextElement = n;
					min = n->NElementsIAmStrongNeighbourOf;
				}
			}
		}
		return nextElement;
	}
};