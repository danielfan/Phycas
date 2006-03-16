#include "phycas/force_include.h"
#include <cassert>
#include <vector>
#include "phycas/trees/tree_counter.hpp"
using std::vector;

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of rooted trees having `m' internal nodes. Calls RecalculateCounts function if the `counts' 
|	array is not valid for rooted trees of `n' taxa.
*/
double TreeCounter::GetRootedCount(
  unsigned n,	/* number of taxa */
  unsigned m)/* number of internal nodes */
	{
	assert(n > 1);
	assert(m < n);
	if (n == 2)
		return 1.0;
	if (n != ntax)
		RecalculateCounts(n);
	return counts[m];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of fully-resolved rooted trees. Calls RecalculateCounts function if the `counts' array is not valid
|	for rooted trees of `n' taxa.
*/
double TreeCounter::GetSaturatedRootedCount(
  unsigned n)/* number of taxa */
	{
	assert(n > 1);
	if (n == 2)
		return 1.0;
	if (n != ntax)
		RecalculateCounts(n);
	return counts[n - 1];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of unrooted trees having `n' taxa and `m' internal nodes. The count returned is just the number of
|	rooted trees for one fewer taxa, so this function simply calls GetRootedCount(n - 1, m). Calls RecalculateCounts 
|	function if the `counts' array is not valid for unrooted trees of `n' taxa.
*/
double TreeCounter::GetUnrootedCount(
  unsigned n,	/* number of taxa */
  unsigned m)/* number of internal nodes */
	{
	assert(n > 2);
	return (n == 3 ? 1.0 : GetRootedCount(n - 1, m));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of fully-resolved unrooted trees. Calls RecalculateCounts function if the `counts' array is not valid
|	for unrooted trees of `n' taxa.
*/
double TreeCounter::GetSaturatedUnrootedCount(
  unsigned n)/* number of taxa */
	{
	assert(n > 2);
	return (n == 3 ? 1.0 : GetSaturatedRootedCount(n - 1));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns total number of rooted trees of all possible resolutions for `n' taxa. Calls RecalculateCounts function if 
|	the `counts' array is not valid for `n' taxa and rooted trees.
*/
double TreeCounter::GetTotalRooted(
  unsigned n)/* number of taxa */
	{
	assert(n > 1);
	if (n == 2)
		return 1.0;
	if (n != ntax)
		RecalculateCounts(n);
	return counts[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns total number of unrooted trees of all possible resolutions for `n' taxa. This number is just the total 
|	number of rooted trees for one fewer taxa, so this function simply calls GetTotalRooted(n - 1). Calls 
|	RecalculateCounts function if the `counts' array is not valid for `n' taxa and unrooted trees.
*/
double TreeCounter::GetTotalUnrooted(
  unsigned n) /* number of taxa */
	{
	assert(n > 2);
	return (n == 3 ? 1.0 : GetTotalRooted(n - 1));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputed counts vector for the new number of taxa supplied using the method outlined by J. Felsenstein in his 2003
|	book and also in Felsenstein (1978; The number of evolutionary trees. Syst. Zool. 27: 27-33) and Felsenstein
|	1981; Syst. Zool. 30: 122).
*/
void TreeCounter::RecalculateCounts(
  unsigned n) /* number of taxa */
	{
	assert(n > 2);
	if (n == ntax)
		return;
	unsigned first_n = 3;
	bool adding_taxa = true;
	if (ntax == 0)
		{
		counts.clear();
		counts.push_back(0.0);
		counts.push_back(1.0);
		}
	else if (ntax < n)
		first_n = ntax + 1;
	else if (ntax > n)
		adding_taxa = false;
	
	if (adding_taxa)
		{
		for (unsigned k = first_n; k <= n; ++k)
			{
			counts.push_back(0.0);
			RecalcNext(k);
			}
		}
	else
		{
		for (unsigned k = ntax; k > n; --k)
			{
			RecalcPrev(k);
			counts.erase(counts.end() - 1);
			}
		}
	ntax = n;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes `counts' array for `n' taxa given that it is already correctly computed for `n' - 1 taxa. Assumes that 
|	size of `counts' array is `n', which allows for first element to store total count and m to range from 1 to `n' - 1.
|	A new element should have been added to the `counts' array before this function was called.
*/
void TreeCounter::RecalcNext(unsigned n)
	{
	assert(counts.size() == n);
	counts[0] = 1.0;
	counts[1] = 1.0;
	double a = 1.0;
	for (unsigned m = 2; m < n; ++m)
		{
		double b = counts[m];
		double c = a*(n + m - 2);
		if (m < n - 1)
			c += b*m;
		counts[m] = c;
		counts[0] += c;
		a = b;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes `counts' array for `n' - 1 taxa given that it is already correctly computed for `n' taxa. Assumes that
|	size of `counts' array is `n'. Last element should be deleted from the `counts' array after this function is called.
*/
void TreeCounter::RecalcPrev(unsigned n)
	{
	assert(counts.size() == n);
	counts[0] = 1.0;
	counts[1] = 1.0;
	for (unsigned m = 2; m < n - 1; ++m)
		{
		double c = counts[m] - (double)(m + n - 2) * counts[m - 1];
		c /= (double) m;
		counts[m] = c;
		counts[0] += c;
		}
	}
