#if !defined (PHO_DISTRIBUTION_MANAGER_H)
#define PHO_DISTRIBUTION_MANAGER_H

#include "ncl/nxs_defs.hpp"
#include "phycas/rand/distribution_description.hpp"
class NxsCommandManager;
class DistributionSettings;
class DistributionDescription;

/*----------------------------------------------------------------------------------------------------------------------
|	class that stores DistributionDescription by user-assigned names
*/
class DistributionManager 
	{
	public :
		void							AddDistribution(const std::string &name, const DistributionDescription &descrip);
		const DistributionDescription * GetDistribution(const std::string &) const;
		VecString						GetReservedNames() const;
		
		void 							SetupDistributionCommand(NxsCommandManager *cmdMgr);
		
	protected:
		DistributionManager();
	
		CmdResult						HandleNewDistribution(DistributionSettings *);
		
		typedef std::map<std::string, DistributionDescription, NxsStringNoCaseLess> DistribMap;
		DistribMap			distributions;	/*map of name to distribution description */
		const VecString		reservedDistributionNames;

		friend DistributionManagerCreator;
	};
		
/*----------------------------------------------------------------------------------------------------------------------
|	returns a vector of words that cannot be used as names for distributions
*/
inline VecString DistributionManager::GetReservedNames() const
	{
	return reservedDistributionNames;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns pointer to a (const) DistributionDescription with the name `n'
*/
inline const DistributionDescription *DistributionManager::GetDistribution(
  const std::string &n) const
	{
	DistribMap::const_iterator dIt = distributions.find(n);
	if (dIt == distributions.end())
		return NULL;
	return &(dIt->second);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	adds the pair {n, descrip}, overwriting if n is already the name of a distribution description
*/
inline void DistributionManager::AddDistribution(
  const std::string &n, 
   const DistributionDescription &descrip)
   {
   distributions[n] = descrip;
   }
	
#endif
