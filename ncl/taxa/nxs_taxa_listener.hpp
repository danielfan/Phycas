#ifndef NCL_NXS_TAXA_LISTENER_H
#define NCL_NXS_TAXA_LISTENER_H

class BaseTaxaManager;
/*----------------------------------------------------------------------------------------------------------------------
|	Interface that is the base for classes that need to be alerted when the taxon status changes
*/
class NxsTaxaListener: boost::noncopyable
	{
	
	public:
	enum TaxaChangeType
			{
			kTaxaCleared,
			kTaxaAdded,
			kTaxaDeleted,
			kTaxaExcluded,
			kTaxaIncluded,
			kTaxaReordered,
			kTaxaMgrDestroyed
			};
		virtual void 	TaxaChanged(BaseTaxaManager *, TaxaChangeType) = 0;
		void 		 	ChangeManager(BaseTaxaManager *);
		virtual void	ManagerIsDying()
							{
							taxaMgr = NULL;
							TaxaChanged(NULL, kTaxaMgrDestroyed);
							}
	
	protected:
		NxsTaxaListener(BaseTaxaManager *t);
		virtual ~NxsTaxaListener();
		
		BaseTaxaManager * taxaMgr; /* alias to the taxaManager*/
		
	};

#endif
