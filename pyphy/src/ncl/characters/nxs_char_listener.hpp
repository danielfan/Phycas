#ifndef CHAR_LISTENER_H
#define CHAR_LISTENER_H

class NxsCharactersManager;
/*----------------------------------------------------------------------------------------------------------------------
|	Interface that is the base for classes that need to be alerted when the taxon status changes
*/
class NxsCharListener: boost::noncopyable
	{
	public:
	enum CharChangeType
			{
			kCharsCleared,
			kCharsAdded,
			kCharsDeleted,
			kCharsExcluded,
			kCharsIncluded,
			kCharsReordered,
			kCharsMgrDestroyed
			};
		virtual void 	CharsChanged(NxsCharactersManager *, CharChangeType) = 0;
		void		 	ChangeManager(NxsCharactersManager *);
		virtual void ManagerIsDying()
					{
					charMgr = NULL;
					CharsChanged(NULL, kCharsMgrDestroyed);
					}
		
	protected:
		NxsCharListener(NxsCharactersManager *);
		virtual ~NxsCharListener();
		
		NxsCharactersManager *charMgr;
	};
	
#endif
