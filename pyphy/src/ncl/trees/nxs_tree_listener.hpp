#ifndef TREE_LISTENER_H
#define TREE_LISTENER_H

class NxsTreesManager;
/*----------------------------------------------------------------------------------------------------------------------
|	Interface that is the base for classes that need to be alerted when the taxon status changes
*/
class NxsTreeListener: boost::noncopyable
	{
	public:
	enum TreeChangeType
			{
			kTreesCleared,
			kTreesAdded,
			kTreesDeleted,
			kTreesReordered,
			kTreesMgrDestroyed
			};
		virtual void 	TreesChanged(NxsTreesManager *, TreeChangeType) = 0;
		void 			ChangeManager(NxsTreesManager *);
		virtual void	ManagerIsDying()
							{
							treeMgr = NULL;
							TreesChanged(NULL, kTreesMgrDestroyed);
							}
		
	protected:
		NxsTreeListener(NxsTreesManager *);
		virtual ~NxsTreeListener();
		
		NxsTreesManager *treeMgr;
	};
#endif
