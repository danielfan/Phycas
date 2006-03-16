#ifndef MANAGER_WITH_LISTENERS_H
#define MANAGER_WITH_LISTENERS_H

template<class ListenerType> class ManagerWithListeners
	{
	protected:
		std::vector<ListenerType *> listeners;
	public:
		void	AddListener(ListenerType *l)
			{
			listeners.push_back(l);
			}
	
		void	RemoveListener(ListenerType *l)
			{
			listeners.erase(std::remove(listeners.begin(), listeners.end(), l));
			}
		virtual ~ManagerWithListeners()
			{
			typedef typename std::vector<ListenerType *>::iterator ListenerIt;
			for (ListenerIt lIt = listeners.begin(); lIt != listeners.end(); ++lIt)
				(*lIt)->ManagerIsDying();
			}
	};


#endif
