/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
