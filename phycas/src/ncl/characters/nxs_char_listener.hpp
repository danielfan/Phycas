/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
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
