#include "phycas/force_include.h"
//  defining PHO_MEMCHK_CPP means that the functions in this file won't use the 
//	the overloaded new and delete.  This prevents odd behaviour resulting when the 
//	memory database storage containers are forced to resize

#if defined (MONITORING_ALLOCATION) && !defined(NDEBUG)
#define PHO_MEMCHK_CPP

bool PhorestMemoryAccountant::recording = false;

/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor does nothing currently.
*/
PhorestMemoryAccountant::PhorestMemoryAccountant() 
	: nBad(0L), 
	nUnknown(0L), 
	nAllocs(0L), 
	numBytesAllocated(0L),
	peakAllocation(0L),
	currentlyAllocated(0L)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing currently.
*/
PhorestMemoryAccountant::~PhorestMemoryAccountant()
	{
	filenames.erase(filenames.begin(), filenames.end());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Fills in fields of temporary PhorestMemoryInfo structure tmp and pushes it onto the mem_info vector.
*/
void PhorestMemoryAccountant::AddMemoryInfo(
  void			*p,		/* pointer to the allocated memory */
  const char	*fName,		/* name of file where allocation occurred */
  unsigned		ln,		/* line number where allocation occurred */
  unsigned long b,		/* number of bytes allocated (defaults to 0L, which means number of bytes is not being tracked) */
  bool 			is_arr)	/* true if the allocated using the new [] operator */
	{
	const string fn(fName);
	if (!PhorestMemoryAccountant::recording) 
		return;
	++nAllocs;
	numBytesAllocated += b;
	currentlyAllocated += b;
	if (currentlyAllocated > peakAllocation)
		peakAllocation = currentlyAllocated;
		
	// Attempt to find fn in the stored vector of file names
	//
	VecFileName::iterator i = find(filenames.begin(), filenames.end(), fn);
	if (i == filenames.end())
		{
		// fn has not previously been encountered
		//
		tmp.filename_index = filenames.size();
		filenames.push_back(fn);
		}
	else
		{
		// fn has not previously been encountered
		//
		tmp.filename_index = (unsigned) (i - filenames.begin());
		}
	
	tmp.ptr = p;
	tmp.flag = is_arr ? (unsigned) PhorestMemoryInfo::Mem_Array : 0U;
	tmp.num_bytes = b;
	tmp.line_number = ln;
	mem_info.push_back(tmp);
	++nBad;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Finds pointer p in mem_info vector and sets it to NULL, marking it as having been deleted.
*/
bool PhorestMemoryAccountant::MarkDeleted(
  void *p,		/* the pointer to be marked */
  bool is_arr)	/*true if the memory is being freed with the delete [] operator */
	{
	if (!recording || p == NULL)
		return false;

	bool found = false;
	bool middleArrayError = false;
	VecMemInfo::iterator i;
	for (i = mem_info.begin(); i != mem_info.end(); ++i)
		{
		if (i->ptr == p)
			{
			if (is_arr)
				{
				if (!(i->flag & PhorestMemoryInfo::Mem_Array))
					i->flag |= PhorestMemoryInfo::Mem_Op_Err_Free_Array;
				}
			else
				{
				if (i->flag & PhorestMemoryInfo::Mem_Array)
					i->flag |= PhorestMemoryInfo::Mem_Op_Err_New_Array;
				}
			i->ptr = NULL;
			--nBad;
			found = true;
			currentlyAllocated -= i->num_bytes;
			break;
			}
		else if (i->ptr < p && ((unsigned long) i->ptr + (unsigned long) i->num_bytes > (unsigned long) p))
			{
			middleArrayError = true;
			i->flag |= PhorestMemoryInfo::Mem_Mid_Array_Delete;
			assert(middleArrayError == false);
			}
		}
	assert(! (found && middleArrayError));	// shouldn't be possible to trip the middle of an array error and later find the mem object
	if (!found)
		++nUnknown;
	return found;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Summarizes memory allocations recorded, displaying total number of allocations, number of allocations currently
|	not deleted, and file name and line number for remaining undeleted elements.
*/
void PhorestMemoryAccountant::Summarize(
  ostream &out)	/* output stream for showing summary */
	{
	out << "\n\nMemory Report" << endl;

	out << "\nTotal allocations: " << nAllocs << " (" << numBytesAllocated << " bytes)" << endl;
	out << "\nLargest Memory Requirement: " << peakAllocation << "  bytes" << endl;
	if (nBad == 0L)
		out << "0 undeleted memory elements!" << endl;
	else if (nBad == 1L)
		out << "There was 1 undeleted memory element" << endl;
	else
		out << "There were " << nBad << " undeleted memory elements" << endl;

	// Compute approximate amount of memory required to just track allocs
	//
	unsigned z = mem_info.size() * sizeof(PhorestMemoryInfo);
	z += sizeof(PhorestMemoryAccountant);
	VecFileName::iterator j;
	for (j = filenames.begin(); j != filenames.end(); ++j)
		z += (*j).length();
	out << "The memory tracker itself required at least " << z << " bytes" << endl;

	VecMemInfo::iterator i;
	for (i = mem_info.begin(); i != mem_info.end(); ++i)
		{
		PhorestMemoryInfo &mi = (*i);
		if (mi.ptr != NULL || mi.flag > 1)
			{
			out << "  ";
			if (mi.num_bytes > 0L)
				out << mi.num_bytes << " bytes ";
			else
				out << "unknown number of bytes ";
			if (mi.ptr != NULL)
				out << "remaining from allocation at ";
			if (mi.flag & PhorestMemoryInfo::Mem_Op_Err_Free_Array)
				out << "allocated with new but freed with delete [].  Allocation at ";
			else if (mi.flag & PhorestMemoryInfo::Mem_Op_Err_New_Array)
				out << "allocated with new [] but freed with delete.  Allocation at ";
			else if (mi.flag & PhorestMemoryInfo::Mem_Mid_Array_Delete)
				out << "allocated but deletion occurred from the middle. Allocation at ";
			string s = filenames[mi.filename_index];
			out << s;
			out << " (" << mi.line_number << ')' << endl;
			}
		}
	}
#endif //if defined (MONITORING_ALLOCATION) && !defined(NDEBUG)


