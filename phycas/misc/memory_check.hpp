#ifndef PHO_MEMORYACCOUNTANT
#	define PHO_MEMORYACCOUNTANT

#	if defined (MONITORING_ALLOCATION) && !defined(NDEBUG)
#		define CREATE_MEMCHK					\
		PhorestMemoryAccountant memAccountant;	\
		memTracker = &memAccountant;			\
		memTracker->StartRecording();
#		define MEMCHK_REPORT(o)					\
		memTracker->StopRecording();			\
		memAccountant.Summarize(o);
#	else
#		define CREATE_MEMCHK
#		define MEMCHK_REPORT(o)
#	endif

#	if defined(MONITORING_ALLOCATION) && !defined(NDEBUG)

	class PhorestMemoryAccountant;

#	if defined (INSTANTIATE_MEMCHK)
		//struct dmalloc_t dmalloc;
		PhorestMemoryAccountant* memTracker = NULL;
#	else
		//extern struct dmalloc_t dmalloc;
		extern PhorestMemoryAccountant* memTracker;
#	endif

	/*----------------------------------------------------------------------------------------------------------------------
	|	Structure used to store accounting information for one individual memory allocation.
	*/
	class PhorestMemoryInfo
		{
		public:
		void		*ptr;
		unsigned 	filename_index;
		unsigned 	line_number;
		unsigned 	num_bytes;
		enum	{	Mem_Array 				= 0x01,	/* was allocated using the new [] operator */
					Mem_Op_Err_Free_Array 	= 0x02, /* was allocated using the new, but freed with delete [] */
					Mem_Op_Err_New_Array 	= 0x04, /* was allocated using the new[], but freed with delete */
					Mem_Mid_Array_Delete	= 0x08  /* delete called with an address in the middle of array */
				};
		unsigned	flag;
		bool		array_allocation;
		PhorestMemoryInfo(bool is_arr = false);
		};
		
	inline PhorestMemoryInfo::PhorestMemoryInfo(
	  bool is_arr)	/*true if the memory was allocated using new [] */
		{
		flag = is_arr ? (unsigned) Mem_Array : 0U;
		ptr = (void *)NULL;
		filename_index = 0L;
		line_number = 0L;
		num_bytes = 0L;
		}

	/*----------------------------------------------------------------------------------------------------------------------
	|	Keeps track of memory allocations and deletions if AddMemoryInfo is called after every call to the new operator and 
	|	MarkDeleted is called for every call to the delete or delete [] operator, which can be done using macros so that 
	|	memory accounting is not done in the release version.
	*/
	class PhorestMemoryAccountant
		{
		typedef vector<PhorestMemoryInfo>	VecMemInfo;
		typedef vector<string>				VecFileName;

		VecFileName			filenames;			/* vector of source file names used to provide filename_index element in the PhorestMemoryInfo struct */
		PhorestMemoryInfo	tmp;				/* workspace used when adding a new mem_info element */
		VecMemInfo			mem_info;			/* vector of allocation records */
		unsigned			nBad;				/* number of remaining allocs that have not yet been deleted */
		unsigned			nUnknown;			/* number of deletes on non-NULL pointers that are not in our database */
		unsigned			nAllocs;			/* total number of allocs that we have caught with our overload of new */
		unsigned long		numBytesAllocated;	/* total number of bytes from allocs that we have caught with our overload of new */
		unsigned long		currentlyAllocated;	/* the number of bytes that are allocated but haven't been freed (only catchs allocations through our new operator) */
		unsigned long		peakAllocation;		/* the most bytes that are ever allocated at one time (only catchs allocations through our new operator) */
		public:
		STATIC_DATA bool			recording;			/* if true, records allocations and deletions; otherwise, ignores them */
			PhorestMemoryAccountant();
			~PhorestMemoryAccountant();

			void StartRecording();
			void StopRecording();

			void AddMemoryInfo(void *p, const char * fn, unsigned ln, unsigned long b = 0L, bool is_arr = false);
			bool MarkDeleted(void *p, bool is_arr = false);

			void Summarize(ostream &out);
		};

	inline void PhorestMemoryAccountant::StartRecording()
		{
		recording = true;
		}

	inline void PhorestMemoryAccountant::StopRecording()
		{
		recording = false;
		}

	inline void *operator new (size_t size, const char *file, int line) throw (std::bad_alloc)
		{
		void *p = malloc(size);
		if (p == NULL)
			throw std::bad_alloc();
		if (memTracker != NULL && PhorestMemoryAccountant::recording)
			memTracker->AddMemoryInfo(p, file, (unsigned)line, size, false);
		return p;
		}

	inline void operator delete (void *p)
		{
		if (p != NULL)
			{
			if (memTracker != NULL && PhorestMemoryAccountant::recording && memTracker->MarkDeleted(p, false))
				free(p);
			else	
				free(p);
			}
		}

#	if !defined (_MSC_VER)
		inline void * operator new [] (size_t size, const char *file, int line)
			{
		 	void *p = malloc (size);
			if (p == NULL)
				throw std::bad_alloc();
			if (memTracker != NULL && PhorestMemoryAccountant::recording)
					memTracker->AddMemoryInfo(p, file, (unsigned) line, size, true);
			return p;
			}

		inline void operator delete [](void *p)
			{
			if (p != NULL)
				{
				if (memTracker != NULL && PhorestMemoryAccountant::recording && (memTracker->MarkDeleted(p, true)))
					free(p);
				else
					free(p);
				}
			}
#	endif

#	if !defined(PHO_MEMCHK_CPP)
			//#define NEW new (dmalloc, __FILE__, __LINE__)
#		define NEW new (__FILE__, __LINE__)
#		define new NEW
#	endif

#	endif // #if defined(MONITORING_ALLOCATION) && !defined(NDEBUG)
#endif // #ifndef PHO_MEMORYACCOUNTANT
