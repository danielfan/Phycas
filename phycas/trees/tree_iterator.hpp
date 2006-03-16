#if !defined (PHO_TREE_ITERATORS_H)
#include "phycas/trees/node_iterator.hpp"
#include "phycas/trees/tree.hpp"

				
inline pre_iterator Tree::pre_begin() 
	{
	return pre_iterator(firstPreorder);
	}
	
inline const_pre_iterator Tree::pre_begin() const 
	{
	return const_pre_iterator(firstPreorder);
	}
	
inline post_iterator Tree::post_begin() 
	{
	return post_iterator(lastPreorder);
	}
	
inline const_post_iterator Tree::post_begin() const 
	{
	return const_post_iterator(lastPreorder);
	}
	
inline pre_iterator Tree::pre_end() 
	{
	return pre_iterator(NULL);
	}
	
inline const_pre_iterator Tree::pre_end() const 
	{
	return const_pre_iterator(NULL);
	}
	
inline post_iterator Tree::post_end() 
	{
	return post_iterator(NULL);
	}
	
inline const_post_iterator	Tree::post_end() const 
	{
	return const_post_iterator(NULL);
	}

inline tip_iterator Tree::tips_begin() 
	{
	return tip_iterator(firstPreorder);
	}
	
inline tip_iterator Tree::tips_end() 
	{
	return tip_iterator(NULL);
	}

inline tip_iterator Tree::leaves_begin() 
	{
	tip_iterator r(firstPreorder);
	++r;
	return r;
	}
	
inline tip_iterator Tree::leaves_end() 
	{
	return tip_iterator(NULL);
	}
	
inline const_tip_iterator	Tree::tips_end() const 
	{
	return const_tip_iterator(NULL);
	}
inline const_tip_iterator	Tree::tips_begin() const 
	{
	return const_tip_iterator(firstPreorder);
	}

inline const_tip_iterator Tree::leaves_begin() const
	{
	tip_iterator r(firstPreorder);
	++r;
	return r;
	}
	
inline const_tip_iterator Tree::leaves_end() const
	{
	return tip_iterator(NULL);
	}


#endif

