#if POLPY_NEWWAY

#if ! defined(NCAT_MOVE_INL)
#define NCAT_MOVE_INL

namespace phycas
{

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
inline void NCatMove::accept()
	{
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns current value of `phi', the probability that an addcat move (as opposed to a delcat move) will be proposed each
|	time update() is called.
*/
inline double NCatMove::getPhi() const
	{
	return phi;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets the value of `phi', the probability that an addcat move (as opposed to a delcat move) will be proposed each
|	time update() is called. Assumes `p' is in the interval (0,1).
*/
inline void NCatMove::setPhi(double new_phi)
	{
	assert(new_phi > 0.0);
	assert(new_phi < 1.0);
	phi = new_phi;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns current value of `ncat_max', the maximum number of rate categories encountered.
*/
inline unsigned NCatMove::getNCatMax() const
	{
	return ncat_max;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets the value of `ncat_max', the maximum number of rate categories encountered. Use this member function with care: if
|	you already know that the tree is equipped for, say, 4 rate categories, then you can call setNCatMax(4), but be aware
|	that this would later result in a crash if in fact the tree was only equipped for 3 rate categories. The safest thing to
|	do is to let `ncat_max' start at its default value of 1.
*/
inline void NCatMove::setNCatMax(unsigned new_ncat_max)
	{
	assert(ncat_max > 0);
	ncat_max = new_ncat_max;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns current value of `lambda', the parameter determining the Poisson prior on the number of categories.
*/
inline double NCatMove::getLambda() const
	{
	return lambda;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets the value of `lambda', the parameter determining the Poisson prior on the number of categories. Assumes 
|	`new_lambda' is greater than zero.
*/
inline void NCatMove::setLambda(double new_lambda)
	{
	assert(new_lambda > 0.0);
	lambda = new_lambda;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns current value of `L', the upper limit of the valid interval for unnormalized relative rate parameters.
*/
inline double NCatMove::getL() const
	{
	return L;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets the value of `L', the upper limit of the valid interval for unnormalized relative rate parameters. Assumes `new_L' 
	is greater than zero. There is normally no reason to set this to anything other than the default value of 1.0.
*/
inline void NCatMove::setL(double new_L)
	{
	assert(new_L > 0.0);
	L = new_L;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns current value of `s', the number of fake relative rates between each real relative rate.
*/
inline unsigned NCatMove::getS() const
	{
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets the value of `s', the number of fake relative rates between each real relative rate. The value of `s' determines 
|	the informativeness of the order statistics prior on relative rates. Larger values of `s' tend to make the differences 
|	between adjacent relative rates more homogeneous. Setting `s' to 0 provides a flat prior (as if relative rates were
|	drawn from a Uniform(0,L) distribution and then ordered from lowest to highest). There is normally no reason to set 
|	this to anything other than the default value of 1. Assumes `new_s' is greater than or equal to zero.
*/
inline void NCatMove::setS(unsigned new_s)
	{
	assert(new_s >= 0);
	s = new_s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns current value of `cat_prob_prior', the prior distribution for unnormalized category probability parameters. This
|	distribution is used to choose new category probabilities in an ADDCAT move.
*/
inline ProbDistShPtr NCatMove::getCatProbPrior() const
	{
	return cat_prob_prior;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets the value of `cat_prob_prior', the prior distribution for unnormalized category probability parameters. This
|	distribution is used to choose new category probabilities in an ADDCAT move. Assumes that `new_cat_prob_prior' is not
|	an empty pointer.
*/
inline void NCatMove::setCatProbPrior(ProbDistShPtr new_cat_prob_prior)
	{
	assert(new_cat_prob_prior); //@POL should we also check to make sure the distribution pointed to is either Gamma or Exponential?
	cat_prob_prior = new_cat_prob_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of `addcat_move_proposed' data member, which is true if the last move proposed was an
|	add-edge move and false if last move proposed was a delete-edge move.
*/
inline bool NCatMove::addCatMoveProposed() const
	{
	return addcat_move_proposed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ln_hastings', which is the natural log of the Hastings ratio for this move and which is 
|	computed in both NCatMove::ProposeAddEdgeMove and NCatMove::ProposeDeleteEdgeMove.
*/
inline double NCatMove::getLnHastingsRatio() const
	{
	return ln_hastings;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ln_jacobian', which is the natural log of the Jacobian for this move and which is computed in
|	both NCatMove::ProposeAddEdgeMove and NCatMove::ProposeDeleteEdgeMove.
*/
inline double NCatMove::getLnJacobian() const
	{
	return ln_jacobian;
	}

} // namespace phycas

#endif

#endif //POLPY_NEWWAY

