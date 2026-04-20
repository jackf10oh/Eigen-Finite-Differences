// RowEvaluator.hpp
//
// linops own internal evaluator similar to the eigen library's. 
//
// JAF 4/18/2026 

#ifndef ROWEVALUATORS_H
#define ROWEVALUATORS_H 

// Forward declare the actual implementation 
template</* TODO what does this actually need passed in......*/>
struct RowEvaluator_impl; 

template<class U>
struct RowEvaluator : public RowEvaluator_impl</* TODO pull eveything impl needs with traits...*/>{}; 

template<class Derived>
struct RowEvaluator_impl<PartialDerivBase<Derived>> : RowEvaluator_impl<Derived>{}; 

template<>
struct RowEvaluator_impl<NthPartialDeriv>
{
  // TODO heavy lifting here.....
}

template<class Operator, class ArgType1, class ArgType2>
struct RowEvaluator_impl< NwiseUnaryOp<Operator, ArgType1, ArgType2> >
{
  // traits or something idk 

  // store left and right sides row evaluator

}

template<class Operator, class ArgType>
struct RowEvaluator_impl< NwiseUnaryOp<Operator, ArgType> >
{
  // traits or something 

  // storeage o
}

#endif // RowEvaluator.hpp 