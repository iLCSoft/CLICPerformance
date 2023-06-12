// This file's extension implies that it's C, but it's really -*- C++ -*-.
// $Id$

#include <cmath>
#include <functional>


#ifndef CXXUTILS_FPCOMPARE_H
#define CXXUTILS_FPCOMPARE_H


// Decide whether we need to use volatile or not.
#if defined(__FLT_EVAL_METHOD__) && \
  (__FLT_EVAL_METHOD__ == 2 || __FLT_EVAL_METHOD__ < 0)
  // __FLT_EVAL_METHOD__ < 0 means unspecified.
  // Be pessimistic in that case.
# define CXXUTILS_FPCOMPARE_VOLATILE volatile
#elif defined(__i386__) && !defined(__SSE2__)
  // On x86, gcc -msse -mfpmath=sse is observed to _not_ generate
  // sse fp instructions, but does set __FLT_EVAL_METHOD__ to 0.
  // -msse2 -mfpmath=sse does seem to work as expected.
  // Special-case this for now; should follow up with a gcc bug report
  // if this still happens in current releases.
# define CXXUTILS_FPCOMPARE_VOLATILE volatile
#else
# define CXXUTILS_FPCOMPARE_VOLATILE
#endif


namespace CxxUtils {
namespace fpcompare {


inline
bool equal (double a, double b)
{
  CXXUTILS_FPCOMPARE_VOLATILE double va = a;
  CXXUTILS_FPCOMPARE_VOLATILE double vb = b;
  return va == vb;
}


inline
bool equal (float a, float b)
{
  CXXUTILS_FPCOMPARE_VOLATILE float va = a;
  CXXUTILS_FPCOMPARE_VOLATILE float vb = b;
  return va == vb;
}


inline
bool greater (double a, double b)
{
  CXXUTILS_FPCOMPARE_VOLATILE double va = a;
  CXXUTILS_FPCOMPARE_VOLATILE double vb = b;
  return va > vb;
}


inline
bool greater (float a, float b)
{
  CXXUTILS_FPCOMPARE_VOLATILE float va = a;
  CXXUTILS_FPCOMPARE_VOLATILE float vb = b;
  return va > vb;
}


inline
bool less (double a, double b)
{
  CXXUTILS_FPCOMPARE_VOLATILE double va = a;
  CXXUTILS_FPCOMPARE_VOLATILE double vb = b;
  return va < vb;
}


inline
bool less (float a, float b)
{
  CXXUTILS_FPCOMPARE_VOLATILE float va = a;
  CXXUTILS_FPCOMPARE_VOLATILE float vb = b;
  return va < vb;
}


inline
bool greater_equal (double a, double b)
{
  CXXUTILS_FPCOMPARE_VOLATILE double va = a;
  CXXUTILS_FPCOMPARE_VOLATILE double vb = b;
  return va >= vb;
}


inline
bool greater_equal (float a, float b)
{
  CXXUTILS_FPCOMPARE_VOLATILE float va = a;
  CXXUTILS_FPCOMPARE_VOLATILE float vb = b;
  return va >= vb;
}


inline
bool less_equal (double a, double b)
{
  CXXUTILS_FPCOMPARE_VOLATILE double va = a;
  CXXUTILS_FPCOMPARE_VOLATILE double vb = b;
  return va <= vb;
}


inline
bool less_equal (float a, float b)
{
  CXXUTILS_FPCOMPARE_VOLATILE float va = a;
  CXXUTILS_FPCOMPARE_VOLATILE float vb = b;
  return va <= vb;
}


} // namespace fpcompare
} // namespace CxxUtils


#endif // not CXXUTILS_FPCOMPARE_H
