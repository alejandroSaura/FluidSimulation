using System.Collections.Generic;
using UnityEngine;
using System;

namespace global
{
   

    public static class GlobalVar
    {
        public const float FLT_EPSILON = 1.192092896e-07F;        /* smallest such that 1.0+FLT_EPSILON != 1.0 */
        public const float FourPi = 4.0f * 3.1415926535897932384626433832795f;
        public const float OneOverFourPi = 1.0f / FourPi;
        public const float Nudge = 1.0f + FLT_EPSILON;
        public static float sAvoidSingularity = Mathf.Pow(float.MinValue, 1.0f / 3.0f);
        public static float sTiny = Mathf.Exp(0.5f * (Mathf.Log(FLT_EPSILON) + Mathf.Log(float.Epsilon)));

        public static float finvsqrtf(float val)
        {
            long i = (long)val;             // Exploit IEEE 754 inner workings.
            i = 0x5f3759df - (i >> 1);          // From Taylor's theorem and IEEE 754 format.
            float y = (float)i;              // Estimate of 1/sqrt(val) close enough for convergence using Newton's method.
            float f = 1.5f;        // Derived from Newton's method.
            float x = val * 0.5f;  // Derived from Newton's method.
            y = y * (f - (x * y * y));        // Newton's method for 1/sqrt(val)
            y = y * (f - (x * y * y));        // Another iteration of Newton's method
            return y;
        }
    }
    
    // add-on to the standard List class to resize the list
    public static class ListExtra
    {
        public static void Resize<T>(List<T> list, int size, T element = default(T))
        {
            int count = list.Count;

            if (size < count)
            {
                list.RemoveRange(size, count - size);
            }
            else if (size > count)
            {
                if (size > list.Capacity)   // Optimization
                    list.Capacity = size;

                int elementsToInsert = size - count;
                for(int i = 0; i < elementsToInsert; ++i)
                {
                    list.Add( (T)Activator.CreateInstance(typeof(T), new object[] {}) );
                }
                //list.AddRange(System.Linq.Enumerable.Repeat(element, size - count));
            }
        }
    }


}

