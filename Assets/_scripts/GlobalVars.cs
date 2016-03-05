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
        public static float sAvoidSingularity = Mathf.Pow(float.Epsilon, 1.0f / 3.0f);
        public static float sTiny = Mathf.Exp(0.5f * (Mathf.Log(FLT_EPSILON) + Mathf.Log(float.Epsilon)));

        public static float finvsqrtf(float x)
        {
            float xhalf = 0.5f * x;
            int i = BitConverter.ToInt32(BitConverter.GetBytes(x), 0);
            i = 0x5f3759df - (i >> 1);
            x = BitConverter.ToSingle(BitConverter.GetBytes(i), 0);
            x = x * (1.5f - xhalf * x * x);
            x = x * (1.5f - xhalf * x * x);
            return x;
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

