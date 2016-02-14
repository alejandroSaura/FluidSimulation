using System.Collections.Generic;

namespace global
{
    public static class GlobalVar
    {
        public const float FLT_EPSILON = 1.192092896e-07F;        /* smallest such that 1.0+FLT_EPSILON != 1.0 */
        public const float Nudge = 1.0f + FLT_EPSILON;

        
    }

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

                list.AddRange(System.Linq.Enumerable.Repeat(element, size - count));
            }
        }
    }


}

