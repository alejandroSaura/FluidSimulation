using UnityEngine;
using System.Collections;

public abstract class IgridItem
{
    Vector3 x;
    Vector3 y;
    Vector3 z;

    public static IgridItem operator *(float operand, IgridItem item)
    {
        item.x *= operand;
        item.y *= operand;
        item.z *= operand;

        return item;
    }
    public static IgridItem operator +(IgridItem item1, IgridItem item2)
    {
        item1.x += item2.x;
        item1.y += item2.y;
        item1.z += item2.z;

        return item1;
    }

}
