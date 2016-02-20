using UnityEngine;
using System.Collections;

public abstract class IMatrix3x3
{
    Vector3 x;
    Vector3 y;
    Vector3 z;

    public static IMatrix3x3 operator *(float operand, IMatrix3x3 item)
    {
        item.x *= operand;
        item.y *= operand;
        item.z *= operand;

        return item;
    }
    public static IMatrix3x3 operator +(IMatrix3x3 item1, IMatrix3x3 item2)
    {
        item1.x += item2.x;
        item1.y += item2.y;
        item1.z += item2.z;

        return item1;
    }
}
