using UnityEngine;
using System.Collections;

public class Vector : IgridItem
{
    public Vector3 v;    

    public Vector()
    {
        v = Vector3.zero;
    }

    public static Vector operator +(Vector item1, Vector item2)
    {
        item1.v.x += item2.v.x;
        item1.v.y += item2.v.y;
        item1.v.z += item2.v.z;

        return item1;
    }
}
