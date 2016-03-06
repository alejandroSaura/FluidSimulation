using UnityEngine;
using System.Collections;
using System;

public class Matrix3x3 : IgridItem
{
    public Vector3 x;
    public Vector3 y;
    public Vector3 z;
    
    public Matrix3x3()
    {
        x = Vector3.zero;
        y = Vector3.zero;
        z = Vector3.zero;
    }

    public static Vector3 operator *(Vector3 vector, Matrix3x3 m)
    {
        return new Vector3(Vector3.Dot(m.x , vector), Vector3.Dot(m.y , vector), Vector3.Dot(m.z , vector));

    }
}
