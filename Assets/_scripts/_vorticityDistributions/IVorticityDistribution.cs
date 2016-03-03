using UnityEngine;
using System.Collections;
using System;

public interface IVorticityDistribution
{
    Vector3 GetDomainSize();
    Vector3 GetCenter();
    void AssignVorticity(ref Vector3 vorticity, Vector3 position, Vector3 vCenter);
}