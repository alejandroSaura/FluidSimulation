using UnityEngine;
using System.Collections;
using System;

// Specify vorticity in the shape of a vortex ring
public class VortexRing : MonoBehaviour, IVorticityDistribution
{
    public float mRadius = 1.5f;
    public float mThickness = 0.5f;
    //public Vector3 mDirection;

    public Vector3 GetCenter()
    {
        return transform.position;
    }

    public VortexRing(float fRadius, float fThickness)
    {
        mRadius = fRadius;
        mThickness = fThickness;
        //mDirection = vDirection;        
    }    

    public Vector3 GetDomainSize()
    {
        float boxSideLength = 2.0f * (mRadius + mThickness); // length of side of virtual cube
        return new Vector3(1.0f, 1.0f, 1.0f) * boxSideLength;
    }

    public void AssignVorticity(ref Vector3 vorticity, Vector3 position, Vector3 vCenter)
    {
        Vector3 mDirection = transform.forward;

        Vector3 vFromCenter = position - vCenter; // displacement from ring center to vorton position
        float tween = Vector3.Dot(vFromCenter, mDirection); // projection of position onto axis
        Vector3 vPtOnLine = vCenter + mDirection * tween; // closest point on axis to vorton position
        Vector3 vRho = position - vPtOnLine; // direction radially outward from annulus core
        float rho = vRho.magnitude; // distance from axis
        float distAlongDir = Vector3.Dot(mDirection, vFromCenter); // distance along axis of vorton position
        float radCore = Mathf.Sqrt(Mathf.Pow(rho - mRadius, 2) + Mathf.Pow(distAlongDir, 2)); // distance from annular core
        if (radCore < mThickness)
        {   // Probe position is inside vortex core.
            float vortProfile = radCore < mThickness ? 0.5f * (Mathf.Cos(Mathf.PI * radCore / mThickness) + 1.0f) : 0.0f;
            float vortPhi = vortProfile;
            Vector3 rhoHat = vRho;                        // direction radially away from annular core
            rhoHat.Normalize();
            Vector3 phiHat = Vector3.Cross(mDirection, rhoHat);         // direction along annular core
            vorticity = vortPhi * phiHat;
        }
        else
        {
            vorticity = Vector3.zero;
        }
    }

    void OnDrawGizmos()
    {
        DebugExtension.DrawCircle(transform.position, transform.forward, Color.white, mRadius);

        DebugExtension.DrawCircle(transform.position + transform.right * mRadius, transform.up, Color.white, mThickness / 2);
        DebugExtension.DrawCircle(transform.position - transform.right * mRadius, transform.up, Color.white, mThickness / 2);

        DebugExtension.DrawCircle(transform.position + transform.up * mRadius, transform.right, Color.white, mThickness / 2);
        DebugExtension.DrawCircle(transform.position - transform.up * mRadius, transform.right, Color.white, mThickness / 2);

        Vector3 pos = transform.position + (transform.up * mRadius + transform.right * mRadius) / Mathf.Sqrt(2);
        DebugExtension.DrawCircle(pos, transform.up - transform.right, Color.white, mThickness / 2);

        pos = transform.position - (transform.up * mRadius + transform.right * mRadius) / Mathf.Sqrt(2);
        DebugExtension.DrawCircle(pos, transform.up - transform.right, Color.white, mThickness / 2);

        pos = transform.position + (transform.up * mRadius - transform.right * mRadius) / Mathf.Sqrt(2);
        DebugExtension.DrawCircle(pos, transform.up + transform.right, Color.white, mThickness / 2);

        pos = transform.position - (transform.up * mRadius - transform.right * mRadius) / Mathf.Sqrt(2);
        DebugExtension.DrawCircle(pos, transform.up + transform.right, Color.white, mThickness / 2);

        DebugExtension.DrawArrow(transform.position, transform.forward, Color.white);

    }


}