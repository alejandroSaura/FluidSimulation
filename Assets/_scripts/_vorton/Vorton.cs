using UnityEngine;
using System.Collections;

// Vortex particle
public class Vorton : IgridItem     
{

    Vector3 mPosition;   ///< Position (in world units) of center of vortex particle
    public Vector3 position
    {
        get { return mPosition; }
        set { mPosition = value; }
    }
    Vector3 mVorticity;   ///< Vorticity of vortex particle
    public Vector3 vorticity
    {
        get { return mVorticity; }
        set { mVorticity = value; }
    }
    float mRadius;  ///< Radius of vortex particle
    public float radius
    {
        get { return mRadius; }
        set { mRadius = value; }
    }
    Vector3 mVelocity;   ///< Velocity of this vorton -- used to cache value obtained during advection, to optimize collision response.
    public Vector3 velocity
    {
        get { return mVelocity; }
        set { mVelocity = value; }
    }


    /*
        Construct a vortex particle
    */
    public Vorton()
    {
        mPosition = new Vector3(0, 0, 0);
        mVorticity = new Vector3(0, 0, 0);
        mRadius = 0;
        mVelocity = new Vector3(0, 0, 0);
    }

    public Vorton( Vector3 vPos , Vector3 vVort , float fRadius = 0.0f )    
    {
        mPosition = vPos;
        mVorticity = vVort;
        mRadius = fRadius;
        mVelocity = new Vector3(0, 0, 0);
    }

    public Vorton( Vorton that )
    {
        mPosition = that.mPosition;
        mVorticity = that.mVorticity;
        mRadius = that.mRadius;    
    }   

    /*
    Compute velocity induced by this vortex element (a vorton)
    vVelocity - (in/out) variable in which to accumulate velocity
    vPosQuery - position where we want to know velocity
    
    mRadius currently serves double-duty for two things
    which should probably be kept separate.
    One is the radius of the finite-size vorton,
    where the vorticity distribution inside the radius
    is finite, to avoid evaluating a singularity.
    The other is the volume of the "infinitesimal"
    volume element, used to compute a contribution
    to a velocity field.

    */
    public void AccumulateVelocity(ref Vector vVelocity , Vector3 vPosQuery )
    {
        if (mRadius == 0) return;

        Vector3 vNeighborToSelf = vPosQuery - mPosition;
        float radius2 = mRadius * mRadius;
        float dist2 = vNeighborToSelf.sqrMagnitude + global.GlobalVar.sAvoidSingularity;
        float oneOverDist = global.GlobalVar.finvsqrtf(dist2);
        Vector3 vNeighborToSelfDir = vNeighborToSelf * oneOverDist;

        /* If the reciprocal law is used everywhere then when 2 vortices get close, they tend to jettison. */
        /* Mitigate this by using a linear law when 2 vortices get close to each other. */
        float distLaw = (dist2 < radius2)
                ?   /* Inside vortex core - linear law */
                (oneOverDist / radius2)
                :   /* Outside vortex core */
                (oneOverDist / dist2);

        // this is the implementation of the formula:
        // "v = (1/4phi) * (vorticity x distance)/distance^3" - velocity from vorticity
        // with a improvement to mitigate the singularity inside the vorton
        vVelocity.v += Vector3.Cross(global.GlobalVar.OneOverFourPi * (8.0f * radius2 * mRadius) * mVorticity, vNeighborToSelf * distLaw);
    }

    /*
    Compute vorticity required to obtain a given velocity, due to a single vorton at a given position.

    This assigns the vorticity
    w = 4 Pi r^2 v / volumeElement

    where

    r is the distance from the vorton (which here is also the radius of the vorton)
    v is the velocity induced by the vorton
    volumeElement is the volume occupied by the vorton
    w_hat is r_hat cross v_hat.

    This assumes v and r are orthogonal, so this is a very special-purpose
    routine.

    This routine also assumes this vorton's position and radius are where they need to be.

    */
    public void AssignByVelocity(Vector3 vQueryPosition, Vector3 velocity )
    {
        Vector3 vPosRelative = vQueryPosition - mPosition;
        float dist = vPosRelative.magnitude;
        mVorticity = Vector3.Cross(global.GlobalVar.FourPi * dist * vPosRelative , velocity / ( 8.0f * mRadius* mRadius * mRadius )) ;
    }



}
