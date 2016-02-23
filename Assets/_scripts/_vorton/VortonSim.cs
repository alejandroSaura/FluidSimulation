using UnityEngine;
using System.Collections.Generic;
using System;

/*
    Dynamic simulation of a fluid, using tiny vortex elements.

    This implements a portion of a fluid simulation, and effectively
    neglects boundary conditions.  This module defers the enforcement
    of boundary conditions to another module. (FluidBodySim)

*/
public class VortonSim
{

    List<Vorton> mVortons;   ///< Dynamic array of tiny vortex elements
    NestedGrid<Vorton> mInfluenceTree;   ///< Influence tree
    UniformGrid<Vector> mVelGrid;   ///< Uniform grid of velocity values
    Vector3 mMinCorner;   ///< Minimal corner of axis-aligned bounding box
    Vector3 mMaxCorner;   ///< Maximal corner of axis-aligned bounding box
    float mViscosity;   ///< Viscosity.  Used to compute viscous diffusion.
    Vector3 mCirculationInitial;   ///< Initial circulation, which should be conserved when viscosity is zero.
    Vector3 mLinearImpulseInitial;   ///< Initial linear impulse, which should be conserved when viscosity is zero.
    Vector3 mAverageVorticity;   ///< Hack, average vorticity used to compute a kind of viscous vortex diffusion.
    float mFluidDensity;   ///< Uniform density of fluid.
    float mMassPerParticle;   ///< Mass of each fluid particle (vorton or tracer).
    List<Particle> mTracers;   ///< Passive tracer particles

    // vorton simulation constructor
    public VortonSim(float viscosity = 0.0f, float density = 1.0f)
    {
        mMinCorner = new Vector3(float.MaxValue, float.MaxValue, float.MaxValue);
        mMaxCorner = (-mMinCorner);
        mViscosity = (viscosity);
        mCirculationInitial = new Vector3(0.0f, 0.0f, 0.0f);
        mLinearImpulseInitial = new Vector3(0.0f, 0.0f, 0.0f);
        mAverageVorticity = new Vector3(0.0f, 0.0f, 0.0f);
        mFluidDensity = (density);
        mMassPerParticle = (0.0f);
    }

    /*
        Initialize a vortex particle fluid simulation

        This method assumes the vortons have been initialized.
        That includes removing any vortons embedded inside
        rigid bodies.
    */
    public void Initialize(uint numTracersPerCellCubeRoot)
    {
        //if (global.GlobalVar.useMultiThreads)
        //{
        //    // Query environment for number of processors on this machine.
        //    int numberOfProcessors = SystemInfo.processorCount;            
        //}       

        //ConservedQuantities(mCirculationInitial, mLinearImpulseInitial);
        //ComputeAverageVorticity();
        //CreateInfluenceTree(); // This is a marginally superfluous call.  We only need the grid geometry to seed passive tracer particles.

        //InitializePassiveTracers(numTracersPerCellCubeRoot); // Create particles

        //{
        //    float domainVolume = mInfluenceTree[0].GetExtent().x * mInfluenceTree[0].GetExtent().y * mInfluenceTree[0].GetExtent().z;
        //    if (0.0f == mInfluenceTree[0].GetExtent().z)
        //    {   // Domain is 2D in XY plane.
        //        domainVolume = mInfluenceTree[0].GetExtent().x * mInfluenceTree[0].GetExtent().y;
        //    }
        //    const float totalMass = domainVolume * mFluidDensity;
        //    const unsigned numTracersPerCell = POW3(numTracersPerCellCubeRoot);
        //    mMassPerParticle = totalMass / float(mInfluenceTree[0].GetGridCapacity() * numTracersPerCell);
        //}
    }

    /*
        Update axis-aligned bounding box corners to include given point

        vMinCorner - minimal corner of axis-aligned bounding box
        vMaxCorner - maximal corner of axis-aligned bounding box
        vPoint - point to include in bounding box

    */
    public void UpdateBoundingBox(ref Vector3 vMinCorner , ref Vector3 vMaxCorner , Vector3 vPoint )
    {
        vMinCorner.x = Mathf.Min(vPoint.x , vMinCorner.x );
        vMinCorner.y = Mathf.Min(vPoint.y , vMinCorner.y );
        vMinCorner.z = Mathf.Min(vPoint.z , vMinCorner.z );
        vMaxCorner.x = Mathf.Max(vPoint.x , vMaxCorner.x );
        vMaxCorner.y = Mathf.Max(vPoint.y , vMaxCorner.y );
        vMaxCorner.z = Mathf.Max(vPoint.z , vMaxCorner.z );
    }

    /*
        Kill the tracer at the given index
    */
    void KillTracer(int iTracer)
    {
        mTracers.RemoveAt(iTracer);        
    }

    void Clear()
    {
        mVortons.Clear();
        mInfluenceTree.Clear();
        mVelGrid.Clear();
        mTracers.Clear();
    }

}
