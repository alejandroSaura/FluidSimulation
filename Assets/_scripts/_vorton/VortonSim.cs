using UnityEngine;
using System.Collections.Generic;
using System;

public class VortonSim : MonoBehaviour {

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

}
