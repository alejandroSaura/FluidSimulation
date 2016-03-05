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
    public bool mUseMultithreads;
    public bool bMultithreads
    {
        set { mUseMultithreads = value; }
    }

    float mViscosity;   ///< Viscosity.  Used to compute viscous diffusion.
    public float viscosity
    {        
        set { mViscosity = value; }
    }
    float mFluidDensity;   ///< Uniform density of fluid.
    public float density
    {
        set { mFluidDensity = value; }
    }
    

    List<Vorton> mVortons;   ///< Dynamic array of tiny vortex elements
    public List<Vorton> GetVortons() { return mVortons; }

    NestedGrid<Vorton> mInfluenceTree;   ///< Influence tree
    UniformGrid<Vector> mVelGrid;   ///< Uniform grid of velocity values

    Vector3 mMinCorner;   ///< Minimal corner of axis-aligned bounding box
    Vector3 mMaxCorner;   ///< Maximal corner of axis-aligned bounding box
    
    Vector3 mCirculationInitial;   ///< Initial circulation, which should be conserved when viscosity is zero.
    Vector3 mLinearImpulseInitial;   ///< Initial linear impulse, which should be conserved when viscosity is zero.
    Vector3 mAverageVorticity;   ///< Hack, average vorticity used to compute a kind of viscous vortex diffusion.
    
    float mMassPerParticle;   ///< Mass of each fluid particle (vorton or tracer).
    List<ParticleSystem.Particle> mTracers;   ///< Passive tracer particles


    // vorton simulation constructor
    public VortonSim(float viscosity = 0.0f, float density = 1.0f, bool useMultithreads = true)
    {
        mUseMultithreads = useMultithreads;

        mMinCorner = new Vector3(float.MaxValue, float.MaxValue, float.MaxValue);
        mMaxCorner = (-mMinCorner);
        mViscosity = (viscosity);
        mCirculationInitial = new Vector3(0.0f, 0.0f, 0.0f);
        mLinearImpulseInitial = new Vector3(0.0f, 0.0f, 0.0f);
        mAverageVorticity = new Vector3(0.0f, 0.0f, 0.0f);
        mFluidDensity = (density);
        mMassPerParticle = (0.0f);

        mVortons = new List<Vorton>();
        mInfluenceTree = new NestedGrid<Vorton>();
        mVelGrid = new UniformGrid<Vector>();
        mTracers = new List<ParticleSystem.Particle>();
    }

    /*
        Initialize a vortex particle fluid simulation

        This method assumes the vortons have been initialized.
        That includes removing any vortons embedded inside
        rigid bodies.
    */
    public void Initialize(List<ParticleSystem.Particle> tracers)
    {
        if (mUseMultithreads)
        {
            // Query environment for number of processors on this machine.
            int numberOfProcessors = SystemInfo.processorCount;
        }

        ConservedQuantities();
        ComputeAverageVorticity();

        mTracers = tracers;
        CreateInfluenceTree(); // Create influence tree taking into consideration the tracers and the Vortons
        // Calculate particles physic properties
        {
            float domainVolume = mInfluenceTree[0].GetExtent().x * mInfluenceTree[0].GetExtent().y * mInfluenceTree[0].GetExtent().z;
            if (0.0f == mInfluenceTree[0].GetExtent().z)
            {   // Domain is 2D in XY plane.
                domainVolume = mInfluenceTree[0].GetExtent().x * mInfluenceTree[0].GetExtent().y;
            }
            float totalMass = domainVolume * mFluidDensity;
            //uint numTracersPerCell = (numTracersPerCellCubeRoot * numTracersPerCellCubeRoot * numTracersPerCellCubeRoot);
            float numTracersPerCell = (
                tracers.Count / 
                (float)(mInfluenceTree[0].GetNumCells(0) * mInfluenceTree[0].GetNumCells(1) * mInfluenceTree[0].GetNumCells(2))
                );
            mMassPerParticle = totalMass / (float)(mInfluenceTree[0].GetGridCapacity() * numTracersPerCell);
        }


    }

    /*
        Create nested grid vorticity influence tree.

        Each layer of this tree represents a simplified, aggregated version of
        all of the information in its "child" layer, where each
        "child" has higher resolution than its "parent".

        Derivation:

        Using conservation properties, I_0 = I_0' , I_1 = I_1' , I_2 = I_2'

        I_0 : wx d = w1x d1 + w2x d2
            : wy d = w1y d1 + w2y d2
            : wz d = w1z d1 + w2z d2

        These 3 are not linearly independent:
        I_1 : ( y wz - z wy ) d = ( y1 wz1 - z1 wy1 ) d1 + ( y2 wz2 - z2 wy2 ) d2
            : ( z wx - x wz ) d = ( z1 wx1 - x1 wz1 ) d1 + ( z2 wx2 - x2 wz2 ) d2
            : ( x wy - y wx ) d = ( x1 wy1 - y1 wx1 ) d1 + ( x2 wy2 - y2 wx2 ) d2

        I_2 : ( x^2 + y^2 + z^2 ) wx d = (x1^2 + y1^2 + z1^2 ) wx1 d1 + ( x2^2 + y2^2 + z2^2 ) wx2 d2
            : ( x^2 + y^2 + z^2 ) wy d = (x1^2 + y1^2 + z1^2 ) wy1 d1 + ( x2^2 + y2^2 + z2^2 ) wy2 d2
            : ( x^2 + y^2 + z^2 ) wz d = (x1^2 + y1^2 + z1^2 ) wz1 d1 + ( x2^2 + y2^2 + z2^2 ) wz2 d2

        Can replace I_2 with its magnitude:
              ( x^2  + y^2  + z^2  ) ( wx^2  + wy^2  + wz^2  )^(1/2) d
            = ( x1^2 + y1^2 + z1^2 ) ( wx1^2 + w1y^2 + w1z^2 )^(1/2) d1
            + ( x2^2 + y2^2 + z2^2 ) ( wx2^2 + w2y^2 + w2z^2 )^(1/2) d2
    */
    void CreateInfluenceTree()
    {
        FindBoundingBox(); // Find axis-aligned bounding box that encloses all vortons and tracers.

        // Create skeletal nested grid for influence tree.
        float numElements;

        int numVortons = mVortons.Count;
        numElements = numVortons; // Default approach. Loss of precission when adding tracers... but fast.

        //float boundingBoxVolume = (mMaxCorner - mMinCorner).x * (mMaxCorner - mMinCorner).y * (mMaxCorner - mMinCorner).z;
        //float vortonVolume = 4*(1.0f/3.0f)*Mathf.PI*
        //    (mVortons[0].radius * mVortons[0].radius * mVortons[0].radius);
        // Accurate approach that maintains precission. Very Slow.
        //numElements = boundingBoxVolume / vortonVolume;         

        {
            UniformGrid<Vorton> ugSkeleton = new UniformGrid<Vorton>();   ///< Uniform grid with the same size & shape as the one holding aggregated information about mVortons.
            ugSkeleton.DefineShape((uint)numElements, mMinCorner, mMaxCorner, true);
            mInfluenceTree.Initialize(ugSkeleton); // Create skeleton of influence tree.
        }

        MakeBaseVortonGrid();

        int numLayers = mInfluenceTree.GetDepth();
        for (uint uParentLayer = 1; uParentLayer < numLayers; ++uParentLayer)
        {   // For each layer in the influence tree...
            AggregateClusters(uParentLayer);
        }
    }


    /*
        Create base layer of vorton influence tree.

        This is the leaf layer, where each grid cell corresponds(on average) to
        a single vorton.Some cells might contain multiple vortons and some zero.

        Each cell effectively has a single "supervorton" which its parent layers
        in the influence tree will in turn aggregate.

        This implementation of gridifying the base layer is NOT suitable
        for Eulerian operations like approximating spatial derivatives

        of vorticity or solving a vector Poisson equation, because this
        routine associates each vortex with a single corner point of the
        grid cell that contains it.  To create a grid for Eulerian calculations,
        each vorton would contribute to all 8 corner points of the grid
        cell that contains it.

        We could rewrite this to suit "Eulerian" operations, in which case
        we would want to omit "size" and "position" since the grid would

        implicitly represent that information.  That concern goes hand-in-hand
        with the method used to compute velocity from vorticity.

        Ultimately we need to make sure theoretically conserved quantities behave as expected.

        This method assumes the influence tree skeleton has already been created,
        and the leaf layer initialized to all "zeros", meaning it contains no vortons.
    */
    public void MakeBaseVortonGrid()
    {
        int numVortons = mVortons.Count;

        UniformGrid<VortonClusterAux> ugAux = new UniformGrid<VortonClusterAux>(mInfluenceTree[0]) ; // Temporary auxilliary information used during aggregation.
        ugAux.Init();

        // Compute preliminary vorticity grid.
        for (int uVorton = 0; uVorton < numVortons; ++uVorton)
        {   // For each vorton in this simulation...            
            Vector3 position = mVortons[uVorton].position;
            uint uOffset = mInfluenceTree[0].OffsetOfPosition(position);
                        
            float vortMag = mVortons[uVorton].vorticity.magnitude;

            mInfluenceTree[0][uOffset].position += mVortons[uVorton].position * vortMag; // Compute weighted position -- to be normalized later.
            mInfluenceTree[0][uOffset].vorticity += mVortons[uVorton].vorticity; // Tally vorticity sum.
            mInfluenceTree[0][uOffset].radius = mVortons[uVorton].radius; // Assign volume element size.
            // OBSOLETE. See comments below: UpdateBoundingBox( rVortonAux.mMinCorner , rVortonAux.mMaxCorner , rVorton.mPosition ) ;
            ugAux[uOffset].mVortNormSum += vortMag; // Accumulate vorticity on the VortonClusterAux
        }

        // Post-process preliminary grid (VortonClusterAux); normalize center-of-vorticity and compute sizes, for each grid cell.
        uint[] num = {
            mInfluenceTree[0].GetNumPoints( 0 ) ,
            mInfluenceTree[0].GetNumPoints( 1 ) ,
            mInfluenceTree[0].GetNumPoints( 2 )
        };
        uint numXY = num[0] * num[1];
        uint[] idx = new uint[3];
        for (idx[2] = 0; idx[2] < num[2]; ++idx[2])
        {
           uint zShift = idx[2] * numXY;
            for (idx[1] = 0; idx[1] < num[1]; ++idx[1])
            {
                uint yzShift = idx[1] * num[0] + zShift;
                for (idx[0] = 0; idx[0] < num[0]; ++idx[0])
                {
                    uint offset = idx[0] + yzShift;
                    VortonClusterAux rVortonAux = ugAux[offset];
                    if (rVortonAux.mVortNormSum != float.Epsilon)
                    {   // This cell contains at least one vorton.
                        // Normalize weighted position sum to obtain center-of-vorticity.
                        mInfluenceTree[0][offset].position /= rVortonAux.mVortNormSum;
                    }
                }
            }
        }

    }


    void AggregateClusters(uint uParentLayer)
    {
        // number of cells in each grid cluster
        uint[] pClusterDims = mInfluenceTree.GetDecimations((int)uParentLayer);

        uint[] numCells = { mInfluenceTree[uParentLayer].GetNumCells(0), mInfluenceTree[uParentLayer].GetNumCells(1), mInfluenceTree[uParentLayer].GetNumCells(2) };
        uint numXY = mInfluenceTree[uParentLayer].GetNumPoints(0) * mInfluenceTree[uParentLayer].GetNumPoints(1);
        // (Since this loop writes to each parent cell, it should readily parallelize without contention.)
        uint[] idxParent = new uint[3];
        for (idxParent[2] = 0; idxParent[2] < numCells[2]; ++idxParent[2])
        {
            uint offsetZ = idxParent[2] * numXY;
            for (idxParent[1] = 0; idxParent[1] < numCells[1]; ++idxParent[1])
            {
                uint offsetYZ = idxParent[1] * mInfluenceTree[uParentLayer].GetNumPoints(0) + offsetZ;
                for (idxParent[0] = 0; idxParent[0] < numCells[0]; ++idxParent[0])
                {   // For each cell in the parent layer...
                    uint offsetXYZ = idxParent[0] + offsetYZ;

                    UniformGrid<Vorton> rChildLayer = mInfluenceTree[uParentLayer - 1];
                    VortonClusterAux vortAux = new VortonClusterAux();
                    uint[] clusterMinIndices = mInfluenceTree.GetChildClusterMinCornerIndex(pClusterDims, idxParent); ;
                    
                    uint[] increment = { 0, 0, 0 };
                    uint numXchild = rChildLayer.GetNumPoints(0);
                    uint numXYchild = numXchild * rChildLayer.GetNumPoints(1);
                    // For each cell of child layer in this grid cluster...
                    for (increment[2] = 0; increment[2] < pClusterDims[2]; ++increment[2])
                    {
                        uint childOffsetZ = (clusterMinIndices[2] + increment[2]) * numXYchild;
                        for (increment[1] = 0; increment[1] < pClusterDims[1]; ++increment[1])
                        {
                            uint childOffsetYZ = (clusterMinIndices[1] + increment[1]) * numXchild + childOffsetZ;
                            for (increment[0] = 0; increment[0] < pClusterDims[0]; ++increment[0])
                            {
                                uint childOffsetXYZ = (clusterMinIndices[0] + increment[0]) + childOffsetYZ;
                                Vorton rVortonChild = rChildLayer[childOffsetXYZ];
                                float vortMag = rVortonChild.vorticity.magnitude;

                                // Aggregate vorton cluster from child layer into parent layer:
                                mInfluenceTree[uParentLayer][offsetXYZ].position += rVortonChild.position * vortMag;
                                mInfluenceTree[uParentLayer][offsetXYZ].vorticity += rVortonChild.vorticity;
                                vortAux.mVortNormSum += vortMag;
                                if (rVortonChild.radius != 0.0f)
                                {
                                    mInfluenceTree[uParentLayer][offsetXYZ].radius = rVortonChild.radius;
                                }
                            }
                        }
                    }
                    // Normalize weighted position sum to obtain center-of-vorticity.
                    // (See analogous code in MakeBaseVortonGrid.)
                    mInfluenceTree[uParentLayer][offsetXYZ].position /= vortAux.mVortNormSum;
                }
            }
        }
    }

    /*
        Computes the total circulation and linear impulse of all vortons in this simulation.

        vCirculation - Total circulation, the volume integral of the vorticity.
        vLinearImpulse - Volume integral of circulation weighted by position.
    */
    void ConservedQuantities()
    {
        // Zero accumulators.
        mCirculationInitial = mLinearImpulseInitial = new Vector3(0.0f, 0.0f, 0.0f);
        int numVortons = mVortons.Count;
        for (int iVorton = 0; iVorton < numVortons; ++iVorton)
        {   // For each vorton in this simulation...
            Vorton Vorton = mVortons[iVorton];
            float volumeElement = (mVortons[iVorton].radius * mVortons[iVorton].radius * mVortons[iVorton].radius) * 8.0f;
            // Accumulate total circulation.
            mCirculationInitial += mVortons[iVorton].vorticity * volumeElement;
            // Accumulate total linear impulse.
            mLinearImpulseInitial += Vector3.Cross(mVortons[iVorton].position, mVortons[iVorton].vorticity * volumeElement);
        }
    }

    /*
        Computes the average vorticity of all the vortons in this simulation.
        This is use to compute a hacky, non-physical aproximation to viscous vortex diffusion.
    */
    void ComputeAverageVorticity()
    {
        mAverageVorticity = Vector3.zero;
        int numVortons = mVortons.Count;
        for (int iVorton = 0; iVorton < numVortons; ++iVorton)
        {   // For each vorton in this simulation...
            mAverageVorticity += mVortons[iVorton].vorticity;
        }
        mAverageVorticity /= (float)numVortons;
    }

    public void AddVorton(Vorton vorton)
    {
        mVortons.Add(vorton);
    }

    public void Update(float timeStep)
    {
        //Debug.Log("Updating Vorton Simulation");
        CreateInfluenceTree();
    }

    public void Clear()
    {
        mVortons.Clear();
        mInfluenceTree.Clear();
        mVelGrid.Clear();
        mTracers.Clear();
    }

    void FindBoundingBox()
    {
        int numVortons = mVortons.Count;
        mMinCorner.x = mMinCorner.y = mMinCorner.z = float.MaxValue;
        mMaxCorner = -mMinCorner;

        for (int iVorton = 0; iVorton < numVortons; ++iVorton)
        {   // For each vorton in this simulation...
            // Find corners of axis-aligned bounding box.
            UpdateBoundingBox(ref mMinCorner, ref mMaxCorner, mVortons[iVorton].position);
        }

        int numTracers = mTracers.Count;
        for (int iTracer = 0; iTracer < numTracers; ++iTracer)
        {   // For each passive tracer particle in this simulation...
            
            // Find corners of axis-aligned bounding box.
            UpdateBoundingBox(ref mMinCorner, ref mMaxCorner, mTracers[iTracer].position);
        }

        // Slightly enlarge bounding box to allow for round-off errors.
        Vector3 extent = (mMaxCorner - mMinCorner) ;
        Vector3 nudge = (extent * global.GlobalVar.FLT_EPSILON) ;
        mMinCorner -= nudge;
        mMaxCorner += nudge;
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

    public void DebugDrawVortons()
    {
        for (int i = 0; i < mVortons.Count; ++i)
        {
            Gizmos.color = new Color(0.8f, 0.8f, 1, 0.5f);
            Gizmos.DrawSphere(mVortons[i].position, mVortons[i].radius);
            Gizmos.color = Color.red;
            Gizmos.DrawLine(mVortons[i].position, mVortons[i].position + mVortons[i].vorticity.normalized * mVortons[i].radius);
        }
    }

    public void DebugDrawInfluenceGrid()
    {        
        if (mInfluenceTree != null && mInfluenceTree.mLayers.Count != 0)
        {
            for (int u = 0; u < mInfluenceTree.mLayers.Count; ++u)
            {
                //if (u != 0) continue;
                UniformGrid<Vorton> grid = mInfluenceTree.mLayers[u];

                //Gizmos.color = new Vector4(1, 1 - 1 / ((float)u + 1), 1 - 1 / ((float)u + 1), 1.1f - 1 / ((float)u + 1));

                Vector3 cellExtent = grid.GetCellExtent();
                Vector3 gridExtent = grid.GetExtent();
                Vector3 numCells = new Vector3(grid.GetNumCells(0), grid.GetNumCells(1), grid.GetNumCells(2));
                Vector3 gridOrigin = grid.GetMinCorner();

                for (int i = 0; i < numCells.x; ++i)
                {
                    for (int j = 0; j < numCells.y; ++j)
                    {
                        for (int k = 0; k < numCells.z; ++k)
                        {
                            if (mInfluenceTree[0][0] == null) break;

                            uint[] indices = { (uint)i, (uint)j, (uint)k };
                            uint offset = mInfluenceTree[(uint)u].OffsetFromIndices(indices);
                            float vorticity = mInfluenceTree[(uint)u][offset].vorticity.magnitude;

                            Gizmos.color = new Vector4(vorticity, u/mInfluenceTree.mLayers.Count, u/mInfluenceTree.mLayers.Count, vorticity);
                            Gizmos.DrawWireCube(gridOrigin + new Vector3(cellExtent.x * i + cellExtent.x / 2, cellExtent.y * j + cellExtent.y / 2, cellExtent.z * k + cellExtent.z / 2), cellExtent);
                            
                        }
                    }
                }
            }

        }
    }
    

}
