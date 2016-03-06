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
    int numberOfProcessors;

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
    
    //TO-DO: make it private
    public List<Vorton> mVortons;   ///< Dynamic array of tiny vortex elements
    public List<Vorton> GetVortons() { return mVortons; }

    NestedGrid<Vorton> mInfluenceTree;   ///< Influence tree
    UniformGrid<Vector> mVelGrid;   ///< Uniform grid of velocity values

    Vector3 mMinCorner;   ///< Minimal corner of axis-aligned bounding box
    Vector3 mMaxCorner;   ///< Maximal corner of axis-aligned bounding box
    
    Vector3 mCirculationInitial;   ///< Initial circulation, which should be conserved when viscosity is zero.
    Vector3 mLinearImpulseInitial;   ///< Initial linear impulse, which should be conserved when viscosity is zero.
    Vector3 mAverageVorticity;   ///< Hack, average vorticity used to compute a kind of viscous vortex diffusion.
    
    float mMassPerParticle;   ///< Mass of each fluid particle (vorton or tracer).
    public List<ParticleSystem.Particle> mTracers;   ///< Passive tracer particles


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
    public void Initialize()
    {
        if (mUseMultithreads)
        {
            // Query environment for number of processors on this machine.
            numberOfProcessors = SystemInfo.processorCount;
        }

        ConservedQuantities();
        ComputeAverageVorticity();

        
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
                mTracers.Count / 
                (float)(mInfluenceTree[0].GetNumCells(0) * mInfluenceTree[0].GetNumCells(1) * mInfluenceTree[0].GetNumCells(2))
                );
            mMassPerParticle = totalMass / (float)(mInfluenceTree[0].GetGridCapacity() * numTracersPerCell);
        }


    }

    public void Update(float timeStep)
    {
        //Debug.Log("Updating Vorton Simulation");
        CreateInfluenceTree();

        ComputeVelocityGrid();

        StretchAndTiltVortons(timeStep);

        DiffuseVorticityPSE(timeStep);

        AdvectVortons(timeStep);

        AdvectTracers(timeStep);
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
        //float vortonVolume = 4 * (1.0f / 3.0f) * Mathf.PI *
        //    (mVortons[0].radius * mVortons[0].radius * mVortons[0].radius);
        ////Accurate approach that maintains precission.Very Slow.
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
        Compute velocity due to vortons, for every point in a uniform grid
        This routine assumes CreateInfluenceTree has already executed.
    */
    void ComputeVelocityGrid()
    {
        mVelGrid.Clear();                                  // Clear any stale velocity information
        mVelGrid.CopyShape(mInfluenceTree[0]);           // Use same shape as base vorticity grid. (Note: could differ if you want.)
        mVelGrid.Init();                                   // Reserve memory for velocity grid.

        uint numZ = mVelGrid.GetNumPoints(2);

        if (mUseMultithreads)
        {
            // Estimate grain size based on size of problem and number of processors.
            //const unsigned grainSize = MAX2(1, numZ / gNumberOfProcessors);
            // Compute velocity grid using multiple threads.
            //parallel_for(tbb::blocked_range<size_t>(0, numZ, grainSize), VortonSim_ComputeVelocityGrid_TBB(this));
        }
        else
        {
            ComputeVelocityGridSlice(0, numZ);
        }
        
    }


    /*
        Compute velocity due to vortons, for a subset of points in a uniform grid
        izStart - starting value for z index
        izEnd - ending value for z index   

        This routine assumes CreateInfluenceTree has already executed,
        and that the velocity grid has been allocated.
    */
    void ComputeVelocityGridSlice(uint izStart, uint izEnd)
    {
        int numLayers = mInfluenceTree.GetDepth();

        Vector3 vMinCorner = mVelGrid.GetMinCorner();
        float nudge = 1.0f - 2.0f * global.GlobalVar.FLT_EPSILON;
        Vector3 vSpacing = mVelGrid.GetCellSpacing() * nudge;
        uint[] dims =   { mVelGrid.GetNumPoints( 0 )
                        , mVelGrid.GetNumPoints( 1 )
                        , mVelGrid.GetNumPoints( 2 ) };
        uint numXY = dims[0] * dims[1];
        uint[] idx = new uint[3];

        for (idx[2] = izStart; idx[2] < izEnd; ++idx[2])
        {   // For subset of z index values...
            Vector3 vPosition;
            // Compute the z-coordinate of the world-space position of this gridpoint.
            vPosition.z = vMinCorner.z + (float)(idx[2]) * vSpacing.z;
            // Precompute the z contribution to the offset into the velocity grid.
            uint offsetZ = idx[2] * numXY;
            for (idx[1] = 0; idx[1] < dims[1]; ++idx[1])
            {   // For every gridpoint along the y-axis...
                // Compute the y-coordinate of the world-space position of this gridpoint.
                vPosition.y = vMinCorner.y + (float)(idx[1]) * vSpacing.y;
                // Precompute the y contribution to the offset into the velocity grid.
                uint offsetYZ = idx[1] * dims[0] + offsetZ;
                for (idx[0] = 0; idx[0] < dims[0]; ++idx[0])
                {   // For every gridpoint along the x-axis...
                    // Compute the x-coordinate of the world-space position of this gridpoint.
                    vPosition.x = vMinCorner.x + (float)(idx[0]) * vSpacing.x;
                    // Compute the offset into the velocity grid.
                    uint offsetXYZ = idx[0] + offsetYZ;

                    // Compute the fluid flow velocity at this gridpoint, due to all vortons.
                    uint[] zeros = { 0 , 0 , 0 } ; // Starter indices for recursive algorithm
                    mVelGrid[ offsetXYZ ] = ComputeVelocity( vPosition , zeros , numLayers - 1  ) ;
                }
            }
        }
    }


    /*
        Compute velocity at a given point in space, due to influence of vortons

        vPosition - point in space whose velocity to evaluate
        indices - indices of cell to visit in the given layer
        iLayer - which layer to process

        returns velocity at vPosition, due to influence of vortons
        This is a recursive algorithm with time complexity O(log(N)). 
        The outermost caller should pass in mInfluenceTree.GetDepth().

    */
    Vector ComputeVelocity(Vector3 vPosition , uint[] indices, int iLayer )
    {
        Vector velocityAccumulator = new Vector();

        UniformGrid<Vorton> rChildLayer = mInfluenceTree[(uint)iLayer - 1];

        uint[] pClusterDims = mInfluenceTree.GetDecimations(iLayer);
        uint[] clusterMinIndices =  mInfluenceTree.GetChildClusterMinCornerIndex(pClusterDims, indices);

        Vector3 vGridMinCorner = rChildLayer.GetMinCorner();
        Vector3 vSpacing = rChildLayer.GetCellSpacing();
        uint[] increment = new uint[3];
        uint numXchild = rChildLayer.GetNumPoints(0);
        uint numXYchild = numXchild * rChildLayer.GetNumPoints(1);

        // The larger this is, the more accurate (and slower) the evaluation.
        // Reasonable values lie in [0.00001,4.0].
        // Setting this to 0 leads to very bad errors, but values greater than (tiny) lead to drastic improvements.
        // Changes in margin have a quantized effect since they effectively indicate how many additional
        // cluster subdivisions to visit.
        float marginFactor = 0.0001f; // 0.4f ; // ship with this number: 0.0001f ; test with 0.4
                                                   // When domain is 2D in XY plane, min.z==max.z so vPos.z test below would fail unless margin.z!=0.
        Vector3 margin = marginFactor * vSpacing + (0.0f == vSpacing.z ? new Vector3(0, 0, float.Epsilon) : Vector3.zero);

        // For each cell of child layer in this grid cluster...
        for (increment[2] = 0; increment[2] < pClusterDims[2]; ++increment[2])
        {
            uint[] idxChild = new uint[3];
            idxChild[2] = clusterMinIndices[2] + increment[2];
            Vector3 vCellMinCorner, vCellMaxCorner;
            vCellMinCorner.z = vGridMinCorner.z + (float)(idxChild[2]) * vSpacing.z;
            vCellMaxCorner.z = vGridMinCorner.z + (float)(idxChild[2] + 1) * vSpacing.z;
            uint offsetZ = idxChild[2] * numXYchild;
            for (increment[1] = 0; increment[1] < pClusterDims[1]; ++increment[1])
            {
                idxChild[1] = clusterMinIndices[1] + increment[1];
                vCellMinCorner.y = vGridMinCorner.y + (float)(idxChild[1]) * vSpacing.y;
                vCellMaxCorner.y = vGridMinCorner.y + (float)(idxChild[1] + 1) * vSpacing.y;
                uint offsetYZ = idxChild[1] * numXchild + offsetZ;
                for (increment[0] = 0; increment[0] < pClusterDims[0]; ++increment[0])
                {
                    idxChild[0] = clusterMinIndices[0] + increment[0];
                    vCellMinCorner.x = vGridMinCorner.x + (float)(idxChild[0]) * vSpacing.x;
                    vCellMaxCorner.x = vGridMinCorner.x + (float)(idxChild[0] + 1) * vSpacing.x;
                    if (
                            (iLayer > 1)
                        && (vPosition.x >= vCellMinCorner.x - margin.x)
                        && (vPosition.y >= vCellMinCorner.y - margin.y)
                        && (vPosition.z >= vCellMinCorner.z - margin.z)
                        && (vPosition.x < vCellMaxCorner.x + margin.x)
                        && (vPosition.y < vCellMaxCorner.y + margin.y)
                        && (vPosition.z < vCellMaxCorner.z + margin.z)
                      )
                    {   // Test position is inside childCell and currentLayer > 0...
                        // Recurse child layer.
                        Vector upVelocicy = ComputeVelocity(vPosition, idxChild, iLayer - 1);
                        velocityAccumulator += upVelocicy;
                    }
                    else
                    {   // Test position is outside childCell, or reached leaf node.
                        //    Compute velocity induced by cell at corner point x.
                        //    Accumulate influence, storing in velocityAccumulator.
                        uint offsetXYZ = idxChild[0] + offsetYZ;
                        
                        // Add velocity due to this vorton to the accumulator
                        rChildLayer[offsetXYZ].AccumulateVelocity(ref velocityAccumulator, vPosition);
                    }
                }
            }
        }


        return velocityAccumulator;
    }


    /*
        Stretch and tilt vortons using velocity field
        timeStep - amount of time by which to advance simulation
        uFrame - frame counter

        see J. T. Beale, A convergent three-dimensional vortex method with
                grid-free stretching, Math. Comp. 46 (1986), 401-24, April.

        This routine assumes CreateInfluenceTree has already executed.

    */
    void StretchAndTiltVortons(float timeStep)
    {
        if ((0.0f == mVelGrid.GetExtent().x)
        || (0.0f == mVelGrid.GetExtent().y)
        || (0.0f == mVelGrid.GetExtent().z))
        {   // Domain is 2D, so stretching & tilting does not occur.
            return;
        }

        // Compute all gradients of all components of velocity.
        UniformGrid<Matrix3x3> velocityJacobianGrid = new UniformGrid<Matrix3x3>(mVelGrid);
        velocityJacobianGrid.Init();

        UniformGridMath.ComputeJacobian(ref velocityJacobianGrid, mVelGrid);

        int numVortons = mVortons.Count;

        for (int offset = 0; offset < numVortons; ++offset)
        {   // For each vorton...
            Matrix3x3 velJac = (Matrix3x3) velocityJacobianGrid.Interpolate(mVortons[offset].position);
            Vector3 stretchTilt = mVortons[offset].vorticity * velJac;    // Usual way to compute stretching & tilting
            mVortons[offset].vorticity += /* fudge factor for stability */ 0.5f * stretchTilt * timeStep;
        }

    }


    /*
        Diffuse vorticity using a particle strength exchange method.

        This routine partitions space into cells using the same grid
        as the "base vorton" grid.  Each vorton gets assigned to the
        cell that contains it.  Then, each vorton exchanges some
        of its vorticity with its neighbors in adjacent cells.

        This routine makes some simplifying assumptions to speed execution:

            -   Distance does not influence the amount of vorticity exchanged,
                except in as much as only vortons within a certain region of
                each other exchange vorticity.  This amounts to saying our kernel,
                eta, is a top-hat function.

            -   Theoretically, if an adjacent cell contains no vortons
                then this simulation should generate vorticity within
                that cell, e.g. by creating a new vorton in the adjacent cell.

            -   This simulation reduces the vorticity of each vorton, alleging
                that this vorticity is dissipated analogously to how energy
                dissipates at Kolmogorov microscales.  This treatment is not
                realistic but it retains qualitative characteristics that we
                want, e.g. that the flow dissipates at a rate related to viscosity.
                Dissipation in real flows is a more complicated phenomenon.

        see Degond & Mas-Gallic (1989): The weighted particle method for
            convection-diffusion equations, part 1: the case of an isotropic viscosity.
            Math. Comput., v. 53, n. 188, pp. 485-507, October.

        timeStep - amount of time by which to advance simulation      

        This routine assumes CreateInfluenceTree has already executed.

    */
    void DiffuseVorticityPSE( float timeStep)
    {
        // Phase 1: Partition vortons

        // Create a spatial partition for the vortons.
        // Each cell contains a dynamic array of integers
        // whose values are offsets into mVortons.
        UniformGrid<IntList> ugVortRef = new UniformGrid<IntList>(mInfluenceTree[0] );
        ugVortRef.Init() ;

        int numVortons = mVortons.Count;

        for(int offset = 0 /* Start at 0th vorton */ ; offset<numVortons ; ++ offset )
        {   // For each vorton...            
            // Insert the vorton's offset into the spatial partition.
            ugVortRef[mVortons[offset].position].list.Add(offset);
        }

        // Phase 2: Exchange vorticity with nearest neighbors

        uint nx = ugVortRef.GetNumPoints( 0 ) ;
        uint nxm1 = nx - 1;
        uint ny = ugVortRef.GetNumPoints( 1 ) ;
        uint nym1 = ny - 1;
        uint nxy = nx * ny;
        uint nz = ugVortRef.GetNumPoints( 2 ) ;
        uint nzm1 = nz - 1;

        uint[] idx = new uint[3];
        for( idx[2] = 0 ; idx[2] < nzm1 ; ++ idx[2] )
        {   // For all points along z except the last...

            uint offsetZ0 = idx[2] * nxy;
            uint offsetZp = (idx[2] + 1) * nxy;

            for( idx[1] = 0 ; idx[1] < nym1 ; ++ idx[1] )
            {   // For all points along y except the last...

                uint offsetY0Z0 = idx[1] * nx + offsetZ0;
                uint offsetYpZ0 = (idx[1] + 1) * nx + offsetZ0;
                uint offsetY0Zp = idx[1] * nx + offsetZp;

                for( idx[0] = 0 ; idx[0] < nxm1 ; ++ idx[0] )
                {   // For all points along x except the last...

                    uint offsetX0Y0Z0 = idx[0] + offsetY0Z0;

                    for( int ivHere = 0; ivHere<ugVortRef[offsetX0Y0Z0].list.Count ; ++ ivHere )
                    {   // For each vorton in this gridcell...

                        int rVortIdxHere  = ugVortRef[offsetX0Y0Z0].list[ivHere] ;
                        Vorton rVortonHere = mVortons[rVortIdxHere] ;
                        Vector3 rVorticityHere  = rVortonHere.vorticity ;

                        // Diffuse vorticity with other vortons in this same cell:
                        for( int ivThere = ivHere + 1; ivThere<ugVortRef[offsetX0Y0Z0].list.Count ; ++ ivThere )
                        {   // For each OTHER vorton within this same cell...

                            int rVortIdxThere = ugVortRef[offsetX0Y0Z0].list[ivThere] ;
                            Vorton rVortonThere    = mVortons[rVortIdxThere] ;
                            Vector3 rVorticityThere = rVortonThere.vorticity ;

                            Vector3 vortDiff = rVorticityHere - rVorticityThere;
                            Vector3 exchange = 2.0f * mViscosity * timeStep * vortDiff;    // Amount of vorticity to exchange between particles.

                            mVortons[rVortIdxHere].vorticity -= exchange ;   // Make "here" vorticity a little closer to "there".
                            mVortons[rVortIdxThere].vorticity += exchange ;   // Make "there" vorticity a little closer to "here".
                        }

                        // Diffuse vorticity with vortons in adjacent cells:
                        {
                            uint offsetXpY0Z0 = idx[0] + 1 + offsetY0Z0; // offset of adjacent cell in +X direction
                            for( int ivThere = 0; ivThere<ugVortRef[offsetXpY0Z0].list.Count ; ++ ivThere )
                            {   // For each vorton in the adjacent cell in +X direction...

                                int rVortIdxThere = ugVortRef[offsetXpY0Z0].list[ivThere] ;
                                Vorton rVortonThere = mVortons[rVortIdxThere] ;
                                Vector3  rVorticityThere = rVortonThere.vorticity ;

                                Vector3 vortDiff = rVorticityHere - rVorticityThere;
                                Vector3 exchange = mViscosity * timeStep * vortDiff;    // Amount of vorticity to exchange between particles.

                                mVortons[rVortIdxHere].vorticity -= exchange ;   // Make "here" vorticity a little closer to "there".
                                mVortons[rVortIdxThere].vorticity += exchange ;   // Make "there" vorticity a little closer to "here".
                            }
                        }

                        {
                            uint offsetX0YpZ0 = idx[0] + offsetYpZ0; // offset of adjacent cell in +Y direction
                            for( int ivThere = 0; ivThere<ugVortRef[offsetX0YpZ0].list.Count ; ++ ivThere )
                            {   // For each vorton in the adjacent cell in +Y direction...
                                int  rVortIdxThere = ugVortRef[offsetX0YpZ0].list[ivThere] ;
                                Vorton rVortonThere = mVortons[rVortIdxThere] ;
                                Vector3 rVorticityThere = rVortonThere.vorticity ;

                                Vector3 vortDiff = rVorticityHere - rVorticityThere;
                                Vector3 exchange = mViscosity * timeStep * vortDiff;    // Amount of vorticity to exchange between particles.

                                mVortons[rVortIdxHere].vorticity -= exchange ;   // Make "here" vorticity a little closer to "there".
                                mVortons[rVortIdxThere].vorticity += exchange ;   // Make "there" vorticity a little closer to "here".
                            }
                        }

                        {
                            uint offsetX0Y0Zp = idx[0] + offsetY0Zp; // offset of adjacent cell in +Z direction
                            for( int ivThere = 0; ivThere<ugVortRef[offsetX0Y0Zp].list.Count ; ++ ivThere )
                            {   // For each vorton in the adjacent cell in +Z direction...
                                int rVortIdxThere = ugVortRef[offsetX0Y0Zp].list[ivThere] ;
                                Vorton rVortonThere = mVortons[rVortIdxThere] ;
                                Vector3 rVorticityThere = rVortonThere.vorticity ;

                                Vector3 vortDiff = rVorticityHere - rVorticityThere;
                                Vector3 exchange = mViscosity * timeStep * vortDiff;    // Amount of vorticity to exchange between particles.

                                mVortons[rVortIdxHere].vorticity -= exchange ;   // Make "here" vorticity a little closer to "there".
                                mVortons[rVortIdxThere].vorticity += exchange ;   // Make "there" vorticity a little closer to "here".
                            }
                        }

                        // Dissipate vorticity.  See notes in header comment.
                        mVortons[rVortIdxHere].vorticity -= mViscosity* timeStep * mVortons[rVortIdxHere].vorticity;   // Reduce "here" vorticity.
                    }
                }
            }
        }
    }


    /*
        Advect vortons using velocity field

        timeStep - amount of time by which to advance simulation        

    */
    void AdvectVortons( float timeStep )
    {
        int numVortons = mVortons.Count;

        for( int offset = 0; offset<numVortons ; ++ offset )
        {   // For each vorton...
            Vorton rVorton = mVortons[offset] ;

            Vector3 velocity = ((Vector) mVelGrid.Interpolate(rVorton.position)).v;
            mVortons[offset].position += velocity * timeStep;
            mVortons[offset].velocity = velocity ;  // Cache this for use in collisions with rigid bodies.
        }
    }

    /*
        Advect passive tracers using velocity field

        timeStep - amount of time by which to advance simulation        

    */
    void AdvectTracers( float timeStep)
    {
        int numTracers = mTracers.Count;

        if (mUseMultithreads)
        {
            // Estimate grain size based on size of problem and number of processors.
            int grainSize = Mathf.Max(1, numTracers / numberOfProcessors);
            // Advect tracers using multiple threads.
            //parallel_for(tbb::blocked_range<size_t>(0, numTracers, grainSize), VortonSim_AdvectTracers_TBB(this, timeStep, uFrame));
        }
        else
        {
            AdvectTracersSlice(timeStep, 0, numTracers);
        }
    
    }

    /*
        Advect (subset of) passive tracers using velocity field

        timeStep - amount of time by which to advance simulation
        itStart - index of first tracer to advect
        itEnd - index of last tracer to advect

    */
    void AdvectTracersSlice( float timeStep , int itStart, int itEnd )
    {
        for( int offset = itStart; offset<itEnd ; ++ offset )
        {   // For each passive tracer in this slice...
            Vector3 velocity = ((Vector)mVelGrid.Interpolate(mTracers[offset].position)).v;

            ParticleSystem.Particle newTracer = new ParticleSystem.Particle();
            newTracer.position = mTracers[offset].position + velocity * timeStep;
            newTracer.velocity = velocity ; // Cache for use in collisions

            newTracer.startColor = new Color32(255, 255, 255, 50);
            newTracer.startSize = 0.5f;
            newTracer.lifetime = 9999;

            mTracers[offset] = newTracer;
        }
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

    public void DebugDrawVelocityGrid()
    {
        if (mVelGrid != null)
        {
            Vector3 cellExtent = mVelGrid.GetCellExtent();
            Vector3 gridExtent = mVelGrid.GetExtent();
            Vector3 numCells = new Vector3(mVelGrid.GetNumCells(0), mVelGrid.GetNumCells(1), mVelGrid.GetNumCells(2));
            Vector3 gridOrigin = mVelGrid.GetMinCorner();

            for (int i = 0; i < numCells.x; ++i)
            {
                for (int j = 0; j < numCells.y; ++j)
                {
                    for (int k = 0; k < numCells.z; ++k)
                    {
                        uint[] indices = { (uint)i, (uint)j, (uint)k };
                        uint offset = mVelGrid.OffsetFromIndices(indices);
                        Vector3 center = mVelGrid.PositionFromIndices(indices) + cellExtent/2;
                        if (mVelGrid[offset].v.magnitude == 0) break;

                        Color color = new Vector4(0.8f, 0.8f, 1, 1);
                        DebugExtension.DrawArrow(center, mVelGrid[offset].v/8, color);
                    }
                }
            }   
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
