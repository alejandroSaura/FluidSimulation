using UnityEngine;
using System.Collections.Generic;

[RequireComponent(typeof(ParticleSystem))]
public class FluidSimulation : MonoBehaviour
{
    public bool debugVortons = true;
    public bool debugInfluenceTree = true;

    public bool useMultiThreads = true;

    public float viscosity = 0.05f;
    public float density = 1.0f;
    
    public float magnitude = 20.0f;
    public int numCellsPerDim = 16;    
    public uint numTracersPer = 3;
    int numVortonsMax;

    public bool BOUNDARY_NO_SLIP_NO_THRU = true;
    // BOUNDARY_RESPECTS_AMBIENT_FLOW only takes effect if BOUNDARY_NO_SLIP_NO_THRU is true
    public bool BOUNDARY_RESPECTS_AMBIENT_FLOW = true;
    // BOUNDARY_AMBIENT_FLOW_OMITS_VORTON_OLD_POSITION only takes effect if BOUNDARY_RESPECTS_AMBIENT_FLOW is true.
    public bool BOUNDARY_AMBIENT_FLOW_OMITS_VORTON_OLD_POSITION = true;

    public bool FLOW_AFFECTS_BODY = true;

    VortonSim mVortonSim;
    //IVorticityDistribution vorticityDistribution;

	void Start ()
    {
        Random.seed = (int)(Time.realtimeSinceStartup*1000);

        numVortonsMax = numCellsPerDim * numCellsPerDim * numCellsPerDim;

        mVortonSim = new VortonSim(viscosity, density, useMultiThreads);
        //mVortonSim.Clear();

        //IVorticityDistribution vorticityDistribution = gameObject.GetComponentInChildren<IVorticityDistribution>();
        //AssignVorticity(2.0f * magnitude, (uint)numVortonsMax, vorticityDistribution);

        // distribute numVortonsMax between all vortices of the simulation.
        GameObject[] vortices = GameObject.FindGameObjectsWithTag("Vortex");
        for(int i = 0; i < vortices.Length; ++i)
        {
            AssignVorticity(2.0f * magnitude, (uint)(numVortonsMax/vortices.Length), vortices[i].GetComponent<IVorticityDistribution>());
        }        

        Initialize(numTracersPer);
    }    

    /*
        Initialize a fluid-and-body simulation

        This routine will remove particles that are embedded
        inside rigid bodies.

        This method assumes the vortons have been initialized,
        i.e. that their initial positions, vorticities and
        radius have all been set.
    */
    void Initialize(uint numTracersPerCellCubeRoot)
    {
        //RemoveEmbeddedParticles();
        mVortonSim.Initialize(numTracersPerCellCubeRoot);
        //RemoveEmbeddedParticles();
    }

    // to handle changes in editor
    void OnValidate()
    {
        if (mVortonSim != null)
        {
            mVortonSim.bMultithreads = useMultiThreads;
            mVortonSim.density = density;
            mVortonSim.viscosity = viscosity;
        }
    }

    /*
        Update the fluid and solve boundary conditions
        timeStep - change in virtual time since last update        
    */
    void Update()
    {
        // Update fluid, temporarily ignoring rigid bodies and boundary conditions.
        mVortonSim.Update(Time.deltaTime);

        // Apply boundary conditions and calculate impulses to apply to rigid bodies.
        //SolveBoundaryConditions();        
    }

    void AssignVorticity(float fMagnitude, uint numVortonsMax, IVorticityDistribution vorticityDistribution)
    {
        Vector3 vDimensions = vorticityDistribution.GetDomainSize(); // length of each side of grid box
        Vector3 vCenter = (vorticityDistribution).GetCenter(); // Center of vorticity distribution
        Vector3 vMin = (vCenter - 0.5f * vDimensions) ; // Minimum corner of box containing vortons
        Vector3 vMax = (vMin + vDimensions) ; // Maximum corner of box containing vortons
        UniformGridGeometry skeleton = new UniformGridGeometry(numVortonsMax, vMin, vMax, true ) ;
        // number of grid cells in each direction of virtual uniform grid
        int[] numCells = {   (int)Mathf.Max( 1 , skeleton.GetNumCells(0))
                         ,   (int)Mathf.Max( 1 , skeleton.GetNumCells(1))
                         ,   (int)Mathf.Max( 1 , skeleton.GetNumCells(2)) };

        // Total number of cells should be as close to numVortonsMax as possible without going over.
        // Worst case allowable difference would be numVortonsMax=7 and numCells in each direction is 1 which yields a ratio of 1/7.
        // But in typical situations, the user would like expect total number of virtual cells to be closer to numVortonsMax than that.
        // E.g. if numVortonsMax=8^3=512 somehow yielded numCells[0]=numCells[1]=numCells[2]=7 then the ratio would be 343/512~=0.67.
        while (numCells[0] * numCells[1] * numCells[2] > numVortonsMax)
        {   // Number of cells is excessive.
            // This can happen when the trial number of cells in any direction is less than 1 -- then the other two will likely be too large.
            numCells[0] = (int)Mathf.Max(1, numCells[0] / 2);
            numCells[1] = (int)Mathf.Max(1, numCells[1] / 2);
            numCells[2] = (int)Mathf.Max(1, numCells[2] / 2);
        }

        float[] oneOverN = { 1.0f / (float)(numCells[0]), 1.0f / (float)(numCells[1]), 1.0f / (float)(numCells[2]) };
        Vector3 gridCellSize = new Vector3(vDimensions.x * oneOverN[0] , vDimensions.y * oneOverN[1], vDimensions.z * oneOverN[2] ) ;
        float vortonRadius = Mathf.Pow(gridCellSize.x * gridCellSize.y * gridCellSize.z, 1.0f / 3.0f) * 0.5f;
        if (0.0f == vDimensions.z)
        {   // z size is zero, so domain is 2D.
            vortonRadius = Mathf.Pow(gridCellSize.x * gridCellSize.y, 0.5f) * 0.5f;
        }
        Vector3 vNoise = (0.0f * gridCellSize) ;

        //-----------------------------------------------------------------------
        // Iterate through each point in a uniform grid.
        // If probe position is inside vortex core, add a vorton there.
        // This loop could be rewritten such that it only visits points inside the core,
        // but this loop structure can readily be reused for a wide variety of configurations.
        Vector3 position = new Vector3(0.0f, 0.0f, 0.0f); // vorton position
        int[] index = new int[3]; // index of each position visited       

        for (index[2] = 0; index[2] < numCells[2]; ++index[2])
        {   // For each z-coordinate...
            position.z = ((float)(index[2]) + 0.25f) * gridCellSize.z + vMin.z;
            for (index[1] = 0; index[1] < numCells[1]; ++index[1])
            {   // For each y-coordinate...
                position.y = ((float)(index[1]) + 0.25f) * gridCellSize.y + vMin.y;
                for (index[0] = 0; index[0] < numCells[0]; ++index[0])
                {   // For each x-coordinate...
                    position.x = ((float)(index[0]) + 0.25f) * gridCellSize.x + vMin.x;
                    position += Random.Range(-1, 1) * (vNoise);
                    Vector3 vorticity = Vector3.zero;
                    vorticityDistribution.AssignVorticity(ref vorticity, position, vCenter);
                    Vorton vorton = new Vorton(position, vorticity * fMagnitude, vortonRadius);
                    if (vorticity.sqrMagnitude > global.GlobalVar.sTiny)
                    {   // Vorticity is significantly non-zero.
                        mVortonSim.AddVorton(vorton);                        
                    }
                }
            }
        }
        //-----------------------------------------------------------------------
    }

    void OnDrawGizmos()
    {
        if (mVortonSim != null)
        {
            if(debugInfluenceTree) mVortonSim.DebugDrawInfluenceGrid();
            if(debugVortons) mVortonSim.DebugDrawVortons();
        }        
    }

}
