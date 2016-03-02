using UnityEngine;
using System.Collections.Generic;

public class FluidSimulation : MonoBehaviour {

    public float viscosity = 0.05f;
    public float density = 1.0f;

    public float fRadius = 1.0f;
    public float fThickness = 1.0f;
    public float fMagnitude = 20.0f;
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

	void Start ()
    {
        numVortonsMax = numCellsPerDim * numCellsPerDim * numCellsPerDim;

        mVortonSim = new VortonSim(viscosity, density);
        mVortonSim.Clear();

        List<Vorton> vortons = mVortonSim.GetVortons();
        //AssignVorticity(vortons, 2.0f * fMagnitude, numVortonsMax, VortexRing(fRadius, fThickness, Vec3(1.0f, 0.0f, 0.0f)));
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

}
