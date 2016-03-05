using UnityEngine;
using System.Collections.Generic;


public class Tracers : MonoBehaviour {

    public uint numberOfTracers = 512;
    UniformGridGeometry grid;

    public List<ParticleSystem.Particle> tracers;
    
    void Awake ()
    {
        Vector3 posMin = transform.position - 0.5f * new Vector3(transform.localScale.x, transform.localScale.y, transform.localScale.z);
        Vector3 posMax = transform.position + 0.5f * new Vector3(transform.localScale.x, transform.localScale.y, transform.localScale.z);

        // Define base grid to instantiate one particle per cell
        grid = new UniformGridGeometry(numberOfTracers, posMin, posMax, false);
        tracers = new List<ParticleSystem.Particle>();              

        CreateTracers();
    }

    void CreateTracers()
    {
        Vector3 numCells = new Vector3(grid.GetNumCells(0), grid.GetNumCells(1), grid.GetNumCells(2));
        Vector3 relativeCenterOfGrid = new Vector3(grid.GetCellExtent().x/2, grid.GetCellExtent().y / 2, grid.GetCellExtent().z / 2);
        for (uint i = 0; i < numCells.x; ++i)
        {
            for (uint j = 0; j < numCells.y; ++j)
            {
                for (uint k = 0; k < numCells.z; ++k)
                {
                    ParticleSystem.Particle newTracer = new ParticleSystem.Particle();
                    uint[] indices = { i, j, k };
                    newTracer.position = (grid.PositionFromIndices(indices) + relativeCenterOfGrid) ;
                    newTracer.startColor = new Color32(255, 255, 255, 255);
                    newTracer.startSize = 0.2f;
                    newTracer.lifetime = 9999;
                    tracers.Add(newTracer);
                }
            }
        }

    }

    void LateUpdate()
    {
        // TO-DO: This will be done in the fluid simulation
        //particleSystem.SetParticles(tracers.ToArray(), tracers.Count);
    }

    void OnDrawGizmos()
    {
        DebugExtension.DebugLocalCube(transform, new Vector3(1, 1, 1));

        if(grid != null)
        {
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
                        Gizmos.color = new Vector4(0.7f, 0.7f, 1, 0.2f);
                        Gizmos.DrawWireCube(gridOrigin + new Vector3(cellExtent.x * i + cellExtent.x / 2, cellExtent.y * j + cellExtent.y / 2, cellExtent.z * k + cellExtent.z / 2), cellExtent);
                    }
                }
            }
        }
    }
	
	
}
