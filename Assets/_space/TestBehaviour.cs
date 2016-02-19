using UnityEngine;
using System.Collections.Generic;


public class TestBehaviour : MonoBehaviour
{
    UniformGrid<int> uniformGrid;

    void Start()
    {        
        uniformGrid = new UniformGrid<int>(4096, new Vector3(1,1,1), new Vector3(5, 5, 5), true);
        uniformGrid.Init();
    }

    // to view the grid.
    void OnDrawGizmos()
    {
        if (uniformGrid != null)
        {
            Gizmos.color = new Vector4(0.5f, 0.5f, 1, 0.5f);

            Vector3 cellExtent = uniformGrid.GetCellExtent();
            Vector3 gridExtent = uniformGrid.GetExtent();
            Vector3 numCells = new Vector3(uniformGrid.GetNumCells(0), uniformGrid.GetNumCells(1), uniformGrid.GetNumCells(2));
            Vector3 gridOrigin = uniformGrid.GetMinCorner();

            for (int i = 0; i < numCells.x; ++i)
            {
                for (int j = 0; j < numCells.y; ++j)
                {
                    for (int k = 0; k < numCells.z; ++k)
                    {
                        Gizmos.DrawWireCube(gridOrigin + new Vector3(cellExtent.x * i, cellExtent.y * j, cellExtent.z * k), cellExtent);
                    }
                }
            }

        }
    }
        
}
