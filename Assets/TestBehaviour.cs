using UnityEngine;
using System.Collections.Generic;


public class TestBehaviour : MonoBehaviour
{
    public bool[] seeLayer;

    NestedGrid<Vorton> nestedGrid;

    void Start()
    {
        Vector m = new Vector();
        m = (Vector)(5 * m);

        nestedGrid = new NestedGrid<Vorton>();
        nestedGrid.Initialize(new UniformGrid<Vorton>(4096, new Vector3(0, 0, 0), new Vector3(5, 5, 5), true));

        seeLayer = new bool[nestedGrid.GetDepth()];
    }

    // to view the grid.
    void OnDrawGizmos()
    {
        if (nestedGrid != null && nestedGrid.mLayers.Count != 0)
        {
            for (int u = 0; u < nestedGrid.mLayers.Count; ++u)
            {
                if (seeLayer[u] == false) continue;

                UniformGrid<Vorton> grid = nestedGrid.mLayers[u];

                Gizmos.color = new Vector4(1, 1 - 1 / ((float)u + 1), 1 - 1 / ((float)u + 1), 1.1f - 1 / ((float)u + 1));

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
                            Gizmos.DrawWireCube(gridOrigin + new Vector3(cellExtent.x * i + cellExtent.x/2, cellExtent.y * j + cellExtent.y / 2, cellExtent.z * k + cellExtent.z / 2), cellExtent);
                        }
                    }
                }
            }

        }
    }

}
