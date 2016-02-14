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


/*
    The shape of this grid is such that 
    the "minimal corner" point resides at indices {0,0,0} 
    and the "maximal corner" point resides at indices { Nx-1,Ny-1,Nz-1}.

    The number of /points/ in each direction i is N_i.
    A cell is defined by the 8 points that lie at its corners.
    The number of /cells/ in each direction i is (N_i-1).

    The size of a side i of each cell is therefore
    s_i = (posMax-posMin)_i / (N_i-1) .
*/

public class UniformGridGeometry
{
    Vector3 mMinCorner; // starting position
    public Vector3 GetMinCorner() { return (mMinCorner); }
    Vector3 mGridExtent; // total grid size
    public Vector3 GetExtent() { return mGridExtent; }
    Vector3 mCellExtent; // cell size
    public Vector3 GetCellExtent() { return mCellExtent; }
    Vector3 mCellsPerExtent; // units relation between cell and total extent
    public Vector3 GetCellsPerExtent() { return mCellsPerExtent; }
    protected int[] mNumPoints = new int[3];
    public uint GetNumCells( uint index ) { return (uint)mNumPoints[index] - 1 ; } // cells per side

    public uint GetGridCapacity() { return (uint)(mNumPoints[0] * mNumPoints[1] * mNumPoints[2]); }



public UniformGridGeometry() // Creates empty Uniform Grid geometry
    {
        mMinCorner = new Vector3(0.0f, 0.0f, 0.0f);
        mGridExtent = new Vector3(0.0f, 0.0f, 0.0f);
        mCellExtent = new Vector3(0.0f, 0.0f, 0.0f);
        mCellsPerExtent = new Vector3(0.0f, 0.0f, 0.0f);
        mNumPoints[0] = mNumPoints[1] = mNumPoints[2] = 0;
    }

    /* 
        Copy constructor
    */
    public UniformGridGeometry (UniformGridGeometry gridToCopy)
    {
        // Do nothing
    }

    /* 
        Copy shape information from another UniformGrid into this one
    */
    void CopyShape(UniformGridGeometry src )
    {
        Decimate(src, 1);
    }

    // Construct a uniform grid that fits the given geometry.
    public UniformGridGeometry(uint uNumElements, Vector3 posMin , Vector3 posMax , bool bPowerOf2)
    {
            DefineShape(uNumElements , posMin, posMax, bPowerOf2 );            
    }

    /*
    Define the shape a uniform grid such that it fits the given geometry.
        
        Params:
        uNumElements - number of elements this container will contain.
        posMin - minimal coordinate of axis-aligned bounding box.
        posMax - maximal coordinate of axis-aligned bounding box.
        bPowerOf2 - whether to make each grid dimension a power of 2.
                Doing so simplifies grid subdivision, if this grid will be used in a hierarchical grid.

        This makes a uniform grid of cells, where each cell is the same size
        and the side of each cell is nearly the same size.  If the cells are
        3-dimensional then that means each cell is a box, nearly a cube.
        The number of dimensions of the region depends on the actual size of
        the region.  If any size component is smaller than a small threshold
        then this class considers that component to be zero, and reduces the
        dimensionality of the region.  For example, if the region size is
        (2,3,0) then this class considers the region to have 2 dimensions
        (x and y) since the z size is zero.  In this example, the cells
        would be nearly square rectangles (instead of boxes).
    */
    void DefineShape(uint uNumElements, Vector3 posMin, Vector3 posMax, bool bPowerOf2)
    {
        mMinCorner = posMin;          
        mGridExtent = (posMax - posMin) * global.GlobalVar.Nudge; // slightly expand size to ensure robust containment even with roundoff

        Vector3 vSizeEffective = mGridExtent;

        int numDims = 3;   // Number of dimensions to region.
        if (0.0f == vSizeEffective.x)
        {   // X size is zero so reduce dimensionality
            vSizeEffective.x = 1.0f; // This component will not contribute to the total region volume/area/length.
            mGridExtent.x = 0.0f;
            --numDims;
        }
        if (0.0f == vSizeEffective.y)
        {   // Y size is zero so reduce dimensionality
            vSizeEffective.y = 1.0f; // This component will not contribute to the total region volume/area/length.
            mGridExtent.y = 0.0f;
            --numDims;
        }
        if (0.0f == vSizeEffective.z)
        {   // Z size is zero so reduce dimensionality
            vSizeEffective.z = 1.0f; // This component will not contribute to the total region volume/area/length.
            mGridExtent.z = 0.0f;
            --numDims;
        }

        // Compute region volume, area or length (depending on dimensionality).
        float volume = vSizeEffective.x * vSizeEffective.y * vSizeEffective.z;
        float cellVolumeCubeRoot = Mathf.Pow(volume / (uNumElements), -1.0f / (numDims));   // Approximate size of each cell in grid.
                                                                                            // Compute number of cells in each direction of uniform grid.
                                                                                            // Choose grid dimensions to fit as well as possible, so that the total number
                                                                                            // of grid cells is nearly the total number of elements in the contents.

        int[] numCells = { (int)Mathf.Max( 1 , (uint)( mGridExtent.x * cellVolumeCubeRoot + 0.5f ) ) ,
                            (int)Mathf.Max( 1 , (uint)( mGridExtent.y * cellVolumeCubeRoot + 0.5f ) ) ,
                            (int)Mathf.Max( 1 , (uint)( mGridExtent.z * cellVolumeCubeRoot + 0.5f ) ) };

        if (bPowerOf2)
        {   // Choose number of gridcells to be powers of 2.
            // This will simplify subdivision in a NestedGrid.
            numCells[0] = Mathf.ClosestPowerOfTwo(numCells[0]);
            numCells[1] = Mathf.ClosestPowerOfTwo(numCells[1]);
            numCells[2] = Mathf.ClosestPowerOfTwo(numCells[2]);
        }
        while (numCells[0] * numCells[1] * numCells[2] >= uNumElements * 8)
        {   // Grid capacity is excessive.
            // This can occur when the trial numCells is below 0.5 in which case the integer arithmetic loses the subtlety.
            numCells[0] = (int)Mathf.Max(1, numCells[0] / 2);
            numCells[1] = (int)Mathf.Max(1, numCells[1] / 2);
            numCells[2] = (int)Mathf.Max(1, numCells[2] / 2);
        }
        mNumPoints[0] = numCells[0] + 1; // Increment to obtain number of points.
        mNumPoints[1] = numCells[1] + 1; // Increment to obtain number of points.
        mNumPoints[2] = numCells[2] + 1; // Increment to obtain number of points.

        PrecomputeSpacing();
    }

    protected void PrecomputeSpacing()
    {
        mCellExtent.x = mGridExtent.x / (float)(GetNumCells(0));
        mCellExtent.y = mGridExtent.y / (float)(GetNumCells(1));
        mCellExtent.z = mGridExtent.z / (float)(GetNumCells(2));
        mCellsPerExtent.x = (float)(GetNumCells(0)) / mGridExtent.x;
        mCellsPerExtent.y = (float)(GetNumCells(1)) / mGridExtent.y;
        if (0.0f == GetExtent().z)
        {   // Avoid divide-by-zero for 2D domains that lie in the XY plane.
            mCellsPerExtent.z = 1.0f / float.MinValue;
        }
        else
        {
            mCellsPerExtent.z = (float)(GetNumCells(2)) / mGridExtent.z;
        }
    }

    /*
        Create a lower-resolution uniform grid based on another
        src - Source uniform grid upon which to base dimensions of this one
        iDecimation - amount by which to reduce the number of grid cells in each dimension. Typically this would be 2.

        The number of cells is decimated.  The number of points is different.
    */
    void Decimate( UniformGridGeometry src , int iDecimation )
    {
        mGridExtent = src.mGridExtent ;
        mMinCorner = src.mMinCorner ;
        mNumPoints[0] = (int)(src.GetNumCells( 0 ) / iDecimation + 1) ;
        mNumPoints[1] = (int)(src.GetNumCells( 1 ) / iDecimation + 1) ;
        mNumPoints[2] = (int)(src.GetNumCells( 2 ) / iDecimation + 1) ;

        if( iDecimation > 1 )
        {   // Decimation could reduce dimension and integer arithmetic could make value be 0, which is useless if src contained any data.
            mNumPoints[0] = Mathf.Max( 2 , mNumPoints[0] ) ;
            mNumPoints[1] = Mathf.Max( 2 , mNumPoints[1] ) ;
            mNumPoints[2] = Mathf.Max( 2 , mNumPoints[2] ) ;
        }

        PrecomputeSpacing();
    }

    // In this region we translate between position, offset and indices
    #region locationTranslators

    /*
        Compute X,Y,Z grid cell indices from offset into contents array.
        indices - Individual X,Y,Z component grid cell indices.
        offset - Offset into mContents.
    */
    public uint[] IndicesFromOffset(uint offset )
    {
        uint[] indices = new uint[3];
        indices[2] = (uint)(offset / ( mNumPoints[0] * mNumPoints[1])) ;
        indices[1] = (uint)(( offset - indices[2] * mNumPoints[0] * mNumPoints[1]) / mNumPoints[0]);
        indices[0] = (uint)(offset - mNumPoints[0] * ( indices[1] + mNumPoints[1] * indices[2] )) ;

        return indices;
    }

    /*
        Compute position of minimal corner of grid cell with given indices
        position - position of minimal corner of grid cell
        indices - grid cell indices.

        Rarely if ever would you want to compute position from indices in this way.
        Typically, this kind of computation occurs inside a triply-nested loop,
        in which case the procedure should compute each component
        separately.  Furthermore, such a routine would cache
        GetCellSpacing instead of computing it each iteration.

    */
    public Vector3 PositionFromIndices(uint[] indices)
    {
        Vector3 vPosition;
        vPosition.x = GetMinCorner().x + (float)( indices[0] ) * mCellExtent.x ;
        vPosition.y = GetMinCorner().y + (float)( indices[1] ) * mCellExtent.y ;
        vPosition.z = GetMinCorner().z + (float)( indices[2] ) * mCellExtent.z ;

        return vPosition;
    }

    /*
        brief Get position of grid cell minimum corner.
        vPos - position of grid cell minimum corner
        offset - offset into contents array

        Each grid cell spans a region (whose size is given by GetCellSpacing)
        starting at a location which this routine returns.  So the grid cell
        with the given offset spans the region from vPos (as this routine
        assigns) to vPos + mCellExtent.

        Derived class provides actual contents array.
    */
    public Vector3 PositionFromOffset(uint offset)
    {
        Vector3 vPos = new Vector3();
        uint[] indices = new uint[3];
        indices = IndicesFromOffset(offset );
        vPos.x = GetMinCorner().x + (float)( indices[0] ) * mCellExtent.x ;
        vPos.y = GetMinCorner().y + (float)( indices[1] ) * mCellExtent.y ;
        vPos.z = GetMinCorner().z + (float)( indices[2] ) * mCellExtent.z ;

        return vPos;
    }


    /*
        Compute indices into contents array of a point at a given position
        vPosition - position of a point.  It must be within the region of this container.
        indices - Indices into contents array of a point at vPosition.

        Derived class defines the actual contents array.
    */
    public uint[] IndicesOfPosition(Vector3 vPosition)
    {
        uint[] indices = new uint[3];
        // Notice the peculiar test here.  vPosition may lie slightly outside of the extent give by vMax.
        Vector3 vPosRel = (vPosition - GetMinCorner()) ;   // position of given point relative to container region
        Vector3 vIdx = new Vector3(vPosRel.x * GetCellsPerExtent().x , vPosRel.y * GetCellsPerExtent().y , vPosRel.z* GetCellsPerExtent().z ) ;
        indices[0] = (uint)(vIdx.x );
        indices[1] = (uint)(vIdx.y );
        indices[2] = (uint)(vIdx.z );
        return indices;
    }

    /*
        Compute offset into contents array of a point at a given position
        vPosition - position of a point.  It must be within the region of this container.
        Returns Offset into contents array of a point at vPosition.

        Derived class defines the actual contents array.
    */
    public uint OffsetOfPosition(Vector3 vPosition)
    {
        uint[] indices = new uint[3];
        indices = IndicesOfPosition(vPosition);
        uint offset = (uint)(indices[0] + mNumPoints[0] * (indices[1] + mNumPoints[1] * indices[2]));
        return offset;
    }   

    #endregion

}

public class UniformGrid<ItemT> : UniformGridGeometry
{
    // Data:
    List<ItemT> mContents = new List<ItemT>();   // 3D array of items.

    /*
        Construct an empty UniformGrid.
    */
    public UniformGrid() : base() { }

    /*
        Construct a uniform grid container that fits the given geometry.
    */
    public UniformGrid(uint uNumElements, Vector3 posMin, Vector3 posMax, bool bPowerOf2) : base(uNumElements, posMin, posMax, bPowerOf2) { }

    /*
        Copy shape from given uniform grid
    */
    public UniformGrid(UniformGridGeometry gridToCopy) : base(gridToCopy ) { }


    // Operator [] overload
    public ItemT this[uint offset]
    {
        get
        {
            return mContents[(int)offset];
        }        
    }
    public ItemT this[Vector3 vPosition]
    {
        get
        {
            return mContents[(int)OffsetOfPosition(vPosition)];
        }
    }

    public void Init()
    {        
        global.ListExtra.Resize<ItemT>(mContents, (int)GetGridCapacity());
    }
}
