using UnityEngine;
using System.Collections.Generic;
using System;



public class NestedGrid<ItemT> where ItemT : IMatrix3x3
{
    public List<UniformGrid<ItemT>> mLayers;   // Dynamic array of UniformGrids
    uint[][] mDecimations;   // Cache of cluster sizes

    // Construc a blank nested uniform grid spatial partition
    public NestedGrid()
    {
        mLayers = new List<UniformGrid<ItemT>>();
    }

    // Construct an unpopulated nested uniform grid, based on a given layer
    public NestedGrid(UniformGrid<ItemT> layer)
    {
        mLayers = new List<UniformGrid<ItemT>>();
        Initialize(layer);
    }

    public void Initialize(UniformGrid<ItemT> srcLayer)
    {
        mLayers.Clear();
        mLayers.TrimExcess();

        int numLayers = PrecomputeNumLayers(srcLayer);
        mLayers.Capacity = numLayers;

        AddLayer(srcLayer, 1);
        int index = 1;
        while (mLayers[index - 1] != null && mLayers[index - 1].GetGridCapacity() > 8) // a cell has 8 corners
        { // layer to decimate has more than 1 cell
            AddLayer(mLayers[index - 1], 2); // Initialize child layer based on decimation of its parent grid.
            ++index;
        }

        PrecomputeDecimations();
    }

    /* 
        Add a layer to the top of the nested grid.
        layerTemplate - UNiformGridGeometry defining child layer.
        iDecimation - amount by which to decimate child layer.

        This facilitates building the tree from leaves to root.
        This method laso allocates memory for the newly added layer,
        and initializes its content with the default constructor.
    */

    public void AddLayer(UniformGridGeometry layerTemplate, int iDecimation)
    {
        mLayers.Add(new UniformGrid<ItemT>());
        int count = mLayers.Count;
        mLayers[count - 1].Decimate(layerTemplate, iDecimation);
        mLayers[count - 1].Init();
    }

    // Return number of layers in tree
    public int GetDepth() { return mLayers.Count; }

    /*
        Get layer (uniformGrid) at specific depth of the tree.
        index - depth of the layer to obtain. 
            0 means the leaf layer and GetDepth()-1 means root layer.
    */
    public UniformGrid<ItemT> this[uint offset] // Operator [] overload
    {
        get
        {
            return mLayers[(int)offset];
        }
    }

    public uint[] GetDecimations(int parentLayerIndex)
    {
        return mDecimations[parentLayerIndex];
    }

    /*
        Get indices of minimal cell in child layer of cluster represented by specified cell in parent layer.

        Each cell in a parent layer represents a grid cluster of typically 8 cells
        in the child layer.  This routine calculates the index of the "minimal"
        cell in the child layer grid cluster, i.e. the cell in the child layer
        which corresponds to minimum corner cell of the grid cluster represented
        by the cell in the parent layer with the specified index.        
    */
    public int[] GetChildClusterMinCornerIndex(int[] decimations, int[] indicesOfParentCell)
    {
        int[] clusterMinIndices = new int[3];
        clusterMinIndices[0] = indicesOfParentCell[0] * decimations[0];
        clusterMinIndices[1] = indicesOfParentCell[1] * decimations[1];
        clusterMinIndices[2] = indicesOfParentCell[2] * decimations[2];

        return clusterMinIndices;
    }

    public void Clear()
    {
        for (int iLayer = 0; iLayer < GetDepth(); ++iLayer)
        {
            mLayers[iLayer].Clear();
        }
        mLayers.Clear();
    }

    private NestedGrid(NestedGrid<ItemT> other) {} // disallow copy construction


    /*
        Precompute the total number of layers this nested grid will contain          
        src - UniformGrid upon which this NestedGrid is based.
    */
    int PrecomputeNumLayers(UniformGrid<ItemT> src)
    {
        int numLayers = 1;    // Tally src layer.
        uint[] numPoints = new uint[3];
        numPoints[0] = src.GetNumPoints(0);
        numPoints[1] = src.GetNumPoints(1);
        numPoints[2] = src.GetNumPoints(2);
        uint size = numPoints[0] * numPoints[1] * numPoints[2];
        while( size > 8 /* a cell has 8 corners */ )
        {   // Layer has more than 1 cell.
            ++ numLayers ;
            // Decimate number of cells (where #cells = #points-1):
            numPoints[0] = (uint)Mathf.Max( (float)( numPoints[0] - 1 ) / 2 , 1 ) + 1 ;
            numPoints[1] = (uint)Mathf.Max((float)( numPoints[1] - 1 ) / 2 , 1 ) + 1 ;
            numPoints[2] = (uint)Mathf.Max((float)( numPoints[2] - 1 ) / 2 , 1 ) + 1 ;
            size = numPoints[0] * numPoints[1] * numPoints[2] ;
        }
        return numLayers ;
     }

    /*
        Precompute decimations for each layer.

        This provides the number of grid cells per cluster
        of a child of each layer.

        The child layer has index one less than the parent layer index.
        That implies there is no such thing as "parent layer 0".
        Layer 0  has no children. That further implies there is no
        meaningful value for decimations at iParentLayer==0.
    */
    void PrecomputeDecimations()
    {
        int numLayers = GetDepth();

        // Precompute decimations for each layer.
        mDecimations = new uint[numLayers][];
        for (int i = 0; i < numLayers; ++i)
        {// initialize the array
            mDecimations[i] = new uint[3];
        }

        for (uint iLayer = 1; iLayer < numLayers; ++iLayer)
        {   // For each parent layer...            
            mDecimations[iLayer] = ComputeDecimations(iLayer);
        }
        // Layer 0 is strictly a child (i.e. has no children), so has no decimations.
        // Assign the values with useless nonsense to make this more obvious.
        mDecimations[0][0] = mDecimations[0][1] = mDecimations[0][2] = 0;
    }


    /*
        Compute decimations, in each direction, for specified parent layer

        decimations - (out) ratio of dimensions between child layer and its parent.
        iParentLayer - index of parent layer.
                            Child has index iParentLayer-1.
                            Layer 0 has no child so providing "0" is invalid.

        This method effectively gives the number of child cells in each
        grid cluster that a parent cell represents.

        Each non-leaf layer in this NestedGrid is a decimation of its child
        layer. Typically that decimation is 2 in each direction, but the
        decimation can also be 1, or, more atypically, any other integer.
        Each child typically has twice as many cells in each direction as
        its parent.

        This assumes each parent has an integer decimation of its child.
    */
    uint[] ComputeDecimations(uint iParentLayer) 
    {
        uint[] decimations = new uint[3];
        UniformGrid<ItemT> parent = (this)[iParentLayer] ;
        UniformGrid<ItemT> child  = (this)[iParentLayer - 1 ] ;
        decimations[0] = child.GetNumCells( 0 ) / parent.GetNumCells( 0 ) ;
        decimations[1] = child.GetNumCells( 1 ) / parent.GetNumCells( 1 ) ;
        decimations[2] = child.GetNumCells( 2 ) / parent.GetNumCells( 2 ) ;

        return decimations;
    }


}
