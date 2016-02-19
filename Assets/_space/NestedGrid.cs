using UnityEngine;
using System.Collections.Generic;



public class NestedGrid<ItemT>
{
    List<UniformGrid<ItemT>> mLayers;   // Dynamic array of UniformGrids
    int[] mDecimations;   // Cache of cluster sizes

    // Construc a blank nested uniform grid spatial partition
    NestedGrid()           
    {
        mDecimations = new int[3];
    }

    // Construct an unpopulated nested uniform grid, based on a given layer
    NestedGrid(UniformGrid<ItemT> layer)
    {
        mDecimations = new int[3];
        Initialize(layer);
    }

    void Initialize(UniformGrid<ItemT> srcLayer)
    {
        mLayers.Clear();
        mLayers.TrimExcess();

        int numLayers = PrecomputeNumLayers(srcLayer);
        mLayers.Capacity = numLayers;

        AddLayer(srcLayer, 1);
        int index = 1;
        while(mLayers[index - 1] != null && mLayers[index-1].GetGridCapacity()>8) // a cell has 8 corners
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

    public void AddLayer( UniformGridGeometry layerTemplate , int iDecimation )
    {
        mLayers.Add(new UniformGrid<ItemT>()) ;
        int count = mLayers.Count;
        mLayers[count-1].Decimate(layerTemplate , iDecimation );
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

    public int GetDecimations(int parentLayerIndex)
    {
        return mDecimations[parentLayerIndex];
    }
    

}
