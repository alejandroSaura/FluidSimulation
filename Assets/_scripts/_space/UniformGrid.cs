using UnityEngine;
using System.Collections.Generic;
using System;

/*
    Templated container for fast spatial lookups and insertions
*/
public class UniformGrid<ItemT> : UniformGridGeometry where ItemT : IgridItem
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
    public UniformGrid(UniformGridGeometry gridToCopy) : base(gridToCopy) { }
    
    public Vector3 GetCellSpacing() { return GetCellExtent(); }

    // Operator [] overload
    public ItemT this[uint offset]
    {
        get
        {
            return mContents[(int)offset];
        }
        set
        {
            mContents[(int)offset] = value;
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

    public override void DefineShape(uint numElements, Vector3 posMin, Vector3 posMax, bool bPowerOf2)
    {
        mContents.Clear();
        mContents.TrimExcess();
        base.DefineShape(numElements, posMin, posMax, bPowerOf2);
    }

    // return number of cells which have been assigned values
    public int Size()
    {
        return mContents.Count;
    }

    public override void Decimate(UniformGridGeometry src, int iDecimation)
    {
        base.Decimate(src, iDecimation);
    }

    /*
        Compute statistics of data in an uniform grid
        min - minimum of all values in grid
        max -maximun of all values in grid
    */
    //public void ComputeStatistics(out ItemT min, out ItemT max)
    //{
    //    max = min = this[0];
    //    uint numCells = GetGridCapacity();
    //    for (uint offset = 0; offset < numCells; ++offset)
    //    {
    //        ItemT rVal = (this)[offset];
    //        min = (ItemT)MIN2<ItemT>(min, rVal);
    //        max = (ItemT)MAX2<ItemT>(max, rVal);
    //    }
    //}

    //T MIN2<T>(T x, T y) where T : IItemT
    //{
    //    if (x.CompareTo(y) < 0)
    //        return x;
    //    else
    //        return y;
    //}
    //T MAX2<T>(T x, T y) where T : IItemT
    //{
    //    if (x.CompareTo(y) > 0)
    //        return x;
    //    else
    //        return y;
    //}


    /*
        Interpolate values from Grid to get value at a given position
        vPosition - position to sample
    */
    public IgridItem Interpolate(Vector3 vPosition)
    {
        uint[] indices = new uint[3];
        indices = IndicesOfPosition(vPosition);
        Vector3 vMinCorner;
        vMinCorner = PositionFromIndices(indices);
        uint offsetX0Y0Z0 = OffsetFromIndices(indices);

        Vector3 vDiff = vPosition - vMinCorner; // Relative location of position within its containing grid cell.
        Vector3 tween = new Vector3(vDiff.x * GetCellsPerExtent().x, vDiff.y * GetCellsPerExtent().y, vDiff.z * GetCellsPerExtent().z);
        Vector3 oneMinusTween = new Vector3(1, 1, 1) - tween;

        uint numXY = GetNumPoints(0) * GetNumPoints(1);

        uint offsetX1Y0Z0 = offsetX0Y0Z0 + 1;
        uint offsetX0Y1Z0 = offsetX0Y0Z0 + GetNumPoints(0);
        uint offsetX1Y1Z0 = offsetX0Y0Z0 + GetNumPoints(0) + 1;
        uint offsetX0Y0Z1 = offsetX0Y0Z0 + numXY;
        uint offsetX1Y0Z1 = offsetX0Y0Z0 + numXY + 1;
        uint offsetX0Y1Z1 = offsetX0Y0Z0 + numXY + GetNumPoints(0);
        uint offsetX1Y1Z1 = offsetX0Y0Z0 + numXY + GetNumPoints(0) + 1;

        IgridItem vResult = ( oneMinusTween.x * oneMinusTween.y * oneMinusTween.z * (this)[offsetX0Y0Z0]
                    + tween.x * oneMinusTween.y * oneMinusTween.z * (this)[offsetX1Y0Z0]
                    + oneMinusTween.x * tween.y * oneMinusTween.z * (this)[offsetX0Y1Z0]
                    + tween.x * tween.y * oneMinusTween.z * (this)[offsetX1Y1Z0]
                    + oneMinusTween.x * oneMinusTween.y * tween.z * (this)[offsetX0Y0Z1]
                    + tween.x * oneMinusTween.y * tween.z * (this)[offsetX1Y0Z1]
                    + oneMinusTween.x * tween.y * tween.z * (this)[offsetX0Y1Z1]
                    + tween.x * tween.y * tween.z * (this)[offsetX1Y1Z1] );

        return vResult;
    }

    // Insert given value into grid at given position
    //public void Insert(Vector3 vPosition, ItemT item)
    //{
    //    uint[] indices = new uint[3]; // Indices of grid cell containing position.
    //    indices = IndicesOfPosition(vPosition);
    //    Vector3 vMinCorner;
    //    vMinCorner = PositionFromIndices(indices);
    //    uint offsetX0Y0Z0 = OffsetFromIndices(indices);
    //    Vector3 vDiff = vPosition - vMinCorner; // Relative location of position within its containing grid cell.
    //    Vector3 tween = new Vector3(vDiff.x * GetCellsPerExtent().x, vDiff.y * GetCellsPerExtent().y, vDiff.z * GetCellsPerExtent().z);
    //    Vector3 oneMinusTween = new Vector3(1.0f, 1.0f, 1.0f) - tween;
    //    uint numXY = GetNumPoints(0) * GetNumPoints(1);
    //    uint offsetX1Y0Z0 = offsetX0Y0Z0 + 1;
    //    uint offsetX0Y1Z0 = offsetX0Y0Z0 + GetNumPoints(0);
    //    uint offsetX1Y1Z0 = offsetX0Y0Z0 + GetNumPoints(0) + 1;
    //    uint offsetX0Y0Z1 = offsetX0Y0Z0 + numXY;
    //    uint offsetX1Y0Z1 = offsetX0Y0Z0 + numXY + 1;
    //    uint offsetX0Y1Z1 = offsetX0Y0Z0 + numXY + GetNumPoints(0);
    //    uint offsetX1Y1Z1 = offsetX0Y0Z0 + numXY + GetNumPoints(0) + 1;
    //    (this)[offsetX0Y0Z0] += oneMinusTween.x * oneMinusTween.y * oneMinusTween.z * item;
    //    (this)[offsetX1Y0Z0] += tween.x * oneMinusTween.y * oneMinusTween.z * item;
    //    (this)[offsetX0Y1Z0] += oneMinusTween.x * tween.y * oneMinusTween.z * item;
    //    (this)[offsetX1Y1Z0] += tween.x * tween.y * oneMinusTween.z * item;
    //    (this)[offsetX0Y0Z1] += oneMinusTween.x * oneMinusTween.y * tween.z * item;
    //    (this)[offsetX1Y0Z1] += tween.x * oneMinusTween.y * tween.z * item;
    //    (this)[offsetX0Y1Z1] += oneMinusTween.x * tween.y * tween.z * item;
    //    (this)[offsetX1Y1Z1] += tween.x * tween.y * tween.z * item;
    //}

    public override void Clear()
    {
        mContents.Clear();
        base.Clear();
    }

}
