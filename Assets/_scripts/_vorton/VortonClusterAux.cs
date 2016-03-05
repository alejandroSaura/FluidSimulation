using UnityEngine;
using System.Collections;

/*
    Auxiliary information used to calculate a representation of a cluster of vortons

    In this fluid simulation, a "vorton cluster" refers to a collection of vortons
    which are all inside the same region.

    Each vorton influences the motion of all other vortons.

    This fluid simulation approximates the influence of vorton clusters
    by representing them in a simpler form, i.e. a form which has
    fewer parameters and less computational complexity than expressing
    the influence of all individual vortons in the cluster.

    This structure is used to store information used to compute parameters
    of that simplified form.  This structure is NOT used to represent
    the simplified form itself, during the influence calculation.
*/

public class VortonClusterAux : IgridItem
{
    public float mVortNormSum = float.Epsilon;

    /*
        Construct auxiliary info for aggregating vortons

        Initial values for mVortNormSum is initially almost 0
        because it accumulates the sum of |vorticity|.
        mVortNormSum is not, however, exactly 0; instead
        it is FLT_MIN (the smallest value a float can
        represent) in order to avoid divide-by-zero.
        The only time when mVortNormSum would be 0
        is when a grid cluster contains no vortons,
        in which case the center-of-vorticity is undefined.
    */
    public VortonClusterAux()
    {
        mVortNormSum = float.Epsilon;
    }

}
