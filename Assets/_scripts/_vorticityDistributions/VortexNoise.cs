using UnityEngine;

// Specify vorticity in the shape of a vortex ring
public class VortexNoise : MonoBehaviour, IVorticityDistribution
{
    public Vector3 mBox = new Vector3(1, 1, 1);
    public Vector3 mAmplitude = new Vector3(1, 1, 1);

    static int counter = 1;

    public Vector3 GetCenter()
    {
        return transform.position;
    }

    public VortexNoise(Vector3 vBox)
    {
        mBox = vBox;
        mAmplitude = new Vector3(1, 1, 1);
    }

    public Vector3 GetDomainSize()
    {
        return mBox;
    }

    public void AssignVorticity(ref Vector3 vorticity, Vector3 position, Vector3 vCenter)
    {       

        //generate rnd number
        Random.seed = System.DateTime.Now.Millisecond + counter;
        float n1 = Random.Range(-0.5f, 0.5f);
        Random.seed = System.DateTime.Now.Millisecond + counter + 10;
        float n2 = Random.Range(-0.5f, 0.5f);
        Random.seed = System.DateTime.Now.Millisecond + counter + 20;
        float n3 = Random.Range(-0.5f, 0.5f);

        vorticity = new Vector3(n1 * mAmplitude.x, n2 * mAmplitude.y, n3 * mAmplitude.z) ;

        ++counter;
    }

    void OnDrawGizmos()
    {
        DebugExtension.DrawLocalCube(transform, mBox);
    }


}