using UnityEngine;
using System.Collections;

public static class UniformGridMath
{
    public static void ComputeJacobian(ref UniformGrid<Matrix3x3> jacobian, UniformGrid<Vector> vec)
    {
       Vector3 spacing = vec.GetCellSpacing();
        // Avoid divide-by-zero when z size is effectively 0 (for 2D domains)
        Vector3 reciprocalSpacing = new Vector3(1.0f / spacing.x , 1.0f / spacing.y, spacing.z > float.Epsilon ? 1.0f / spacing.z : 0.0f ) ;
        Vector3 halfReciprocalSpacing = (0.5f * reciprocalSpacing) ;

        uint[] dims = { vec.GetNumPoints(0), vec.GetNumPoints(1), vec.GetNumPoints(2) };
        uint[] dimsMinus1 = { vec.GetNumPoints(0) - 1, vec.GetNumPoints(1) - 1, vec.GetNumPoints(2) - 1 };
        uint numXY = dims[0] * dims[1];
        uint[] index = new uint[3];

        
        // Compute derivatives for interior (i.e. away from boundaries).
        for (index[2] = 1; index[2] < dimsMinus1[2]; ++index[2])
        {
            //ASSIGN_Z_OFFSETS
            uint offsetZM = numXY * (index[2] - 1);
            uint offsetZ0 = numXY * index[2];
            uint offsetZP = numXY * (index[2] + 1); 

            for (index[1] = 1; index[1] < dimsMinus1[1]; ++index[1])
            {
                //ASSIGN_YZ_OFFSETS                                                   
                uint offsetYMZ0 = dims[0] * (index[1] - 1) + offsetZ0;
                uint offsetY0Z0 = dims[0] * index[1] + offsetZ0;
                uint offsetYPZ0 = dims[0] * (index[1] + 1) + offsetZ0;
                uint offsetY0ZM = dims[0] * index[1] + offsetZM;
                uint offsetY0ZP = dims[0] * index[1] + offsetZP;

                for (index[0] = 1; index[0] < dimsMinus1[0]; ++index[0])
                {
                    //ASSIGN_XYZ_OFFSETS                                      
                    uint offsetX0Y0Z0 = index[0] + offsetY0Z0;
                    uint offsetXMY0Z0 = index[0] - 1 + offsetY0Z0;
                    uint offsetXPY0Z0 = index[0] + 1 + offsetY0Z0;
                    uint offsetX0YMZ0 = index[0] + offsetYMZ0;
                    uint offsetX0YPZ0 = index[0] + offsetYPZ0;
                    uint offsetX0Y0ZM = index[0] + offsetY0ZM;
                    uint offsetX0Y0ZP = index[0] + offsetY0ZP;

                    Matrix3x3 rMatrix = jacobian[offsetX0Y0Z0];
                    /* Compute d/dx */
                    rMatrix.x = (vec[offsetXPY0Z0].v - vec[offsetXMY0Z0].v) * halfReciprocalSpacing.x;
                    /* Compute d/dy */
                    rMatrix.y = (vec[offsetX0YPZ0].v - vec[offsetX0YMZ0].v) * halfReciprocalSpacing.y;
                    /* Compute d/dz */
                    rMatrix.z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0ZM].v) * halfReciprocalSpacing.z;
                }
            }
        }

        #region boundariesDerivatives
        // Compute derivatives for boundaries: 6 faces of box. ------------------------------------------------------------
        // In some situations, these macros compute extraneous data.
        // A tiny bit more efficiency could be squeezed from this routine,
        // but it turns out to be well under 1% of the total expense.
                

        // Compute derivatives for -X boundary.
        index[0] = 0;
        for (index[2] = 0; index[2] < dims[2]; ++index[2])
        {
            //ASSIGN_Z_OFFSETS
            uint offsetZM = numXY * (index[2] - 1);
            uint offsetZ0 = numXY * index[2];
            uint offsetZP = numXY * (index[2] + 1);

            for (index[1] = 0; index[1] < dims[1]; ++index[1])
            {
                //ASSIGN_YZ_OFFSETS                                                   
                uint offsetYMZ0 = dims[0] * (index[1] - 1) + offsetZ0;
                uint offsetY0Z0 = dims[0] * index[1] + offsetZ0;
                uint offsetYPZ0 = dims[0] * (index[1] + 1) + offsetZ0;
                uint offsetY0ZM = dims[0] * index[1] + offsetZM;
                uint offsetY0ZP = dims[0] * index[1] + offsetZP;

                {
                    //ASSIGN_XYZ_OFFSETS                                      
                    uint offsetX0Y0Z0 = index[0] + offsetY0Z0;
                    uint offsetXMY0Z0 = index[0] - 1 + offsetY0Z0;
                    uint offsetXPY0Z0 = index[0] + 1 + offsetY0Z0;
                    uint offsetX0YMZ0 = index[0] + offsetYMZ0;
                    uint offsetX0YPZ0 = index[0] + offsetYPZ0;
                    uint offsetX0Y0ZM = index[0] + offsetY0ZM;
                    uint offsetX0Y0ZP = index[0] + offsetY0ZP;

                    //COMPUTE_FINITE_DIFF                                                                                                             
                    if (index[0] == 0)
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.x; }   
                    else if (index[0] == dimsMinus1[0])
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetX0Y0Z0].v - vec[offsetXMY0Z0].v) * reciprocalSpacing.x; }   
                    else
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetXMY0Z0].v) * halfReciprocalSpacing.x; }
                       
                    if (index[1] == 0)
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.y; }   
                    else if (index[1] == dimsMinus1[1])
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0Y0Z0].v - vec[offsetX0YMZ0].v) * reciprocalSpacing.y; }   
                    else
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0YMZ0].v) * halfReciprocalSpacing.y; } 
                      
                    if (index[2] == 0)
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.z; }   
                    else if (index[2] == dimsMinus1[2])
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0Z0].v - vec[offsetX0Y0ZM].v) * reciprocalSpacing.z; }   
                    else
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0ZM].v) * halfReciprocalSpacing.z; }

                }
            }
        }

        // Compute derivatives for -Y boundary.
        index[1] = 0;
        for (index[2] = 0; index[2] < dims[2]; ++index[2])
        {
            //ASSIGN_Z_OFFSETS
            uint offsetZM = numXY * (index[2] - 1);
            uint offsetZ0 = numXY * index[2];
            uint offsetZP = numXY * (index[2] + 1);

            {
                //ASSIGN_YZ_OFFSETS                                                   
                uint offsetYMZ0 = dims[0] * (index[1] - 1) + offsetZ0;
                uint offsetY0Z0 = dims[0] * index[1] + offsetZ0;
                uint offsetYPZ0 = dims[0] * (index[1] + 1) + offsetZ0;
                uint offsetY0ZM = dims[0] * index[1] + offsetZM;
                uint offsetY0ZP = dims[0] * index[1] + offsetZP;

                for (index[0] = 0; index[0] < dims[0]; ++index[0])
                {
                    //ASSIGN_XYZ_OFFSETS                                      
                    uint offsetX0Y0Z0 = index[0] + offsetY0Z0;
                    uint offsetXMY0Z0 = index[0] - 1 + offsetY0Z0;
                    uint offsetXPY0Z0 = index[0] + 1 + offsetY0Z0;
                    uint offsetX0YMZ0 = index[0] + offsetYMZ0;
                    uint offsetX0YPZ0 = index[0] + offsetYPZ0;
                    uint offsetX0Y0ZM = index[0] + offsetY0ZM;
                    uint offsetX0Y0ZP = index[0] + offsetY0ZP;

                    //COMPUTE_FINITE_DIFF                                                                                                             
                    if (index[0] == 0)
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.x; }
                    else if (index[0] == dimsMinus1[0])
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetX0Y0Z0].v - vec[offsetXMY0Z0].v) * reciprocalSpacing.x; }
                    else
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetXMY0Z0].v) * halfReciprocalSpacing.x; }

                    if (index[1] == 0)
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.y; }
                    else if (index[1] == dimsMinus1[1])
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0Y0Z0].v - vec[offsetX0YMZ0].v) * reciprocalSpacing.y; }
                    else
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0YMZ0].v) * halfReciprocalSpacing.y; }

                    if (index[2] == 0)
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.z; }
                    else if (index[2] == dimsMinus1[2])
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0Z0].v - vec[offsetX0Y0ZM].v) * reciprocalSpacing.z; }
                    else
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0ZM].v) * halfReciprocalSpacing.z; }
                }
            }
        }

        // Compute derivatives for -Z boundary.
        index[2] = 0;
        {
            //ASSIGN_Z_OFFSETS
            uint offsetZM = numXY * (index[2] - 1);
            uint offsetZ0 = numXY * index[2];
            uint offsetZP = numXY * (index[2] + 1);

            for (index[1] = 0; index[1] < dims[1]; ++index[1])
            {
                //ASSIGN_YZ_OFFSETS                                                   
                uint offsetYMZ0 = dims[0] * (index[1] - 1) + offsetZ0;
                uint offsetY0Z0 = dims[0] * index[1] + offsetZ0;
                uint offsetYPZ0 = dims[0] * (index[1] + 1) + offsetZ0;
                uint offsetY0ZM = dims[0] * index[1] + offsetZM;
                uint offsetY0ZP = dims[0] * index[1] + offsetZP;

                for (index[0] = 0; index[0] < dims[0]; ++index[0])
                {
                    //ASSIGN_XYZ_OFFSETS                                      
                    uint offsetX0Y0Z0 = index[0] + offsetY0Z0;
                    uint offsetXMY0Z0 = index[0] - 1 + offsetY0Z0;
                    uint offsetXPY0Z0 = index[0] + 1 + offsetY0Z0;
                    uint offsetX0YMZ0 = index[0] + offsetYMZ0;
                    uint offsetX0YPZ0 = index[0] + offsetYPZ0;
                    uint offsetX0Y0ZM = index[0] + offsetY0ZM;
                    uint offsetX0Y0ZP = index[0] + offsetY0ZP;

                    //COMPUTE_FINITE_DIFF                                                                                                             
                    if (index[0] == 0)
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.x; }
                    else if (index[0] == dimsMinus1[0])
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetX0Y0Z0].v - vec[offsetXMY0Z0].v) * reciprocalSpacing.x; }
                    else
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetXMY0Z0].v) * halfReciprocalSpacing.x; }

                    if (index[1] == 0)
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.y; }
                    else if (index[1] == dimsMinus1[1])
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0Y0Z0].v - vec[offsetX0YMZ0].v) * reciprocalSpacing.y; }
                    else
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0YMZ0].v) * halfReciprocalSpacing.y; }

                    if (index[2] == 0)
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.z; }
                    else if (index[2] == dimsMinus1[2])
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0Z0].v - vec[offsetX0Y0ZM].v) * reciprocalSpacing.z; }
                    else
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0ZM].v) * halfReciprocalSpacing.z; }
                }
            }
        }

        // Compute derivatives for +X boundary.
        index[0] = dimsMinus1[0];
        for (index[2] = 0; index[2] < dims[2]; ++index[2])
        {
            //ASSIGN_Z_OFFSETS
            uint offsetZM = numXY * (index[2] - 1);
            uint offsetZ0 = numXY * index[2];
            uint offsetZP = numXY * (index[2] + 1);

            for (index[1] = 0; index[1] < dims[1]; ++index[1])
            {
                //ASSIGN_YZ_OFFSETS                                                   
                uint offsetYMZ0 = dims[0] * (index[1] - 1) + offsetZ0;
                uint offsetY0Z0 = dims[0] * index[1] + offsetZ0;
                uint offsetYPZ0 = dims[0] * (index[1] + 1) + offsetZ0;
                uint offsetY0ZM = dims[0] * index[1] + offsetZM;
                uint offsetY0ZP = dims[0] * index[1] + offsetZP;

                {
                    //ASSIGN_XYZ_OFFSETS                                      
                    uint offsetX0Y0Z0 = index[0] + offsetY0Z0;
                    uint offsetXMY0Z0 = index[0] - 1 + offsetY0Z0;
                    uint offsetXPY0Z0 = index[0] + 1 + offsetY0Z0;
                    uint offsetX0YMZ0 = index[0] + offsetYMZ0;
                    uint offsetX0YPZ0 = index[0] + offsetYPZ0;
                    uint offsetX0Y0ZM = index[0] + offsetY0ZM;
                    uint offsetX0Y0ZP = index[0] + offsetY0ZP;

                    //COMPUTE_FINITE_DIFF                                                                                                             
                    if (index[0] == 0)
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.x; }
                    else if (index[0] == dimsMinus1[0])
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetX0Y0Z0].v - vec[offsetXMY0Z0].v) * reciprocalSpacing.x; }
                    else
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetXMY0Z0].v) * halfReciprocalSpacing.x; }

                    if (index[1] == 0)
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.y; }
                    else if (index[1] == dimsMinus1[1])
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0Y0Z0].v - vec[offsetX0YMZ0].v) * reciprocalSpacing.y; }
                    else
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0YMZ0].v) * halfReciprocalSpacing.y; }

                    if (index[2] == 0)
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.z; }
                    else if (index[2] == dimsMinus1[2])
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0Z0].v - vec[offsetX0Y0ZM].v) * reciprocalSpacing.z; }
                    else
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0ZM].v) * halfReciprocalSpacing.z; }
                }
            }
        }


        // Compute derivatives for +Y boundary.
        index[1] = dimsMinus1[1];
        for (index[2] = 0; index[2] < dims[2]; ++index[2])
        {
            //ASSIGN_Z_OFFSETS
            uint offsetZM = numXY * (index[2] - 1);
            uint offsetZ0 = numXY * index[2];
            uint offsetZP = numXY * (index[2] + 1);

            {
                //ASSIGN_YZ_OFFSETS                                                   
                uint offsetYMZ0 = dims[0] * (index[1] - 1) + offsetZ0;
                uint offsetY0Z0 = dims[0] * index[1] + offsetZ0;
                uint offsetYPZ0 = dims[0] * (index[1] + 1) + offsetZ0;
                uint offsetY0ZM = dims[0] * index[1] + offsetZM;
                uint offsetY0ZP = dims[0] * index[1] + offsetZP;

                for (index[0] = 0; index[0] < dims[0]; ++index[0])
                {
                    //ASSIGN_XYZ_OFFSETS                                      
                    uint offsetX0Y0Z0 = index[0] + offsetY0Z0;
                    uint offsetXMY0Z0 = index[0] - 1 + offsetY0Z0;
                    uint offsetXPY0Z0 = index[0] + 1 + offsetY0Z0;
                    uint offsetX0YMZ0 = index[0] + offsetYMZ0;
                    uint offsetX0YPZ0 = index[0] + offsetYPZ0;
                    uint offsetX0Y0ZM = index[0] + offsetY0ZM;
                    uint offsetX0Y0ZP = index[0] + offsetY0ZP;

                    //COMPUTE_FINITE_DIFF                                                                                                             
                    if (index[0] == 0)
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.x; }
                    else if (index[0] == dimsMinus1[0])
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetX0Y0Z0].v - vec[offsetXMY0Z0].v) * reciprocalSpacing.x; }
                    else
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetXMY0Z0].v) * halfReciprocalSpacing.x; }

                    if (index[1] == 0)
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.y; }
                    else if (index[1] == dimsMinus1[1])
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0Y0Z0].v - vec[offsetX0YMZ0].v) * reciprocalSpacing.y; }
                    else
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0YMZ0].v) * halfReciprocalSpacing.y; }

                    if (index[2] == 0)
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.z; }
                    else if (index[2] == dimsMinus1[2])
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0Z0].v - vec[offsetX0Y0ZM].v) * reciprocalSpacing.z; }
                    else
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0ZM].v) * halfReciprocalSpacing.z; }
                }
            }
        }

        // Compute derivatives for +Z boundary.
        index[2] = dimsMinus1[2];
        {
            //ASSIGN_Z_OFFSETS
            uint offsetZM = numXY * (index[2] - 1);
            uint offsetZ0 = numXY * index[2];
            uint offsetZP = numXY * (index[2] + 1);

            for (index[1] = 0; index[1] < dims[1]; ++index[1])
            {
                //ASSIGN_YZ_OFFSETS                                                   
                uint offsetYMZ0 = dims[0] * (index[1] - 1) + offsetZ0;
                uint offsetY0Z0 = dims[0] * index[1] + offsetZ0;
                uint offsetYPZ0 = dims[0] * (index[1] + 1) + offsetZ0;
                uint offsetY0ZM = dims[0] * index[1] + offsetZM;
                uint offsetY0ZP = dims[0] * index[1] + offsetZP;

                for (index[0] = 0; index[0] < dims[0]; ++index[0])
                {
                    //ASSIGN_XYZ_OFFSETS                                      
                    uint offsetX0Y0Z0 = index[0] + offsetY0Z0;
                    uint offsetXMY0Z0 = index[0] - 1 + offsetY0Z0;
                    uint offsetXPY0Z0 = index[0] + 1 + offsetY0Z0;
                    uint offsetX0YMZ0 = index[0] + offsetYMZ0;
                    uint offsetX0YPZ0 = index[0] + offsetYPZ0;
                    uint offsetX0Y0ZM = index[0] + offsetY0ZM;
                    uint offsetX0Y0ZP = index[0] + offsetY0ZP;

                    //COMPUTE_FINITE_DIFF                                                                                                             
                    if (index[0] == 0)
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.x; }
                    else if (index[0] == dimsMinus1[0])
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetX0Y0Z0].v - vec[offsetXMY0Z0].v) * reciprocalSpacing.x; }
                    else
                    { jacobian[offsetX0Y0Z0].x = (vec[offsetXPY0Z0].v - vec[offsetXMY0Z0].v) * halfReciprocalSpacing.x; }

                    if (index[1] == 0)
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.y; }
                    else if (index[1] == dimsMinus1[1])
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0Y0Z0].v - vec[offsetX0YMZ0].v) * reciprocalSpacing.y; }
                    else
                    { jacobian[offsetX0Y0Z0].y = (vec[offsetX0YPZ0].v - vec[offsetX0YMZ0].v) * halfReciprocalSpacing.y; }

                    if (index[2] == 0)
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0Z0].v) * reciprocalSpacing.z; }
                    else if (index[2] == dimsMinus1[2])
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0Z0].v - vec[offsetX0Y0ZM].v) * reciprocalSpacing.z; }
                    else
                    { jacobian[offsetX0Y0Z0].z = (vec[offsetX0Y0ZP].v - vec[offsetX0Y0ZM].v) * halfReciprocalSpacing.z; }
                }
            }
            #endregion
        }








    }
}
