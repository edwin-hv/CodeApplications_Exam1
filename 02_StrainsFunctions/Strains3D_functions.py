import numpy as np
import sympy as sp

def transform_strain_tensor(strain_tensor, direction_cosines):
    """
    Transform a strain tensor given in certain reference system to other strain tensor, 
    in other system of reference given for direction cosines.

    Parameters
    ----------
    -  strain_tensor(np.array) : array with the original strain tensor
    -  direction_cosines(np.array) : array with the direction cosines which describe de position 
                                            of the new coordinate system regarding the original cosine system 
                                            defined as the text book

    Returns
    -------
     - transformed_strain_tensor(np.array) : stress tensor in the new coordinate system
    """
    # Extract rows as vectors
    l1_m1_n1 = direction_cosines[0, :].T  # First row
    l2_m2_n2 = direction_cosines[1, :].T  # Second row
    l3_m3_n3 = direction_cosines[2, :].T  # Third row
    
    x_prime_plane_strains = direction_cosines @ strain_tensor @  l1_m1_n1 
    y_prime_plane_strains = direction_cosines @ strain_tensor @  l2_m2_n2 
    z_prime_plane_strains = direction_cosines @ strain_tensor @  l3_m3_n3 
    print("x'strains are:")
    display( x_prime_plane_strains)
    print("y'strains are:")
    display( y_prime_plane_strains)
    print("z'strains are:")
    display( z_prime_plane_strains)
    
    Ex_prime  = x_prime_plane_strains[0]
    O5_Gxy_prime = x_prime_plane_strains[1]
    O5_Gxz_prime = x_prime_plane_strains[2]
    Ey_prime  = y_prime_plane_strains[1]
    O5_Gyz_prime = y_prime_plane_strains[2]
    Ez_prime  = z_prime_plane_strains[2]
  
    

    transformed_strain_tensor = np.array([[Ex_prime, O5_Gxy_prime, O5_Gxz_prime],
                                           [O5_Gxy_prime, Ey_prime, O5_Gyz_prime],
                                           [O5_Gxz_prime, O5_Gyz_prime, Ez_prime]])
 
    return transformed_strain_tensor


def get_strain_invariants(strain_tensor):
    """
    compute strain invariants given a strain tensor.

    Parameters
    ----------
    -  strain_tensor(np.array) : array with the strain tensor

    Returns
    -------
     - invariants(list) : with I1, I2, I3 
    """
    Ex  = strain_tensor[0][0]
    Gxy = 2 * strain_tensor[0][1]
    Gxz = 2 * strain_tensor[0][2]
    Ey  = strain_tensor[1][1]
    Gyz = 2 * strain_tensor[1][2]
    Ez  = strain_tensor[2][2]
    
    I1 = Ex + Ey + Ez
    I2 = Ex*Ey + Ey*Ez  + Ex*Ez  - 0.25*((Gxy**2) + (Gyz**2) + (Gxz**2))
    I3 = np.linalg.det(strain_tensor)
    print("I1:")
    display(I1)
    print("I2:")
    display(I2)
    print("I3:")
    display(I3)

    invariants = [I1, I2,I3 ]

    return invariants

def get_principal_strains(invariants):
    """
    compute principal strains E1,E2,E3 from the invariants belonging to a certain strain state
    Parameters
    ----------
    -  invariants(list) : with I1, I2, I3 ordered 

    Returns
    -------
     - principal_strains(list) : with S1, S2, S3 
    """

    I1 = invariants[0]
    I2 = invariants[1]
    I3 = invariants[2]
    print(I1,I2,I3)

    # a = 1  # coeficiente de x^3 
    # b = -I1  # coeficiente de x^2
    # c =  I2  # coeficiente de x
    # d = -I3  # término constante

    coeficientes = [1, -I1, I2, -I3 ]
    pincipal_strains = np.roots(coeficientes)

    return  pincipal_strains


def get_strains_from_stresses(stress_tensor,E,G,v):
    """
    compute strain tensor for a given stress tensor by using constitutive relations
    Parameters
    ----------
    -  stress_tensor(np.array) : stress tensor given for the point
    -  E(int) : Elasticity modulous of the material
    -  G(int) : Shear modulous of the material
    -  v(float) : Pisson ratio
    
    Returns
    -------
     - strains_tensor(np.array) : strain tensor for the point
    """
    
    
    Sx  = stress_tensor[0][0]
    Txy = stress_tensor[0][1]
    Txz = stress_tensor[0][2]
    Sy  = stress_tensor[1][1]
    Tyz = stress_tensor[1][2]
    Sz  = stress_tensor[2][2]

    Ex = (1/E) *( Sx - v*(Sy + Sz))
    Ey = (1/E) *( Sy - v*(Sx + Sz))
    Ez = (1/E) *( Sz - v*(Sx + Sy))
    Gxy = Txy/G
    Gyz = Tyz/G
    Gxz = Txz/G
    
    
    strain_tensor =  np.array([[Ex      , 0.5*Gxy , 0.5*Gxz  ],
                               [0.5*Gxy , Ey      , 0.5*Gyz  ],
                               [0.5*Gxz , 0.5*Gyz , Ez       ]]) 
    return(strain_tensor)