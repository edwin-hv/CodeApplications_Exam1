import numpy as np
import sympy as sp

def transform_stress_tensor(stress_tensor, direction_cosines):
    """
    Transform a stress tensor given in certain reference system to other stress tensor, 
    in other system of reference given for direction cosines.

    Parameters
    ----------
    -  stress_tensor(np.array) : array with the original stress tensor
    -  direction_cosines(np.array) : array with the direction cosines which describe de position 
                                            of the new coordinate system regarding the original cosine system 
                                            defined as the text book

    Returns
    -------
     - transformed_stress_tensor(np.array) : stress tensor in the new coordinate system
    """
    # Extract rows as vectors
    l1_m1_n1 = direction_cosines[0, :].T  # First row
    l2_m2_n2 = direction_cosines[1, :].T  # Second row
    l3_m3_n3 = direction_cosines[2, :].T  # Third row
    
    x_prime_plane_streses = direction_cosines @ stress_tensor @  l1_m1_n1 
    y_prime_plane_streses = direction_cosines @ stress_tensor @  l2_m2_n2 
    z_prime_plane_streses = direction_cosines @ stress_tensor @  l3_m3_n3 
    print("x'stresses are:")
    display( x_prime_plane_streses)
    print("y'stresses are:")
    display( y_prime_plane_streses)
    print("z'stresses are:")
    display( z_prime_plane_streses)
    
    Sx_prime  = x_prime_plane_streses[0]
    Txy_prime = x_prime_plane_streses[1]
    Txz_prime = x_prime_plane_streses[2]
    Sy_prime  = y_prime_plane_streses[1]
    Tyz_prime = y_prime_plane_streses[2]
    Sz_prime  = z_prime_plane_streses[2]
  
    

    transformed_stress_tensor = np.array([[Sx_prime, Txy_prime, Txz_prime],
                                           [Txy_prime, Sy_prime, Tyz_prime],
                                           [Txz_prime, Tyz_prime, Sz_prime]])
 
    return transformed_stress_tensor


def get_stress_invariant(stress_tensor):
    """
    compute stress invariant given a stress tensor.

    Parameters
    ----------
    -  stress_tensor(np.array) : array with the original stress tensor

    Returns
    -------
     - invariants(list) : with I1, I2, I3 ordered from biggest to smallest
    """
    Sx  = stress_tensor[0][0]
    Txy = stress_tensor[0][1]
    Txz = stress_tensor[0][2]
    Sy  = stress_tensor[1][1]
    Tyz = stress_tensor[1][2]
    Sz  = stress_tensor[2][2]
    
    I1 = Sx + Sy + Sz
    I2 = Sx*Sy + Sy*Sz  + Sx*Sz - (Txy**2) - (Tyz**2) - (Txz**2)
    I3 = np.linalg.det(stress_tensor)
    print("I1:")
    display(I1)
    print("I2:")
    display(I2)
    print("I3:")
    display(I3)

    invariants = [I1, I2,I3 ]

    return invariants


def get_principal_stresses(invariants):
    """
    compute principal stresses S1,S2,S3 from the invariants belonging to a certain stress state
    Parameters
    ----------
    -  invariants(list) : with I1, I2, I3 ordered 

    Returns
    -------
     - principal_stresses(list) : with S1, S2, S3 ordered from biggest to smallest
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
    pincipal_stresses = np.roots(coeficientes)

    return  pincipal_stresses


def get_principal_directions(stress_tensor,principal_stress):
    """
    compute principal direction for a princiapl stress of a certain stresses state.

    Parameters
    ----------
    -  invariants(list) : with I1, I2, I3 ordered 

    Returns
    -------
     - principal_stresses(list) : with S1, S2, S3 ordered from biggest to smallest
    """

    Sx = stress_tensor[0][0]
    Sy = stress_tensor[1][1]
    Sz = stress_tensor[2][2]
    
    Sx_modified = Sx - principal_stress
    Sy_modified = Sy - principal_stress
    Sz_modified = Sz - principal_stress
    
    stress_tensor_modified = np.copy(stress_tensor)
    np.fill_diagonal( stress_tensor_modified , [ Sx_modified , Sy_modified,   Sz_modified ])

    
    # Calcular la descomposición en valores singulares (SVD)
    U, S, Vt = np.linalg.svd( stress_tensor_modified)
 
    # El rango de A es el número de valores singulares no nulos
    rank = np.sum(S > 1e-10)  # Ajusta el umbral según sea necesario
    
    # El número de soluciones libres es 3 - rango
    n = stress_tensor_modified.shape[1] - rank
    
    if n > 0:
        # Para obtener la solución general del sistema homogéneo
        # Las columnas de Vt correspondientes a los valores singulares cero son las soluciones
        null_space = Vt[rank:].T  # Espacio nulo
        print("El espacio nulo (soluciones) es:")
        print(null_space)
    else:
        print("El sistema no tiene soluciones no triviales.")

    return  


