from loadmodules import *

def calculate_value_relative_to_vector(loaded_snap, value, vector):
    if value not in loaded_snap.data.keys() or len(loaded_snap.data[value].shape) <= 1:
        print("cannot compute ", value, " relative to the motion axis")
        return
    else:
        print("computing ", value, "relative to ", vector)
        vector_size = np.sqrt((vector ** 2).sum())
        loaded_snap.data[value+"_v"] = (loaded_snap.data[value] * vector).sum(axis=1) / vector_size
        loaded_snap.data[value+"_u"] = np.sqrt((loaded_snap.data[value] ** 2).sum(axis=1) -
                                               loaded_snap.data[value+"_v"]**2)

        return loaded_snap

def project_vector(v,r):
    dist = np.sqrt((r*r).sum(axis=1))
    return ((r*v).sum(axis=1)) / dist