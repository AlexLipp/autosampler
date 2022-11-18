import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from auto_catchments import fast_delete, process_topo


def snap_to_drainage(model_grid, sample_sites, threshold):
    """Moves sample sites to nearest channel node with area greater than threshold
    Args:
        model_grid (RasterModelGrid): Model-grid with flowrouted across it
        (see auto_catchments.process_topo

        sample_sites (Nx2 float array): Coordinates [x,y] (in units of the base DEM, e.g., km) of sampled
        localities. These will not necessarily lie on a modelled channel in the DEM.

        threshold (float): Drainage area (in units of base DEM) threshold, above which we define
        a node to be a `channel'. Nodes with area less than this cannot be snapped to.

    Returns:

        snapped (Nx2 float array): Coordinates [x,y] of the snapped, drainage aligned sample sites
        (in units of the base DEM, e.g., km).
    """

    channels = model_grid.at_node["drainage_area"] > threshold
    channels_y_ind, channels_x_ind = np.unravel_index(
        np.ravel(np.where(channels)), model_grid.shape
    )
    channels_xy = np.array([channels_x_ind * model_grid.dx, channels_y_ind * model_grid.dy])

    snapped = np.zeros(sample_sites.shape)
    for i in range(sample_sites.shape[0]):
        sample = sample_sites[i, :]
        diff = channels_xy - sample[:, np.newaxis]
        closest_channel_coords = channels_xy[:, np.argmin(np.sum(diff**2, axis=0))]
        snapped[i, :] = closest_channel_coords
    return snapped


def coords_to_node_ids(model_grid, coordinates):
    """Converts true coordinates into node IDs
    args:
        model_grid (RasterModelGrid): Model-grid with flowrouted across it
        (see auto_catchments.process_topo

        coordinates (Nx2 float array): Coordinates [x,y] for which nodes are being found

    returns:
        nodes: corresponding node IDs for each coordinate pair"""

    n_samples = coordinates.shape[0]
    nodes = np.zeros(n_samples)
    for i in np.arange(n_samples):
        x, y = coordinates[i, 0], coordinates[i, 1]
        x_ind = int(np.floor(x / model_grid.dx))
        y_ind = int(np.floor(y / model_grid.dy))
        nodes[i] = np.ravel_multi_index((y_ind, x_ind), model_grid.shape)
    return nodes


def get_subcatchments(model_grid, target_nodes, sample_names):
    """Subdivides a drainage network into subcatchments based on provided sample nodes

    args:
        model_grid (RasterModelGrid): Model-grid with flowrouted across it
        (see auto_catchments.process_topo

        target_nodes (int array): IDs of the sample nodes

        sample_names (str array): The sample name for each sample site.
        This *MUST* have the same order as target_nodes

    Returns:
        sample_catchments (nested dict{{}}): Returns a nested dictionary structure
        describing the information about the sub-catchment. The key 'areaID' is the
        integer to which each sub-catchment is assigned going from 1 (upstream) to N
        (most downstream). sample_catchments[areaID] returns a dictionary which stores:
        1) the SampleName for that sample site (using key "SampleName"); 2) the nodeIDs of
        the catchment nodes (using key "catchment_nodes") and 3) the nodeID of the sample
        site (using key "node")"""

    ordered_nodes = model_grid.at_node["flow__upstream_node_order"]
    receiver_at_node = model_grid.at_node["flow__receiver_node"]
    uV = ordered_nodes.copy()  # unvisited nodes
    sample_catchments = {}
    initial_len = len(uV)
    area_ID = 1
    # Iterate through nodes from upstream to downstream
    for i in np.arange(initial_len - 1, -1, -1):
        node = uV[i]  # Node in network
        if node in target_nodes:
            sample_name = sample_names[np.ravel(np.where(target_nodes == node))][0]
            unvis_up_nodes = [node]
            # Loop through unvisited nodes upstream
            for new_up_node in uV[i + 1 :]:
                # If this node drains to a node in our list we add it to the list
                # Hence we progressively go upstream building the subcatchment
                if receiver_at_node[new_up_node] in unvis_up_nodes:
                    unvis_up_nodes.append(new_up_node)
                # When we find a node not in the subcatchment we stop as we have reached a drainage divide
                else:
                    break
            sample_catchments[area_ID] = {
                "SampleName": sample_name,
                "catchment_nodes": unvis_up_nodes,
                "node": node,
            }  # Append catchment to dictionary
            uV = fast_delete(uV, unvis_up_nodes)  # Remove visited nodes from array
            area_ID += 1  # Update new node ID
    return sample_catchments


def get_subcatchment_map(model_grid, sample_catchment_dict):
    """Generates a map (2D numpy array) of the identified sub-catchments
    args:
        model_grid (RasterModelGrid): Model-grid with flowrouted across it
        (see auto_catchments.process_topo

        sample_catchment_dict (dict): Dictionary defining sub-catchments
        (see `get_subcatchments`)

    Returns:
        out (NxM int array): Map of the areaIDs for each sub-catchment.
        See output from `save_outputs` and `get_subcatchments` for how these relate
        to the sample sites.
    """

    out = np.zeros(model_grid.shape).flatten()
    for label, catchment_dict in sample_catchment_dict.items():
        out[catchment_dict["catchment_nodes"]] = label
    out = out.reshape(model_grid.shape)
    return out


def load_sample_data(path_to_file):
    """Loads sample site data into memory:
    args:
        path_to_file (str): Path to file which contains the coordinates of
        sample sites with corresponding sample names.
        This file must have format:

        SampleName | x coordinate | y coordinate
        -----------------------------------------
           string  |    float     |     float"""

    data = pd.read_csv(path_to_file).to_numpy()
    sample_sites = data[:, 1:].astype(float)
    sample_names = data[:, 0].astype(str)
    return sample_names, sample_sites


def node_to_coord(model_grid, node):
    """Converts node IDs into true coordinates
    args:
        model_grid (RasterModelGrid): Model-grid with flowrouted across it
        (see auto_catchments.process_topo

        node (int)): corresponding node IDs for each coordinate pair

    returns:
        coordinates (tuple): Coordinates (y,x) for the input node ID"""
    y, x = np.unravel_index(node, model_grid.shape)
    coordindates = y * model_grid.dy, x * model_grid.dx
    return coordindates


def save_outputs(model_grid, sample_catchment_dict, catchment_map):
    """Save outputs from sample-processing to disk in appropriate format

    args:
        model_grid (RasterModelGrid): Model-grid with flowrouted across it
        (see auto_catchments.process_topo

        sample_catchment_dict (dict): Dictionary defining sub-catchments
        (see `get_subcatchments`)

        catchment_map (NxM int array): Map of the areaIDs for each sub-catchment.
        See output from `save_outputs` and `get_subcatchments` for how these relate
        to the sample sites.

    Returns:
        None:

    Generates two output files.
    1. "fitted_localities.csv" is a table of the coordinates of the *fitted* sample sites,
        alongside the associated sample name, and the AreaID of the sub-catchment. This areaID
        can be used to interpret/interoperate with the 'area_IDs.asc' output.
    2. "area_IDs.asc" is a map of the delineated sub-catchments. The area_IDs map onto the
        sample sites given in "fitted_localities.csv". This is an ESRI ASCII raster file format
        appropriate for use in most GIS software.
    """

    node_IDs = []
    sample_names = []
    xs = []
    ys = []
    for label, catchment_dict in sample_catchment_dict.items():
        node_IDs.append(label)  # node ID
        sample_names.append(catchment_dict["SampleName"])
        x, y = node_to_coord(model_grid, catchment_dict["node"])
        xs.append(x)
        ys.append(y)

    outdf = pd.DataFrame({"Catchment ID": node_IDs, "SampleName": sample_names, "x": xs, "y": ys})
    outdf.to_csv("fitted_localities.csv", index=False)

    if os.path.exists("area_IDs.asc"):  # Allows over-writing of .asc files
        os.remove("area_IDs.asc")
    _ = model_grid.add_field("area_IDs", catchment_map)
    model_grid.save("area_IDs.asc", names="area_IDs")


def main():
    path_to_topo = sys.argv[1]
    path_to_samples = sys.argv[2]
    catchment_threshold = float(sys.argv[3])
    print("Loading in topographic data...")
    mg = process_topo(path_to_topo)
    sample_names, sample_sites = load_sample_data(path_to_samples)
    snapped_localities = snap_to_drainage(mg, sample_sites, catchment_threshold)
    target_nodes = coords_to_node_ids(mg, snapped_localities)
    sub_catchment_dict = get_subcatchments(mg, target_nodes, sample_names)
    catchment_map = get_subcatchment_map(mg, sub_catchment_dict)
    save_outputs(mg, sub_catchment_dict, catchment_map)
    plt.title("Subcatchments: Samples aligned to drainage")
    plt.imshow(catchment_map, origin="lower", cmap="nipy_spectral")
    cb = plt.colorbar()
    cb.set_label("Catchment ID")
    plt.show()


if __name__ == "__main__":
    main()
