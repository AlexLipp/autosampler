import os
import sys
from datetime import datetime

import landlab
import matplotlib.pyplot as plt
import numpy as np
from landlab.components import FlowAccumulator, SinkFillerBarnes
from matplotlib.colors import LogNorm


# TODO: Add doc strings to everything
# TODO: Generate a readME


def fast_delete(x, elements):
    """Removes `elements` from numpy array by value faster than list.remove
    Args:
        x (array): The array to remove elements from.
        elements (array like): The elements to remove from x
    Returns:
        (array): The array with elements removed
    """
    indices = np.ravel([np.where(x == e) for e in elements])
    return np.delete(x, indices)


def get_sample_nodes_by_area(model_grid, target_area):
    """Finds sample sites which best divide DEM into ~equal sub-catchments
    Args:
        model_grid (RasterModelGrid): The landscape to partition. This is
        a LandLab RasterModelGrid object that must have flow routed across it using D8 method.
        It is recommended that sinks are filled too to allow for continuous flow paths.

        target_area: The area (in units of the model grid NOT number of nodes) which basins should be
        larger than. Note that this is a *minimum* value and basins may be this size (or greater)
        but no smaller.

    Returns:
        {int: list of ints}: Dictionary, the keys of which are the node IDs of the identified sample sites.
        The items for each key are the IDs of the nodes in the subcatchment delineated by that sample.
        Node IDs can be turned into coordinates using `process_output_dict` or `np.unravel_index`.
    """

    print("~~~~~~~ Beginning Calculation ~~~~~~~")
    # Node array contining downstream-to-upstream ordered list of node
    ordered_nodes = model_grid.at_node["flow__upstream_node_order"]
    receiver_at_node = model_grid.at_node["flow__receiver_node"]
    cell_area = model_grid.dx * model_grid.dy
    nodes_per_samp = target_area / cell_area

    uV = ordered_nodes.copy()  # unvisited nodes
    print("Removing catchments smaller than target area")
    # Remove nodes from uV if they are within catchments smaller than target (as they will never be visited)
    num_sinks = len(np.where(model_grid.at_node["flow__sink_flag"])[0])
    counter = 0
    for sink in np.where(model_grid.at_node["flow__sink_flag"])[0]:
        if counter % 1000 == 0:
            print(
                "\t Processed sink",
                counter,
                "of",
                num_sinks,
                ",",
                datetime.now().strftime("%H:%M:%S"),
            )
        if model_grid.at_node["drainage_area"][sink] < target_area:
            i = np.where(ordered_nodes == sink)[0][0]
            bad_nodes = [sink]
            for up_node in ordered_nodes[i + 1 :]:
                if receiver_at_node[up_node] in bad_nodes:
                    bad_nodes.append(up_node)
                else:
                    uV = fast_delete(uV, bad_nodes)
                    break
        counter += 1
    sample_nodes = {}
    counter = 0
    print("Target area = ", nodes_per_samp)
    print("Looping through all unvisited nodes upstream to downstream")
    initial_len = len(uV)
    # Iterate through nodes from upstream to downstream
    for i in np.arange(initial_len - 1, -1, -1):
        if counter % 1000 == 0:
            print(
                "\t Processing node",
                counter,
                "of",
                initial_len,
                ",",
                datetime.now().strftime("%H:%M:%S"),
            )
        node = uV[i]  # Node in network
        # Initiate list of unvisited nodes that are in the upstream catchment of node in question
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
        # If number of nodes in new subcatchment greater than threshold we add it to output
        if len(unvis_up_nodes) > nodes_per_samp:
            print("\t * Found a sample locality *")
            sample_nodes[
                node
            ] = unvis_up_nodes  # Add node to list with corresponding catchment
            uV = fast_delete(
                uV, unvis_up_nodes
            )  # Remove the new catchment  from array of unvisited nodes
        counter += 1
    print("Found", len(sample_nodes.keys()), "sample localities")
    area_sizes = [len(areas) * cell_area for _, areas in sample_nodes.items()]
    mean, std = np.mean(area_sizes), np.std(area_sizes)
    print("Average area per basin = ", mean, "+/-", std)
    print("~~~~~~~ Finished Calculation ~~~~~~~")
    return sample_nodes


def viz_drainage_area(model_grid):
    """ "Visualises drainage area logarithmically.
    Args:
        model_grid (RasterModelGrid): LandLab grid with drainage routed across it.

    Returns:
        None

    Produces instance of matplotlib.plt"""

    plt.imshow(
        model_grid.at_node["drainage_area"].reshape(model_grid.shape),
        norm=LogNorm(),
        origin="lower",
    )
    cb = plt.colorbar()
    cb.set_label("Drainage Area")
    plt.xlabel("x")
    plt.ylabel("y")


def process_output_dict(node_catchment_dict, model_grid):
    """Reformats the output dictionary.
    Args:
        node_catchment_dict ({int: list of ints}): Output from `get_sample_nodes_by_area`
        model_grid (RasterModelGrid): LandLab grid with drainage routed across it.
    Returns:
        np.array(N,3): 3 column table output with cols: Sample ID, x-coordinate, y-coordinate. Has
        number of rows equal to number of allocated samples.

        np.array(nx,ny): Map of the identified sub-catchments, with ID corresponding to the Sample ID in
        associated table. 2D array with same dimensions as model grid/DEM. Areas not covered by sample sites
        given the NaN value of -999.

    """

    out_area = np.zeros(model_grid.shape).flatten() - 999
    N = 1
    Ns, xs, ys = [], [], []
    for node, upst_nodes in node_catchment_dict.items():
        x, y = np.unravel_index(node, model_grid.shape)
        Ns += [N]
        xs += [x * model_grid.dx]
        ys += [y * model_grid.dy]
        out_area[upst_nodes] = N
        N += 1
    out_area = out_area.reshape(model_grid.shape)
    return (np.array([Ns, xs, ys]).T, out_area)


def save_results(locs, areas, model_grid):
    """Saves output as files with appropriate names.
    Args:
    locs (np.array(N,3)): 3 column table output with cols: Sample ID, x-coordinate, y-coordinate. Has
            number of rows equal to number of allocated samples. See `process_output_dict`

    areas (np.array(nx,ny)): Map of the identified sub-catchments, with ID corresponding to the Sample ID in
        associated table. 2D array with same dimensions as model grid/DEM. Areas not covered by sample sites
        given the NaN value of -999. See `process_output_dict`

    model_grid (RasterModelGrid): LandLab grid with drainage routed across it.

    Returns:
        None

    Produces two files:
        1. "sample_sites.csv" a file containing the sample site localities, given in `locs` input
        2. "area_IDs.asc" a map of the delineated sub-catchments, given in `areas` input. This is an
        ESRI ASCII raster file format appropriate for use in most GIS software.

    """
    np.savetxt(
        "sample_sites.csv", X=locs, delimiter=",", header="Area ID, x, y", comments=""
    )
    if os.path.exists("area_IDs.asc"):  # Allows over-writing of .asc files
        os.remove("area_IDs.asc")
    _ = model_grid.add_field("area_IDs", areas)
    model_grid.save("area_IDs.asc", names="area_IDs")


def viz_results(locs, areas, model_grid):
    """Visuaslises identified sample localities and associated sub-catchments.
    Args:
        locs (np.array(N,3)): 3 column table output with cols: Sample ID, x-coordinate, y-coordinate. Has
                number of rows equal to number of allocated samples. See `process_output_dict`

        areas (np.array(nx,ny)): Map of the identified sub-catchments, with ID corresponding to the Sample ID in
            associated table. 2D array with same dimensions as model grid/DEM. Areas not covered by sample sites
            given the NaN value of -999. See `process_output_dict`

        model_grid (RasterModelGrid): LandLab grid with drainage routed across it.

    Returns:
        None

    Produces instance of matplotlib.plt
    """
    plt.figure(figsize=(12, 5))
    areas[areas < 0] = np.nan
    plt.subplot(1, 2, 1)
    plt.imshow(areas, origin="lower", cmap="nipy_spectral")
    cb = plt.colorbar()
    cb.set_label("Area ID")
    plt.title("Sample areas")
    plt.scatter(x=locs[:, 2], y=locs[:, 1], c="black", marker="x", s=50)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.subplot(1, 2, 2)
    viz_drainage_area(model_grid=model_grid)
    plt.scatter(x=locs[:, 2], y=locs[:, 1], c="red", marker="x", s=50)
    plt.tight_layout()


def process_topo(path):
    """Turns a topographic data file (as .asc) into a LandLab model grid with drainage routed across it. 
    
    Args: 
        path (string): Path to the ESRI ASCII file. 
    Returns: 
        RasterModelGrid: An initialised LandLab grid with sinks filled and drainage routed across it. 
    """

    ascii_data = landlab.io.esri_ascii.read_esri_ascii(path)
    model_grid = ascii_data[0]
    model_grid.add_field("topographic__elevation", ascii_data[1])
    print("Filling sinks (can be slow)")
    sb = SinkFillerBarnes(model_grid, ignore_overfill=True)
    sb.run_one_step()
    print("Running flow-routing")
    frr = FlowAccumulator(model_grid, flow_director="FlowDirectorD8")
    frr.run_one_step()
    return model_grid


def main():
    print("Loading in topographic data...")
    path_to_topo = sys.argv[1]
    area_per_basin = float(sys.argv[2])
    mg = process_topo(path_to_topo)
    viz_drainage_area(mg)
    plt.show()
    sample_nodes_catchments = get_sample_nodes_by_area(mg, area_per_basin)
    localities, node_map = process_output_dict(sample_nodes_catchments, mg)
    save_results(localities, node_map, mg)
    viz_results(localities, node_map, mg)
    plt.show()


if __name__ == "__main__":
    main()
