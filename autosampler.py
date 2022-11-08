import sys

import landlab
import matplotlib.pyplot as plt
import numpy as np
from landlab.components import FlowAccumulator
from matplotlib.colors import LogNorm

# Load in topography
# Route drainage across it

# TODO: Try with real DEMs 
# TODO: Add doc strings to everything 
# TODO: Delete "area_IDs.asc" if present
# TODO: Generate a readME 

def get_sample_nodes(model_grid, num_samps, increment=0.01):
    print("Setting up")

    # Node array containing downstream-to-upstream ordered list of node
    ordered_nodes = model_grid.at_node["flow__upstream_node_order"]
    total_area = len(ordered_nodes)
    A_per_samp = total_area / num_samps

    node_area_dict = {}
    for node in ordered_nodes:
        node_area_dict[node] = landlab.utils.get_watershed_nodes(model_grid, node)

    samps_allocated = 0
    iterations = 0
    print("Starting iterations")
    while samps_allocated < num_samps:

        V = []  # visited nodes
        uV = ordered_nodes.tolist()  # unvisited otes
        sample_nodes = {}
        samps_allocated = 0

        print("Trying with target area = ", A_per_samp)
        for i in np.arange(len(ordered_nodes) - 1, -1, -1):
            # Iterate from upstream to downstream
            node = ordered_nodes[i]  # Node in network
            all_upst_nodes = node_area_dict[node]  # Get all the upstream nodes
            unvisited_upst_nodes = all_upst_nodes[
                np.isin(all_upst_nodes, uV)
            ]  # Identify those which are 'new'
            if len(unvisited_upst_nodes) > A_per_samp:
                print("Found a sample locality:", node)
                sample_nodes[node] = unvisited_upst_nodes  # Add node to list
                for x in unvisited_upst_nodes:  # Remove added nodes to unvisited list
                    uV.remove(x)
                samps_allocated += 1
                if samps_allocated == num_samps:
                    break
        iterations += 1
        print("# Finished iteration", iterations, "#")
        if samps_allocated < num_samps:
            print("Only allocated:", samps_allocated, "of", num_samps)
            print("Trying again with smaller target area...")
            print("___________________________")
        else:
            print("All samples allocated")
            print("Finishing...")
        A_per_samp = A_per_samp * (
            1 - increment
        )  # initial target number of nodes per samp
    return sample_nodes


def viz_drainage_area(model_grid):
    plt.imshow(
        model_grid.at_node["drainage_area"].reshape(model_grid.shape),
        norm=LogNorm(),
        origin="lower",
    )
    cb = plt.colorbar()
    cb.set_label("Drainage Area (# of nodes)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


def process_output_dict(node_catchment_dict, model_grid):
    out_area = np.zeros(model_grid.shape).flatten() - 999
    N = 1
    Ns, xs, ys = [], [], []
    for node, upst_nodes in node_catchment_dict.items():
        x, y = np.unravel_index(node, model_grid.shape)
        Ns += [N]
        xs += [x]
        ys += [y]
        out_area[upst_nodes] = N
        N += 1
    out_area = out_area.reshape(model_grid.shape)
    return (np.array([Ns, xs, ys]).T, out_area)


def save_results(locs, areas, model_grid):
    np.savetxt(
        "sample_sites.csv", X=locs, delimiter=",", header="Area ID, x, y", comments=""
    )
    _ = model_grid.add_field("area_IDs", areas)
    model_grid.save("area_IDs.asc", names="area_IDs")


def viz_results(locs, areas):
    areas[areas < 0] = np.nan
    plt.imshow(areas)
    cb = plt.colorbar()
    cb.set_label("Area ID")
    plt.scatter(
        x=locs[:, 2], y=locs[:, 1], c="black", marker="x", s=50, facecolor="white"
    )
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


def main():
    print("Loading in topographic data...")
    path_to_topo = sys.argv[1]
    number_of_samps = int(sys.argv[2])
    ascii_data = landlab.io.esri_ascii.read_esri_ascii(path_to_topo)
    mg = ascii_data[0]
    mg.add_field("topographic__elevation", ascii_data[1])
    print("Running flow-routing")
    frr = FlowAccumulator(mg, flow_director="FlowDirectorD8")
    frr.run_one_step()
    viz_drainage_area(mg)
    sample_nodes_catchments = get_sample_nodes(mg, number_of_samps, increment=0.1)
    localities, node_map = process_output_dict(sample_nodes_catchments, mg)
    save_results(localities, node_map, mg)
    viz_results(localities, node_map)


if __name__ == "__main__":
    main()
