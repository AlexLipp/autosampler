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
    indices = np.ravel([np.where(x == e) for e in elements])
    return np.delete(x, indices)


def get_sample_nodes_by_area(model_grid, target_area):

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
    samps_allocated = 0
    counter = 0
    print("Target area = ", nodes_per_samp)
    print("Looping through all unvisited nodes upstream to downstream")
    initial_len = len(uV)
    visited = 0
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
        unvis_up_nodes = [
            node
        ]  # List of unvisited nodes that are in the catchment of node in question
        # Loop through unvisited nodes upstream
        for new_up_node in uV[i + 1 :]:
            if receiver_at_node[new_up_node] in unvis_up_nodes:
                unvis_up_nodes.append(new_up_node)
            else:
                break
        if len(unvis_up_nodes) > nodes_per_samp:
            print("\t * Found a sample locality *")
            sample_nodes[node] = unvis_up_nodes  # Add node to list
            uV = fast_delete(uV, unvis_up_nodes)  # Remove visited nodes from list
            samps_allocated += 1
        counter += 1
    print("Found", samps_allocated, "sample localities")
    area_sizes = [len(areas) * cell_area for _, areas in sample_nodes.items()]
    mean, std = np.mean(area_sizes), np.std(area_sizes)
    print("Average area per basin = ", mean, "+/-", std)
    print("~~~~~~~ Finished Calculation ~~~~~~~")
    return sample_nodes


def viz_drainage_area(model_grid):
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
    if os.path.exists("area_IDs.asc"):  # Allows over-writing of .asc files
        os.remove("area_IDs.asc")
    _ = model_grid.add_field("area_IDs", areas)
    model_grid.save("area_IDs.asc", names="area_IDs")


def viz_results(locs, areas, model_grid):
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
