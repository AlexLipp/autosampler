import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
from auto_catchments import fast_delete, process_topo


def snap_to_drainage(model_grid, sample_sites, threshold):

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


def coords_to_node_ids(model_grid, sample_sites):
    n_samples = sample_sites.shape[0]
    target_nodes = np.zeros(n_samples)
    for i in np.arange(n_samples):
        x, y = sample_sites[i, 0], sample_sites[i, 1]
        x_ind = int(np.floor(x / model_grid.dx))
        y_ind = int(np.floor(y / model_grid.dy))
        target_nodes[i] = np.ravel_multi_index((y_ind, x_ind), model_grid.shape)
    return target_nodes


def get_subcatchments(model_grid, target_nodes, sample_names):
    # TODO Allow labels to be a string by default
    ordered_nodes = model_grid.at_node["flow__upstream_node_order"]
    receiver_at_node = model_grid.at_node["flow__receiver_node"]
    uV = ordered_nodes.copy()  # unvisited nodes
    sample_catchments = {}
    initial_len = len(uV)
    label = 1
    # Iterate through nodes from upstream to downstream
    for i in np.arange(initial_len - 1, -1, -1):
        node = uV[i]  # Node in network
        if node in target_nodes:
            sample_name = sample_names[np.ravel(np.where(target_nodes==node))][0]
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
            sample_catchments[label] = {"SampleName": sample_name,"catchment_nodes": unvis_up_nodes, "node": node}  # Append catchment to dictionary
            uV = fast_delete(uV, unvis_up_nodes)  # Remove visited nodes from array
            label+=1 # Update new node ID
    return sample_catchments


def get_subcatchment_map(model_grid, sample_catchment_dict):
    out = np.zeros(model_grid.shape).flatten()
    for label, catchment_dict in sample_catchment_dict.items():
        out[catchment_dict["catchment_nodes"]] = label
    return out.reshape(model_grid.shape)


def load_sample_data(path_to_file):
    data = pd.read_csv(path_to_file).to_numpy()
    sample_sites = data[:, 1:].astype(float)
    sample_names = data[:, 0].astype(str)
    return sample_names, sample_sites

def node_to_coord(model_grid,node):
    y,x=np.unravel_index(node,model_grid.shape)
    return(y*model_grid.dy,x*model_grid.dx)

def save_outputs(model_grid,sample_dict,catchment_map):
    node_IDs=[]
    sample_names = []
    xs = []
    ys = []
    for label,catchment_dict in sample_dict.items():
        node_IDs.append(label) # node ID
        sample_names.append(catchment_dict["SampleName"])
        x,y = node_to_coord(model_grid,catchment_dict["node"])
        xs.append(x)
        ys.append(y)

    outdf = pd.DataFrame({"Catchment ID": node_IDs,"SampleName":sample_names,"x":xs,"y":ys})
    outdf.to_csv("fitted_localities.csv",index=False)


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
    sub_catchment_dict = get_subcatchments(mg, target_nodes,sample_names)
    catchment_map = get_subcatchment_map(mg, sub_catchment_dict)
    save_outputs(mg,sub_catchment_dict,catchment_map)

    plt.title("Subcatchments: Samples aligned to drainage")
    plt.imshow(catchment_map, origin="lower", cmap="nipy_spectral")
    cb = plt.colorbar()
    cb.set_label("Catchment ID")
    plt.show()


if __name__ == "__main__":
    main()
