"""Tools to draw figures from the MaBoSS results."""

import operator
import sys
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import networkx as nx
import matplotlib as mpl

def persistent_color(palette, state):

    reordered = " -- ".join(sorted(state.split(" -- ")))
    if state in palette:
        return palette[state]
    if reordered in palette:
        return palette[reordered]
    if state == "Others":
        color = "lightgray"
    else:
        count = len(palette)
        if "Others" not in palette:
            count += 1
        color = "C%d" % (count % 10)
    palette[state] = color
    palette[reordered] = color
    return color

def register_states_for_color(palette, collection):
    [persistent_color(palette, state) for state in collection]

def make_plot_trajectory(time_table, ax, palette, legend=True, error_table=None):
    register_states_for_color(palette, time_table.columns.values)
    color_list = [persistent_color(palette, idx) \
                    for idx in time_table.columns.values]
    time_table.plot(ax=ax, color=color_list, legend=legend, yerr=error_table)
    if legend:
        plt.legend(loc=(1.1, 0))


def plot_node_prob(time_table, ax, palette, legend=True, error_table=None):
    """Plot the probability of each node being up over time."""
    
    register_states_for_color(palette, time_table.columns.values)
    color_list = [persistent_color(palette, idx) \
                    for idx in time_table.columns.values]
    if len(time_table.columns) == 0:
        print("No data to plot !", file=sys.stderr)
        
    else:        
        if ax is None:
            _, ax = plt.subplots(1,1)
        
    
        time_table.plot(ax=ax, color=color_list, legend=legend, yerr=error_table)
        if not np.any(time_table):
            plt.ylim(0, 1)
        if legend:
            plt.legend(loc='upper right')
        return ax


def plot_piechart(plot_table, ax, palette, embed_labels=False, autopct=4, \
                    prob_cutoff=0.01, legend=True, nil_label=None):
    plot_line = plot_table.iloc[-1].rename("")  # Takes the last time point

    others = plot_line[plot_line <= prob_cutoff].sum()

    plot_line = plot_line[plot_line > prob_cutoff]
    plot_line.sort_values(ascending=False, inplace=True)

    if others:
        plot_line.at["Others"] = others

    register_states_for_color(palette, plot_line.index)

    plotting_labels = []
    color_list = []
    for value_index, value in enumerate(plot_line):
        state = plot_line.index.values[value_index]
        color_list.append(persistent_color(palette, state))
        if embed_labels and value >= 0.1:
            if state == "<nil>" and nil_label is not None:
                plotting_labels.append(nil_label)
            else:
                plotting_labels.append(state)
        else:
            plotting_labels.append("")

    opts = {}
    if autopct:
        cutoff = autopct if type(autopct) is not bool else 4
        opts.update(autopct=lambda p: '%1.1f%%' % p if p >= cutoff else "")
    else:
        opts.update(labeldistance=0.4)

    if others:
        plot_line = plot_line.rename({"Others": "Others (%1.2f%%)" % (others*100)})

    ax.pie(plot_line, labels=plotting_labels, radius=1.2,
           startangle=90, colors=color_list, **opts)
    ax.axis('equal')
    if legend:
        if nil_label is not None:
            ax.legend([nil_label if label == "<nil>" else label for label in plot_line.index.values], loc=(0.9, 0.2))
        else:
            ax.legend(plot_line.index.values, loc=(0.9, 0.2))


def plot_fix_point(table, ax, palette):
    """Plot a piechart representing the fixed point probability."""
    palette['no_fp'] = '#121212'
    prob_list = []
    color_list = []
    labels = []
    for i in range(len(table)):
        prob = table['Proba'][i]
        state = table['State'][i]
        prob_list.append(prob)
        color_list.append(persistent_color(palette, state))
        labels.append('FP '+str(i+1))
    prob_ns = 1 - sum(prob_list)
    if prob_ns > 0.01:
        prob_list.append(prob_ns)
        color_list.append(persistent_color(palette, "no_fp"))
        labels.append('no_fp')
    ax.pie(prob_list, labels=labels, colors=color_list)

def plot_observed_graph(table, ax, prune=True):

    if prune:
        pruned_table = table.copy()
        for i in range(len(table)):
            if table.iloc[i,:].sum() == 0 and table.iloc[:,i].sum() == 0:
                pruned_table.drop(table.columns[i], axis=1, inplace=True)
                pruned_table.drop(table.index[i], axis=0, inplace=True)
        table = pruned_table

    G = nx.from_pandas_adjacency(table,  create_using=nx.DiGraph())

    edge_colors = [edge['weight'] for _, edge in G.edges.items()]

    pos = nx.drawing.nx_agraph.graphviz_layout(G, prog="neato")
    # fig = plt.figure(figsize=(12,6), dpi=100)

    nx.draw(G, pos, with_labels=True, edgelist=[], node_size=0, ax=ax)
    cmap = mpl.colors.ListedColormap(plt.cm.Blues(np.linspace(0.2, 1, 100)))

    edges = nx.draw_networkx_edges(
        G,
        pos,
        arrowstyle="->",
        arrowsize=10,
        edgelist=G.edges().keys(),
        edge_color=edge_colors,
        edge_cmap=cmap,
        width=2,
        connectionstyle='arc3, rad = 0.1',
        ax=ax
    )

    pc = mpl.collections.PatchCollection(edges, cmap=cmap)
    pc.set_array(edge_colors)
    plt.colorbar(pc, ax=ax)

