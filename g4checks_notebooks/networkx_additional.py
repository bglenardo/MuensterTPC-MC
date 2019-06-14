import networkx as nx
from networkx.convert import _prep_create_using
from networkx.utils import not_implemented_for

# Additional functions for networkx drawing
##################################################
### Additional functions for networkx plotting ###
##################################################

#---------- Define graph from dataframe -------------
def nx_from_pandas_dataframe(df, source, target, edge_attr=None,
        create_using=None):
    """Return a graph from Pandas DataFrame.

    The Pandas DataFrame should contain at least two columns of node names and
    zero or more columns of node attributes. Each row will be processed as one
    edge instance.

    Note: This function iterates over DataFrame.values, which is not
    guaranteed to retain the data type across columns in the row. This is only
    a problem if your row is entirely numeric and a mix of ints and floats. In
    that case, all values will be returned as floats. See the
    DataFrame.iterrows documentation for an example.

    Parameters
    ----------
    df : Pandas DataFrame
        An edge list representation of a graph

    source : str or int
        A valid column name (string or iteger) for the source nodes (for the
        directed case).

    target : str or int
        A valid column name (string or iteger) for the target nodes (for the
        directed case).

    edge_attr : str or int, iterable, True
        A valid column name (str or integer) or list of column names that will
        be used to retrieve items from the row and add them to the graph as edge
        attributes. If `True`, all of the remaining columns will be added.

    create_using : NetworkX graph
        Use specified graph for result.  The default is Graph()

    See Also
    --------
    to_pandas_dataframe

    Examples
    --------
    Simple integer weights on edges:

    >>> import pandas as pd
    >>> import numpy as np
    >>> r = np.random.RandomState(seed=5)
    >>> ints = r.random_integers(1, 10, size=(3,2))
    >>> a = ['A', 'B', 'C']
    >>> b = ['D', 'A', 'E']
    >>> df = pd.DataFrame(ints, columns=['weight', 'cost'])
    >>> df[0] = a
    >>> df['b'] = b
    >>> df
       weight  cost  0  b
    0       4     7  A  D
    1       7     1  B  A
    2      10     9  C  E
    >>> G=nx.from_pandas_dataframe(df, 0, 'b', ['weight', 'cost'])
    >>> G['E']['C']['weight']
    10
    >>> G['E']['C']['cost']
    9
    """

    g = _prep_create_using(create_using)

    # Index of source and target
    src_i = df.columns.get_loc(source)
    tar_i = df.columns.get_loc(target)
    if edge_attr:
        # If all additional columns requested, build up a list of tuples
        # [(name, index),...]
        if edge_attr is True:
            # Create a list of all columns indices, ignore nodes
            edge_i = []
            for i, col in enumerate(df.columns):
                if col is not source and col is not target:
                    edge_i.append((col, i))
        # If a list or tuple of name is requested
        elif isinstance(edge_attr, (list, tuple)):
            edge_i = [(i, df.columns.get_loc(i)) for i in edge_attr]
        # If a string or int is passed
        else:
            edge_i = [(edge_attr, df.columns.get_loc(edge_attr)),]

        # Iteration on values returns the rows as Numpy arrays
        for row in df.values:
            g.add_edge(row[src_i], row[tar_i], {i:row[j] for i, j in edge_i})
    
    # If no column names are given, then just return the edges.
    else:
        for row in df.values:
            if row[src_i]<row[tar_i]:
                g.add_edge(row[src_i], row[tar_i])
            else:
                g.add_edge(row[tar_i], row[src_i])

    return g


#---------- Define hierarchic positions -------------
def nx_hierarchy_pos(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5, 
                  pos = None, parent = None):
    '''If there is a cycle that is reachable from root, then this will see infinite recursion.
       G: the graph
       root: the root node of current branch
       width: horizontal space allocated for this branch - avoids overlap with other branches
       vert_gap: gap between levels of hierarchy
       vert_loc: vertical location of root
       xcenter: horizontal location of root
       pos: a dict saying where all nodes go if they have been assigned
       parent: parent of this branch.'''
    if pos == None:
        pos = {root:(xcenter,vert_loc)}
    else:
        pos[root] = (xcenter, vert_loc)
    neighbors = list(G.neighbors(root))
    if parent != None:   #this should be removed for directed graphs.
        neighbors.remove(parent)  #if directed, then parent not in neighbors.
    if len(neighbors)!=0:
        dx = width/len(neighbors) 
        nextx = xcenter - width/2 - dx/2
        for neighbor in neighbors:
            nextx += dx
            pos = nx_hierarchy_pos(G,neighbor, width = dx, vert_gap = vert_gap, 
                                   vert_loc = vert_loc-vert_gap, xcenter=nextx, pos=pos, 
                                   parent = root)
    return pos


#---------- Draw additional labels on nodes -------------
def nx_draw_labels(G, pos, yoffset=0.035,
                         labels=None,
                         font_size=12,
                         font_color='k',
                         font_family='sans-serif',
                         font_weight='normal',
                         alpha=1.0,
                         bbox=None,
                         ax=None,
                         **kwds):
    """Draw node labels on the graph G.

    Parameters
    ----------
    G : graph
       A networkx graph

    pos : dictionary
       A dictionary with nodes as keys and positions as values.
       Positions should be sequences of length 2.

    labels : dictionary, optional (default=None)
       Node labels in a dictionary keyed by node of text labels

    font_size : int
       Font size for text labels (default=12)

    font_color : string
       Font color string (default='k' black)

    font_family : string
       Font family (default='sans-serif')

    font_weight : string
       Font weight (default='normal')

    alpha : float
       The text transparency (default=1.0)

    ax : Matplotlib Axes object, optional
       Draw the graph in the specified Matplotlib axes.

    Returns
    -------
    dict
        `dict` of labels keyed on the nodes

    Examples
    --------
    >>> G=nx.dodecahedral_graph()
    >>> labels=nx.draw_networkx_labels(G,pos=nx.spring_layout(G))

    Also see the NetworkX drawing examples at
    http://networkx.github.io/documentation/latest/gallery.html


    See Also
    --------
    draw()
    draw_networkx()
    draw_networkx_nodes()
    draw_networkx_edges()
    draw_networkx_edge_labels()
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cbook as cb
    except ImportError:
        raise ImportError("Matplotlib required for draw()")
    except RuntimeError:
        print("Matplotlib unable to open display")
        raise

    if ax is None:
        ax = plt.gca()

    if labels is None:
        labels = dict((n, n) for n in G.nodes())

    # set optional alignment
    horizontalalignment = kwds.get('horizontalalignment', 'center')
    verticalalignment = kwds.get('verticalalignment', 'center')
    
    if bbox is None:
            bbox = dict(boxstyle='round',
                        ec=(1.0, 1.0, 1.0),
                        fc=(1.0, 1.0, 1.0),
                        )

    text_items = {}  # there is no text collection so we'll fake one
    i=0
    for n, label in labels.items():
        (x, y) = pos[n]
        #print((x,y))
        if not cb.is_string_like(label):
            label = str(label)  # this will cause "1" and 1 to be labeled the same
        t = ax.text(x, y+yoffset,
                  label,
                  size=font_size,
                  color=font_color[i],
                  family=font_family,
                  weight=font_weight,
                  horizontalalignment=horizontalalignment,
                  verticalalignment=verticalalignment,
                  transform=ax.transData,
                  bbox=bbox,
                  clip_on=True,
                  )
        text_items[n] = t
        i=i+1

    return text_items