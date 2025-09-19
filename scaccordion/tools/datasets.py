import importlib.resources
#from .utils import parse_CrossTalkeR
import networkx as nx

# def load_peng2019():
#     """

#     Return a dict with the LR pairs from the PDAC patients Described in Peng 2019
    
#     """
#     # This is a stream-like object. If you want the actual info, call
#     # stream.read()
#     stream = pkg_resources.resource_filename("scaccordion", 'data/peng_ctker.Rds')
#     return parse_CrossTalkeR(stream)



def load_peng2019_metadata():
    """

    Return a dict with the LR pairs from the PDAC patients Described in Peng 2019
    
    """
    # This is a stream-like object. If you want the actual info, call
    # stream.read()
    import pandas as pd
    stream = importlib.resources.files("scaccordion.data")/ 'Peng_PDAC_metaacc.h5ad'
    return pd.read_hdf(stream)


def make_gl_ciclic_graph():
   import pandas as pd
   import networkx as nx
   stream = importlib.resources.files("scaccordion.data")/ 'glexample.hdf'
   return (nx.from_pandas_edgelist(pd.read_hdf(stream),create_using=nx.DiGraph, edge_attr=True))