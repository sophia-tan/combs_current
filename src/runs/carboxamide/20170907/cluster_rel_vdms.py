import sys                                                                   
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs

picklepath = '/Users/npolizzi/Projects/combs/results/carboxamide/rel_vdms/20170830/'
combs.cluster.Interactamer.get_centroids(picklepath, radius=0.2)
